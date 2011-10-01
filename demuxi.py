#!/usr/bin/env python
# encoding: utf-8
"""
linker.py

Copyright (c) 2009-2010 Brant C. Faircloth.  All rights reserved.

Provides:
    - parsing and error correction of hierarchically tagged
        next generation sequence reads
        
Requires:
    - Python > 2.6.x
    - MySQL
    - MySQLdb
    - Biopython
    - progress.py
    - numpy

Usage:

    `python linker.py --configuration=linker-py.conf`

"""

import pdb

import os
import sys
import re
import time
import numpy
import string
import cPickle
import sqlite3
import argparse
import progress
import itertools
import ConfigParser

from multiprocessing import Process, Queue, JoinableQueue, cpu_count

from tools.sequence.fasta import FastaQualityReader
from demuxi.lib import FullPaths, ListQueue, Tagged, Parameters

from Bio import pairwise2
from Bio.SeqIO import QualityIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import SingleLetterAlphabet

def motd():
    '''Startup info'''
    motd = '''
    ##############################################################
    #                     linker.py                              #
    # Provides:                                                  #
    #   - parsing and error correction of hierarchically tagged  #
    #     next generation sequence reads                         #
    #                                                            #
    #                                                            #
    # Copyright (c) 2009-2010 Brant C. Faircloth                 #
    # 621 Charles E. Young Drive                                 #
    # University of California, Los Angeles, 90095, USA          #
    ##############################################################\n
    '''
    print motd

def revComp(seq):
    '''Return reverse complement of seq'''
    bases = string.maketrans('AGCTagct','TCGAtcga')
    # translate it, reverse, return
    return seq.translate(bases)[::-1]

def revCompTags(tags):
    '''Return the reverse complements of a tag dictionary'''
    revTags = {}
    for tag in tags:
        revTags[revComp(tag)] = tags[tag]
    return revTags

def trim(sequence, left=None, right=None):
    '''Trim a given sequence given left and right offsets'''
    if left and right:
        sequence = sequence[left:right]
    elif left:
        sequence = sequence[left:]
    elif right:
        sequence = sequence[:right]
    return sequence

def matches(tag, seq_match_span, tag_match_span, allowed_errors):
    '''Determine the gap/error counts for a particular match'''
    # deal with case where tag match might be perfect, but extremely gappy, 
    # e.g. ACGTCGTGCGGA-------------------------ATC
    if tag_match_span.count('-') > allowed_errors or \
        seq_match_span.count('-') > allowed_errors:
        return 0, 0
    else:
        #pdb.set_trace()
        seq_array   = numpy.array(list(seq_match_span))
        tag_array   = numpy.array(list(tag_match_span))
        matches     = sum(seq_array == tag_array)
        error       = sum(seq_array != tag_array) + (len(tag) - \
            len(tag_match_span.replace('-','')))
        # I didn't like the way that the original method at
        # http://github.com/chapmanb/bcbb/tree/master treats gaps 
        return matches, error

def smithWaterman(seq, tags, allowed_errors):
    '''Smith-Waterman alignment method for aligning tags with their respective
    sequences.  Only called when regular expression matching patterns fail.
    Inspired by http://github.com/chapmanb/bcbb/tree/master'''
    high_score = {'tag':None, 'seq_match':None, 'mid_match':None, 'score':None, 
        'start':None, 'end':None, 'matches':None, 'errors':allowed_errors}
    for tag in tags:
        seq_match, tag_match, score, start, end = pairwise2.align.localms(seq, 
        tag, 5.0, -4.0, -9.0, -0.5, one_alignment_only=True)[0]
        seq_match_span  = seq_match[start:end]
        tag_match_span  = tag_match[start:end]
        match, errors   = matches(tag, seq_match_span, tag_match_span, allowed_errors)
        if match >= len(tag)-allowed_errors and match > high_score['matches'] \
            and errors <= high_score['errors']:
            high_score['tag'] = tag
            high_score['seq_match'] = seq_match
            high_score['tag_match'] = tag_match
            high_score['score'] = score
            high_score['start'] = start
            high_score['end'] = end
            high_score['matches'] = match
            high_score['seq_match_span'] = seq_match_span
            high_score['errors'] = errors
    if high_score['matches']:
        return high_score['tag'], high_score['matches'], \
        high_score['seq_match'], high_score['seq_match_span'], \
        high_score['start'], high_score['end']
    else:
        return None

def qualTrimming(sequence, min_score=10):
    '''Remove ambiguous bases from 5' and 3' sequence ends'''
    s = str(sequence.seq)
    sl = list(s)
    for q in enumerate(sequence.letter_annotations["phred_quality"]):
        if q[1] < min_score:
            sl[q[0]] = 'N'
    s = ''.join(sl)
    # find runs of ambiguous bases at 5' and 3' ends
    left_re, right_re = re.compile('^N+'),re.compile('N+$')
    left_trim, right_trim = re.search(left_re, s), re.search(right_re, s)
    if left_trim:
        left_trim = left_trim.end()
    if right_trim:
        right_trim = right_trim.end()
    return trim(sequence, left_trim, right_trim)

def mid_trim(tagged, tags, max_gap_char, mid_len, fuzzy, errors):
    """Remove the MID tag from the sequence read"""
    #if sequence.id == 'MID_No_Error_ATACGACGTA':
    #    pdb.set_trace()
    mid = leftLinker(tagged.read.sequence, tags, max_gap_char, mid_len, fuzzy,
        errors, gaps = True)
    if mid:
        tagged.mid, tagged.mid_match, tagged_mid_type = mid[0],mid[1],mid[4]
        tagged.read = tagged.read.slice(mid[3],len(tagged.read.sequence), False)
        return tagged
    else:
        return None

def SWMatchPos(seq_match_span, start, stop):
    # slice faster than ''.startswith()
    if seq_match_span[0] == '-':
        start = start + seq_match_span.count('-')
    else:
        stop = stop - seq_match_span.count('-')
    return start, stop

def leftLinker(s, tags, max_gap_char, mid_len, fuzzy, errors, gaps=False):
    '''Matching methods for left linker - regex first, followed by fuzzy (SW)
    alignment, if the option is passed'''
    for tag in tags:
        if gaps:
            regex = re.compile(('^%s') % (tag))
        else:
            regex = re.compile(('^[acgtnACGTN]{0,%s}%s') % (max_gap_char, tag))
        match = regex.search(s)
        if match:
            m_type = 'regex'
            start, stop = match.start(), match.end()
            # by default, this is true
            seq_match = tag
            break
    #if s == 'ACCTCGTGCGGAATCGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAG':
    #    pdb.set_trace()
    if not match and fuzzy:
        match = smithWaterman(s[:max_gap_char + mid_len], tags, errors)
        # we can trim w/o regex
        if match:
            m_type = 'fuzzy'
            tag = match[0]
            seq_match = match[3]
            start, stop = SWMatchPos(match[3],match[4], match[5])
    if match:
        return tag, m_type, start, stop, seq_match
    else:
        return None

def rightLinker(s, tags, max_gap_char, mid_len, fuzzy, errors, gaps=False):
    '''Mathing methods for right linker - regex first, followed by fuzzy (SW)
    alignment, if the option is passed'''
    revtags = revCompTags(tags)
    for tag in revtags:
        if gaps:
            regex = re.compile(('%s$') % (tag))
        else:
            regex = re.compile(('%s[acgtnACGTN]{0,%s}$') % (tag, max_gap_char))
        match = regex.search(s)
        if match:
            m_type = 'regex'
            start, stop = match.start(), match.end()
            # by default, this is true
            seq_match = tag
            break
    if not match and fuzzy:
        match = smithWaterman(s[-(mid_len + max_gap_char):], revtags, errors)
        # we can trim w/o regex
        if match:
            m_type = 'fuzzy'
            tag = match[0]
            seq_match = match[3]
            start, stop = SWMatchPos(match[3],match[4], match[5])
    if match:
        return revComp(tag), m_type, start, stop, seq_match
    else:
        return None

def linkerTrim(sequence, tags, max_gap_char, mid_len, fuzzy, errors):
    '''Use regular expression and (optionally) fuzzy string matching
    to locate and trim linkers from sequences'''
    m_type  = False
    left    = leftLinker(sequence, tags, max_gap_char, mid_len, fuzzy, errors)
    right   = rightLinker(sequence, tags, max_gap_char, mid_len, fuzzy, errors)
    s = str(sequence.seq)
    if left and right and left[0] == right[0]:
        # we can have lots of conditional matches here
        if left[2] <= max_gap_char and right[2] >= (len(s) - (len(right[0]) +\
        max_gap_char)):
            trimmed = trim(sequence, left[3], right[2])
            # left and right are identical so largely pass back the left
            # info... except for m_type which can be a combination
            tag, m_type, seq_match = left[0], left[1]+'-'+right[1]+'-both', \
            left[4]
        else:
            pass
    elif left and right and left[0] != right[0]:
        # flag
        if left[2] <= max_gap_char and right[2] >= (len(s) - (len(right[0]) +\
        max_gap_char)):
            trimmed = None
            tag, m_type, seq_match = None, 'tag-mismatch', None
    elif left:
        if left[2] <= max_gap_char:
            trimmed = trim(sequence, left[3])
            tag, m_type, seq_match = left[0], left[1]+'-left', left[4]
        else:
            # flag
            pass
    elif right:
        if right[2] >= (len(s) - (len(right[0]) + max_gap_char)):
            trimmed = trim(sequence, None, right[2])
            tag, m_type, seq_match = right[0], right[1]+'-right', right[4]
        else:
            # flag
            pass
    if m_type:
        try:
            return tag, trimmed, seq_match, tags[tag], m_type
        except:
            return tag, trimmed, seq_match, None, m_type
    else:
        return None

def reverse(items, null=False):
    '''build a reverse dictionary from a list of tuples'''
    l = []
    if null:
        items += ((None, None),)
    for i in items:
        t = (i[1],i[0])
        l.append(t)
    return dict(l)

def create_db_and_new_tables(db_name):
    conn = sqlite3.connect(db_name)
    cur = conn.cursor()
    cur.execute("PRAGMA foreign_keys = ON")
    try:
        cur.execute('''CREATE TABLE tags (
            id integer PRIMARY KEY AUTOINCREMENT,
            name text,
            mid text,
            mid_seq text,
            mid_match text,
            mid_method text,
            linker text,
            linker_seq text,
            linker_match text,
            linker_method text,
            cluster text,
            concat_seq text,
            concat_match text,
            concat_method text,
            n_count integer,
            untrimmed_len integer,
            seq_trimmed text,
            trimmed_len text,
            record blob)
            ''')
        cur.execute('''CREATE TABLE sequence (
            id INTEGER,
            untrimmed_len, integer,
            trimmed_len integer,
            seq_trimmed text,
            FOREIGN KEY(id) REFERENCES tags(id) DEFERRABLE INITIALLY
            DEFERRED)''')
        cur.execute("CREATE INDEX idx_sequence_cluster on tags(cluster)")
    except sqlite3.OperationalError, e:
        #pdb.set_trace()
        if 'already exists' in e[0]:
            answer = raw_input("\n\tDatabase already exists.  Overwrite [Y/n]? ")
            #pdb.set_trace()
            if answer == "Y" or answer == "YES":
                os.remove(db_name)
                conn, cur = create_db_and_new_tables(db_name)
            else:
                sys.exit()
        else:
            raise sqlite3.OperationalError, e
    return conn, cur

def concatCheck(sequence, all_tags, all_tags_regex, reverse_linkers, **kwargs):
    '''Check screened sequence for the presence of concatemers by scanning 
    for all possible tags - after the 5' and 3' tags have been removed'''
    s = str(sequence.seq)
    m_type = None
    # do either/or to try and keep speed up, somewhat
    #if not kwargs['fuzzy']:
    #pdb.set_trace()
    for tag in all_tags_regex:
        match = re.search(tag, s)
        if match:
            tag = tag.pattern
            m_type = 'regex-concat'
            seq_match = tag
            break
    if not match and ['fuzzy']:
    #else:
        match = smithWaterman(s, all_tags, 1)
        # we can trim w/o regex
        if match:
            tag = match[0]
            m_type = 'fuzzy-concat'
            seq_match = match[3]
    if m_type:
        return tag, m_type, seq_match
    else:
        return None, None, None

def get_sequence_count(input):
    '''Determine the number of sequence reads in the input'''
    handle = open(input, 'rU')
    lines = handle.read().count('>')
    handle.close()
    return lines

def singleproc(job, results, params):
    for sequence in job:
        # for now, we'll keep this here
        tagged = Tagged(sequence)
        if params.qual_trim:
            #pdb.set_trace()
            tagged.read = tagged.read.trim(params.min_qual, False)
        if params.mid_trim:
            tagged = mid_trim(tagged, params.tags, params.mid_gap, \
                    params.mid_len, params.fuzzy, params.allowed_errors)
            #if mid:
            #    tag.mid, tag.sequence, tag.seq_match, tag.m_type,
            #        tag.reverse_mid, tags =
            #        mid[0], mid[1], mid[2], mid[3],
            #        params.reverse_mid[seqRecord.mid], params.tags[seqRecord.mid]

        if tagged:
            results.put(tagged)
            pdb.set_trace()
    return results

def multiproc(jobs, results, params):
    """locate linker sequences in a read, returning a record object"""
    while True:
        job = jobs.get()
        if job is None:
            break
        _ = singleproc(job, results, params)
    """#pdb.set_trace()
    tags = params.tags
    if params.midTrim:
        # search on 5' (left) end for MID
        mid = midTrim(seqRecord.sequence, params.tags, params.midGap, params.mid_len, params.fuzzy, params.allowed_errors)
        if mid:
            # if MID, search for exact matches (for and revcomp) on Linker
            # provided no exact matches, use fuzzy matching (Smith-Waterman) +
            # error correction to find Linker
            seqRecord.mid           = mid[0]
            seqRecord.sequence      = mid[1]
            seqRecord.seq_match     = mid[2]
            seqRecord.m_type        = mid[3]
            seqRecord.reverse_mid   = params.reverse_mid[seqRecord.mid]
            tags                    = params.tags[seqRecord.mid]
    #pdb.set_trace()
    if params.linkerTrim:
        linker = linkerTrim(seqRecord.sequence, tags, params.linkerGap,
            params.mid_len, params.fuzzy, params.allowed_errors)
        if linker:
            if linker[0]:
                seqRecord.l_tag             = linker[0]
                seqRecord.sequence          = linker[1]
                seqRecord.l_seq_match       = linker[2]
                seqRecord.l_critter         = linker[3]
                seqRecord.l_m_type          = linker[4]
                seqRecord.reverse_linker    = params.reverse_linkers[seqRecord.l_tag]
            # deal with tag-mismatch
            if not linker[0] and linker[4]:
                seqRecord.l_m_type          = linker[4]
    # check for concatemers
    if params.concat:
        if l_trimmed and len(l_trimmed.seq) > 0:
            concat_tag, concat_type, concat_seq_match = concatCheck(l_trimmed, 
                all_tags, all_tags_regex, reverse_linkers, fuzzy=params.fuzzy)
        else:
            concat_tag, concat_type, concat_seq_match = None, None, None
    else:
        concat_tag, concat_type, concat_seq_match = None, None, None
    # pickle the sequence record, so we can store it as a BLOB in MySQL, we
    # can thus resurrect it as a sequence object when we need it next.
    sequence_pickle = cPickle.dumps(seqRecord.sequence,1)
    cur.execute('''INSERT INTO sequence (name, mid, mid_seq, mid_match, 
        mid_method, linker, linker_seq, linker_match, linker_method, cluster, 
        concat_seq, concat_match, concat_method, n_count, untrimmed_len, 
        seq_trimmed, trimmed_len, record) 
        VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)''', 
        (seqRecord.sequence.id, seqRecord.reverse_mid, seqRecord.mid, \
        seqRecord.seq_match, seqRecord.m_type, seqRecord.reverse_linker, \
        seqRecord.l_tag, seqRecord.l_seq_match, seqRecord.l_m_type, \
        seqRecord.l_critter, seqRecord.concat_tag, \
        seqRecord.concat_seq_match, seqRecord.concat_type, seqRecord.nCount, \
        len(seqRecord.unmod.seq), seqRecord.sequence.seq, \
        len(seqRecord.sequence.seq), sequence_pickle))
    #pdb.set_trace()
    cur.close()
    conn.commit()
    # keep our connection load low
    conn.close()
    return
    """

def get_args():
    """get arguments (config file location)"""
    parser = argparse.ArgumentParser(description = "demuxi.py:  sequence " + \
        "demultiplexing for hierarchically-tagged samples")
    parser.add_argument('config', help="The input configuration file",
            action=FullPaths)
    return parser.parse_args()

def split_reads_into_groups(fasta, qual, num_reads, num_procs):
    reads = FastaQualityReader(fasta, qual)
    job_size = num_reads/num_procs
    print "Parsing reads into groups of {} reads".format(job_size)
    i = iter(reads)
    chunk = list(itertools.islice(i,job_size))
    while chunk:
        yield chunk
        chunk = list(itertools.islice(i, job_size))

def main():
    '''Main loop'''
    start_time = time.time()
    motd()
    args = get_args()
    print 'Started: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(start_time))
    # build our configuration object w/ input params
    conf = ConfigParser.ConfigParser()
    conf.read(args.config)
    params = Parameters(conf)
    # create the db and tables, returninn connection
    # and cursor
    conn, cur = create_db_and_new_tables(params.db)
    # get read count of input
    num_reads = get_sequence_count(params.fasta)
    if params.num_procs > 1:
        # split reads into generator objects based on equal
        # split across cores
        work = split_reads_into_groups(params.fasta, params.quality,
                num_reads, params.num_procs)
    else:
        work = FastaQualityReader(params.fasta, params.quality)
    # MULTICORE
    if params.num_procs > 1:
        jobs = Queue()
        results = JoinableQueue()
        # We're stacking groups of jobs on the work
        # Queue, conceivably to save the overhead of
        # placing them on there one-by-one.
        for unit in work:
            jobs.put(unit)
        # setup the processes for the jobs
        print "Starting {} workers".format(params.num_procs)
        # start the worker processes
        [Process(target = multiproc, args=(jobs, results, params)).start()
            for i in xrange(params.num_procs)]
        # we're putting single results on the results Queue so
        # that the db can (in theory) consume them at
        # a rather consistent rate rather than in spurts
        #for unit in xrange(num_reads):
        for unit in xrange(18):
            #enter_to_db(results.get())
            r = results.get()
            print r.unmod
            results.task_done()
        # make sure we put None at end of Queue
        # in an amount equiv. to num_procs
        for unit in xrange(params.num_procs):
            jobs.put(None)
        # join the results, so that they can finish
        results.join()
        # close up our queues
        jobs.close()
        results.close()
    # SINGLECORE
    else:
        results = ListQueue()
        singleproc(work, results, params)
        for r in results:
            print r.unmod

    
    print '\n'
    cur.close()
    conn.close()
    end_time = time.time()
    print 'Ended: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(end_time))
    print '\nTime for execution: ', (end_time - start_time)/60, 'minutes'

if __name__ == '__main__':
    main()
