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

def matches(tag, seq_match_span, tag_match_span, allowed_errors):
    '''Determine the gap/error counts for a particular match'''
    # deal with case where tag match might be perfect, but extremely gappy, 
    # e.g. ACGTCGTGCGGA-------------------------ATC
    if tag_match_span.count('-') > allowed_errors or \
        seq_match_span.count('-') > allowed_errors:
        return 0, 0
    else:
        seq_array   = numpy.array(list(seq_match_span))
        tag_array   = numpy.array(list(tag_match_span))
        matches     = sum(seq_array == tag_array)
        error       = sum(seq_array != tag_array) + (len(tag) - \
            len(tag_match_span.replace('-','')))
        return matches, error

def align(seq, tags, allowed_errors):
    '''Alignment method for aligning tags with their respective
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

def mid_trim(tagged, tags, max_gap_char, mid_len, fuzzy, errors):
    """Remove the MID tag from the sequence read"""
    #if sequence.id == 'MID_No_Error_ATACGACGTA':
    #    pdb.set_trace()
    mid = leftLinker(tagged.read.sequence, tags, max_gap_char, mid_len, fuzzy,
        errors, gaps = True)
    if mid:
        tagged.mid, tagged.m_type, tagged.seq_match = mid[0],mid[1],mid[4]
        tagged.read = tagged.read.slice(mid[3],len(tagged.read.sequence), False)
    return tagged

def get_align_match_position(seq_match_span, start, stop):
    # slice faster than ''.startswith()
    if seq_match_span[0] == '-':
        start = start + seq_match_span.count('-')
    else:
        stop = stop - seq_match_span.count('-')
    return start, stop

def leftLinker(s, tags, max_gap_char, tag_len, fuzzy, errors, gaps=False):
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
    if not match and fuzzy:
        match = align(s[:max_gap_char + tag_len], tags, errors)
        # we can trim w/o regex
        if match:
            m_type = 'fuzzy'
            tag = match[0]
            seq_match = match[3]
            start, stop = get_align_match_position(match[3],match[4], match[5])
    if match:
        return tag, m_type, start, stop, seq_match
    else:
        return None

def rightLinker(s, tags, max_gap_char, tag_len, fuzzy, errors, gaps=False):
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
        match = align(s[-(tag_len + max_gap_char):], revtags, errors)
        # we can trim w/o regex
        if match:
            m_type = 'fuzzy'
            tag = match[0]
            seq_match = match[3]
            start, stop = get_align_match_position(match[3],match[4], match[5])
    if match:
        return revComp(tag), m_type, start, stop, seq_match
    else:
        return None

def both_tags_within_gaps(sequence, left, right, max_gap):
    if left[2] <= max_gap and \
                right[2] >= (len(sequence) - (len(right[0]) + max_gap)):
        return True

def linkerTrim(tagged, tags, max_gap_char, tag_len, fuzzy, errors):
    '''Use regular expression and (optionally) fuzzy string matching
    to locate and trim linkers from sequences'''

    left = leftLinker(tagged.read.sequence, tags, max_gap_char, tag_len, fuzzy, errors)
    right = rightLinker(tagged.read.sequence, tags, max_gap_char, tag_len, fuzzy, errors)
    
    # we can have 3 types of matches - tags on left and right sides,
    # tags on left side only, tags on right side only, mismatching tags
    # and no tags at all
    if left is not None \
            and right is not None \
            and left[0] == right[0] and \
            both_tags_within_gaps(tagged.read.sequence, left, right, max_gap_char):
        # trim the read
        tagged.read = tagged.read.slice(left[3], right[2], False)
        # left and right are identical so largely pass back the left
        # info... except for m_type which can be a combination
        tagged.l_tag, tagged.l_seq_match = left[0], left[4]
        tagged.l_m_type = "{}-{}-both".format(left[1], right[1])
        tagged.l_critter = tags[tagged.l_tag]

    elif left and right \
            and left[0] != right[0] \
            and both_tags_within_gaps(tagged.read.sequence, left, right, max_gap_char):
        # these are no good.  check for within gaps
        # to make sure it's not a spurious match
        tagged.trimmed = None
        tagged.l_tag, tagged.l_seq_match = None, None
        tagged.l_m_type = "tag-mismatch"
        tagged.l_critter = None

    elif left and left[2] <= max_gap_char:
        tagged.read = tagged.read.slice(left[3], len(tagged.read.sequence), False)
        tagged.l_tag, tagged.l_seq_match = left[0], left[4]
        tagged.l_seq_match = "{}-left".format(left[1])
        tagged.l_critter = tags[tagged.l_tag]

    elif right and right[2] >= (len(s) - (len(right[0]) + max_gap_char)):
        tagged.read = tagged.read.slice(0, right[2], False)
        tagged.l_tag, tagged.l_seq_match = right[0], right[4]
        tagged.l_m_type = "{}-right".format(right[1])
        tagged.l_critter = tags[tagged.l_tag]

    else:
        trimmed = None
        tagged.l_tag, tagged.l_m_type, tagged.l_seq_match = None, None, None
        tagged.l_critter = None
    
    return tagged

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
        match = align(s, all_tags, 1)
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
    return sum([1 for line in open(input, 'rU') if line.startswith('>')])

def singleproc(job, results, params):
    for sequence in job:
        # set tags = empty
        tags = None
        # for now, we'll keep this here
        tagged = Tagged(sequence)
        if params.qual_trim:
            #pdb.set_trace()
            tagged.read = tagged.read.trim(params.min_qual, False)
        if params.mid_trim:
            tagged = mid_trim(tagged, params.tags, params.mid_gap, \
                    params.mid_len, params.fuzzy, params.allowed_errors)
            if tagged.mid:
                tagged.reverse_mid = params.reverse_mid[tagged.mid]
        if params.linker_trim:
            # if we never trimmed the MID, then the tags are in 
            # the param object
            if not params.mid_trim:
                tags = params.tags
            # if we're tagging hierarchically, then we should get a
            # mid.
            elif tagged.mid is not None:
                # reduce out list of possible tags to only those that
                # that go w/ this particular MID
                tags = params.tags[tagged.mid]
            else:
                tags = None
            if tags:
                linker = linkerTrim(tagged, tags, params.linker_gap,
                            params.linker_len, params.fuzzy, params.allowed_errors)
                if tagged.l_tag:
                    tagged.reverse_linker = params.reverse_linkers[tagged.l_tag]
        if tagged:
            results.put(tagged)
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
            o = results.get()
            print o.mid, o.l_tag
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
        for o in results:
            print o.read.identifier, o.m_type, o.l_m_type
    
    print '\n'
    cur.close()
    conn.close()
    end_time = time.time()
    print 'Ended: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(end_time))
    print '\nTime for execution: ', (end_time - start_time)/60, 'minutes'

if __name__ == '__main__':
    main()
