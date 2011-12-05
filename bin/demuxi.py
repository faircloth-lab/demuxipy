#!/usr/bin/env python
# encoding: utf-8

"""
File: demuxi.py
Author: Brant Faircloth

Created by Brant Faircloth on 05 October 2011 14:10 PDT (-0700)
Copyright (c) 2011 Brant C. Faircloth. All rights reserved.

Description: 

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
import itertools
import ConfigParser

from multiprocessing import Process, Queue, JoinableQueue, cpu_count

from seqtools.sequence.fastq import FastqReader
from seqtools.sequence.fasta import FastaReader
from seqtools.sequence.fasta import FastaQualityReader
from seqtools.sequence.transform import DNA_reverse_complement

from demuxipy import db
from demuxipy import pairwise2
from demuxipy.lib import FullPaths, ListQueue, Tagged, Parameters


def motd():
    """Print a welcome message """
    motd = """
    ###############################################################
    #                       demuxi.py                             #
    #                                                             #
    # Demultiplexing of hierarchically tagged, massively parallel #
    # DNA sequence reads                                          #
    #                                                             #
    #                                                             #
    # Copyright (c) 2009-2011 Brant C. Faircloth                  #
    #                                                             #
    #                                                             #
    # Ecology and Evolutionary Biology                            #
    # 621 Charles E. Young Drive                                  #
    # University of California, Los Angeles, 90095, USA           #
    ###############################################################\n
    """
    print motd

def matches(tag, seq_match_span, tag_match_span, allowed_errors):
    """Determine the gap/error counts for a particular match"""
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
    """Alignment method for aligning tags with their respective
    sequences.  Only called when regular expression matching patterns fail.
    Inspired by http://github.com/chapmanb/bcbb/tree/master"""
    high_score = {'tag':None, 'seq_match':None, 'mid_match':None, 'score':None, 
        'start':None, 'end':None, 'matches':None, 'errors':allowed_errors}
    for tag in tags:
        result = pairwise2.align.localms(seq,
            tag, 5.0, -4.0, -9.0, -0.5, 
            one_alignment_only=True)
        if result:
            seq_match, tag_match, score, start, end = result[0]
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

def get_align_match_position(seq_match_span, start, stop):
    # slice faster than ''.startswith()
    if seq_match_span[0] == '-':
        start = start + seq_match_span.count('-')
    else:
        stop = stop - seq_match_span.count('-')
    return start, stop

def find_left_linker(s, tag_regexes, tag_strings, max_gap_char, tag_len, fuzzy, errors):
    """Matching methods for left linker - regex first, followed by fuzzy (SW)
    alignment, if the option is passed"""
    match = None
    for regex in tag_regexes:
        match = regex.search(s)
        if match is not None:
            m_type = 'regex'
            start, stop = match.start(), match.end()
            # by default, this is true
            tag_matched = regex.pattern.split('}')[1]
            seq_matched = s[start:stop]
            break
    if match is None and fuzzy:
        match = align(s[:max_gap_char + tag_len], tag_strings, errors)
        # we can trim w/o regex
        if match:
            m_type = 'fuzzy'
            tag_matched = match[0]
            seq_matched = match[3]
            start, stop = get_align_match_position(match[3],match[4], match[5])
    if match:
        return tag_matched, m_type, start, stop, seq_matched
    else:
        return None

def find_right_linker(s, tag_regexes, tag_strings, max_gap_char, tag_len,
        fuzzy, errors, tagged):
    """Matching methods for right linker - regex first, followed by fuzzy (SW)
    alignment, if the option is passed"""
    #pdb.set_trace()
    match = None
    for regex in tag_regexes:
        match = regex.search(s)
        if match is not None:
            m_type = 'regex'
            start, stop = match.start(), match.end()
            # by default, this is true
            tag_matched = regex.pattern.split('[')[0]
            seq_matched = s[start:stop]
            break
    if match is None and fuzzy:
        match = align(s[-(tag_len + max_gap_char):], tag_strings, errors)
        # we can trim w/o regex
        if match:
            start_of_slice = len(s) - (tag_len + max_gap_char)
            m_type = 'fuzzy'
            tag_matched = match[0]
            seq_matched = match[3]
            start, stop = get_align_match_position(match[3],match[4], match[5])
            start, stop = start + start_of_slice, stop + start_of_slice
    if match:
        return DNA_reverse_complement(tag_matched), m_type, start, stop, DNA_reverse_complement(seq_matched)
    else:
        return None

def mid_trim(tagged, params):
    """Remove the MID tag from the sequence read"""
    #if sequence.id == 'MID_No_Error_ATACGACGTA':
    #    pdb.set_trace()
    mid = find_left_linker(tagged.read.sequence,
            params.sequence_tags.mids['forward_regex'], 
            params.sequence_tags.mids['forward_string'],
            params.sequence_tags.mid_gap,
            params.sequence_tags.mid_len,
            params.mid_fuzzy,
            params.mid_allowed_errors)
    if mid:
        tagged.mid, tagged.m_type, tagged.seq_match = mid[0],mid[1],mid[4]
        tagged.read = tagged.read.slice(mid[3],len(tagged.read.sequence), False)
        tagged.mid_name = params.sequence_tags.reverse_mid_lookup[tagged.mid]
    return tagged

def find_and_trim_linkers(tagged, params):
    """Use regular expression and (optionally) fuzzy string matching
    to locate and trim linkers from sequences"""
    #if '>MID15_NoError_SimpleX1_NoError_FandR' in tagged.read.identifier:
    #    pdb.set_trace()
    #pdb.set_trace()
    left = find_left_linker(tagged.read.sequence,
            params.sequence_tags.linkers[str(tagged.mid)]['forward_regex'],
            params.sequence_tags.linkers[str(tagged.mid)]['forward_string'],
            params.sequence_tags.linker_gap,
            params.sequence_tags.linker_len,
            params.linker_fuzzy,
            params.linker_allowed_errors)

    right = find_right_linker(tagged.read.sequence,
            params.sequence_tags.linkers[str(tagged.mid)]['reverse_regex'],
            params.sequence_tags.linkers[str(tagged.mid)]['reverse_string'],
            params.sequence_tags.linker_gap,
            params.sequence_tags.linker_len,
            params.linker_fuzzy,
            params.linker_allowed_errors, tagged)

    max_gap_char = params.sequence_tags.linker_gap

    # we can have 5 types of matches - tags on left and right sides,
    # tags on left side only, tags on right side only, mismatching tags
    # and no tags at all.  We take care of matching position (i,e. we
    # want only matches at the ends) by regex and slicing in the
    # search methods above.  If for some reason you want to check
    # for concatemers, then turn that function on.

    if left is not None \
            and right is not None \
            and left[0] == right[0]:
        # trim the read
        tagged.read = tagged.read.slice(
                left[3], right[2],False
            )
        # left and right are identical so largely pass back the left
        # info... except for m_type which can be a combination
        tagged.l_tag, tagged.l_seq_match = left[0], left[4]
        tagged.l_m_type = "{}-{}-both".format(left[1], right[1])

    elif left and right \
            and left[0] != right[0]:
        # these are no good.  check for within gaps
        # to make sure it's not a spurious match
        #pdb.set_trace()
        tagged.trimmed = None
        tagged.l_tag, tagged.l_seq_match = None, None
        tagged.l_m_type = "tag-mismatch"
        tagged.l_critter = None
        tagged.l_name = None

    elif left and left[2] <= max_gap_char:
        tagged.read = tagged.read.slice(left[3], len(tagged.read.sequence), False)
        tagged.l_tag, tagged.l_seq_match = left[0], left[4]
        tagged.l_m_type = "{}-left".format(left[1])

    elif right and right[2] >= (len(tagged.read.sequence) - (len(right[0]) + max_gap_char)):
        tagged.read = tagged.read.slice(0, right[2], False)
        tagged.l_tag, tagged.l_seq_match = right[0], right[4]
        tagged.l_m_type = "{}-right".format(right[1])

    else:
        trimmed = None
        tagged.l_tag, tagged.l_m_type, tagged.l_seq_match = None, None, None
        tagged.l_critter = None
        tagged.l_name = None
    
    if tagged.l_tag:
        tagged.l_critter = params.sequence_tags.cluster_map[str(tagged.mid)][str(tagged.l_tag)]
        tagged.l_name = params.sequence_tags.reverse_linker_lookup[tagged.l_tag]

    return tagged

def concat_check(tagged, params):
    """Check screened sequence for the presence of concatemers by scanning 
    for all possible tags - after the 5' and 3' tags have been removed"""
    s = tagged.read.sequence
    m_type = None
    #pdb.set_trace()
    for tag in params.sequence_tags.all_tags[str(tagged.mid)]['regex']:
        match = tag.search(s)
        if match:
            tagged.concat_tag= tag.pattern
            tagged.concat_m_type = "regex-concat"
            tagged.concat_seq_match = tagged.read.sequence[match.start():match.end()]
            break
    if match is None and params.concat_fuzzy:
        match = align(s,
                params.sequence_tags.all_tags[str(tagged.mid)]['string'], 
                params.concat_allowed_errors
            )
        if match:
            tagged.concat_tag = match[0]
            tagged.concat_m_type = "fuzzy-concat"
            tagged.concat_seq_match = match[3]
    return tagged

def progress(count, interval, big_interval):
    """give a rudimentary indication of progress"""
    if count % big_interval == 0:
        sys.stdout.write('%')
        sys.stdout.flush()
    elif count % interval == 0:
        sys.stdout.write('.')
        sys.stdout.flush()

def singleproc(job, results, params, interval = 1000, big_interval = 10000):
    count = 0
    for sequence in job:
        # for now, we'll keep this here
        tagged = Tagged(sequence)
        # trim
        if params.qual_trim:
            #pdb.set_trace()
            tagged.read = tagged.read.trim(params.min_qual, False)

        # check for MIDs
        if params.mid_trim:
            tagged = mid_trim(tagged, params)

        # check for linkers
        if (params.linker_trim and tagged.mid and params.search == 'MidLinkerGroups') or \
                (params.linker_trim and params.search == 'LinkerGroups'):
            tagged = find_and_trim_linkers(tagged, params)

        # check for concatemers
        if (params.concat_check and len(tagged.read.sequence) > 0) and \
                ((tagged.mid and tagged.l_tag and params.search == 'MidLinkerGroups') or \
                (tagged.mid and params.search == 'MidGroups') or \
                (tagged.l_tag and params.search == 'LinkerGroups')):
            tagged = concat_check(tagged, params)

        count += 1
        progress(count, interval, big_interval)
        results.put(tagged)
    return results

def multiproc(jobs, results, params):
    """locate linker sequences in a read, returning a record object"""
    while True:
        job = jobs.get()
        if job is None:
            break
        _ = singleproc(job, results, params)

def get_args():
    """get arguments (config file location)"""
    parser = argparse.ArgumentParser(description = "demuxi.py:  sequence " + \
        "demultiplexing for hierarchically-tagged samples")
    parser.add_argument('config', help="The input configuration file",
            action=FullPaths)
    return parser.parse_args()

def get_sequence_count(input, kind):
    """Determine the number of sequence reads in the input"""
    if kind == 'fasta':
        return sum([1 for line in open(input, 'rU') if line.startswith('>')])
    elif kind == 'fastq':
        return sum([1 for line in open(input, 'rU') if line.startswith('@')])

def split_fasta_reads_into_groups(reads, num_reads, num_procs):
    job_size = num_reads/num_procs
    print "Parsing reads into groups of {} reads".format(job_size)
    i = iter(reads)
    chunk = list(itertools.islice(i,job_size))
    while chunk:
        yield chunk
        chunk = list(itertools.islice(i, job_size))

def get_work(params):
    if params.fasta and params.quality:
        reads = FastaQualityReader(params.fasta, params.quality)
        # get read count of input
        num_reads = get_sequence_count(params.fasta, 'fasta')
        if params.num_procs > 1:
            # split reads into generator objects based on equal
            # split across cores
            work = split_fasta_reads_into_groups(reads,
                    num_reads, params.num_procs)
        else:
            work = FastaQualityReader(params.fasta, params.quality)
    elif params.fasta and params.quality == None:
        reads = FastaReader(params.fasta)
        # get read count of input
        num_reads = get_sequence_count(params.fasta, 'fasta')
        if params.num_procs > 1:
            # split reads into generator objects based on equal
            # split across cores
            work = split_fasta_reads_into_groups(reads,
                    num_reads, params.num_procs)
        else:
            work = FastaReader(params.fasta)
    elif params.fastq:
        reads = FastqReader(params.fastq)
        # get read count of input
        num_reads = get_sequence_count(params.fasta, 'fastq')
        if params.num_procs > 1:
            # split reads into generator objects based on equal
            # split across cores
            work = split_fasta_reads_into_groups(reads,
                    num_reads, params.num_procs)
        else:
            work = FastqReader(params.fasta, params.quality)
    return num_reads, work

def main():
    """Main loop"""
    start_time = time.time()
    motd()
    args = get_args()
    print 'Started: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(start_time))
    # build our configuration object w/ input params
    conf = ConfigParser.ConfigParser()
    conf.read(args.config)
    params = Parameters(conf)
    #pdb.set_trace()
    # create the db and tables, returninn connection
    # and cursor
    conn, cur = db.create_db_and_new_tables(params.db)
    # get num reads and split up work
    num_reads, work = get_work(params)
    # give some indication of progress for longer runs
    if num_reads > 999:
        sys.stdout.write('Running')

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
        for unit in xrange(num_reads):
            #enter_to_db(results.get())
            db.insert_record_to_db(cur, results.get())
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
        # fake a multiprocessing queue, so stacking and accessing results
        # is identical.
        results = ListQueue()
        singleproc(work, results, params)
        for tagged in results:
            db.insert_record_to_db(cur, tagged)
    conn.commit()
    cur.close()
    conn.close()
    end_time = time.time()
    pretty_end_time = time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(end_time))
    print "\nEnded: {} (run time {} minutes)".format(pretty_end_time,
            round((end_time - start_time)/60, 3))

if __name__ == '__main__':
    main()
