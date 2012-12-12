#!/usr/bin/env python
# encoding: utf-8

"""
File: demuxi.py
Author: Brant Faircloth

Created by Brant Faircloth on 05 October 2011 14:10 PDT (-0700)
Copyright (c) 2011 Brant C. Faircloth. All rights reserved.

Description: 

"""

#import os
import sys
#import re
import gzip
import time
import numpy
#import string
#import cPickle
#import sqlite3
import argparse
import itertools
import ConfigParser

from multiprocessing import Process, Queue, JoinableQueue

#from seqtools.sequence.fastq import FastqReader
from seqtools.sequence.fasta import FastaQualityReader
from seqtools.sequence.transform import DNA_reverse_complement

from demuxipy import db
from demuxipy import pairwise2
from demuxipy.lib import FullPaths, ListQueue, Tagged, Parameters

import pdb


def motd():
    """Print a welcome message """
    motd = """
    ###############################################################
    #                       demuxi.py                             #
    #                                                             #
    # Demultiplexing of hierarchically tagged, massively parallel #
    # DNA sequence reads from 454 data.                           #
    #                                                             #
    #                                                             #
    # Copyright (c) 2009-2012 Brant C. Faircloth                  #
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
        #pdb.set_trace()
        try:
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
        except IndexError:
            pass
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

def find_left_tag(s, tag_regexes, tag_strings, max_gap_char, tag_len, fuzzy, errors):
    """Matching methods for left linker - regex first, followed by fuzzy (SW)
    alignment, if the option is passed"""
    #pdb.set_trace()
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

def find_right_tag(s, tag_regexes, tag_strings, max_gap_char, tag_len,
        fuzzy, errors, tagged, revcomp = True):
    """Matching methods for right linker - regex first, followed by fuzzy (SW)
    alignment, if the option is passed"""
    #if 'MID15_NoError_SimpleX1_NoError_F_NEQ_R' in tagged.read.identifier:
    #    pdb.set_trace()

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
            # correct match_position
            start_of_slice = len(s) - (tag_len + max_gap_char)
            m_type = 'fuzzy'
            tag_matched = match[0]
            seq_matched = match[3]
            start, stop = get_align_match_position(match[3], match[4], match[5])
            start, stop = start + start_of_slice, stop + start_of_slice
    if match and revcomp:
        return DNA_reverse_complement(tag_matched), m_type, start, stop, DNA_reverse_complement(seq_matched)
    elif match and not revcomp:
        return tag_matched, m_type, start, stop, seq_matched
    else:
        return None

def trim_one(tagged, regexes, strings, buff, length, fuzzy, errors, trim = 0):
    """Remove the MID tag from the sequence read"""
    #if sequence.id == 'MID_No_Error_ATACGACGTA':
    #    pdb.set_trace()
    #pdb.set_trace()
    mid = find_left_tag(tagged.read.sequence,
                regexes,
                strings,
                buff,
                length,
                fuzzy,
                errors
            )
    if mid:
        target, match_type, match = mid[0], mid[1], mid[4]
        #tagged.mid, tagged.m_type, tagged.seq_match = mid[0],mid[1],mid[4]
        tagged.read = tagged.read.slice(mid[3] + trim, len(tagged.read.sequence), False)
        #tagged.mid_name = params.sequence_tags.reverse_mid_lookup[tagged.mid]
        return tagged, target, match_type, match
    else:
        return tagged, None, None, None

def trim_two(tagged, fregex, fstring, rregex, rstring, buff,
        length, fuzzy, errors, trim = 0, revcomp = True):
    """Use regular expression and (optionally) fuzzy string matching
    to locate and trim linkers from sequences"""

    left = find_left_tag(tagged.read.sequence,
                fregex,
                fstring,
                buff,
                length,
                fuzzy,
                errors
            )
    
    right = find_right_tag(tagged.read.sequence,
                rregex,
                rstring,
                buff,
                length,
                fuzzy,
                errors,
                tagged,
                revcomp
            )

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
        tagged.read = tagged.read.slice(left[3], right[2], False)
        # left and right are identical so largely pass back the left
        # info... except for m_type which can be a combination
        target, match = left[0], left[4]
        match_type = "{}-{}-both".format(left[1], right[1])

    elif left is not None \
            and right is not None\
            and left[0] != right[0]:
        # these are no good.  check for within gaps
        # to make sure it's not a spurious match
        #pdb.set_trace()
        tagged.trimmed = None
        target, match = None, None
        match_type = "tag-mismatch"

    elif right is None and left and left[2] <= buff:
        tagged.read = tagged.read.slice(left[3], len(tagged.read.sequence), False)
        target, match = left[0], left[4]
        match_type = "{}-left".format(left[1])

    elif left is None and right and right[2] >= (len(tagged.read.sequence) - (len(right[0]) + buff)):
        tagged.read = tagged.read.slice(0, right[2], False)
        target, match = right[0], right[4]
        match_type = "{}-right".format(right[1])

    else:
        target, match_type, match = None, None, None

    return tagged, target, match_type, match


def concat_check(tagged, params):
    """Check screened sequence for the presence of concatemers by scanning 
    for all possible tags - after the 5' and 3' tags have been removed"""
    s = tagged.read.sequence
    m_type = None
    #pdb.set_trace()
    for tag in params.sequence_tags.all_tags[str(tagged.outer_seq)]['regex']:
        match = tag.search(s)
        if match:
            tagged.concat_seq = tag.pattern
            tagged.concat_type = "regex-concat"
            tagged.concat_match = tagged.read.sequence[match.start():match.end()]
            break
    if match is None and params.concat_fuzzy:
        match = align(s,
                params.sequence_tags.all_tags[str(tagged.outer_seq)]['string'], 
                params.concat_allowed_errors
            )
        if match:
            tagged.concat_seq = match[0]
            tagged.concat_type = "fuzzy-concat"
            tagged.concat_match = match[3]
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
        #pdb.set_trace()
        # for now, we'll keep this here
        tagged = Tagged(sequence)
        # trim
        if params.qual_trim:
            #pdb.set_trace()
            tagged.read = tagged.read.trim(params.min_qual, False)
        # check for Outers:
        if (params.search == 'OuterGroups' or params.search == 'OuterInnerGroups'):
            assert params.outer, "Search != True for Outer tags"
            if params.outer_type.lower() == 'single':
                result = trim_one(
                    tagged,
                    params.sequence_tags.outers['forward_regex'],
                    params.sequence_tags.outers['forward_string'],
                    params.sequence_tags.outer_gap,
                    params.sequence_tags.outer_len,
                    params.outer_fuzzy,
                    params.outer_errors
                )
            elif params.outer_type.lower() == 'both':
                if params.outer_orientation.lower() == 'reverse':
                    revcomp = True
                else:
                    revcomp = False
                result = trim_two(
                    tagged,
                    params.sequence_tags.outers['forward_regex'],
                    params.sequence_tags.outers['forward_string'],
                    params.sequence_tags.outers['reverse_regex'],
                    params.sequence_tags.outers['reverse_string'],
                    params.sequence_tags.outer_gap,
                    params.sequence_tags.outer_len,
                    params.outer_fuzzy,
                    params.outer_errors,
                    revcomp
                )
        tagged, tagged.outer_seq, tagged.outer_type, tagged.outer_match = result
        if tagged.outer_seq:
            tagged.outer_name = params.sequence_tags.reverse_outer_lookup[tagged.outer_seq]
        # check for Inners
        if (tagged.outer_seq and (params.search == 'OuterInnerGroups' or params.search == 'InnerGroups')):
            assert params.inner, "Search != for Inner tags."
            if params.inner_type.lower() == 'single':
                result = trim_one(
                    tagged,
                    params.sequence_tags.inners[tagged.outer_seq]['forward_regex'],
                    params.sequence_tags.inners[tagged.outer_seq]['forward_string'],
                    params.sequence_tags.inners[tagged.outer_seq]['reverse_regex'],
                    params.sequence_tags.inners[tagged.outer_seq]['reverse_string'],
                    params.sequence_tags.inner_gap,
                    params.sequence_tags.inner_len,
                    params.inner_fuzzy,
                    params.inner_errors
                )
            elif params.inner_type.lower() == 'both':
                if params.outer_orientation.lower() == 'reverse':
                    revcomp = True
                else:
                    revcomp = False
                result = trim_two(
                    tagged,
                    params.sequence_tags.inners[tagged.outer_seq]['forward_regex'],
                    params.sequence_tags.inners[tagged.outer_seq]['forward_string'],
                    params.sequence_tags.inners[tagged.outer_seq]['reverse_regex'],
                    params.sequence_tags.inners[tagged.outer_seq]['reverse_string'],
                    params.sequence_tags.inner_gap,
                    params.sequence_tags.inner_len,
                    params.inner_fuzzy,
                    params.inner_errors,
                    revcomp
                )
            tagged, tagged.inner_seq, tagged.inner_type, tagged.inner_match = result
            if tagged.inner_seq:
                tagged.inner_name = params.sequence_tags.reverse_inner_lookup[tagged.inner_seq]

        # lookup cluster name; should => None, None is no outers or inners
        tagged.cluster = params.sequence_tags.cluster_map[str(tagged.outer_seq)][str(tagged.inner_seq)]

        # check for concatemers
        if (params.concat_check and len(tagged.read.sequence) > 0) and \
                ((tagged.outer_seq and tagged.inner_seq and params.search ==
                    'OuterInnerGroups') or \
                (tagged.outer_seq and params.search == 'OuterGroups') or \
                (tagged.inner_seq and params.search == 'InnerGroups')):
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
    elif kind == 'fastq' and input.endswith('gz'):
        return sum([1 for line in gzip.open(input, 'rb')]) / 4
    else:
        return sum([1 for line in open(input, 'rU')]) / 4


def split_fasta_reads_into_groups(reads, num_reads, num_procs):
    job_size = num_reads / num_procs
    sys.stdout.write("Parsing reads into groups of {} reads\n".format(job_size))
    sys.stdout.flush()
    i = iter(reads)
    chunk = list(itertools.islice(i, job_size))
    while chunk:
        yield chunk
        chunk = list(itertools.islice(i, job_size))


def imerge(a, b):
    for i, j in itertools.izip(a, b):
        yield i, j


def get_work(params):
    if params.fasta and params.quality:
        reads = FastaQualityReader(params.fasta, params.quality)
        # get read count of input
        num_reads = get_sequence_count(params.fasta, 'fasta')
        if params.num_procs > 1:
            # split reads into generator objects based on equal
            # split across cores
            work = split_fasta_reads_into_groups(
                    reads,
                    num_reads,
                    params.num_procs
                )
        else:
            work = FastaQualityReader(params.fasta, params.quality)
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
    # create the db and tables, returning connection
    # and cursor
    conn, cur = db.create_db_and_new_tables(params.db)
    # get num reads and split up work
    num_reads, work = get_work(params)
    #pdb.set_trace()
    # MULTICORE
    if params.multiprocessing and params.num_procs > 1:
        jobs = Queue()
        results = JoinableQueue()
        # We're stacking groups of jobs on the work
        # Queue, conceivably to save the overhead of
        # placing them on there one-by-one.
        for unit in work:
            jobs.put(unit)
        # setup the processes for the jobs
        sys.stdout.write("Starting {} workers\n".format(params.num_procs))
        sys.stdout.flush()
        sys.stdout.write('Running')
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
            pdb.set_trace()
    conn.commit()
    cur.close()
    conn.close()
    end_time = time.time()
    pretty_end_time = time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(end_time))
    print "\nEnded: {} (run time {} minutes)".format(pretty_end_time,
            round((end_time - start_time)/60, 3))

if __name__ == '__main__':
    main()
