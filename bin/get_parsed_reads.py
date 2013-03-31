#!/usr/bin/env python
# encoding: utf-8
"""
get_parsed_reads.py

Created by Brant Faircloth on 2009-12-17.
Copyright (c) 2009 Brant Faircloth. All rights reserved.

USAGE:  python get_parsed_reads.py --configuration=db.conf --all --trimmed
"""

import os
import sys
import time
import numpy
import errno
import argparse
import ConfigParser
from collections import defaultdict
from seqtools.sequence.fasta import FastaQualityReader, FastaWriter

import pdb


class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))


def is_dir(dirname):
    if not os.path.isdir(dirname):
        msg = "{0} is not a directory".format(dirname)
        raise argparse.ArgumentTypeError(msg)
    else:
        return dirname


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Extract demuxipy reads from the database""")
    parser.add_argument(
            "fasta",
            action=FullPaths,
            help="""The input fasta file"""
        )
    parser.add_argument(
            "quality",
            action=FullPaths,
            help="""The input quality file"""
        )
    parser.add_argument(
            "output",
            action=FullPaths,
            help="""A directory to hold the output"""
        )
    parser.add_argument(
            "--cluster",
            type=str,
            default='all',
            help="""Grab sequence for a particular cluster"""
        )
    parser.add_argument(
            "--remap",
            action=FullPaths,
            type=str,
            default=None,
            help="""Given an input Conf file, remap cluster names to something else"""
        )
    parser.add_argument(
            "--remap-section",
            dest="remap_section",
            type=str,
            default=None,
            help="""Use the specified section within the Conf file"""
        )
    parser.add_argument(
            "--header",
            choices=['cluster', 'mid', 'name', 'id'],
            nargs='+',
            help="""The elements to add to the fasta header""",
        )
    parser.add_argument(
            "--filter-type",
            dest="filter_type",
            choices=['regex', 'fuzzy', 'both'],
            default='regex',
            help="""Filter the sequences based on match-types""",
        )
    parser.add_argument(
            "--filter-length",
            dest="filter_length",
            type=int,
            default=0,
            help="""Filter the sequence based on a minimum length"""
        )
    args = parser.parse_args()
    if args.remap:
        assert args.remap_section, parser.error("If you are remapping with a Conf file, you must pass a Section name.")
    return args


class Header:
    def __init__(self, h):
        h = h.replace('>', 'name=')
        dh = dict((record.split('=') for record in h.split(' ')))
        for k, v in dh.iteritems():
            setattr(self, k, v)


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def get_remap_dictionary(remap, section):
    config = ConfigParser.ConfigParser()
    config.read(remap)
    return dict(config.items(section))


def write_sequences(record, header, output, sample_map, count):
    if sample_map is not None:
        header.name = sample_map[header.cluster.lower()]
    else:
        header.name = header.cluster
    record.identifier += ' name={}'.format(header.name)
    # create the cluster-specific output directory if not exists
    outdir = os.path.join(output, header.name)
    mkdir_p(outdir)
    outf = FastaWriter(
            os.path.join(outdir, "{}.fasta".format(header.name)),
            os.path.join(outdir, "{}.qual".format(header.name)),
            mode='a'
        )
    if count != 0 and count % 1000 == 0:
        sys.stdout.write('.')
        sys.stdout.flush()
    outf.write(record)
    outf.close()
    count += 1
    return count, header


def filter_match_type(args, header):
    filtered = True
    if args.filter_type == 'regex':
        regex_status = ['regex', 'regex-left', 'regex-right', 'regex-regex-both']
        if header.outer in regex_status and header.inner in regex_status:
            filtered = False
    elif args.filter_type == 'fuzzy':
        fuzzy_status = ['fuzzy', 'fuzzy-left', 'fuzzy-right', 'fuzzy-fuzzy-both']
        if header.outer in fuzzy_status and header.inner in fuzzy_status:
            filtered = False
    elif args.filter_type == 'both':
        both_status = [
            'regex',
            'regex-left',
            'regex-right',
            'regex-regex-both',
            'fuzzy',
            'fuzzy-left',
            'fuzzy-right',
            'fuzzy-fuzzy-both',
            'regex-fuzzy-both',
            'fuzzy-regex-both',
            ]
        if header.outer in both_status and header.inner in both_status:
            filtered = False
    return filtered


def filter_length(args, header):
    filtered = True
    if int(header.length) >= args.filter_length:
        filtered = False
    else:
        filtered = True
    return filtered


def get_stats(values):
    stats = (
            len(values),
            numpy.sum(values),
            numpy.mean(values),
            numpy.median(values),
            1.96 * (numpy.std(values, ddof=1) / numpy.sqrt(len(values))),
            numpy.nanmin(values),
            numpy.nanmax(values)
        )
    return stats


def print_stats(stats):
    print "n:{0[0]:.>18d}\nsum:{0[1]:.>16,.2f}\nmean:{0[2]:.>15.2f}\nmedian:{0[3]:.>13.2f}\n95CI{0[4]:.>16.2f}\nmin:{0[5]:.>16.2f}\nmax:{0[6]:.>16.2f}".format(stats)


def get_read_length_summary_stats(summary):
    print "Read length summary stats"
    cl = [numpy.mean(v['length']) for k, v in summary.iteritems()]
    stats = get_stats(cl)
    print_stats(stats)


def get_read_count_summary_stats(summary):
    print "Read count summary stats"
    cc = [numpy.sum(v['count']) for k, v in summary.iteritems()]
    stats = get_stats(cc)
    print_stats(stats)


def get_read_summary_per_sample(summary):
    print "cluster\tn\tbp\tmean\tmedian\t95CI\tmin\tmax"
    for k, v in summary.iteritems():
        #pdb.set_trace()
        stats = get_stats(v['length'])
        print "{0}\t{1[0]:d}\t{1[1]:.0f}\t{1[2]:.2f}\t{1[3]:.2f}\t{1[4]:.2f}\t{1[5]:.0f}\t{1[6]:.0f}".format(
            k,
            stats,
            )


def main():
    start_time = time.time()
    args = get_args()
    if args.remap:
        sample_map = get_remap_dictionary(args.remap, args.remap_section)
    else:
        sample_map = None
    print 'Started: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(start_time))
    try:
        os.makedirs(args.output)
    except OSError as exc:
        if exc.errno == errno.EEXIST:
            raise OSError("The output directory exists.  Please remove it manually.")
    # get all sequences where cluster != null
    # put them in monolithic file
    # change the sequence header
    sys.stdout.write("Processing (1 dot / 1000 sequences)")
    sys.stdout.flush()
    count = 0
    summary = defaultdict(lambda: defaultdict(list))
    if args.cluster == 'all':
        for record in FastaQualityReader(args.fasta, args.quality):
            # parse fasta header to get info we need
            header = Header(record.identifier)
            match_filter = filter_match_type(args, header)
            length_filter = filter_length(args, header)
            if not match_filter and not length_filter:
                count, header = write_sequences(record, header, args.output, sample_map, count)
                summary[header.name]['count'].append(1)
                summary[header.name]['length'].append(int(header.length))
    print "\n\n"
    get_read_length_summary_stats(summary)
    print "\n"
    get_read_count_summary_stats(summary)
    print "\n"
    get_read_summary_per_sample(summary)
    end_time = time.time()
    print "\n\n"
    print 'Ended: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(end_time))
    print 'Time for execution: ', (end_time - start_time) / 60, 'minutes'


if __name__ == '__main__':
    main()

