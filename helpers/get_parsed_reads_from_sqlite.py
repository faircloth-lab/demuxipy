"""
File: get_parsed_reads_from_sqlite.py
Author: Brant Faircloth

Created by Brant Faircloth on 04 December 2011 13:12 PST (-0800)
Copyright (c) 2011 Brant C. Faircloth. All rights reserved.

Description: 

"""

import os
import sys
import sqlite3
import argparse
import ConfigParser

from seqtools.sequence import fasta
from seqtools.fs.args import is_dir, FullPaths

try:
    import cPickle as pickle
except:
    print "Using pickle instead of cPickle"
    import pickle

import pdb

def get_args():
    """get CLI arguments"""
    parser = argparse.ArgumentParser(description="""Get reads from demuxipy
            database""")
    parser.add_argument('config', help = "The demuxipy configuration file", 
            action = FullPaths)
    parser.add_argument('outdir', help = """The ouput directory in which to store
        fasta files""", type = is_dir, action = FullPaths)
    return parser.parse_args()

def get_groups(cur):
    cur.execute("SELECT DISTINCT(cluster) FROM tags")
    return cur.fetchall()

def get_and_write_seqs(cur, group, fw):
    if group == "NoTag":
        cur.execute("""SELECT record FROM sequence, tags where sequence.id = tags.id
            and tags.cluster IS NULL""")
    else:
        cur.execute("""SELECT record FROM sequence, tags where sequence.id = tags.id
            and tags.cluster = ?""", (group,))
    #pdb.set_trace()
    for row in cur:
        seq = pickle.loads(str(row[0]))
        fw.write(seq)
        #pdb.set_trace()


def main():
    """ """
    args = get_args()
    conf = ConfigParser.ConfigParser()
    conf.read(args.config)
    db = conf.get('Database','DATABASE')
    conn = sqlite3.connect(db)
    cur = conn.cursor()
    groups = get_groups(cur)
    for group in groups:
        if group[0] is None:
            group = "NoTag"
        else:
            group = group[0]
        outfile = "{}.fasta".format(group)
        outpath = os.path.join(args.outdir, outfile)
        fw = fasta.FastaWriter(outpath) 
        get_and_write_seqs(cur, group, fw)
        fw.close()
        # get reads in group from db

if __name__ == '__main__':
    main()
