

import os
import sys
import numpy
import oursql
import argparse
import ConfigParser
import cPickle as cpickle

from collections import defaultdict
from tools.fs.args import is_dir, is_sheet
from tools.sequence import fasta

import pdb


def get_args():
    parser = argparse.ArgumentParser(description='Match UCE probes to assembled contigs and store the data')
    parser.add_argument('configuration', help='The configuration file for the db')
    parser.add_argument('wildcard', help='An SQL wildcard expression to limit resuts by cluster')
    parser.add_argument('output',help='The output directory in which to store fasta file', type = is_dir)
    parser.add_argument('--sample_map', help = 'A xlsx or csv file mapping plate'+ 
            ' location to sample name', type = is_sheet)
    parser.add_argument('--quality', help = 'Output quality values as *.qual'+
            ' files', action='store_true')
    parser.add_argument('--mira', help = 'Format filenames for MIRA', action = 'store_true')
    return parser.parse_args()

def get_db_reads(cur, wildcard, qual):
    if not qual:
        query = '''SELECT name, mid, linker, cluster, seq_trimmed FROM sequence WHERE
            cluster LIKE "{}"'''.format(wildcard)
        cur.execute(query)
    else:
        query = '''SELECT name, mid, linker, cluster, record FROM sequence WHERE
            cluster LIKE "{}"'''.format(wildcard)
        cur.execute(query)

    return cur.fetchall()

def main():
    """"""
    args = get_args()
    conf = ConfigParser.ConfigParser()
    conf.read(args.configuration)
    conn = oursql.connect(user = conf.get('Database','USER'), 
          passwd = conf.get('Database','PASSWORD'), 
          db = conf.get('Database','DATABASE')
          )
    cur = conn.cursor()
    reads = defaultdict(list) 
    for read in get_db_reads(cur, args.wildcard, args.quality):
        reads[read[3]].append(read)
    #pdb.set_trace()
    for k,v in reads.iteritems():
        if not args.mira:
            fname = "{}.fasta".format(k)
            qname = "{}.qual".format(k)
        else:
            fname = "{}_in.454.fasta".format(k)
            qname = "{}_in.454.fasta.qual".format(k)
        if not args.quality:
            outf = fasta.FastaWriter(os.path.join(args.output, fname))
        else:
            outf = fasta.FastaWriter(os.path.join(args.output, fname),
                    os.path.join(args.output, qname))
        for read in v:
            seq = fasta.FastaSequence()
            seq.identifier = "{}|{}|{}|{}".format(read[0],read[1],read[2],k)
            if not args.quality:
                seq.sequence = read[4]
            else:
                #pdb.set_trace()
                record = cpickle.loads(read[4])
                seq.sequence = record.seq
                seq.quality = \
                    numpy.array(record.letter_annotations['phred_quality'])
            outf.write(seq)
        outf.close()
    pdb.set_trace()
    

if __name__ == '__main__':
    main()
