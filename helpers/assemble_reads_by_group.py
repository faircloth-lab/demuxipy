
import os
import sys
import glob
import shutil
import argparse


from collections import defaultdict
from tools.fs.args import is_dir, FullPaths

import pdb

def get_args():
    parser = argparse.ArgumentParser(description='Match UCE probes to assembled contigs and store the data')
    parser.add_argument('input', help = 'The input directory containing the'+
        ' fasta files', type = is_dir, action = FullPaths)
    parser.add_argument('output',help='The output directory in which to store'+
        ' the assembly files', type = is_dir, action = FullPaths)
    return parser.parse_args()

def create_dir(output):
    if not os.path.exists(output):
        os.mkdir(output)
        # make a directory for only the fasta reads
        outfasta = os.path.join(output, 'fasta')
        os.mkdir(outfasta)
    else:
        print "Path exists.  Remove manually."
        sys.exit()
    return outfasta

def main():
    args = get_args()
    outfasta = create_dir(args.output)
    # change to input directory - MIRA must be run from fasta
    # file location
    if os.getcwd() != args.input:
        print "Changing working directory to ", args.input
        os.chdir(args.input)
    for fasta in glob.glob(os.path.join(args.input, '*.fasta')):
        reads = sum([1 for line in open(fasta, 'rU').readlines() if \
            line.startswith('>')])
        #pdb.set_trace()
        basename = os.path.basename(fasta)
        dirname = basename.split('.')[0].rstrip('_in')
        cwd = os.path.dirname(fasta)
        if reads > 1:
            print "Assembling and symlinking", dirname
            # run MIRA and junk stuff that goes to stdout
            sysstring = "mira --project={} --job=denovo,est,accurate,454" + \
                " 454_SETTINGS -notraceinfo > /dev/null 2>&1"
            sysstring = sysstring.format(dirname)
            os.system(sysstring)
            outdirname = "{}_assembly".format(dirname)
            outdirpth = os.path.join(args.output, outdirname)
            #pdb.set_trace()
            shutil.move(os.path.join(cwd, outdirname), outdirpth)
            # link the unpadded fasta files to the outfasta
            fasta = os.path.join(outdirpth,
                    "{0}_d_results/{0}_out.unpadded.fasta".format(dirname))
            outfastaname = os.path.join(outfasta, dirname + '.fasta')
            os.symlink(fasta, outfastaname)
        else:
            print "There is only 1 read for {}".format(basename)
            print "Symlinking", dirname
            outfastaname = os.path.join(outfasta, dirname + '.fasta')
            # just symlink the unpadded fasta to the same spot
            os.symlink(fasta, outfastaname)


if __name__ == '__main__':
    main()
