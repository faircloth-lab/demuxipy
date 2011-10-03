"""
File: lib.py
Author: Brant Faircloth

Created by Brant Faircloth on 01 October 2011 11:10 PDT (-0700)
Copyright (c) 2011 Brant C. Faircloth. All rights reserved.

Description:  common files for demuxi.py

"""
import os
import re
import sys
import argparse
from collections import defaultdict
from multiprocessing import cpu_count
from tools.sequence.fasta import FastaSequence
from tools.sequence.transform import DNA_reverse_complement

import pdb

class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))

class ListQueue(list):
    def __init__(self):
        list.__init__(self)
    
    def put(self, item):
        """append an item to the list"""
        self.append(item)
    
    def get(self):
        """return an item from the list"""
        return self.pop()

    def __repr__(self):
        return "ListQueue"

class Parameters():
    '''linkers.py run parameters'''
    def __init__(self, conf):
        self.conf = conf
        #try:
        self.fasta        = self.conf.get('Sequence','fasta')
        self.quality      = self.conf.get('Sequence','quality')
        #except:
        #    self.fastq       = self.conf.get('Input','fastq')
        self.db               = self.conf.get('Database','DATABASE')
        self.qual_trim        = self.conf.getboolean('Quality', 'QualTrim')
        self.min_qual         = self.conf.getint('Quality', 'MinQualScore')
        self.mid_trim         = self.conf.getboolean('Mid (Outer) Tags','MidTrim')
        self.mid_fuzzy        = self.conf.getboolean('Mid (Outer) Tags','MidFuzzyMatching')
        self.mid_allowed_errors = self.conf.getint('Mid (Outer) Tags','MidAllowedErrors')
        self.mid_gap          = self.conf.getint('Mid (Outer) Tags','MidGap')
        self.linker_trim      = self.conf.getboolean('Linker (Inner) Tags', 'LinkerTrim')
        self.linker_gap       = self.conf.getint('Linker (Inner) Tags','LinkerGap')
        self.linker_fuzzy     = self.conf.getboolean('Linker (Inner) Tags','LinkerFuzzyMatching')
        self.linker_allowed_errors = self.conf.getint('Linker (Inner) Tags','LinkerAllowedErrors')
        self.concat_check    = self.conf.getboolean('Concatemers','ConcatemerChecking')
        self.concat_fuzzy    = self.conf.getboolean('Concatemers','ConcatemerFuzzyMatching')
        self.concat_errors   = self.conf.getboolean('Concatemers','ConcatemerAllowedErrors')
        self.search          = self.conf.get('Search','SearchFor')
        self.mids            = None
        self.reverse_mid     = None
        self.linkers         = None
        self.reverse_linkers = None
        self.tags            = None
        self.all_tags        = None
        self.all_tags_regex  = None
        # compute # cores for computation; leave 1 for db and 1 for sys
        if conf.get('Multiprocessing','processors').lower() == 'auto' and cpu_count > 2:
            self.num_procs = cpu_count() - 1
        elif conf.get('Multiprocessing','processors').lower() != 'auto' and \
            cpu_count >= conf.getint('Multiprocessing','processors'):
                self.num_procs = conf.getint('Multiprocessing','processors')
        else:
            self.num_procs = 1
        self._setup()
    
    def __repr__(self):
        return '''<linkers.parameters run values>'''
    
    def _build_both_tag_libraries(self, all_mids, all_linkers):
        '''Create a tag-library from the mids and the linkers which allows us to 
        track which organisms go with which MID+linker combo'''
        self.tags = defaultdict(lambda : defaultdict(str))
        self.mids = defaultdict(str)
        self.linkers = defaultdict(str)
        for c in self.clust:
            pdb.set_trace()
            if self.mid_trim and self.linker_trim:
                m,l = c[0].replace(' ','').split(',')
                org = c[1]
                self.mids[m] = all_mids[m]
                self.linkers[l] = all_linkers[l]
                self.tags[all_mids[m]][all_linkers[l]] = org
        self.mids_f_regex = _regex_builder(self.mids.values(), self.mid_gap)
        self.linkers_f_regex = _regex_builder(self.linkers.values(), 
                self.linker_gap)
        self.linkers_r_regex = _regex_builder(self.linkers.values(),
                self.linker_gap, reverse = True)

    def _build_linker_tag_library(self, all_linkers)
        for c in self.clust:
            l = c[0]
            org = c[1]
            self.linkers[l] = all_linkers[l]
            self.tags[self.linkers[l]] = org
        self.linkers_f_regex = _regex_builder(self.linkers.values(), 
                self.linker_gap)
        self.linkers_r_regex = _regex_builder(self.linkers.values(),
                self.linker_gap, reverse = True)

    def _build_mid_tag_library(self, all_mids)
        for c in self.clust:
            l = c[0]
            org = c[1]
            self.mids[m] all_mids[m]
            self.tags[self.mids[l]] = org
        self.mids_f_regex = _regex_builder(self.mids.values(), self.mid_gap)

    def _setup(self):
        if self.mid_trim and self.linker_trim and self.search == 'MidLinkerGroups':
            # get the list of taxa to tags
            self.clust = self.conf.items(self.search)
            # we need a list of all MIDs possible and their length
            all_mids, self.mid_len = self._mid()
            # we need a list of all linkers possible and their length
            all_linkers, self.linker_len = self._linkers()
            # reduce MIDS and linkers to only those we're using
            self._build_both_tag_libraries(all_mids, all_linkers)

        elif self.mid_trim and not self.linker_trim and self.search == "MidGroup":
            # get the list of taxa to tags
            self.clust = self.conf.items(self.search)
            # we need a list of all MIDs possible and their length
            all_mids, self.mid_len = self._mid()
            # reduce MIDS and linkers to only those we're using
            self._build_mid_tag_library(all_mids)

        elif not self.mid_trim and self.linker_trim and self.search == "LinkerGroup":
            # get the list of taxa to tags
            self.clust = self.conf.items(self.search)
            # we need a list of all linkers possible and their length
            all_linkers, self.linker_len = self._linkers()
            # reduce MIDS and linkers to only those we're using
            self._build_linker_tag_libraries(all_linkers)
            
        # do we check for concatemers?
        if self.concat_check:
            self._allPossibleTags()

    def _regex_builder(self, tags, gap, reverse = True):
        if not reverse:
            return [re.compile('^[acgtnACGTN]{{0,{}}}{}'.format(gap, seq)) for name, seq in
                tags.iteritems()]
        else:
            return [re.compile('{}[acgtnACGTN]{{0,{}}}$'.format(DNA_reverse_complement(seq), 
                gap)) for name, seq in tags.iteritems()]

    
    def _linkers(self):
        linkers         = dict(self.conf.items('Linker'))
        lset = set([len(l) for l in linkers.values()])
        assert len(lset) == 1, "Your linker sequences are difference lengths"
        linker_len = lset.pop()
        return linkers, linker_len
    
    def _mid(self):
        mids = dict(self.conf.items('Mid'))
        mset = set([len(m) for m in mids.values()])
        assert len(mset) == 1, "Your MID sequences are different lengths"
        mid_len = mset.pop()
        return mids, mid_len

    def _allPossibleTags(self):
        '''Create regular expressions for the forward and reverse complements
        of all of the tags sequences used in a run'''
        # at = all tags; rat = reverse complement all tags
        self.all_tags = []
        self.all_tags_regex = []
        for c in self.clust:
            if self.mids and self.linkers:
                m,l = c[0].replace(' ','').split(',')
            elif not self.mids and self.linkers:
                l = c[0]
            elif self.mids and not self.linkers:
                pass
            self.all_tags.append(self.linkers[l])
            self.all_tags_regex.append(re.compile('%s' % self.linkers[l]))
            self.all_tags.append(DNA_reverse_complement(self.linkers[l]))
            self.all_tags_regex.append(re.compile('%s' % DNA_reverse_complement(self.linkers[l])))

class Tagged():
    '''Trimming, tag, and sequence data for individual reads'''
    def __init__(self, sequence):
        # super(Params, self).__init__()
        assert isinstance(sequence,FastaSequence), \
            'The Record class must be instantiated with a FastaSequence object'
        self.read               = sequence # a biopython sequence object
        self.mid                = None
        self.mid_seq            = None
        self.reverse_mid        = None
        self.seq_match          = None
        self.m_type             = None
        self.l_seq              = None
        self.l_tag              = None
        self.l_seq_match        = None
        self.l_critter          = None
        self.l_m_type           = None
        self.reverse_linker     = None
        self.concat_tag         = None
        self.concat_m_type      = None
        self.concat_seq_type    = None
        self.concat_count       = None
        self.concat_seq_match   = None
    
    #def __repr__(self):
    #    return '''<linkers.record for %s>''' % self.identifier

def reverse(items, null=False):
    '''build a reverse dictionary from a list of tuples'''
    l = []
    if null:
        items += ((None, None),)
    for i in items:
        t = (i[1],i[0])
        l.append(t)
    return dict(l)
