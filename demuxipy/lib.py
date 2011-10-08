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
import ConfigParser
from collections import defaultdict
from multiprocessing import cpu_count
from seqtools.sequence.fasta import FastaSequence
from seqtools.sequence.transform import DNA_reverse_complement

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
        try:
            self.fasta        = os.path.abspath(os.path.expanduser(
                    self.conf.get('Sequence','fasta').strip("'")))
            self.quality      = os.path.abspath(os.path.expanduser( \
                    self.conf.get('Sequence','quality').strip("'")))

        except ConfigParser.NoOptionError:
            self.fastq       = self.conf.get('Sequence','fastq').strip("'")
        except ConfigParser.NoOptionError:
            print "Cannot find valid sequence files in [Sequence] section of {}".format(self.conf)
        self.db               = self.conf.get('Database','DATABASE')
        self.qual_trim        = self.conf.getboolean('Quality', 'QualTrim')
        self.min_qual         = self.conf.getint('Quality', 'MinQualScore')
        self.outer         = self.conf.getboolean('Mid (Outer) Tags','MidTrim')
        self.outer_fuzzy        = self.conf.getboolean('Mid (Outer) Tags','MidFuzzyMatching')
        self.outer_errors = self.conf.getint('Mid (Outer) Tags','MidAllowedErrors')
        self.inner      = self.conf.getboolean('Linker (Inner) Tags', 'LinkerTrim')
        self.inner_fuzzy     = self.conf.getboolean('Linker (Inner) Tags','LinkerFuzzyMatching')
        self.inner_errors = self.conf.getint('Linker (Inner) Tags','LinkerAllowedErrors')
        self.concat_check    = self.conf.getboolean('Concatemers','ConcatemerChecking')
        self.concat_fuzzy    = self.conf.getboolean('Concatemers','ConcatemerFuzzyMatching')
        self.concat_allowed_errors   = self.conf.getboolean('Concatemers','ConcatemerAllowedErrors')
        self.search          = self.conf.get('Search','SearchFor')
        all_outer             = self._get_all_outer()
        all_inner          = self._get_all_inner()
        self.sequence_tags   = SequenceTags(all_outer, all_inner, self.search,
                self.conf.items(self.search), self.conf.getint('Mid (Outer) Tags','MidGap'), 
                self.conf.getint('Linker (Inner) Tags','LinkerGap'),
                self.concat_check)
        self.multiprocessing = conf.get('Multiprocessing', 'Multiprocessing') 
        # compute # cores for computation; leave 1 for db and 1 for sys
        if self.multiprocessing == True:
            if conf.get('Multiprocessing','processors').lower() == 'auto' and cpu_count > 2:
                self.num_procs = cpu_count() - 1
            elif conf.get('Multiprocessing','processors').lower() != 'auto' and \
                    cpu_count >= conf.getint('Multiprocessing','processors'):
                self.num_procs = conf.getint('Multiprocessing','processors')
        else:
            self.num_procs = 1

    def __repr__(self):
        return '''<linkers.parameters run values>'''

    def _get_all_outer(self):
        # if only linkers, you don't need MIDs
        if self.search == 'MidGroups' or 'MidLinkerGroups':
            return dict(self.conf.items('Mids'))
        else:
            return None

    def _get_all_inner(self):
        # if only linkers, you don't need MIDs
        if self.search == 'LinkerGroups' or 'MidLinkerGroups':
            return dict(self.conf.items('Linkers'))
        else:
            return None

class SequenceTags():
    """ """
    def __init__(self, outers, inners, search, group, outer_gap, inner_gap,
            concat):
        self.outers = None
        self.inners = None
        self.cluster_map = None
        self.all_tags = None
        self.outer_gap = outer_gap
        self.inner_gap = inner_gap
        if outers:
            # map mid sequences to names
            self.reverse_outer_lookup = self._reverse_dict(outers)
            self.outer_len = self._get_length(outers)
        if inners:
            # map linker sequences to names
            self.reverse_inner_lookup = self._reverse_dict(inners)
            self.inner_len = self._get_length(inners)
        # pare down the list of linkers and MIDS to those we've used
        self._generate_clusters_and_get_cluster_tags(outers, inners, search,
                group)
        # do we check for concatemers?
        if concat:
            self._all_possible_tags(search)

    def _get_length(self, tags):
        #linkers         = dict(self.conf.items('Linker'))
        lset = set([len(l) for l in tags.values()])
        assert len(lset) == 1, "Your {} sequences are difference lengths".format(name)
        return lset.pop()

    def _build_regex(self, tags, gap, rev = False):
        if not rev:
            return [re.compile('^[acgtnACGTN]{{0,{}}}{}'.format(gap, seq)) for seq in
                    tags]
        else:
            return [re.compile('{}[acgtnACGTN]{{0,{}}}$'.format(seq, gap)) 
                    for seq in tags]

    def _reverse_dict(self, d):
        return {v:k for k,v in d.iteritems()}

    def _generate_clusters_and_get_cluster_tags(self, all_outers, all_inners, search,
            group):

        self.cluster_map = defaultdict(lambda : defaultdict(str))

        if search == 'MidLinkerGroups':
            self.outers = defaultdict(list)
            self.inners = defaultdict(lambda : defaultdict(list))
            for row in group:
                m,l = row[0].replace(' ','').split(',')
                org = row[1]
                self.outers['forward_string'].append(all_outers[m])
                #self.linkers['string'].appendd(all_linkers[l])
                self.inners[all_outers[m]]['forward_string'].append(all_inners[l])
                self.inners[all_outers[m]]['reverse_string'].append(DNA_reverse_complement(all_inners[l]))
                j = "{},{}".format(all_outers[m],all_inners[l])
                self.cluster_map[all_outers[m]][all_inners[l]] = org
            
            self.outers['forward_regex'] = \
                    self._build_regex(self.outers['forward_string'], 
                    self.outer_gap)
            for m in self.inners:
                self.inners[m]['forward_regex'] = \
                    self._build_regex(self.inners[m]['forward_string'], 
                    self.inner_gap)
                self.inners[m]['reverse_regex'] = \
                    self._build_regex(self.inners[m]['reverse_string'], 
                    self.inner_gap, rev = True)

        elif search == 'MidGroups':
            for row in group:
                self.outers = defaultdict(list)
                m,l = row[0].replace(' ','').split(',')
                org = row[1]
                self.outers['forward_string'].append(all_outers[m])
                #self.linkers['string'].append(all_linkers[l])
                self.cluster_map[all_outers[m]]['None'] = org
            self.outers['forward_regex'] = \
                    self._build_regex(self.outers['forward_string'], 
                    self.outer_gap)

        elif search == 'LinkerGroups':
            self.inners = defaultdict(lambda : defaultdict(list))
            for row in group:
                m,l = row[0].replace(' ','').split(',')
                org = row[1]
                self.linkers[str(self.outers)]['forward_string'].append(all_inners[l])
                self.linkers[str(self.outers)]['reverse_string'].append(DNA_reverse_complement(all_inners[l]))
                self.cluster_map['None'][all_inners[m]] = org
            for m in self.inners:
                self.inners[m]['forward_regex'] = \
                    self._build_regex(self.inners[m]['forward_string'], 
                    self.inner_gap)
                self.inners[m]['reverse_regex'] = \
                    self._build_regex(self.inners[m]['reverse_string'],
                    self.inner_gap, rev = True)
        #pdb.set_trace()

    def _all_possible_tags(self, search):
        '''Create regular expressions for the forward and reverse complements
        of all of the tags sequences used in a run'''
        # at = all tags; rat = reverse complement all tags
        self.all_tags = defaultdict(lambda : defaultdict())
        for m in self.inners:
            self.all_tags[m]['string'] = self.inners[m]['forward_string'] + self.inners[m]['reverse_string']
            self.all_tags[m]['regex'] = [re.compile(t) for t in
                    self.all_tags[m]['string']]

class Tagged():
    '''Trimming, tag, and sequence data for individual reads'''
    def __init__(self, sequence):
        # super(Params, self).__init__()
        assert isinstance(sequence,FastaSequence), \
            'The Record class must be instantiated with a FastaSequence object'
        self.read               = sequence # a biopython sequence object
        self.outer                = None
        self.outer_name           = None
        self.outer_seq            = None
        self.outer_match          = None
        self.outer_type         = None
        self.inner_name         = None
        self.inner_seq          = None
        self.inner_match        = None
        self.inner_type         = None
        self.cluster            = None
        self.concat_seq         = None
        self.concat_type        = None
        self.concat_match       = None
    
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
