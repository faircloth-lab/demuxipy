"""
File: lib.py
Author: Brant Faircloth

Created by Brant Faircloth on 01 October 2011 11:10 PDT (-0700)
Copyright (c) 2011 Brant C. Faircloth. All rights reserved.

Description:  common files for demuxi.py

"""
import os
import sys
import argparse
from multiprocessing import cpu_count
from tools.sequence.fasta import FastaSequence

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
        self.fasta        = self.conf.get('Input','fasta')
        self.quality      = self.conf.get('Input','quality')
        #except:
        #    self.fastq       = self.conf.get('Input','fastq')
        self.db              = self.conf.get('Database','DATABASE')
        self.qual_trim       = self.conf.getboolean('Steps', 'QualTrim')
        self.min_qual         = self.conf.getint('GeneralParameters', 'MinQualScore')
        self.mid_trim         = self.conf.getboolean('Steps','MidTrim')
        self.mid_gap          = self.conf.getint('GeneralParameters','MidGap')
        self.linker_trim      = self.conf.getboolean('Steps', 'LinkerTrim')
        self.linker_gap       = self.conf.getint('GeneralParameters','LinkerGap')
        self.concat          = self.conf.getboolean('GeneralParameters','CheckForConcatemers')
        self.fuzzy           = self.conf.getboolean('GeneralParameters','FuzzyMatching')
        self.allowed_errors  = self.conf.getint('GeneralParameters','AllowedErrors')
        self.mids            = None
        self.reverse_mid     = None
        self.linkers         = None
        self.reverse_linkers = None
        self.clust           = None
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
    
    def _setup(self):
        if self.mid_trim and self.linker_trim:
            self._mid()
            self._linkers()
            self.clust       = self.conf.items('MidLinkerGroups')
            self._tagLibrary()

        elif self.mid_trim and not self.linker_trim:
            self_mid()
            self.clust       = self.conf.items('MidGroup')
            self._tagLibrary()
            
        elif not self.mid_trim and self.linker_trim:
            self._linkers()
            self.clust       = self.conf.items('LinkerGroups')
            self._tagLibrary()
            
        # do we check for concatemers?
        if self.concat:
            self._allPossibleTags()
    
    def _linkers(self):
        self.linkers         = dict(self.conf.items('Linker'))
        self.reverse_linkers = reverse(self.conf.items('Linker'), True)
        lset = set([len(l) for l in self.linkers.values()])
        assert len(lset) == 1, "Your linker sequences are difference lengths"
        self.linker_len = lset.pop()
    
    def _mid(self):
        self.mids            = dict(self.conf.items('Mid'))
        self.reverse_mid     = reverse(self.conf.items('Mid'), True)
        mset = set([len(m) for m in self.mids.values()])
        assert len(mset) == 1, "Your MID sequences are different lengths"
        self.mid_len = mset.pop()

    
    def _tagLibrary(self):
        '''Create a tag-library from the mids and the linkers which allows us to 
        track which organisms go with which MID+linker combo'''
        self.tags = {}
        for c in self.clust:
            #pdb.set_trace()
            if self.mids and self.linkers:
                m,l = c[0].replace(' ','').split(',')
                org = c[1]
                if self.mids[m] not in self.tags.keys():
                    self.tags[self.mids[m]] = {self.linkers[l]:org}
                else:
                    self.tags[self.mids[m]][self.linkers[l]] = org
            
            elif not self.mids and self.linkers:
                l = c[0]
                org = c[1]
                self.tags[self.linkers[l]] = org
            
            elif self.mids and not self.linker:
                l = c[0]
                org = c[1]
                self.tags[self.mids[l]] = org
    
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
            self.all_tags.append(revComp(self.linkers[l]))
            self.all_tags_regex.append(re.compile('%s' % revComp(self.linkers[l])))

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
        self.concat_type        = None
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
