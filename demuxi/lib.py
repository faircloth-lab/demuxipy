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

class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))

class Parameters():
    '''linkers.py run parameters'''
    def __init__(self, conf):
        self.conf            = conf
        self.db              = self.conf.get('Database','DATABASE')
        self.qualTrim        = self.conf.getboolean('Steps', 'QualTrim')
        self.minQual         = self.conf.getint('GeneralParameters', 'MinQualScore')
        self.midTrim         = self.conf.getboolean('Steps','MidTrim')
        self.midGap          = self.conf.getint('GeneralParameters','MidGap')
        self.linkerTrim      = self.conf.getboolean('Steps', 'LinkerTrim')
        self.linkerGap       = self.conf.getint('GeneralParameters','LinkerGap')
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
        self._setup()
    
    def __repr__(self):
        return '''<linkers.parameters run values>'''
    
    def _setup(self):
        if self.midTrim and self.linkerTrim:
            self._mid()
            self._linkers()
            self.clust       = self.conf.items('MidLinkerGroups')
            self._tagLibrary()

        elif self.midTrim and not self.LinkerTrim:
            self_mid()
            self.clust       = self.conf.items('MidGroup')
            self._tagLibrary()
            
        elif not self.midTrim and self.linkerTrim:
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

class Record():
    '''Trimming, tag, and sequence data for individual reads'''
    def __init__(self, sequence):
        # super(Params, self).__init__()
        assert isinstance(sequence,SeqRecord), \
            'The Record class must be instantiated with a BioPython Seq object'
        self.unmod              = sequence # a biopython sequence object
        self.sequence           = None
        self.nCount             = None
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
    
    def __repr__(self):
        return '''<linkers.record for %s>''' % self.unmod.id

def reverse(items, null=False):
    '''build a reverse dictionary from a list of tuples'''
    l = []
    if null:
        items += ((None, None),)
    for i in items:
        t = (i[1],i[0])
        l.append(t)
    return dict(l)
