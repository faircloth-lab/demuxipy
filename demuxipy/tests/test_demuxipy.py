"""
File: test_demuxipy.py
Author: Brant Faircloth

Created by Brant Faircloth on 09 October 2011 13:10 PDT (-0700)
Copyright (c) 2011 Brant C. Faircloth. All rights reserved.

Description: Tests for demuxipy/lib.py

"""

import os
import unittest
import ConfigParser
from demuxipy import lib

import pdb

class TestParamsInstance(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

class TestSequenceTagsInstance(unittest.TestCase):
    def setUp(self):
        conf = ConfigParser.ConfigParser()
        conf.read('./test-data/demuxi-test.conf')
        self.p = lib.Parameters(conf)

    def _combinatorial_tester(self, o, t, e):
        self.p.inner_orientation = o
        self.p.inner_type = t
        test = self.p._get_sequence_tags(self.p._get_all_outer(), 
                self.p._get_all_inner())
        assert test.cluster_map.keys() == ['None']
        assert len(test.cluster_map['None']) == len(e)
        for k,v in e.iteritems():
            assert k in test.cluster_map['None'].keys()
            assert test.cluster_map['None'][k] == v

    def test_combinatorial_cluster_single_reverse(self):
        """
        duck = CGTCGTGCGGAATC - read - GATTCGCCAGCAGC

        """
        self.p.search = 'InnerCombinatorial'
        expected = {
                'GCTGCTGGCGAATC,AGTTCCGCAGCACG': 'goose',
                'CGTCGTGCGGAATC,AGTTCCGCAGCACG': 'brant',
                'CGTCGTGCGGAATC,GATTCGCCAGCAGC': 'duck',
            }
        self._combinatorial_tester('Reverse','Single', expected)

    def test_combinatorial_cluster_single(self):
        """
        duck = CGTCGTGCGGAATC - read - GCTGCTGGCGAATC

        """
        self.p.search = 'InnerCombinatorial'
        expected = {
                'GCTGCTGGCGAATC,CGTGCTGCGGAACT': 'goose',
                'CGTCGTGCGGAATC,CGTGCTGCGGAACT': 'brant',
                'CGTCGTGCGGAATC,GCTGCTGGCGAATC': 'duck',
            }
        self._combinatorial_tester('Forward','Single', expected)

    def test_combinatorial_cluster_both_reverse(self):
        """
        duck = CGTCGTGCGGAATC - read - GATTCGCCAGCAGC
        duck = GCAGCACGCCTTAG - read - CTAAGCGGTCGTCG

        """
        self.p.search = 'InnerCombinatorial'
        expected = {
                'GCTGCTGGCGAATC,AGTTCCGCAGCACG': 'goose',
                'CGACGACCGCTTAG,TCAAGGCGTCGTGC': 'goose',
                'GCAGCACGCCTTAG,TCAAGGCGTCGTGC': 'brant',
                'CGTCGTGCGGAATC,AGTTCCGCAGCACG': 'brant',
                'CGTCGTGCGGAATC,GATTCGCCAGCAGC': 'duck',
                'GCAGCACGCCTTAG,CTAAGCGGTCGTCG': 'duck'
            }
        self._combinatorial_tester('Reverse','Both', expected)

    def test_combinatorial_cluster_both_forward(self):
        """
        duck = CGTCGTGCGGAATC - read - GCTGCTGGCGAATC
        duck = GCAGCACGCCTTAG - read - CGACGACCGCTTAG
        """
        self.p.search = 'InnerCombinatorial'
        expected = {
                'GCTGCTGGCGAATC,CGTGCTGCGGAACT': 'goose',
                'CGACGACCGCTTAG,GCACGACGCCTTGA': 'goose',
                'CGTCGTGCGGAATC,GCTGCTGGCGAATC': 'duck',
                'GCAGCACGCCTTAG,CGACGACCGCTTAG': 'duck',
                'CGTCGTGCGGAATC,CGTGCTGCGGAACT': 'brant',
                'GCAGCACGCCTTAG,GCACGACGCCTTGA': 'brant'
            }
        self._combinatorial_tester('Forward','Both', expected)














    def tearDown(self):
        del(self.p)




if __name__ == '__main__':
    unittest.main()


