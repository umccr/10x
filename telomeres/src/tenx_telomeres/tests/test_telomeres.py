#!/usr/bin/env python

import logging
import unittest

from tenx_telomeres.hg38_synthetic_telomeres import find_N_boundaries

# Set log level
loglevel = logging.INFO
logging.basicConfig(level=loglevel)
log = logging.getLogger(__name__)

class TestStringMethods(unittest.TestCase):

    def setUp(self):
                      # 0               16                33          47 
        self.src_seq = 'NNNNNNNNNNNNNNNNTAACCCTAACCCTAACCCNNNNNNNNNNNNN'
        self.dst_seq = 'ACCCTAACCCTAACCCTAACCCTAACCCTAACCCNNNNNNNNNNNNN'
        self.tst_seq = None

    def test_N_boundary_positions(self):
        boundaries = find_N_boundaries(self.src_seq)
        self.assertTupleEqual(boundaries, (16, 33))

    def test_elongate_forward_sequence(self):        
        #self.assertEqual(self.tst_seq, self.dst_seq)
        pass

    def test_elongate_reverse_sequence(self):
        pass

    def test_elongate_sequence(self):
        pass

if __name__ == '__main__':
    unittest.main()