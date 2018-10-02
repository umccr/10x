#!/usr/bin/env python

import logging
import unittest

from tenx_telomeres.hg38_synthetic_telomeres import find_N_boundaries, elongate_forward_sequence, elongate_reverse_sequence

# Set log level
loglevel = logging.INFO
logging.basicConfig(level=loglevel)
log = logging.getLogger(__name__)

class TestStringMethods(unittest.TestCase):

    def setUp(self):
                      # 0               16                33          47
        self.src_seq = 'NNNNNNNNNNNNNNNNTAACCCTAACCCTAACCCNNNNNNNNNNNNN'
        self.fwd_seq = 'ACCCTAACCCTAACCCTAACCCTAACCCTAACCCNNNNNNNNNNNNN'
        self.rev_seq = 'NNNNNNNNNNNNNNNNTAACCCTAACCCTAACCCTAACCCTAACCCT'
        self.no_patt = 'NNNNNNNNNNNNNNNNATCATAAtAaaaaccCCANNNNNNNNNNNNN'
        self.kmer_k = 6

    def test_N_boundary_positions(self):
        boundaries = find_N_boundaries(self.src_seq)
        self.assertTupleEqual(boundaries, (16, 33))

    def test_elongate_forward_sequence(self):
        tst_seq = elongate_forward_sequence(self.src_seq)

        self.assertEqual(len(tst_seq), len(self.fwd_seq))
        self.assertEqual(tst_seq, self.fwd_seq)

    def test_elongate_reverse_sequence(self):
        tst_seq = elongate_reverse_sequence(self.src_seq)

        #self.assertEqual(len(tst_seq), len(self.rev_seq))
        self.assertEqual(tst_seq, self.rev_seq)

if __name__ == '__main__':
    unittest.main()