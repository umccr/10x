#!/usr/bin/env python

import math
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
        self.no_patt = 'NNNNNNNNNNNNNNNNATCATAAtAaaaaccCCANNNNNNNNNNNNN'
        self.tst_seq = None
        self.kmer_k = 6

    def test_N_boundary_positions(self):
        boundaries = find_N_boundaries(self.src_seq)
        self.assertTupleEqual(boundaries, (16, 33)) #XXX: review off-by ones

    def test_elongate_forward_sequence(self):
        pos = 0
        boundaries = find_N_boundaries(self.src_seq)
        boundary = boundaries[0]
        boundary_r = boundaries[1]

        kmer_seq = self.src_seq[boundary:boundary + self.kmer_k]

        # How many chunks to elongate and remainder
        chunks = len(self.src_seq[0:boundary]) % self.kmer_k
        chunks_r = len(self.src_seq[0:boundary]) / self.kmer_k

        # Capture remainder of the pattern to fit in sequence
        kmer_seq_r = kmer_seq[math.floor(chunks_r):]
        tmp_seq = kmer_seq_r
        
        # Build forward sequence
        for cnk in range(0, chunks - 2):
            tmp_seq = tmp_seq + kmer_seq

        # Attach inner pattern
        tst_seq = tmp_seq + self.src_seq[boundary:boundary_r] + self.src_seq[boundary_r:]
        tst_seq

        self.assertEqual(len(tst_seq), len(self.dst_seq))
        self.assertEqual(tst_seq, self.dst_seq)

    def test_elongate_reverse_sequence(self):
        pass

if __name__ == '__main__':
    unittest.main()