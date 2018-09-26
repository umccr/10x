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

        kmer_seq = self.src_seq[boundary:boundary + self.kmer_k]
        chunks = len(self.src_seq[0:boundary]) % self.kmer_k

        for pos, chunk in zip(range(0, boundary, self.kmer_k), range(1, chunks)):
            posi = boundary - self.kmer_k * chunk
            log.info(self.src_seq[boundary:posi])
            #log.info(self.tst_seq)
        
        #self.assertEqual(self.tst_seq, self.dst_seq)

    def test_elongate_reverse_sequence(self):
        pass

if __name__ == '__main__':
    unittest.main()