#!/usr/bin/env python

import logging
import pymer
import unittest

from tenx_telomeres.hg38_synthetic_telomeres import TELO_HEXAMERS, HUMAN_TELOMERE, find_N_boundaries, elongate_forward_sequence, elongate_reverse_sequence, determine_hexamers, build_hexamer_table

# Set log level
loglevel = logging.INFO
logging.basicConfig(level=loglevel)
log = logging.getLogger(__name__)

class TestStringMethods(unittest.TestCase):

    def setUp(self):
                      # 0               16               33          46
        self.src_seq = 'NNNNNNNNNNNNNNNNTAACCCTAACCCTAACCCNNNNNNNNNNNNN'
        self.fwd_seq = 'ACCCTAACCCTAACCCTAACCCTAACCCTAACCCNNNNNNNNNNNNN'
        self.rev_seq = 'NNNNNNNNNNNNNNNNTAACCCTAACCCTAACCCTAACCCTAACCCT'
        self.fix_seq = 'NNNNNNNNNNTAACCCTAACCCTAACCCTAACCCNNNNNNNNNNNNN'
        self.rix_seq = 'NNNNNNNNNNNNNNNNTAACCCTAACCCTAACCCTAACCCNNNNNNN'
#                       'NNNNNNNNNNNNNNNNTAACCCTAACCCTAACCCCCCTAACCCTAANNNNNNNN'
        self.no_patt = 'NNNNNNNNNNNNNNNNATCATAAtAaaaaccCCANNNNNNNNNNNNN'
        self.mdl_nnn = 'NNNNNNNNNNNNNNNNCACACANNNNNNCACACANNNNNNNNNNNNN'
        self.kmer_k = 6

    def test_N_boundary_positions(self):
        boundaries = find_N_boundaries(self.src_seq)
        self.assertTupleEqual(boundaries, (16, 33))
        self.assertEqual(self.src_seq[16], 'T')
        self.assertEqual(self.src_seq[33], 'C')

        boundaries = find_N_boundaries(self.fwd_seq)
        self.assertTupleEqual(boundaries, (0, 33))

        boundaries = find_N_boundaries(self.rev_seq)
        self.assertTupleEqual(boundaries, (16, 46))
        self.assertEqual(self.rev_seq[46], 'T')
        self.assertEqual(self.rev_seq[33], 'C')


    def test_elongate_forward_sequence(self):
        boundaries = find_N_boundaries(self.src_seq)
        hexamer_table = build_hexamer_table()

        hexamer_pair = determine_hexamers(self.src_seq, boundaries, hexamer_table)
        tst_seq = elongate_forward_sequence(self.src_seq, hexamer_pair[0], "naive_mode")

        self.assertEqual(len(tst_seq), len(self.fwd_seq))
        self.assertEqual(tst_seq, self.fwd_seq)

    def test_elongate_reverse_sequence(self):
        boundaries = find_N_boundaries(self.src_seq)
        hexamer_table = build_hexamer_table()

        hexamer_pair = determine_hexamers(self.src_seq, boundaries, hexamer_table)
        tst_seq = elongate_reverse_sequence(self.src_seq, hexamer_pair[1], "naive_mode")

        self.assertEqual(len(tst_seq), len(self.rev_seq))
        self.assertEqual(tst_seq, self.rev_seq)

    def test_elongate_forward_sequence_kmer_mode(self):
        boundaries = find_N_boundaries(self.src_seq)
        hexamer_table = build_hexamer_table()

        hexamer_pair = determine_hexamers(self.src_seq, boundaries, hexamer_table)
        tst_seq = elongate_forward_sequence(self.src_seq, hexamer_pair[0], "kmer_mode")

        self.assertEqual(len(tst_seq), len(self.fwd_seq))
        self.assertEqual(tst_seq, self.fwd_seq)
 
    def test_elongate_reverse_sequence_kmer_mode(self):
        boundaries = find_N_boundaries(self.src_seq)
        hexamer_table = build_hexamer_table()

        hexamer_pair = determine_hexamers(self.src_seq, boundaries, hexamer_table)
        tst_seq = elongate_reverse_sequence(self.src_seq, hexamer_pair[1], "kmer_mode")

        self.assertEqual(len(tst_seq), len(self.rev_seq))
        self.assertEqual(tst_seq, self.rev_seq)

    def test_elongate_forward_sequence_fixed_length_mode_partial(self):
        tst_seq = elongate_forward_sequence(self.src_seq, HUMAN_TELOMERE, "fixed_length", 6)

        self.assertEqual(len(tst_seq), len(self.fwd_seq))
        self.assertEqual(tst_seq, self.fix_seq)

    def test_elongate_forward_sequence_fixed_length_mode_full(self):
        tst_seq = elongate_forward_sequence(self.src_seq, HUMAN_TELOMERE, "fixed_length", 16)

        self.assertEqual(len(tst_seq), len(self.fwd_seq))
        self.assertEqual(tst_seq, self.fwd_seq)

    def test_elongate_reverse_sequence_fixed_length_mode_partial(self):
        tst_seq = elongate_reverse_sequence(self.src_seq, 'TAACCC', "fixed_length", 6)

        self.assertEqual(len(tst_seq), len(self.rix_seq))
        self.assertEqual(tst_seq, self.rix_seq)

    def test_elongate_reverse_sequence_fixed_length_mode_full(self):
        tst_seq = elongate_reverse_sequence(self.src_seq, 'TAACCC', "fixed_length", 16)

        self.assertEqual(len(tst_seq), len(self.rix_seq))
        self.assertEqual(tst_seq, self.rev_seq)


    def test_determine_hexamer(self):
        boundaries = find_N_boundaries(self.src_seq)
        hexamer_table = build_hexamer_table()

        # Will find known pattern CCCTAA in sequence
        hexamer_pair = determine_hexamers(self.src_seq, boundaries, hexamer_table)
        self.assertEqual('TAACCC', hexamer_pair[0])
        self.assertEqual('TAACCC', hexamer_pair[1])

        # Should not find any pattern in sequence
        hexamer_pair = determine_hexamers(self.no_patt, boundaries, hexamer_table)
        self.assertNotIn('TAACCC', hexamer_pair)
        self.assertEqual(hexamer_pair, [None, None])

if __name__ == '__main__':
    unittest.main()