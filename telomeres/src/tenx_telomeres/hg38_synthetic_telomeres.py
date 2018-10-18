#!/usr/bin/env python

import math
import gzip
from typing import List
from pathlib import Path
from collections import defaultdict, deque
from itertools import islice, tee
from pymer import ExactKmerCounter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

import logging

# Set log level
loglevel = logging.INFO
logging.basicConfig(level=loglevel)
log = logging.getLogger(__name__)


## Global
FA_IDX = "hg38.fa.idx"
O_OFFSET = 1000

# Telomeric hexamer
KMER_K = 6

# Human telomeric hexamers and complementary sequences
HUMAN_TELOMERE = 'TTAGGG'
TELO_HEXAMERS = defaultdict(list)

# Seed hexamers with all possible orientations
TELO_HEXAMERS[HUMAN_TELOMERE] = [HUMAN_TELOMERE, str(Seq(HUMAN_TELOMERE, generic_dna).complement()), 
                                                 str(Seq(HUMAN_TELOMERE, generic_dna).reverse_complement())]


def find_N_boundaries(seq: str):
    ''' Returns all N-boundaries in a sequence via tuple: (first, second)
    '''
    pos = first = second = 0

    # first N stretch
    for base in seq:
        if 'N' in base:
            pos = pos + 1
        else:
            first = pos
            break

    base = None
    pos = 0

    # last N stretch
    for base in reversed(seq):
        if 'N' in base:
            pos = pos + 1
        else:
            second = len(seq) - pos - 1
            break

    return (first, second)

# XXX: Generalize/merge both elongate functions
# Elongate forward and backward N's, respecting telomeric patterns
def elongate_forward_sequence(seq: str, kmer: str, mode: str):
    # Determine N boundaries in the sequence
    boundary, boundary_r = find_N_boundaries(seq)

    # K-mer telomeric sequence right after the N boundary
    kmer_seq = seq[boundary:boundary + KMER_K]

    # How many chunks to elongate and remainder
    chunks = len(seq[0:boundary]) % KMER_K
    chunks_r = len(seq[0:boundary]) / KMER_K

    if mode == "naive_mode":
        # Just capture remainder of the pattern from the boundary to fit in sequence
        kmer_seq_r = kmer_seq[math.floor(chunks_r):]
        tmp_seq = kmer_seq_r
    elif mode == "kmer_mode":
        # XXX: fairly blunt kmer/pattern transition here
        if kmer[0] is not None:
            tmp_seq = kmer[0]
            kmer_seq = kmer[0]
        else: # just leave N's as they are since no suitable telomeric kmer was found
            tmp_seq = 'N' * KMER_K

    # Build forward sequence
    for _ in range(0, chunks - 2): # XXX 2?
        tmp_seq = tmp_seq + kmer_seq

    # Attach inner pattern
    tmp_seq = tmp_seq + seq[boundary:boundary_r] + seq[boundary_r:]

    return tmp_seq

def elongate_reverse_sequence(seq: str, kmer: str, mode: str):
    # Determine N boundaries in the sequence
    boundary, boundary_r = find_N_boundaries(seq)

    # K-mer telomeric sequence right before the N boundary on the reverse side
    kmer_seq = seq[boundary_r - KMER_K:boundary_r]

    # How many chunks to elongate and remainder
    chunks = len(seq[boundary_r:]) % KMER_K
    chunks_r = len(seq[boundary_r:]) / KMER_K

    # Start with bases present in the forward side
    tmp_seq = seq[0:boundary]

    # Attach inner pattern
    tmp_seq = tmp_seq + seq[boundary:boundary_r]

    if mode == "naive_mode":
        # Just capture remainder of the pattern from the boundary to fit in sequence
        kmer_seq_r = kmer_seq[math.floor(chunks_r):]
        tmp_seq = kmer_seq_r
    elif mode == "kmer_mode":
        # XXX: fairly blunt kmer/pattern transition here
        if kmer[1] is not None:
            tmp_seq = kmer[1]
            kmer_seq = kmer[1]
        else: # just leave N's as they are since no suitable telomeric kmer was found
            tmp_seq = 'N' * KMER_K

    # Build reverse sequence
    for _ in range(0, chunks):
        tmp_seq = tmp_seq + kmer_seq

    # Capture remainder of the pattern to fit in sequence
    kmer_seq_r = kmer_seq[0:math.floor(chunks_r)]
    tmp_seq = tmp_seq + kmer_seq_r

    return tmp_seq

def build_hexamer_table():
    '''
    Builds a table containing hexamers and all its possible rotations.
    
    Useful to determine boundary conditions between N-regions and telomeric
    repeats on the reference genome(s).
    '''
    hexamer_table = defaultdict(list)
    rotated = []

    # Rotate the telomeric pattern to match boundaries
    for pattern in TELO_HEXAMERS[HUMAN_TELOMERE]:
        dq = deque(pattern)
        for rot in range(1, len(pattern)):
            dq.rotate(rot)
            rotated.append(''.join(dq))

        hexamer_table[pattern] = rotated

    return hexamer_table


def determine_hexamers(seq: str, boundaries: tuple, hexamer_table: dict):
    '''
    Takes the sequence seq and finds telomeric kmers in it. Keeps the count
    of found kmers in each direction.
    '''
    fwd_boundary, rev_boundary = boundaries
    fwd_detected, rev_detected = defaultdict(), defaultdict()

    kc_fwd = ExactKmerCounter(KMER_K)
    kc_rev = ExactKmerCounter(KMER_K) 
    
    kc_fwd.consume(str(seq[fwd_boundary:fwd_boundary + O_OFFSET]))
    kc_rev.consume(str(seq[rev_boundary - O_OFFSET:rev_boundary]))

    for _, v in hexamer_table.items():
        for telo in v:
            if kc_fwd[telo] != 0:
                fwd_detected[telo] = kc_fwd[telo]
            if kc_rev[telo] != 0:
                rev_detected[telo] = kc_rev[telo]

    # find most frequent one and use it as best choice for elongate later
    total_detected = [None, None]

    # ValueError: max() arg is an empty sequence
    try:
        total_detected[0] = max(fwd_detected)
    except ValueError:
        pass

    try:
        total_detected[1] = max(rev_detected)
    except ValueError:
        pass

    return total_detected

def fasta_idx(filename):
    ''' Indexes a fasta filename, since SeqIO.to_dict is not very efficient for
        big files, see: https://github.com/biopython/biopython/pull/959 and
        related issues.
    '''
    with gzip.open(filename, 'wt') as hg38_idx:
        SeqIO.index_db(filename, hg38_idx, 'fasta')


#def main(genome_build='data/external/hg38.fa.gz'):
def main(genome_build='data/external/chr11.fa.gz'):

    final_seq = None
    hexamer_table = build_hexamer_table()

    with gzip.open(genome_build, "rt") as hg38_fa:
        record_dict = SeqIO.to_dict(SeqIO.parse(hg38_fa, "fasta"))
        for _, chrom_attrs in record_dict.items():
            sequence = chrom_attrs.seq
            seq_id = chrom_attrs.id
            detected_hexamer_pair = [None, None]

            # Discard _KI_random and _alt assemblies (filter "_"). Also disregard chrM
            # since there are no biologically relevant telomeres there (circular sequence).
            if "_" not in seq_id:
                if "chrM" not in seq_id:
                    fwd_boundary, rev_boundary = find_N_boundaries(sequence)
                    detected_hexamer_pair = determine_hexamers(sequence, (fwd_boundary, rev_boundary), hexamer_table)

                    print("{}\t{}:\t\t{}\t...\t{}\t...\t{}\t{}".format(seq_id.split(':')[0],
                                                                      (fwd_boundary, rev_boundary),
                                                                      sequence[fwd_boundary - 3:fwd_boundary + KMER_K + 10],
                                                                      sequence[rev_boundary - KMER_K - 10:rev_boundary + 4],
                                                                      len(sequence), detected_hexamer_pair))

            # Finally, build the synthetically elongated hg38 build

            final_seq = elongate_forward_sequence(sequence, detected_hexamer_pair, "kmer_mode")
            final_seq = elongate_reverse_sequence(final_seq, detected_hexamer_pair, "kmer_mode")

    final_seq = Seq(final_seq, generic_dna) # convert str to BioPython Seq
    with open("hg38_elongated_telomeres.fasta", "w") as output_handle:
        SeqIO.write(final_seq, output_handle, "fasta")

if __name__ == "__main__":
    main()