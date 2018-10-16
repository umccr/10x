#!/usr/bin/env python

import math
import gzip
from typing import List
from pathlib import Path
from collections import defaultdict, deque
from itertools import islice, tee
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

# XXX: Try to generalize/merge both elongate functions
# Elongate forward and backward N's, respecting telomeric patterns
def elongate_forward_sequence(seq):
    # Determine N boundaries in the sequence
    boundary, boundary_r = find_N_boundaries(seq)

    # K-mer telomeric sequence right after the N boundary
    kmer_seq = seq[boundary:boundary + KMER_K]

    # How many chunks to elongate and remainder
    chunks = len(seq[0:boundary]) % KMER_K
    chunks_r = len(seq[0:boundary]) / KMER_K

    # Capture remainder of the pattern to fit in sequence
    kmer_seq_r = kmer_seq[math.floor(chunks_r):]
    tmp_seq = kmer_seq_r
    
    # Build forward sequence
    for _ in range(0, chunks - 2): # XXX 2?
        tmp_seq = tmp_seq + kmer_seq

    # Attach inner pattern
    tmp_seq = tmp_seq + seq[boundary:boundary_r] + seq[boundary_r:]

    return tmp_seq

def elongate_reverse_sequence(seq):
    # Determine N boundaries in the sequence
    boundary, boundary_r = find_N_boundaries(seq)

    # K-mer telomeric sequence right before the N boundary
    kmer_seq = seq[boundary_r - KMER_K:boundary_r]

    # How many chunks to elongate and remainder
    chunks = len(seq[boundary_r:]) % KMER_K
    chunks_r = len(seq[boundary_r:]) / KMER_K

    # Start with the N boundary
    tst_seq = seq[0:boundary]

    # Attach inner pattern
    tst_seq = tst_seq + seq[boundary:boundary_r]

    # Build reverse sequence
    for _ in range(0, chunks):
        tst_seq = tst_seq + kmer_seq

    # Capture remainder of the pattern to fit in sequence
    kmer_seq_r = kmer_seq[0:math.floor(chunks_r)]
    tst_seq = tst_seq + kmer_seq_r

    return tst_seq


def determine_hexamer(seq: str, boundaries: tuple):
    ''' 
    Builds a table containing hexamers and all its possible rotations.
    
    Useful to determine boundary conditions between N-regions and telomeric
    repeats on the reference genome(s).

    Also takes the sequence seq and tries to find which hexamer pattern it has
    '''
    hexamer_table = defaultdict(list)
    fwd_boundary, rev_boundary = boundaries
    rotated = []

    # Rotate the telomeric pattern to match boundaries
    for pattern in TELO_HEXAMERS[HUMAN_TELOMERE]:
        dq = deque(pattern)
        for rot in range(1, len(pattern)):
            dq.rotate(rot)
            rotated.append(''.join(dq))

        hexamer_table[pattern] = rotated

    # Scan for known telomeric sequences around boundaries
    # XXX: Should detect forward *and* reverse sequences, not halt at the first where it finds them!
    for _, v in hexamer_table.items():
        for kmer in v:
            if kmer in str.upper(str(seq[0:fwd_boundary + O_OFFSET])):        # fwd
                return kmer, "fwd"
            elif kmer in str.upper(str(seq[rev_boundary - O_OFFSET:])):       # rev
                return kmer, "rev"
            else:
                return None, None

    # Should never end up here (tm)
    raise ValueError("Error in parsing telomeric subsequences")

def fasta_idx(filename):
    ''' Indexes a fasta filename, since SeqIO.to_dict is not very efficient for
        big files, see: https://github.com/biopython/biopython/pull/959 and
        related issues.
    '''
    with gzip.open(filename, 'wt') as hg38_idx:
        SeqIO.index_db(filename, hg38_idx, 'fasta')


#def main(genome_build='../../data/external/hg38.fa.gz'):
def main(genome_build='../../data/external/chr11.fa.gz'):
    with gzip.open(genome_build, "rt") as hg38_fa:
        record_dict = SeqIO.to_dict(SeqIO.parse(hg38_fa, "fasta"))
        for _, chrom_attrs in record_dict.items():
            sequence = chrom_attrs.seq
            seq_id = chrom_attrs.id
            detected_hexamer = None

            # Discard _KI_random and _alt assemblies, disregard chrM too
            # since there are no relevant telomeres there (circular sequence)
            if "_" not in seq_id:
                if "chrM" not in seq_id:
                    fwd_boundary, rev_boundary = find_N_boundaries(sequence)
                    detected_hexamer = determine_hexamer(sequence, (fwd_boundary, rev_boundary))

                    #chr1    (10000, 248946421):             taaccctaaccctaaccctaaccctaaccctaaccc    ...     ttagggttagggttagggttaagggttagggttagg    ...     248956422       (TTAGGG, rev)
                    print("{}\t{}:\t\t{}\t...\t{}\t...\t{}\t{}".format(seq_id.split(':')[0],
                                                                      (fwd_boundary, rev_boundary),
                                                                      sequence[fwd_boundary:fwd_boundary + KMER_K + 30],
                                                                      sequence[rev_boundary - KMER_K - 30:rev_boundary],
                                                                      len(sequence), detected_hexamer))

            #final_seq = elongate_forward_sequence(sequence)
            #final_seq = elongate_reverse_sequence(final_seq)

#        with open("hg38_elongated_telomeres.fasta", "w") as output_handle:
#            SeqIO.write(new_hg38, output_handle, "fasta")

if __name__ == "__main__":
    main()