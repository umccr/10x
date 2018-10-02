#!/usr/bin/env python

import math
import gzip
from typing import List
from pathlib import Path
from collections import defaultdict, deque
from itertools import islice, tee
from Bio import SeqIO
from Bio.Seq import MutableSeq, Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

import logging

# Telomeric hexamer
kmer_k = 6

# Human telomeric hexamers
pattern1 = 'ccctaa'
pattern2 = 'ttaggg'

def find_N_boundaries(seq: str):
    ''' Returns all N-boundaries in a sequence via tuple: (first, second)
    '''
    first_found = False
    pos = first = second = 0

    for base in seq:
        if 'N' in base:
            # first N boundary found
            pos = pos+1
        elif not first_found:
            first_found = True
            first = pos
            pos = pos+1
        elif first_found:
            second = pos
            pos = pos+1

    return (first, second)

def elongate_forward_sequence(seq):
    # Determine N boundaries in the sequence
    boundary, boundary_r = find_N_boundaries(seq)

    # K-mer telomeric sequence right after the N boundary
    kmer_seq = seq[boundary:boundary + kmer_k]

    # How many chunks to elongate and remainder
    chunks = len(seq[0:boundary]) % kmer_k
    chunks_r = len(seq[0:boundary]) / kmer_k

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
    kmer_seq = seq[boundary_r - kmer_k:boundary_r]

    # How many chunks to elongate and remainder
    chunks = len(seq[boundary_r:]) % kmer_k
    chunks_r = len(seq[boundary_r:]) / kmer_k

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



def determine_hexamer(hexamer: str, hextable: defaultdict):
    ''' Takes the sequence seq and tries to find which hexamer pattern it has
    '''
    for k, v in hextable.items():
        for kmer in v:
            if str(hexamer) in str(kmer):
                return k
    return None

def build_hexamer_table(hexamer: str, hexamer_table: dict):
    ''' Builds a table containing hexamers and all its possible rotations.
        Useful to determine boundary conditions between N-regions and telomeric
        repeats on the reference genome(s).
    '''
    # Seed table with first non-rotated "canonical" hexamer
    hexamer_table[hexamer] = hexamer

    dq = deque(hexamer)
    rotated = []
    for rot in range(1, len(hexamer)):
        dq.rotate(rot)
        rotated.append(''.join(dq))
    
    hexamer_table[hexamer] = rotated

    return hexamer_table

def main(genome_build='data/processed/hg38_synthetic/new_hg38.fa.gz'):
    new_seq  = ""
    new_hg38 = [] # just a list of records
    hextable = defaultdict(list)

    hextable = build_hexamer_table(pattern1, hextable)
    hextable = build_hexamer_table(pattern2, hextable)

    cnt = 0

    with gzip.open(genome_build, "rt") as hg38_fa:
        record_dict = SeqIO.to_dict(SeqIO.parse(hg38_fa, "fasta"))
        for _, chrom_attrs in record_dict.items():
            # Original and new_seq (mutable) sequence
            sequence = chrom_attrs.seq
            seq_id = chrom_attrs.id

            chrom_length = len(sequence)
            new_seq = sequence.tomutable()

            # XXX: Continue with tests, fixing this code for good
            N_repeats_pos = find_N_boundaries(sequence)

            if N_repeats_pos == 0: # sequence with all N's
                pass
            else:
                hexamer = sequence[N_repeats_pos:N_repeats_pos+6]
                #XXX: Make sure this works in reverse?
                detected_hexamer = determine_hexamer(new_seq[N_repeats_pos:N_repeats_pos+len(hexamer)], hextable)
                print("Detected hexamer: {}".format(detected_hexamer))

                if detected_hexamer is None:
                    print("Cannot detect telomeric hexamer within sequence, skipping")
                    continue
                else:
                    hexamer = detected_hexamer

                print("{} ... {} ... {}".format(sequence[0:10],                                 # first chars of sequence
                                                sequence[N_repeats_pos-10:N_repeats_pos+10],    # Chars around the (filled?) N boundary
                                                sequence[chrom_length-20:chrom_length]))        # last chars of sequence

                pos = N_repeats_pos
                pos2 = N_repeats_pos

                # XXX: Why 4? Re-visit invariant/thinking behind this
                while (pos > 4 or pos2 > 4):
                    pos = pos - len(hexamer)
                    new_seq[pos:pos2] = hexamer
                    pos2 = pos2 - len(hexamer)
                    cnt = cnt+1

                print("{} ... {} ... {}\t\t\t{}".format(new_seq[0:10],
                                                        new_seq[N_repeats_pos-10:N_repeats_pos+10],
                                                        new_seq[chrom_length-20:chrom_length],
                                                        [N_repeats_pos, pos, pos2, chrom_length, cnt])) # "pointer" information


                new_hg38.append(SeqRecord(new_seq.toseq(), id=seq_id, name=seq_id, description=seq_id))

        with open("hg38_elongated_telomeres.fasta", "w") as output_handle:
            SeqIO.write(new_hg38, output_handle, "fasta")
        
if __name__ == "__main__":
    main()