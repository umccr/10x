#!/usr/bin/env python

import re
import sys
import click
import gzip
from typing import List
from pathlib import Path
from collections import defaultdict, deque
from itertools import islice, tee
from Bio import SeqIO
from Bio.Seq import MutableSeq
from Bio.Alphabet import generic_dna

import logging

def find_N_boundary(seq: str):
    ''' Go through the N's until we find some base
    '''
    pos = 0

    for base in seq:
        if 'N' in base:
            pos = pos+1
        else:
            break

    return pos

def determine_hexamer(seq: str):
    if 'ccctaa' in seq:
        return 'ccctaa'
    elif 'ttaggg' in seq:
        return 'ttaggg'

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

def main(genome_build='../../data/processed/hg38_synthetic/new_hg38.fa.gz'):
    hexamer_table = defaultdict(list)

    hexamer_table = build_hexamer_table('ccctaa', hexamer_table)
    hexamer_table = build_hexamer_table('ttaggg', hexamer_table)

    new_seq  = ""

    with gzip.open(genome_build, "rt") as hg38_fa:
        record_dict = SeqIO.to_dict(SeqIO.parse(hg38_fa, "fasta"))
        for _, chrom_attrs in record_dict.items():
            # Original and new (mutable) sequence
            sequence = chrom_attrs.seq
            new_seq = sequence.tomutable()

            N_repeats_pos = find_N_boundary(sequence)
            if N_repeats_pos == 0:
                pass
            else:
                hexamer = sequence[N_repeats_pos:N_repeats_pos+6]
                
                print("{} ... {}".format(sequence[0:10],
                                         sequence[N_repeats_pos-10:N_repeats_pos+10]))

                pos = N_repeats_pos
                pos2 = N_repeats_pos
                while (pos > 0):
                    pos = pos - len(hexamer)
                    new_seq[pos:pos2] = hexamer
                    pos2 = pos2 - len(hexamer)

                print("{} ... {}".format(new_seq[0:10],
                                         new_seq[N_repeats_pos-10:N_repeats_pos+10]))

        
if __name__ == "__main__":
    main()