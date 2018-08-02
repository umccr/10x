#!/usr/bin/env python

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


def main(genome_build='../../data/processed/hg38_synthetic/new_hg38.fa.gz'):
    new_seq  = ""
    cnt = 0

    with gzip.open(genome_build, "rt") as hg38_fa:
        record_dict = SeqIO.to_dict(SeqIO.parse(hg38_fa, "fasta"))
        for _, chrom_attrs in record_dict.items():
            # Original and new (mutable) sequence
            sequence = chrom_attrs.seq
            chrom_length = len(sequence)
            new_seq = sequence.tomutable()

            N_repeats_pos = find_N_boundary(sequence)
            if N_repeats_pos == 0:
                pass
            else:
                hexamer = sequence[N_repeats_pos:N_repeats_pos+6]
                
                print("{} ... {} ... {}".format(sequence[0:10],
                                                sequence[N_repeats_pos-10:N_repeats_pos+10],
                                                sequence[chrom_length-20:chrom_length]))

                pos = N_repeats_pos
                pos2 = N_repeats_pos
                while (pos > 4 or pos2 > 4):
                    pos = pos - len(hexamer)
                    new_seq[pos:pos2] = hexamer
                    pos2 = pos2 - len(hexamer)
                    cnt = cnt+1

                print("{} ... {} ... {}\t\t\t{}".format(new_seq[0:10],
                                                new_seq[N_repeats_pos-10:N_repeats_pos+10],
                                                new_seq[chrom_length-20:chrom_length],
                                                [N_repeats_pos, pos, pos2, chrom_length, cnt]))
        print (cnt)
        
if __name__ == "__main__":
    main()