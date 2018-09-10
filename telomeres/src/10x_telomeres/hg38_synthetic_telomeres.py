#!/usr/bin/env python

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

pattern1 = 'ccctaa'
pattern2 = 'ttaggg'

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

            N_repeats_pos = find_N_boundary(sequence)

            if N_repeats_pos == 0: # sequence with all N's
                pass
            else:
                hexamer = sequence[N_repeats_pos:N_repeats_pos+6]
                #XXX: Make sure this works in reverse?
                detected_hexamer = determine_hexamer(new_seq[N_repeats_pos:N_repeats_pos+len(hexamer)], hextable)
                print("Detected hexamer: {}".format(detected_hexamer))

                if detected_hexamer is None:
                    print("Cannot detect telomeric hexamer within sequence")
                else:
                    hexamer = detected_hexamer

                print("{} ... {} ... {}".format(sequence[0:10],                                 # first chars of sequence
                                                sequence[N_repeats_pos-10:N_repeats_pos+10],    # Chars around the (filled?) N boundary
                                                sequence[chrom_length-20:chrom_length]))        # last chars of sequence

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
                                                        [N_repeats_pos, pos, pos2, chrom_length, cnt])) # "pointer" information


                new_hg38.append(SeqRecord(new_seq.toseq(), id=seq_id, name=seq_id, description=seq_id))

        with open("hg38_elongated_telomeres.fasta", "w") as output_handle:
            SeqIO.write(new_hg38, output_handle, "fasta")
        
if __name__ == "__main__":
    main()