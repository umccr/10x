#!/usr/bin/env python

import sys
import click
import gzip
from typing import List
from pathlib import Path
from collections import defaultdict, deque
from itertools import islice, tee
from Bio import SeqIO
import logging

# Logging boilerplate
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)

motif_size = 6
oliver_offset = 500 # "Can you calculate the number of hexamers around your established transition coordinates for, say, 500bp in each direction?"


def assess_repeats(seq: str):
    seq_s = str(seq).lower()
    
    pattern1 = 'ccctaa'
    pattern3 = 'gggatt'
    pattern2 = 'ttaggg'
    pattern4 = 'aatccc'
    pattern5 = 'taaccc'

    hits = 0

    for kmer in window(seq_s, len(pattern1)):
        if pattern1 in kmer:
            hits = hits + 1
        elif pattern2 in kmer:
            hits = hits + 1
        elif str(pattern3) in kmer:
            hits = hits + 1
        elif str(pattern4) in kmer:
            hits = hits + 1
        elif str(pattern5) in kmer:
            hits = hits + 1
        elif str(reversed(pattern1)) in kmer:
            hits = hits + 1
        elif str(reversed(pattern2)) in kmer:
            hits = hits + 1
        elif str(reversed(pattern3)) in kmer:
            hits = hits + 1
        elif str(reversed(pattern4)) in kmer:
            hits = hits + 1
        elif str(reversed(pattern5)) in kmer:
            hits = hits + 1

    return hits

# https://stackoverflow.com/questions/6822725/rolling-or-sliding-window-iterator
def consume(iterator, n):
    "Advance the iterator n-steps ahead. If n is none, consume entirely."
    # Use functions that consume iterators at C speed.
    if n is None:
        # feed the entire iterator into a zero-length deque
        deque(iterator, maxlen=0)
    else:
        # advance to the empty slice starting at position n
        next(islice(iterator, n, n), None)

def window(iterable, n=motif_size):
    "s -> (s0, ...,s(n-1)), (s1, ...,sn), (s2, ..., s(n+1)), ..."
    iters = tee(iterable, n)

    for i, it in enumerate(iters):
        consume(it, i)
    return zip(*iters)

def rev_string(s: str): 
    return s[::-1]

def scan_record(record: SeqIO, direction: str):
    sequence = record.seq
    chrom_length = len(sequence)
    hits = 0
    pos = 0 
    telomere_found = False # flags first occurrence/hit

    # chromosome, pos_start, pos_end
    bedrow = []
    bedrow.append(record.name)

    if 'reverse' in direction:
        sequence = rev_string(sequence)

    # Go through the N's until we find some base
    for base in sequence:
        if 'N' in base:
            pos = pos+1
        else:
            break
    
    # Fix reverse coordinates
    if 'reverse' in direction:
        pos = chrom_length - pos

    # start
    bedrow.append(pos)

    # Test the kmers for patterns, store the position and count hits
    for hexamer in window(sequence[pos:pos+oliver_offset]):
        hits = assess_repeats(''.join(hexamer)) # fuse tuple back to string
        pos = pos+1
    
    # Fix reverse coordinates
    if 'reverse' in direction:
        pos = chrom_length

    # end
    bedrow.append(pos)

    return bedrow, hits

@click.command()
@click.argument('genome_build', type=click.Path(exists=True))
@click.option('--only-chr', 'only_chr', help='Only processes the specified chromosome')
def main(genome_build='data/external/hg38.fa.gz', only_chr=None):
    bedfile = defaultdict(list)
    hits = 0

    with gzip.open(genome_build, "rt") as hg38_fa:
        record_dict = SeqIO.to_dict(SeqIO.parse(hg38_fa, "fasta"))
        for chrom_name, chrom_attrs in record_dict.items():
            if "_" not in chrom_name: #avoid extra assemblies
                if "chrM" not in chrom_name: #skip Mitochondrial DNA (circular so no point to search for telomeres)
                    #bedrow, hits = scan_record(chrom_attrs, 'forward')
                    bedrow, hits = scan_record(chrom_attrs, 'reverse')
                    print("{bedrow} ## hits: {hits}".format(bedrow=bedrow,
                                                               hits=hits))
                    bedfile[chrom_name].append(bedrow)

if __name__ == "__main__":
    main()