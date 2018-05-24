#!/usr/bin/env python

import sys
import click
import gzip
from pathlib import Path
from collections import defaultdict, deque
from itertools import islice, tee
from Bio import SeqIO
import logging

# Logging boilerplate
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)

allowed_errs = 1 # How many errors in the telomeric pattern are allowed?
length = 1000 # How much raw sequence to display?
parse_threshold = 1000000 # How much raw sequence to parse through?

# XXX: Do not use coordinates, but (BioPython) iterators instead (.next() and so on)
# XXX: Use brentp's cruzdb to compare sequences side by side?
# XXX: K-mer analysis: Tweak length variable and resulting bedfile according to that metric.
# XXX: chrM    -499    501 # after applying Oliver offsets... treat as offlier?
# XXX: Apply Vlad advise, do enumerate instead of iter sliding window.
# XXX: Coordinate scanning per chromosome, i.e: ./telomere_coords.py hg38.gz --chrom='chr11', useful for debugging purposes

def assess_repeats(seq, start, end):
    seq_s = str(seq).lower()
    # pattern = __get_pattern__(seq_s) #XXX: only takes one pattern, check telomerecat source more
    # Search for both patterns instead
    pattern1 = 'ccctaa'
    pattern3 = 'gggatt'
    
    pattern2 = 'ttaggg'
    pattern4 = 'aatccc'

    # as observed in chr5
    pattern5 = 'taaccc'

    seq_s = seq_s[start:end]
    
    hits = 0
    partial_hits = 0

    for kmer in window(seq_s, len(pattern1)):
        kmer_s = ''.join(kmer) # fuse tuple back to string
        if pattern1 in kmer_s:
            hits = hits + 1
        elif pattern2 in kmer_s:
            hits = hits + 1
        elif str(pattern3) in kmer_s:
            hits = hits + 1
        elif str(pattern4) in kmer_s:
            hits = hits + 1
        elif str(pattern5) in kmer_s:
            hits = hits + 1
        elif str(reversed(pattern2)) in kmer_s:
            hits = hits + 1
        elif str(reversed(pattern1)) in kmer_s:
            hits = hits + 1
        elif str(reversed(pattern3)) in kmer_s:
            hits = hits + 1
        elif str(reversed(pattern4)) in kmer_s:
            hits = hits + 1
        elif str(reversed(pattern5)) in kmer_s:
            hits = hits + 1
        else:
            if distance(pattern1, seq) == allowed_errs:
                partial_hits = partial_hits + 1
            if distance(pattern2, seq) == allowed_errs:
                partial_hits = partial_hits + 1

            # XXX: Overriding the above logic with just the distance
            partial_hits = distance(pattern1, kmer_s)

    return hits, partial_hits

# Borrowed from telomerecat's logic:
# https://github.com/jhrf/telomerecat/blob/884c7f830eb2639155113df2c6a7ea4f1154b5fc/telomerecat/telbam2length.py
def __get_pattern__(seq):
    cta, tag = "ccctaa", "ttaggg"
    pattern = None
    if cta in seq or tag in seq:
        if seq.count(cta) > seq.count(tag):
            pattern = cta
        else:
            pattern = tag
    return pattern

def distance(a,b):
    count = 0
    for i in range (0,len(a)):
        if a[i] != b[i]:
            count += 1
    return count

def consume(iterator, n):
    "Advance the iterator n-steps ahead. If n is none, consume entirely."
    # Use functions that consume iterators at C speed.
    if n is None:
        # feed the entire iterator into a zero-length deque
        deque(iterator, maxlen=0)
    else:
        # advance to the empty slice starting at position n
        next(islice(iterator, n, n), None)

def window(iterable, n=6):
    "s -> (s0, ...,s(n-1)), (s1, ...,sn), (s2, ..., s(n+1)), ..."
    iters = tee(iterable, n)

    for i, it in enumerate(iters):
        consume(it, i)
    return zip(*iters)

def scan_record(record, direction):
    coordinates = 0
    sequence = record.seq
    oliver_offset = 500 # "Can you calculate the number of hexamers around your established transition coordinates for, say, 500bp in each direction?"
    chrom_length = len(sequence)

    for char in sequence:
    for char in reversed(seq):

        coordinates = coordinates + 1
        if(char == 'N'):
            continue
        else:
            start = coordinates-oliver_offset
            end = coordinates+length

            # XXX: Reverse
            #            start = chrom_length-coordinates-length-oliver_offset
            #            end = chrom_length-coordinates

            
            bedfile[record.id].append(start)
            hits, partial_hits = assess_repeats(seq, start, end)
            
            log.info("{direction}: Telomere for chrom {chrom} from coord {pos} of {chr_length}".format(direction='forward', chrom=record.id, 
                                                                                                       pos=start, chr_length=chrom_length))
                                                                                                       # XXX: Reverse: format(chrom=record.id, pos=chrom_length-coordinates, chr_length=chrom_length))
            log.info("Sequenc: {seq}".format(seq=seq[start:end]))
            log.info("HexCnt for {chrom}: {mers}, {phits}".format(chrom=record.id, mers=hits, phits=partial_hits))

            coordinates = 0
            hits = 0
            break

    return bedline

@click.command()
@click.argument('genome_build', type=click.Path(exists=True))
@click.option('--only-chr', 'only_chr', help='Only processes the specified chromosome')
def main(genome_build='data/external/hg38.fa.gz', only_chr=None):
    bedfile = defaultdict(list)

    with gzip.open(genome_build, "rt") as hg38_fa:
        record_dict = SeqIO.to_dict(SeqIO.parse(hg38_fa, "fasta"))
        for record in record_dict:
            if "_" not in record.id: #avoid extra assemblies
                if "chrM" not in record.id: #skip Mitochondrial DNA (circular so no point to search for telomeres)
                    bedrow = scan_record(record, 'forward')
                    bedfile[record.id].append(bedrow)

        with open(Path('data/processed/telomere_coords.bed'), "w+") as tel:
            for k, v in bedfile.items():
                tel.write("{}\t{}\t{}\n".format(k, v[0], v[0]+length)) # Forward
                tel.write("{}\t{}\t{}\n".format(k, v[1], v[1]+length)) # Reverse

if __name__ == "__main__":
    """ This is executed when run from the command line """
    main()