#!/usr/bin/env python

import sys
import gzip
from collections import defaultdict, deque
from itertools import islice, tee
from Bio import SeqIO

length = 1000 # How much raw sequence to display?
parse_threshold = 1000000 # How much raw sequence to parse through?

# XXX: Do not use coordinates, but (BioPython) iterators instead (.next() and so on)
# XXX: Use brentp's cruzdb to compare sequences side by side?
# XXX: K-mer analysis: Tweak length variable and resulting bedfile according to that metric.
coordinates = 0
bedfile = defaultdict(list)


def assess_repeats(seq, start, end):
    seq_s = str(seq).lower()
    pattern = __get_pattern__(seq_s)
    seq_s = seq_s[start:end]
    
    hits = 0

    for kmer in window(seq_s, len(pattern)):
        kmer_s = ''.join(kmer)
        if pattern in kmer_s:
            hits = hits + 1
    
    return hits, pattern

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


with gzip.open("data/hg38.fa.gz", "rt") as hg38_fa:
    for record in SeqIO.parse(hg38_fa, "fasta"):
        chrom_length = len(record.seq)
        if "_" not in record.id: #avoid extra assemblies
            seq = record.seq

#XXX: Remove code duplication here
            for char in seq:
                coordinates = coordinates + 1
                if(char == 'N'):
                    continue
                else:
                    start = coordinates
                    end = coordinates+length
                    
                    bedfile[record.id].append(start)
                    hits, pattern = assess_repeats(seq, start, end)
                    
                    print("Forward: Telomere for chrom {chrom} from coord {pos} of {chr_length}".format(chrom=record.id, pos=start, chr_length=chrom_length))
                    print("Sequenc: {seq}".format(seq=seq[start:end]))
                    print("Telpatt: {patt}".format(patt=pattern))
                    print("HexCnt : {mers}".format(mers=hits))

                    coordinates = 0
                    break

            # in reverse now
            for char in reversed(seq):
                coordinates = coordinates + 1
                if(char == 'N'):
                    continue
                else:
                    start = chrom_length-coordinates-length
                    end = chrom_length-coordinates

                    bedfile[record.id].append(start)
                    hits, pattern = assess_repeats(seq, start, end)
                    
                    print("Reverse: Telomere for chrom {chrom} from coord {pos} of {chr_length}".format(chrom=record.id, pos=chrom_length-coordinates, chr_length=chrom_length))
                    print("Sequenc: {seq}".format(seq=seq[start:end]))
                    print("Telpatt: {patt}".format(patt=pattern))
                    print("HexCnt : {mers}".format(mers=hits))
                    
                    coordinates = 0
                    break

print(bedfile)
for k, v in bedfile.items():
    print("{}\t{}\t{}\t".format(k, v[0], v[0]+length))
    print("{}\t{}\t{}\t".format(k, v[1], v[1]+length))