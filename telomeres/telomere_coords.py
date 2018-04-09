#!/usr/bin/env python

import gzip
import sys
from collections import defaultdict, deque
from itertools import islice, tee
from Bio import SeqIO

length = 300 # How much raw sequence to display?
parse_threshold = 100000 # How much raw sequence to parse through?

# XXX: Do not use coordinates, but (BioPython) iterators instead (.next() and so on)
# XXX: Use brentp's cruzdb to compare sequences side by side
# XXX: Count/plot histogram to see TTAGGG or AATCCC pattern frequency over coordinate spans. Tweak length variable according to that metric.
coordinates = 0
bedfile = defaultdict(list)


def assess_repeats(seq, start, end):
    pattern = __get_pattern__(seq)
    hexamers = 0
    if pattern in window(seq):
        hexamers = hexamers + 1

    return hexamers


# Borrowed from telomerecat's logic:
# https://github.com/jhrf/telomerecat/blob/884c7f830eb2639155113df2c6a7ea4f1154b5fc/telomerecat/telbam2length.py
def __get_pattern__(seq):
    cta, tag = "CCCTAA", "TTAGGG"
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
    # Could use enumerate(islice(iters, 1, None), 1) to avoid consume(it, 0), but that's
    # slower for larger window sizes, while saving only small fixed "noop" cost
    for i, it in enumerate(iters):
        consume(it, i)
    return zip(*iters)


with gzip.open("data/hg38.fa.gz", "rt") as hg38_fa:
    for record in SeqIO.parse(hg38_fa, "fasta"):
        chrom_length = len(record.seq)
        if "_" not in record.id: #avoid extra assemblies
            seq = record.seq
#            for char in seq[:parse_threshold]:
            for char in seq:
                coordinates = coordinates + 1
                if(char == 'N'):
                    continue
                else:
                    start = coordinates
                    end = coordinates+length
                    
                    bedfile[record.id].append(start)
                    hexa = assess_repeats(seq, start, end)
                    
                    print("Forward: Telomere for chrom {chrom} from coord {pos} of {chr_length}".format(chrom=record.id, pos=start, chr_length=chrom_length))
                    print("Sequenc: {seq}".format(seq=seq[start:end]))
                    print("HexCnt : {mers}".format(mers=hexa))

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
                    hexa = assess_repeats(seq, start, end)
                    
                    print("Reverse: Telomere for chrom {chrom} from coord {pos} of {chr_length}".format(chrom=record.id, pos=chrom_length-coordinates, chr_length=chrom_length))
                    print("Sequenc: {seq}".format(seq=seq[start:end]))
                    print("HexCnt : {mers}".format(mers=hexa))
                    
                    coordinates = 0
                    break


print(bedfile)
for k, v in bedfile.items():
    print("{}\t{}\t{}\t".format(k, v[0], v[0]+length))
    print("{}\t{}\t{}\t".format(k, v[1], v[1]+length))