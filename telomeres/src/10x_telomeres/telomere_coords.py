#!/usr/bin/env python

import sys
import gzip
from collections import defaultdict, deque
from itertools import islice, tee
from Bio import SeqIO

allowed_errs = 1 # How many errors in the telomeric pattern are allowed?
length = 1000 # How much raw sequence to display?
oliver_offset = 500 # "Can you calculate the number of hexamers around your established transition coordinates for, say, 500bp in each direction?"
parse_threshold = 1000000 # How much raw sequence to parse through?

# XXX: Do not use coordinates, but (BioPython) iterators instead (.next() and so on)
# XXX: Use brentp's cruzdb to compare sequences side by side?
# XXX: K-mer analysis: Tweak length variable and resulting bedfile according to that metric.
# XXX: chrM    -499    501 # after applying Oliver offsets... treat as offlier?
# XXX: Apply Vlad advise, do enumerate instead of iter sliding window.
coordinates = 0
bedfile = defaultdict(list)

def assess_repeats(seq, start, end):
    seq_s = str(seq).lower()
    # pattern = __get_pattern__(seq_s) #XXX: only takes one pattern, check telomerecat source more
    # Search for both patterns instead
    pattern1 = 'ccctaa'
    pattern3 = 'gggatt'
    
    pattern2 = 'ttaggg'
    pattern4 = 'aatccc'
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
        elif str(reversed(pattern2)) in kmer_s:
            hits = hits + 1
        elif str(reversed(pattern1)) in kmer_s:
            hits = hits + 1
        elif str(reversed(pattern3)) in kmer_s:
            hits = hits + 1
        elif str(reversed(pattern4)) in kmer_s:
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


#with gzip.open("data/hg38.fa.gz", "rt") as hg38_fa:
with gzip.open("data/GRCh37.fa.gz", "rt") as hg38_fa:
    for record in SeqIO.parse(hg38_fa, "fasta"):
        chrom_length = len(record.seq)
        if "_" not in record.id: #avoid extra assemblies
            if "chrM" not in record.id: #skip Mitochondrial DNA
                seq = record.seq

    #XXX: Remove code duplication here
                for char in seq:
                    coordinates = coordinates + 1
                    if(char == 'N'):
                        continue
                    else:
                        start = coordinates-oliver_offset
                        end = coordinates+length
                        
                        bedfile[record.id].append(start)
                        hits, partial_hits = assess_repeats(seq, start, end)
                        
                        print("Forward: Telomere for chrom {chrom} from coord {pos} of {chr_length}".format(chrom=record.id, pos=start, chr_length=chrom_length))
                        print("Sequenc: {seq}".format(seq=seq[start:end]))
                        print("HexCnt for {chrom}: {mers}, {phits}".format(chrom=record.id, mers=hits, phits=partial_hits))

                        coordinates = 0
                        hits = 0
                        break

                # in reverse now
                for char in reversed(seq):

                    coordinates = coordinates + 1
                    if(char == 'N'):
                        continue
                    else:
                        start = chrom_length-coordinates-length-oliver_offset
                        end = chrom_length-coordinates

                        bedfile[record.id].append(start)
                        hits, partial_hits = assess_repeats(seq, start, end)
                        
                        print("Reverse: Telomere for chrom {chrom} from coord {pos} of {chr_length}".format(chrom=record.id, pos=chrom_length-coordinates, chr_length=chrom_length))
                        print("Sequenc: {seq}".format(seq=seq[start:end]))
                        print("HexCnt for {chrom}: {mers}, {phits}".format(chrom=record.id, mers=hits, phits=partial_hits))
                        
                        coordinates = 0
                        hits = 0
                        break

print(bedfile)
for k, v in bedfile.items():
    print("{}\t{}\t{}\t".format(k, v[0], v[0]+length)) # Forward
    print("{}\t{}\t{}\t".format(k, v[1], v[1]+length)) # Reverse
