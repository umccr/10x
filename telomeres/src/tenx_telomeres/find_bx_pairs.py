#!/usr/bin/env python
import sys
import pysam
from collections import defaultdict

CHR_FOCUS = 'chr5'
linked_reads = defaultdict(list)
samfile = pysam.AlignmentFile("telomeres/data/interim/chr5_hg38_elongated.bam", "rb")

'''
Use BX tag to match linked reads: https://support.10xgenomics.com/genome-exome/software/pipelines/latest/output/bam
Collects all the possible spans and makes a histogram
'''
# XXX: Filter by a reasonable MAPQ value
# XXX: Filter for max 10X molecule length
for read in samfile.fetch(CHR_FOCUS, 10, 20000):
    start = read.reference_start    # current read coordinates
    end = read.next_reference_start # mate read coordinates
    if read.has_tag('BX'):
        if read.next_reference_name == CHR_FOCUS:
            if end - start < 100: # Arbitrary, just don't want close linked reads
                break
            else:
                linked_reads[read.get_tag('BX')].append((start, end))

samfile.close()

# Construct a CSV with the information gathered
for tag, coordinates in linked_reads.items():
    print("{}, {}".format(tag, coordinates))