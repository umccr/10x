#!/usr/bin/env python
"""
Use BX tag to match linked reads: https://support.10xgenomics.com/genome-exome/software/pipelines/latest/output/bam
Collects all the possible spans and makes a histogram
"""
import sys
import pysam
from collections import defaultdict

CHR_FOCUS = 'chr5'
linked_reads = defaultdict(list)
samfile = pysam.AlignmentFile("telomeres/data/interim/{}_hg38_elongated.bam".format(CHR_FOCUS), "rb")

# XXX: Filter by a reasonable MAPQ value?
for read in samfile.fetch(CHR_FOCUS, 1, 120000):
    current = read.reference_start    # current read coordinates
    #mate = read.next_reference_start # mate read coordinates
    mate_chrom = read.next_reference_name
    if read.has_tag('BX'):
        if read.next_reference_name == CHR_FOCUS:
           linked_reads[read.get_tag('BX')].append(current)
        else:
           linked_reads_other_chroms[read.get_tag]
           # XXX: Rethink data structure for other chroms 

samfile.close()

# Construct a CSV with the information gathered
print("barcode,10X_molecule_length,reads_per_BX_barcode,mate_is_in_chrom")
for tag, coordinates in linked_reads.items():
    print("{}, {}, {}".format(tag, max(coordinates) - min(coordinates), len(coordinates)))
