#!/usr/bin/env python
"""
Use BX tag to match linked reads: https://support.10xgenomics.com/genome-exome/software/pipelines/latest/output/bam
Collects all the possible spans and makes a histogram
"""
import sys
import pysam
from collections import defaultdict

CHR_FOCUS = 'chr5'
linked_reads = defaultdict(defaultdict(list).copy)
samfile = pysam.AlignmentFile("telomeres/data/interim/{}_hg38_elongated.bam".format(CHR_FOCUS), "rb")

# XXX: Filter by a reasonable MAPQ value?
for read in samfile.fetch(CHR_FOCUS, 1, 120000):
    current = read.reference_start    # current read coordinates
    #mate = read.next_reference_start # mate read coordinates
    mate_chrom = read.next_reference_name
    if read.has_tag('BX'):
        if mate_chrom == CHR_FOCUS:
           linked_reads[read.get_tag('BX')]['reads'].append(current)
        elif mate_chrom is not None: # XXX: What does it mean None in this context? No mate? Orphans?
           linked_reads[read.get_tag('BX')]['mates'].append(mate_chrom)
        #else:
        #   linked_reads[read.get_tag('BX')]['orphans'].append(current)

samfile.close()

# XXX: pandas.from_dict instead of loop
# Construct a CSV with the information gathered
print("barcode,10X_molecule_length,reads_per_BX_barcode,total_mates_in_other_chroms")
for tag, lr in linked_reads.items():
    if len(lr['reads']) > 0:
        print("{}, {}, {}, {}".format(tag, max(lr['reads']) - min(lr['reads']), len(lr['reads']), len(lr['mates'])))