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
for read in samfile.fetch(CHR_FOCUS, 1, 120000):
    current = read.reference_start    # current read coordinates
    mate = read.next_reference_start # mate read coordinates
    if read.has_tag('BX'):
        if read.next_reference_name == CHR_FOCUS:
#            if mate - current < 100: # Arbitrary, just don't want close linked reads
#                break
#            else:
           linked_reads[read.get_tag('BX')].append((current, mate))

samfile.close()

# Construct a CSV with the information gathered
print("barcode,10X_molecule_length")
for tag, coordinates in linked_reads.items():
    # Make sure that we follow the 10X molecule read chain all the way up to get the total linked read length,
    # thanks @alhsu! :)
    print("{}, {}".format(tag, coordinates[-1][1] - coordinates[0][0]))
