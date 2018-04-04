#!/usr/bin/env python3
import sys
import os
from numpy import mean
from os.path import basename, dirname, join
import pysam

bam = sys.argv[1]
bamf_bx_out = join(dirname(bam), 'with_bx', basename(bam))
bamf_hqual_out = join(dirname(bam), 'lng_hqual', basename(bam))
bamf_bx_hqual_out = join(dirname(bam), 'with_bx_lng_hqual', basename(bam))

# print('Loading whitelisted barcodes')
# with open('/g/data3/gx8/projects/Saveliev_10X/NA12878-10x/insertions/unmapped_or_mate_is_unmapped/4M-with-alts-february-2016.txt') as f:
#     barcodes = set(l.strip() for l in f.readlines())
# print(f'Read {len(barcodes)} barcodes: ')
# print('...')
# print()

total_i = 0
bx_i = 0
lng_i = 0
hqual_i = 0
lng_hqual_i = 0
bx_lng_hqual_i = 0
paired_i = 0
paired_lng_hqual_i = 0
paired_bx_lng_hqual_i = 0
reads_without_quality = 0

bamf = pysam.AlignmentFile(bam, "rb")
bamf_bx_f = pysam.AlignmentFile(bamf_bx_out, "wb", template=bamf)
bamf_hqual_f = pysam.AlignmentFile(bamf_hqual_out, "wb", template=bamf)
bamf_bx_hqual_f = pysam.AlignmentFile(bamf_bx_hqual_out, "wb", template=bamf)

to_skip = set()

for read in bamf:
    total_i += 1
    bx = read.has_tag('BX')
    lng = read.query_length >= 125
    if not read.query_qualities:
        reads_without_quality += 1
        continue

    hqual = all(q >= 10 for q in read.query_qualities) and mean(read.query_qualities) >= 25

    if bx:
        bx_i += 1
        bamf_bx_f.write(read)
    if lng:
        lng_i += 1
    if hqual:
        hqual_i += 1
    if lng and hqual:
        lng_hqual_i += 1
        if bx:
            bx_lng_hqual_i += 1
    if read.is_paired:
        paired_i += 1
        if lng and hqual:
            paired_lng_hqual_i += 1
            bamf_hqual_f.write(read)
            if bx:
                paired_bx_lng_hqual_i += 1
                bamf_bx_hqual_f.write(read)
    # bx = read.get_tag('BX')[:-2]
    # if bx in barcodes:
        # with_correct_bx_i += 1

    # TODO: try to reconstruct the original BX?

bamf.close()
bamf_bx_f.close()
bamf_hqual_f.close()
bamf_bx_hqual_f.close()

# with open(filtered_fq1, 'w') as fq1_o, open(filtered_fq2, 'w') as fq2_o:
#     for rec1, rec2 in zip(fq1_i, fq2_i):
#         i += 1
#         if not _failed(rec1) and not _failed(rec2):
#             SeqIO.write(rec1, fq1_o, 'fastq')
#             SeqIO.write(rec2, fq2_o, 'fastq')
#             written_i += 1
#         if i % 100000 == 0:
#             print(f'Processed {i}, written {written_i}')

print('Done. Stats:')
print('')
print(f'reads_without_quality: {reads_without_quality}')
print(f'Total: {total_i}')
print(f'bx_i: {bx_i}')
print(f'lng_i: {lng_i}')
print(f'hqual_i: {hqual_i}')
print(f'lng_hqual_i: {lng_hqual_i}')
print(f'bx_lng_hqual_i: {bx_lng_hqual_i}')
print(f'paired_i: {paired_i}')
print(f'paired_lng_hqual_i: {paired_lng_hqual_i}')
print(f'paired_bx_lng_hqual_i: {paired_bx_lng_hqual_i}')

