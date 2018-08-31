#!/usr/bin/env python3
import sys
import os
from numpy import mean
from os.path import basename, dirname, join
import pysam

bam = sys.argv[1]
bamf_hqual_out = sys.argv[2] if len(sys.argv) > 2 else join(dirname(bam), 'not_unmapped_or_mate', basename(bam))

bamf = pysam.AlignmentFile(bam, "rb")
bamf_hqual_f = pysam.AlignmentFile(bamf_hqual_out, "wb", template=bamf)

for read in bamf:
    if read.is_paired and (not read.is_unmapped or not read.mate_is_unmapped):
        bamf_hqual_f.write(read)

bamf.close()
bamf_hqual_f.close()

print('Done. Stats:')

# Using this command instead:
# sambamba view -f bam -F "unmapped or mate_is_unmapped" -t 20 diploid_tumor-ready.bam | samtools sort -n -@ 10 > diploid_tumor-unmapped_or_mate.namesorted.bam
