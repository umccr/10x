#!/usr/bin/env bash

set -x
set -e

PREFIX=$1
VCF=$2
BAM=$3

# convert BAM file to the compact fragment file format containing only haplotype-relevant information. This is a necessary precursor step to running HapCUT2
extractHAIRS --10X 1 --bam ${BAM} --VCF ${VCF} --out ${PREFIX}unlinked_fragment_file

# link fragments into barcoded molecules
LinkFragments.py --bam ${BAM} --VCF ${VCF} --fragments ${PREFIX}unlinked_fragment_file --out ${PREFIX}linked_fragment_file

# use HAPCUT2 to assemble fragment file into haplotype blocks
HAPCUT2 --nf 1 --fragments ${PREFIX}linked_fragment_file --VCF ${VCF} --output ${PREFIX}haplotype_output_file

# generate phased VCF with populated GT tags
fgbio HapCutToVcf -v ${VCF} -i ${PREFIX}haplotype_output_file -o ${VCF}.phased.vcf

set +x
set +e