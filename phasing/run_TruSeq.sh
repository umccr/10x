#!/usr/bin/env bash

set -x
set -e

PREFIX=$1
VCF=$2
BAM=$3

# convert BAM file to the compact fragment file format containing only haplotype-relevant information. This is a necessary precursor step to running HapCUT2
extractHAIRS  --bam ${BAM} --VCF ${VCF} --out ${PREFIX}fragment_file

# use HAPCUT2 to assemble fragment file into haplotype blocks
HAPCUT2  --fragments ${PREFIX}fragment_file --VCF ${VCF} --output ${PREFIX}haplotype_output_file

# generate phased VCF with populated GT tags
fgbio HapCutToVcf -v ${VCF} -i ${PREFIX}haplotype_output_file -o ${PREFIX}phased.vcf

set +x
set +e