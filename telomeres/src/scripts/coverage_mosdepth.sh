#!/bin/bash

HG_BUILD="hg38"

#SBATCH -n 1
#SBATCH -J mosdepth_coverage 
#SBATCH -p vccc 
#SBATCH --mem=16G
#SBATCH -o %J.err
#SBATCH -e %J.out
#SBATCH --time=10:00:00

MOSDEPTH=mosdepth

# Truseq
$MOSDEPTH --by ../bed/GRC38h_real_coords.bed ../final/$HG_BUILD/mosdepth/COLO829BL-$HG_BUILD ../final/$HG_BUILD/COLO829BL-$HG_BUILD-ready.bam
$MOSDEPTH --by ../bed/GRC38h_real_coords.bed ../final/$HG_BUILD/mosdepth/COLO829T-$HG_BUILD ../final/$HG_BUILD/COLO829T-$HG_BUILD-ready.bam
# 10X
$MOSDEPTH --by ../bed/GRC38h_real_coords.bed ../final/$HG_BUILD/mosdepth/COLO829BL-10X-$HG_BUILD ../final/$HG_BUILD/COLO829BL-10X-$HG_BUILD-ready.bam 
$MOSDEPTH --by ../bed/GRC38h_real_coords.bed ../final/$HG_BUILD/mosdepth/COLO829T-10X-$HG_BUILD ../final/$HG_BUILD/COLO829T-10X-$HG_BUILD-ready.bam
