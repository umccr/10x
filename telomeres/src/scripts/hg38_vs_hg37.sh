#!/bin/bash

#SBATCH -n 1
#SBATCH -J hg38_telomeres
#SBATCH -p vccc
#SBATCH --mem=8G
#SBATCH -o %J.err
#SBATCH -e %J.out
#SBATCH --time=180:00:00

SAMBAMBA=/data/cephfs/punim0010/extras/sambamba
### hg38
HG_BUILD="hg38"

# Truseq
#$SAMBAMBA slice ../bams/$HG_BUILD/COLO829BL-ready.bam -L ../bed/GRC38h_real_coords.bed -o ../final/$HG_BUILD/COLO829BL-$HG_BUILD-ready.bam # Normal
#$SAMBAMBA slice ../bams/$HG_BUILD/COLO829T-ready.bam -L ../bed/GRC38h_real_coords.bed -o ../final/$HG_BUILD/COLO829T-$HG_BUILD-ready.bam  # Tumor
# 10X
#$SAMBAMBA slice ../bams/$HG_BUILD/Colo829Bl_10x_EMA.bam -L ../bed/GRC38h_real_coords.bed -o ../final/$HG_BUILD/COLO829BL-10X-$HG_BUILD-ready.bam # Normal
#$SAMBAMBA slice ../bams/$HG_BUILD/Colo829_100pc_10x.bam -L ../bed/GRC38h_real_coords.bed -o ../final/$HG_BUILD/COLO829T-10X-$HG_BUILD-ready.bam  # Tumor

samtools view ../final/$HG_BUILD/COLO829BL-$HG_BUILD-ready.bam | cut -f 3 | sort | uniq -c | sort -rn > ../final/$HG_BUILD/COLO829BL-$HG_BUILD.hist
samtools view ../final/$HG_BUILD/COLO829T-$HG_BUILD-ready.bam | cut -f 3 | sort | uniq -c | sort -rn > ../final/$HG_BUILD/COLO829T-$HG_BUILD.hist

samtools view ../final/$HG_BUILD/COLO829BL-10X-$HG_BUILD-ready.bam | cut -f 3 | sort | uniq -c | sort -rn > ../final/$HG_BUILD/COLO829BL-10X-$HG_BUILD.hist
samtools view ../final/$HG_BUILD/COLO829T-10X-$HG_BUILD-ready.bam | cut -f 3 | sort | uniq -c | sort -rn > ../final/$HG_BUILD/COLO829T-10X-$HG_BUILD.hist

### hg37
#HG_BUILD="hg37"
