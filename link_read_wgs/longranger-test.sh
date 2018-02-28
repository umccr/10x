#!/bin/bash
#SBATCH -p vccc
#SBATCH --time=2-00:00
#SBATCH -o out.longranger-test-%j
#SBATCH -e err.longranger-test-%j
#SBATCH --mem=64G

WORK_DIR=/data/cephfs/punim0010/projects/10X_WGS-test

/data/projects/punim0010/opt/longranger-2.1.3/longranger wgs \
  --id=NA12878-10X-WFU \
  --fastqs=${WORK_DIR}/FASTQ/HNKKVCCXX \
  --reference=/data/projects/punim0010/opt/refdata-longranger-hg19-2.1.0 \
  --jobmode=/data/projects/punim0010/opt/cellranger-2.1.0/slurm.template \
  --mempercore=12 \
  --maxjobs=64 \
  --jobinterval=200

