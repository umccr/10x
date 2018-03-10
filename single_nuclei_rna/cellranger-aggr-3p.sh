#!/bin/bash
#SBATCH -p vccc
#SBATCH --time=2-00:00
##SBATCH --mem=256G
#SBATCH -o out.10x-cellranger-aggr-%j
#SBATCH -e err.10x-cellranger-aggr-%j
#SBATCH --nodes=1
#SBATCH --mem=384G

WORK_DIR=/data/cephfs/punim0010/projects/10X_scRNA-20180226
CSV=${WORK_DIR}/aggr_samples.3p.csv


cellranger aggr \
  --id=final_20180226_3prime \
  --csv=${CSV} \
  --normalize=mapped \
  --jobmode=local \
  --localcores=32 \
  --localmem=320


