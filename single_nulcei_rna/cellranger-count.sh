#!/bin/bash
#SBATCH -p vccc
#SBATCH --time=2-00:00
#SBATCH --array=1-12
#SBATCH --cpus-per-task=32
#SBATCH --mem=384G
#SBATCH -o out.10x-cellranger-%A-%a
#SBATCH -e err.10x-cellranger-%A-%a


# Array job is limited to 12 because samples
# 13_nc_5prime_CUP_1184_B-cell
# 14_nc_5prime_CUP_1184_T-cell
# are VDJ libraries that requires 2x150 reads whereas 26+98 were sequenced.

WORK_DIR=/data/cephfs/punim0010/projects/10X_scRNA-20180226
#FASTQ_DIR=/data/cephfs/punim0010/data/FASTQ/180226_A00130_0038_BH3JKFDMXX/Chromium_20180226
#REF=/data/projects/punim0010/opt/refdata-cellranger-hg19-1.2.0
FASTQ_DIR=${WORK_DIR}/FASTQ
REF=/data/projects/punim0010/opt/refdata-cellranger-hg19-1.2.0-premrna
SAMPLES_FILE=${WORK_DIR}/samples.txt

SAMPLE=`head -n $SLURM_ARRAY_TASK_ID $SAMPLES_FILE | tail -n 1`

if [[ $SAMPLE == *5prime* ]]; then
    CHEMISTRY=fiveprime
elif [[ $SAMPLE == *3prime* ]]; then
    CHEMISTRY=SC3Pv2
fi

echo "[`date`] $SAMPLE"

cellranger count \
  --chemistry=${CHEMISTRY} \
  --id=20180226_${SAMPLE} \
  --fastqs=${FASTQ_DIR} \
  --sample=${SAMPLE} \
  --transcriptome=${REF} \
  --jobmode=local \
  --localcores=32 \
  --localmem=320

echo "[`date`] DONE - $SAMPLE"
