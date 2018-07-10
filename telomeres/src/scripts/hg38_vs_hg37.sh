#!/bin/bash

#SBATCH -n 1
#SBATCH -J hg38_telomeres
#SBATCH -p vccc
#SBATCH --mem=8G
#SBATCH -o %J.err
#SBATCH -e %J.out
#SBATCH --time=180:00:00

SAMBAMBA=/data/cephfs/punim0010/extras/sambamba
#conda activate telomeres

### hg38
HG_BUILD="hg38"
BEDFILE="../bed/hg38_igv_manual.bed"

# Truseq
#$SAMBAMBA slice ../bams/$HG_BUILD/COLO829BL-ready.bam -L $BEDFILE -o ../final/$HG_BUILD/COLO829BL-$HG_BUILD-ready.bam # Normal
#$SAMBAMBA slice ../bams/$HG_BUILD/COLO829T-ready.bam -L $BEDFILE -o ../final/$HG_BUILD/COLO829T-$HG_BUILD-ready.bam  # Tumor
# 10X
#$SAMBAMBA slice ../bams/$HG_BUILD/Colo829Bl_10x_EMA.bam -L $BEDFILE -o ../final/$HG_BUILD/COLO829BL-10X-$HG_BUILD-ready.bam # Normal
#$SAMBAMBA slice ../bams/$HG_BUILD/Colo829_100pc_10x.bam -L $BEDFILE -o ../final/$HG_BUILD/COLO829T-10X-$HG_BUILD-ready.bam  # Tumor

# XXX: MAPQ == 0 means multimapper? Wisdom from Arthur Hsu
#samtools view ../final/$HG_BUILD/COLO829BL-$HG_BUILD-ready.bam | awk '$5 == 0' | cut -f 3,5 | sort | uniq -c | sort -rn

#samtools sort -n ../final/$HG_BUILD/COLO829BL-$HG_BUILD-ready.bam > ../final/$HG_BUILD/COLO829BL-$HG_BUILD.sorted.bam
#samtools sort -n ../final/$HG_BUILD/COLO829T-$HG_BUILD-ready.bam > ../final/$HG_BUILD/COLO829T-$HG_BUILD.sorted.bam
#samtools sort -n ../final/$HG_BUILD/COLO829BL-10X-$HG_BUILD-ready.bam > ../final/$HG_BUILD/COLO829BL-10X-$HG_BUILD.sorted.bam
#samtools sort -n ../final/$HG_BUILD/COLO829T-10X-$HG_BUILD-ready.bam > ../final/$HG_BUILD/COLO829T-10X-$HG_BUILD.sorted.bam

#module load parallel
#parallel samtools index ::: ../final/$HG_BUILD/*.bam
ls ../final/$HG_BUILD/*.bam | xargs -n1 -P5 samtools index

#samtools view ../final/$HG_BUILD/COLO829BL-$HG_BUILD.sorted.bam | cut -f 1,3,5  | sort | uniq -c | sort -rn | awk 'BEGIN{print "chromosome,mmap_count"}; {print $2","$1}' > ../final/$HG_BUILD/COLO829BL-$HG_BUILD.hist.csv
#samtools view ../final/$HG_BUILD/COLO829T-$HG_BUILD.sorted.bam | cut -f 1,3,5 | sort | uniq -c | sort -rn | awk 'BEGIN{print "chromosome,mmap_count"}; {print $2","$1}' > ../final/$HG_BUILD/COLO829T-$HG_BUILD.hist.csv
#samtools view ../final/$HG_BUILD/COLO829BL-10X.sorted.bam | cut -f 1,3,5 | sort | uniq -c | sort -rn | awk 'BEGIN{print "chromosome,mmap_count"}; {print $2","$1}' > ../final/$HG_BUILD/COLO829BL-10X-$HG_BUILD.hist.csv
#samtools view ../final/$HG_BUILD/COLO829T-10X.sorted.bam | cut -f 1,3,5 | sort | uniq -c | sort -rn | awk 'BEGIN{print "chromosome,mmap_count"}; {print $2","$1}' > ../final/$HG_BUILD/COLO829T-10X-$HG_BUILD.hist.csv

samtools idxstats ../final/$HG_BUILD/COLO829BL-$HG_BUILD-ready.bam | awk 'BEGIN{print "chromosome,mapped,unmapped"}; {print $1","$3","$4}' | grep -v "*,0" | grep -v "_alt" | grep -v decoy | grep -v HLA | grep -v _KI | grep -v GL | grep -v EBV > ../final/$HG_BUILD/COLO829BL-sliced-idxstats-$HG_BUILD.hist.csv
samtools idxstats ../final/$HG_BUILD/COLO829T-$HG_BUILD-ready.bam | awk 'BEGIN{print "chromosome,mapped,unmapped"}; {print $1","$3","$4}' | grep -v "*,0" | grep -v "_alt" | grep -v decoy | grep -v HLA | grep -v _KI | grep -v GL | grep -v EBV > ../final/$HG_BUILD/COLO829T-sliced-idxstats-$HG_BUILD.hist.csv
samtools idxstats ../final/$HG_BUILD/COLO829BL-10X-$HG_BUILD-ready.bam | awk 'BEGIN{print "chromosome,mapped,unmapped"}; {print $1","$3","$4}' | grep -v "*,0" | grep -v "_alt" | grep -v decoy | grep -v HLA | grep -v _KI | grep -v GL | grep -v EBV > ../final/$HG_BUILD/COLO829BL-10X-sliced-idxstats-$HG_BUILD.hist.csv
samtools idxstats ../final/$HG_BUILD/COLO829T-10X-$HG_BUILD-ready.bam | awk 'BEGIN{print "chromosome,mapped,unmapped"}; {print $1","$3","$4}' | grep -v "*,0" | grep -v "_alt" | grep -v decoy | grep -v HLA | grep -v _KI | grep -v GL | grep -v EBV > ../final/$HG_BUILD/COLO829T-10X-sliced-idxstats-$HG_BUILD.hist.csv


samtools flagstat ../final/$HG_BUILD/COLO829BL-$HG_BUILD-ready.bam > ../final/$HG_BUILD/COLO829BL-sliced-flagstat-$HG_BUILD.hist.txt
samtools flagstat ../final/$HG_BUILD/COLO829T-$HG_BUILD-ready.bam > ../final/$HG_BUILD/COLO829T-sliced-flagstat-$HG_BUILD.hist.txt
samtools flagstat ../final/$HG_BUILD/COLO829BL-10X-$HG_BUILD-ready.bam > ../final/$HG_BUILD/COLO829BL-10X-sliced-flagstat-$HG_BUILD.hist.txt
samtools flagstat ../final/$HG_BUILD/COLO829T-10X-$HG_BUILD-ready.bam > ../final/$HG_BUILD/COLO829T-10X-sliced-flagstat-$HG_BUILD.hist.txt
