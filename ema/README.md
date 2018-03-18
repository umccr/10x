## Running EMA pipeline on NA12878 WGS 10x data

Based on the notebook https://github.com/arshajii/ema-paper-data/blob/master/experiments.ipynb

### Environment and data

On Spartan:

```
cd /data/cephfs/punim0010/extras/vlad/synced/umccr/10x/ema
```

Load conda environment (conda install -c bioconda seqtk parallel bwa picard -y):
```
# spartan:  export PATH=/data/cephfs/punim0010/extras/10x/miniconda/bin:$PATH ; source activate 10x
# raijin:   export PATH=/g/data3/gx8/extras/10x/miniconda/bin:$PATH
```

Install EMA

```
git clone --recursive https://github.com/arshajii/ema
cd ema
make
cp ema ~/bin/ema
```

Activate an interactive session

```
sinteractive --time=80:00:00 --nodes=1 --ntasks=2 -p vccc --mem=30G -J ema
```


### Prep subsampled data

Subset FastQ

```
cd 10k
seqtk sample -s100 ../ori_ungz_fq/NA12878_WGS_v2_S1_L001_R1_001.fastq 10000 > subset_L001_R1_001.fastq &
seqtk sample -s100 ../ori_ungz_fq/NA12878_WGS_v2_S1_L001_R2_001.fastq 10000 > subset_L001_R2_001.fastq &
```

Renaming to make `ema_wrapper.sh` happy

```
mv subset_L001_R1_001.fastq subset_L001_R1.fastq
mv subset_L001_R2_001.fastq subset_L001_R2.fastq
```

EMA preprocessing

```
cat subset_L001_R*.fastq | ema count -1 - -w ../4M-with-alts-february-2016.txt -o counts_file
ema preproc -1 subset_L001_R1.fastq -2 subset_L001_R2.fastq -w ../4M-with-alts-february-2016.txt -c counts_file -n 2
```

EMA sort buckets

```
parallel -j 2 --xapply "ema sort -1 {1} -2 {2}" \
	::: bucket0*/*1.preproc.fastq \
	::: bucket0*/*2.preproc.fastq
```

Run alignment (spartan)

```
EMAPATH=/home/vlad/bin/ema \
PICARDPATH=/data/cephfs/punim0010/extras/10x/miniconda/envs/10x/share/picard-2.17.11-0/picard.jar \
/data/cephfs/punim0010/extras/vlad/synced/umccr/10x/ema/ema_wrapper.sh \
-r /data/cephfs/punim0010/extras/10x/refdata-GRCh37/fasta/genome.fa \
-R '@RG\tID:NA12878_WGS_v2:LibraryNotSpecified:1:unknown_fc:0\tSM:NA12878_WGS_v2' \
-t 2
```

Run alignment (raijin)

```
EMAPATH=~/bin/ema \
PICARDPATH=/g/data3/gx8/extras/10x/miniconda/share/picard-2.17.11-0/picard.jar \
/g/data3/gx8/extras/10x/ema/util/ema_wrapper.sh \
-r ref/GRCh37.fa \
-R '@RG\tID:NA12878_WGS_v2:LibraryNotSpecified:1:unknown_fc:0\tSM:NA12878_WGS_v2' \
-t 2
```

### Full data (spartan)

```
mkdir full
cd full
ln -s ../ori_ungz_fq/NA12878_WGS_v2_S1_L001_R1_001.fastq NA12878_WGS_v2_S1_L001_R1.fastq
ln -s ../ori_ungz_fq/NA12878_WGS_v2_S1_L001_R2_001.fastq NA12878_WGS_v2_S1_L001_R2.fastq
cat *.fastq | ema count -1 - -w ../4M-with-alts-february-2016.txt -o counts_file

THREADS=29
ema preproc -1 NA12878_WGS_v2_S1_L001_R1.fastq -2 NA12878_WGS_v2_S1_L001_R2.fastq -w ../4M-with-alts-february-2016.txt -c counts_file -n $THREADS
parallel -j $THREADS --xapply "echo ema sort -1 {1} -2 {2} ; ema sort -1 {1} -2 {2}" \
	::: bucket0*/*1.preproc.fastq \
	::: bucket0*/*2.preproc.fastq

EMAPATH=/home/vlad/bin/ema \
PICARDPATH=/data/cephfs/punim0010/extras/10x/miniconda/envs/10x/share/picard-2.17.11-0/picard.jar \
/data/cephfs/punim0010/extras/vlad/synced/umccr/10x/ema/ema_wrapper.sh \
-r /data/cephfs/punim0010/extras/10x/refdata-GRCh37/fasta/genome.fa \
-R '@RG\tID:NA12878_WGS_v2:LibraryNotSpecified:1:unknown_fc:0\tSM:NA12878_WGS_v2' \
-t $THREADS
```

### Full data (raijin)

```
cd full

#!/bin/bash
#PBS -P gx8
#PBS -q normalsp
#PBS -l walltime=96:00:00
#PBS -l mem=128GB
#PBS -l ncpus=31
#PBS -l software=ema
#PBS -l wd

export PATH=/g/data3/gx8/extras/10x/miniconda/bin:$PATH

barcodes=/g/data3/gx8/projects/Saveliev_10X/NA12878-10x/ema/full/4M-with-alts-february-2016.txt

cat *.fastq | ema count -1 - -w $barcodes -o counts_file
THREADS=30
ema preproc -1 NA12878_WGS_v2_S1_L001_R1.fastq -2 NA12878_WGS_v2_S1_L001_R2.fastq -w $barcodes -c counts_file -n $THREADS
parallel -j $THREADS --xapply "echo ema sort -1 {1} -2 {2} ; ema sort -1 {1} -2 {2}" \
	::: bucket0*/*1.preproc.fastq \
	::: bucket0*/*2.preproc.fastq

EMAPATH=~/bin/ema \
PICARDPATH=/g/data3/gx8/extras/10x/miniconda/share/picard-2.17.11-0/picard.jar \
/g/data3/gx8/extras/10x/ema/util/ema_wrapper.sh \
-r ref/GRCh37.fa \
-R '@RG\tID:NA12878_WGS_v2:LibraryNotSpecified:1:unknown_fc:0\tSM:NA12878_WGS_v2' \
-t $THREADS
```





## New version

Install on Spartan:
- git pull
- module load `libzip/1.1.2-GCC-6.2.0` (Spartan) or `gcc/6.2.0` (Raijin)
- make

To use on Raijin, load libzip/gcc6 first:

```
module load gcc/6.2.0
export PATH=/g/data3/gx8/extras/10x/miniconda/bin:$PATH
```

Don't use on Spartan - massive IO issues. But if you need, run:

```
module load libzip/1.1.2-GCC-6.2.0
export PATH=/data/cephfs/punim0010/extras/10x/miniconda/bin:$PATH
```

Enter interactive session:

```
qsub -I -P gx8 -q normalsp -l walltime=96:00:00 -l ncpus=32 -l wd -l mem=128G
```

#### Preprocess

```
cd /data/cephfs/punim0010/extras/vlad/synced/umccr/10x/ema/full

# Make interleaved
paste <(cat NA12878_WGS_v2_S1_L001_R1_001.fastq | paste - - - -) <(cat NA12878_WGS_v2_S1_L001_R2_001.fastq | paste - - - -) | tr "\t" "\n" > NA12878_WGS_v2_S1_L001.fastq

# Count
ema count -w ../4M-with-alts-february-2016.txt -o NA12878_WGS_v2_S1_L001 < NA12878_WGS_v2_S1_L001.fastq 2>NA12878_WGS_v2_S1_L001.log

# Preproc
THREADS=29  # Spartan
THREADS=30  # Raijin

ema preproc -w ../4M-with-alts-february-2016.txt -n 500 -t $THREADS -o output_dir *.ema-ncnt < NA12878_WGS_v2_S1_L001.fastq 2>&1 | tee preproc.log
```

#### Align

```
THREADS=30  # Raijin
ema align -t $THREADS -d -r  -s  | samtools sort -@ 4 -O bam -l 0 -m 4G -o NA12878_WGS_v2.bam -

parallel --bar -j8 "ema align -R $'@RG\tID:NA12878_10x_EMA\tSM:NA12878_10x_EMA' -t 4 -d -r /g/data3/gx8/projects/Saveliev_10X/NA12878-10x/ema/ref/GRCh37.fa -s {} | samtools sort -@ 4 -O bam -l 0 -m 4G -o {}.bam -" ::: output_dir/ema-bin-???

bwa mem -p -t 32 -M -R "@RG\tID:NA12878_10x_EMA\tSM:NA12878_10x_EMA" /g/data3/gx8/projects/Saveliev_10X/NA12878-10x/ema/ref/GRCh37.fa output_dir/ema-bin-nobc |\
  samtools sort -@ 4 -O bam -l 0 -m 4G -o output_dir/ema-bin-nobc.bam
```














