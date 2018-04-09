## Running EMA pipeline on NA12878 WGS 10x data

Based on the notebook https://github.com/arshajii/ema-paper-data/blob/master/experiments.ipynb

### Install EMA

```
module load gcc/6.2.0  # Raijin
git clone --recursive https://github.com/arshajii/ema
cd ema
make
cp ema ~/bin/ema
```

### Environment and data

On Spartan:

```
cd /data/cephfs/punim0010/extras/vlad/synced/umccr/10x/ema
```

Load conda environment (conda install -c bioconda seqtk parallel bwa picard -y):
```
# spartan:  export PATH=/data/cephfs/punim0010/extras/10x/miniconda/bin:$PATH
# raijin:   export PATH=/g/data3/gx8/extras/10x/miniconda/bin:$PATH
```

Activate an interactive session

```
sinteractive --time=80:00:00 --nodes=1 --cpus-per-task=32 -p vccc --mem=256G -J ema
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

### Preprocess

```
cd /data/cephfs/punim0010/extras/vlad/synced/umccr/10x/ema/full

# Make interleaved
paste <(cat NA12878_WGS_v2_S1_L001_R1_001.fastq | paste - - - -) <(cat NA12878_WGS_v2_S1_L001_R2_001.fastq | paste - - - -) | tr "\t" "\n" > NA12878_WGS_v2_S1_L001.fastq

# Count
ema count -w ../4M-with-alts-february-2016.txt -o NA12878_WGS_v2_S1_L001 < NA12878_WGS_v2_S1_L001.fastq 2>NA12878_WGS_v2_S1_L001.log

# Preproc
THREADS=29  # Spartan
THREADS=32  # Raijin

ema preproc -w ../4M-with-alts-february-2016.txt -n 500 -t $THREADS -o output_dir *.ema-ncnt < NA12878_WGS_v2_S1_L001.fastq 2>&1 | tee preproc.log
```

#### Align

```
THREADS=30  # Raijin

parallel -j8 "ema align -R $'@RG\tID:NA12878_10x_EMA\tSM:NA12878_10x_EMA' -t 4 -d -r /g/data3/gx8/projects/Saveliev_10X/NA12878-10x/ema/ref/GRCh37.fa -s {} | samtools sort -@ 4 -O bam -l 0 -m 4G -o {}.bam -" ::: output_dir_2/ema-bin-???

bwa mem -p -t 32 -M -R "@RG\tID:NA12878_10x_EMA\tSM:NA12878_10x_EMA" /g/data3/gx8/projects/Saveliev_10X/NA12878-10x/ema/ref/GRCh37.fa output_dir/ema-bin-nobc |\
  samtools sort -@ 4 -O bam -l 0 -m 4G -o output_dir/ema-bin-nobc.bam

sambamba markdup -t 32 -p -l 0 output_dir/ema-bin-nobc.bam output_dir/ema-bin-nobc-dupsmarked.bam

sambamba merge -t 32 -p ema_final.bam output_dir/*.bam

samtools stats -@ 32 ema_final.bam > ema_final.stats.txt
```


#### Full pipeline for COLO829

```
module load gcc/6.2.0
export PATH=/g/data3/gx8/extras/10x/miniconda/bin:/home/563/vs2870/bin:$PATH

cd /g/data3/gx8/projects/Saveliev_10X/COLO829-10x/ema
ln -s /g/data3/gx8/projects/Hsu_10X_WGS/FASTQ-UMCCR/COLO829-80pc/Colo829_80pc_S1 ori_fq

for fp in ori_fq/*_R1_001.fastq.gz
do
	paste <(pigz -c -d $fp | paste - - - -) <(pigz -c -d ${fp/_R1_/_R2_} | paste - - - -) | awk -F '\t' 'length($2) >= 40' | tr "\t" "\n" >> Colo829_80pc.fastq
done

cat Colo829_80pc.fastq | paste - - - - | awk 'length($2) >= 30 && length($2)' | sed 's/\t/\n/g' | ema count -w 4M-with-alts-february-2016.txt -o Colo829_80pc 2> ema_count.log

ema count -w 4M-with-alts-february-2016.txt -o Colo829_80pc < Colo829_80pc.fastq 2> ema_count.log

ema preproc -w 4M-with-alts-february-2016.txt -n 500 -t 32 -o ema_work *.ema-ncnt < Colo829_80pc.fastq 2>&1 | tee ema_preproc.log

parallel -j8 "ema align -R $'@RG\tID:Colo829_80pc_EMA\tSM:Colo829_80pc_EMA' -t 4 -d -r ref/GRCh37.fa -s {} | samtools sort -@ 4 -O bam -l 0 -m 4G -o {}.bam -" ::: ema_work/ema-bin-???

bwa mem -p -t 32 -M -R "@RG\tID:Colo829_80pc_EMA\tSM:Colo829_80pc_EMA" ref/GRCh37.fa ema_work/ema-bin-nobc |\
  samtools sort -@ 4 -O bam -l 0 -m 4G -o ema_work/ema-bin-nobc.bam

sambamba markdup -t 32 -p -l 0 ema_work/ema-bin-nobc.bam ema_work/ema-bin-nobc-dupsmarked.bam && rm ema_work/ema-bin-nobc.bam

sambamba merge -t 32 -p Colo829_80pc_EMA.bam ema_work/*.bam

samtools stats -@ 32 Colo829_80pc_EMA.bam > Colo829_80pc_EMA.stats.txt
```


#### Subsampling:

```
for fq in ../ori_fq/* ; do
	gunzip -c $fq | head -n700 | gzip -c > ori_fq/$(basename $fq)
done
```


### Not merged version (LATEST)

Create a file `interleave_fq.sh`:

```
#paste <(pigz -c -d $1 | paste - - - - | awk -F '\t' 'length($2) >= 40') <(pigz -c -d ${1/_R1_/_R2_} | paste - - - - | awk -F '\t' 'length($2) >= 40') | tr '\t' '\n'
paste <(pigz -c -d $1 | paste - - - -) <(pigz -c -d ${1/_R1_/_R2_} | paste - - - -) | tr '\t' '\n'
```

```
echo 'paste <(pigz -c -d $1 | paste - - - -) <(pigz -c -d ${1/_R1_/_R2_} | paste - - - -) | tr "\\t" "\\n"' > interleave_fq.sh
```

```
module load gcc/6.2.0
export PATH=/g/data3/gx8/extras/10x/miniconda/bin:/home/563/vs2870/bin:$PATH

SAMPLE=Colo829Bl_10x_80pc

date
parallel -j32 "bash interleave_fq.sh {} | ema count -w 4M-with-alts-february-2016.txt -o {/.} 2>{/.}.log" ::: ori_fq/*_R1_*.gz

date
ls ori_fq/*R1*.gz | xargs -I '{}' bash interleave_fq.sh '{}' | ema preproc -w 4M-with-alts-february-2016.txt -n 500 -t 32 -o ema_work *.ema-ncnt 2>&1 | tee ema_preproc.log

date
parallel -j8 "ema align -R $'@RG\tID:${SAMPLE}_EMA\tSM:${SAMPLE}_EMA' -t 4 -d -r ref/GRCh37.fa -s {} | samtools sort -@ 4 -O bam -l 0 -m 4G -o {}.bam -" ::: ema_work/ema-bin-???

date
bwa mem -p -t 32 -M -R "@RG\tID:${SAMPLE}_EMA\tSM:${SAMPLE}_EMA" ref/GRCh37.fa ema_work/ema-bin-nobc |\
  samtools sort -@ 4 -O bam -l 0 -m 4G -o ema_work/ema-bin-nobc.bam

date
sambamba markdup -t 32 -p -l 0 ema_work/ema-bin-nobc.bam ema_work/ema-bin-nobc-dupsmarked.bam && rm ema_work/ema-bin-nobc.bam

date
sambamba merge -t 32 -p ${SAMPLE}_EMA.bam ema_work/*.bam

date
samtools stats -@ 32 ${SAMPLE}_EMA.bam > ${SAMPLE}_EMA.stats.txt

date
```

<!-- for fq in ../ori_fq/* ; do
	gunzip -c $fq | head -n400000 | gzip -c > ori_fq/$(basename $fq)
done
parallel -j32 "bash interleave_fq.sh {} | ema count -w 4M-with-alts-february-2016.txt -o {/.}" ::: ori_fq/*_R1_*.gz
 -->

### Extract BAMs to explore challenging regions

```
for s in bwa-sort ema_final longranger_decoy_pos_sorted longranger_pos_sorted minimap2-sort
do
	sambamba slice $s.bam 22:42522501-42526883 > $s.CYP2D6.bam
	sambamba index $s.CYP2D6.bam
done

for s in bwa-sort ema_final longranger_decoy_pos_sorted longranger_pos_sorted minimap2-sort
do
	scp -r spa:/data/cephfs/punim0010/projects/Saveliev_10X/NA12878-10x/bcbio_grch37/bams/$s.CYP2D6.bam .
	scp -r spa:/data/cephfs/punim0010/projects/Saveliev_10X/NA12878-10x/bcbio_grch37/bams/$s.CYP2D6.bam.bai .
done


regions.bed:
1	104090000	104310000	AMY1/2
6	31940000	32010000	C4A/B
22	42521302	42549900	CY2D6/7

for s in bwa-sort ema_final longranger_decoy_pos_sorted longranger_pos_sorted minimap2-sort
do
	~/bin/sambamba slice $s.bam -L regions.bed > $s.regions.bam
	sambamba index $s.regions.bam
done

for s in bwa-sort ema_final longranger_decoy_pos_sorted longranger_pos_sorted minimap2-sort
do
	scp -r spa:/data/cephfs/punim0010/projects/Saveliev_10X/NA12878-10x/bcbio_grch37/bams/$s.regions.bam .
	scp -r spa:/data/cephfs/punim0010/projects/Saveliev_10X/NA12878-10x/bcbio_grch37/bams/$s.regions.bam.bai .
done
```

### Run bcbio for QC and variant calling on Raijin:

`cd /g/data3/gx8/projects/Saveliev_10X/COLO829-10x/bcbio_10x/work` ad use this [yaml](bcbio_colo829_raijin_10x.yaml)

Separately produce BWA alignments `cd /g/data3/gx8/projects/Saveliev_10X/COLO829-10x/bcbio_bwa` with [yaml](bcbio_colo829_raijin_bwa.yaml)






















