## Running EMA pipeline on NA12878 WGS 10x data

Based on the notebook https://github.com/arshajii/ema-paper-data/blob/master/experiments.ipynb

### Environment and data

On Spartan:

```
cd /data/cephfs/punim0010/extras/vlad/synced/umccr/10x/ema
```

Install conda environment

```
conda install -c bioconda seqtk -y
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
sinteractive --time=80:00:00 --nodes=30 --ntasks=1 -p vccc --mem=30G -J ema
```

Subset FastQ

```
seqtk sample -s100 ori_ungz_fq/NA12878_WGS_v2_S1_L001_R1_001.fastq 10000 > subset10k_L001_R1_001.fastq &
seqtk sample -s100 ori_ungz_fq/NA12878_WGS_v2_S1_L001_R2_001.fastq 10000 > subset10k_L001_R2_001.fastq &
```

Preprocessing

```
cat 40k_L001_R*_001.fastq | ema count -1 - -w barcode_whitelist.txt -o counts_file
cat 40k_L001_R*_001.fastq | ema preproc -1 - -w barcode_whitelist.txt -c counts_file -n 2
```

Sort buckets

```
ema sort -1 bucket000/*_1.fastq -2 bucket000/*_2.fastq
ema sort -1 bucket001/*_1.fastq -2 bucket001/*_2.fastq
```

Run alignment

```
EMAPATH=/home/vlad/bin/ema \
PICARDPATH=/data/cephfs/punim0010/extras/10x/miniconda/envs/10x/share/picard-2.17.11-0/picard.jar \
./ema_wrapper.sh \
-r /data/projects/punim0010/local/share/bcbio/genomes/Hsapiens/hg19/seq/hg19.fa \
-R $'@RG\tID:foo\tSM:bar' \
-t 2
```
