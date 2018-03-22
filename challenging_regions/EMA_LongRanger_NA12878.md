We want to find regions in NA12787 that are covered with GiaB confidence intervals, captured well with 10x aligners but missed with BWA. In order to do that, first quantizing BAMs with mosdepth:

```
export MOSDEPTH_Q0=LOW
export MOSDEPTH_Q1=MEDIAN
export MOSDEPTH_Q2=CALLABLE
mosdepth -t 30 --no-per-base --quantize 0:5:10: mosdepth/bwa_quantized bwa-sort.bam
mosdepth -t 30 --no-per-base --quantize 0:5:10: mosdepth/ema_quantized ema_final.bam
mosdepth -t 30 --no-per-base --quantize 0:5:10: mosdepth/longranger_quantized NA12878_WGS_v2_phased_possorted_grch37.bam
```

Then extracting LOW for BWA and CALLABLE for Longranger and EMA:

```
zgrep LOW mosdepth/bwa_quantized.quantized.bed.gz > regions/bwa_low.bed
zgrep CALLABLE mosdepth/ema_quantized.quantized.bed.gz  > regions/ema_callable.bed
zgrep CALLABLE mosdepth/longranger_quantized.quantized.bed.gz > regions/longranger_callable.bed
```

Intersecting with confidence intervals, and feeding resulting regions into `eval_vcf`.

```
bedsize bwa_low.bed					#  263306035
bedsize bwa_callable.bed			# 2826062821
bedsize ema_callable.bed			# 2822631804		
bedsize longranger_callable.bed		# 2859373367			
bedsize truth_regions.bed			# 2575064465

bedops -i longranger_callable.bed bwa_low.bed truth_regions.bed > lr_clbl__bwa_low__truth.bed    # 0
bedops -i ema_callable.bed bwa_low.bed truth_regions.bed > ema_clbl__bwa_low__truth.bed          # 943

# Slop the regions a bit:
bedtools slop -b 500 -i ema_clbl__bwa_low__truth.bed -g /g/data/gx8/local/development/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa.fai > ema_clbl__bwa_low__truth.p500.bed

~/bin/sambamba slice ../bwa-sort.bam -L ema_clbl__bwa_low__truth.p500.bed > bwa__ema_clbl__bwa_low__truth.bam
~/bin/sambamba slice ../ema_final.bam -L ema_clbl__bwa_low__truth.p500.bed > ema__ema_clbl__bwa_low__truth.bam
~/bin/sambamba slice ../NA12878_WGS_v2_phased_possorted_grch37.bam -L ema_clbl__bwa_low__truth.p500.bed > lr__ema_clbl__bwa_low__truth.bam
parallel "samtools sort {} > {.}.sort.bam && mv {.}.sort.bam {}" ::: *.bam

# LOCAL
scp rjn:/g/data3/gx8/projects/Saveliev_10X/NA12878-10x/challenging/regions/bwa__ema_clbl__bwa_low__truth.bam .
scp rjn:/g/data3/gx8/projects/Saveliev_10X/NA12878-10x/challenging/regions/ema__ema_clbl__bwa_low__truth.bam .
scp rjn:/g/data3/gx8/projects/Saveliev_10X/NA12878-10x/challenging/regions/lr__ema_clbl__bwa_low__truth.bam .
scp rjn:/g/data3/gx8/projects/Saveliev_10X/NA12878-10x/challenging/regions/ema_clbl__bwa_low__truth.bed .
parallel "samtools index" ::: *.bam
```


```
REG="7:5989546-6069747"
mkdir $REG
~/bin/sambamba slice ../bwa-sort.bam  $REG                              | samtools sort -O bam -o "$REG/bwa.bam"
~/bin/sambamba slice ../ema_final.bam $REG                              | samtools sort -O bam -o "$REG/ema.bam"
~/bin/sambamba slice ../NA12878_WGS_v2_phased_possorted_grch37.bam $REG | samtools sort -O bam -o "$REG/lr.bam"

# LOCAL
REG="7:5989546-6069747"
mkdir $REG
scp rjn:/g/data3/gx8/projects/Saveliev_10X/NA12878-10x/challenging/regions/$REG/bwa.bam "./$REG"
scp rjn:/g/data3/gx8/projects/Saveliev_10X/NA12878-10x/challenging/regions/$REG/ema.bam "./$REG"
scp rjn:/g/data3/gx8/projects/Saveliev_10X/NA12878-10x/challenging/regions/$REG/lr.bam "./$REG"
parallel "samtools index" ::: ./$REG/*.bam

```

