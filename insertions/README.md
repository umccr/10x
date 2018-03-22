See [notebook](10x_novel_insertions.ipynb)

Misc commands:

```
bedtools bamtofastq -i COLO829-10x.sorted.bam -fq COLO829-10x.R1.fq -fq2 COLO829-10x.R2.fq 2>COLO829-10x.unmpaired

bxtools split NA12878_WGS.sorted.bam -a NA12878_WGS -m 100 | sort -k 2 -r -n > NA12878_WGS_bx_counts.txt

sambamba view -f bam -F "unmapped" -t 10 NA12878_WGS.sorted.bam > ../unmapped/NA12878_WGS.sorted.bam

bxtools split COLO829-10x_phased_possorted_bam.bam -a bx_split/COLO829-10x -m 1000 > bx_split/COLO829-10x_bx_counts.txt
```

In `/data/cephfs/punim0010/projects/Saveliev_10X/NA12878-10x/insertions/unmapped_or_mate_is_unmapped/NA12878_WGS_bx_counts.txt`: 
Barcodes: 1,460,274
Reads: 19,613,313

EMA count command on unmapped_or_mate_unmapped reads:
R1: :: Reads with OK barcode: 12,195 out of 5,232,233
R2: :: Reads with OK barcode: 23,052 out of 5,232,233


TODO: check why barcodes are bad? Howe to extract on good barcodes for denovo assembly?

---------

Supernova failed with:

```
Running onfinish handler...

[error] The fraction of input reads having valid barcodes is 1.68 pct, whereas the ideal is at least 80 pct.  This condition could have multiple causes including wrong library type, failed library construction and low sequence quality on the barcode bases.  This could have a severe effect on assembly performance, and Supernova has not been tested on data with these properties, so execution will be terminated.

The following warning(s) were issued prior to encountering an error:


We observe only 36.67 pct of bases on read two with quality scores at least 30. Ideally, we expect at least 75 pct. Data quality issues of this type are difficult to diagnose, but might be caused by defects in sequencing reagents or sequencing instrument condition. Use of low quality reads usually reduces assembly quality.
```

Using atropos to get reads with high quality:

```
#!/bin/bash

FQ1=$1
FQ2=${FQ1/R1/R2}

BASE=$(basename $INP .fq)
OUTDIR=trimmed
mkdir $OUTDIR
fq1out=${OUTDIR}/${BASE}.R1.fq
fq2out=${OUTDIR}/${BASE}.R2.fq

atropos trim \
  --quality-base 33 --format fastq \
  --overlap 8 --no-default-adapters \
  -pe1 ${FQ1} \
  -pe2 ${FQ2} \
  --threads 8 \
  -o ${fq1out} \
  -p ${fq2out} \
  --report-file atropos-report_${BASE}.yaml \
  --report-formats yaml --stats both \
  --sample-id ${BASE} \
  --quality-cutoff=5 --minimum-length=33 --nextseq-trim=25
```

And EMA's preproc to extract reads with correct Barcodes:

```
parallel -j1 "bash interleave_fq.sh {} | ema count -w 4M-with-alts-february-2016.txt -o {/.}" ::: $1 ../atropos/NA12878_WGS.R1.fq

ls ../atropos/NA12878_WGS.R1.fq | xargs -I '{}' bash interleave_fq.sh '{}' | ema preproc -w 4M-with-alts-february-2016.txt -n 1 -t 1 -o ema_work *.ema-ncnt 2>&1 | tee ema_preproc.log
```






