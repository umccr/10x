# 10x
10X dumping ground

## FASTQ Generation - bcl2fastq

[Arthur]: I followed what 10X recommended to generate FASTQ files:

```
# recommended bcl2fastq options (for SC3prime?) by 10x-genomics
/usr/local/bin/bcl2fastq \
        --use-bases-mask Y26n*,I8,Y* \
        --create-fastq-for-index-reads \
        --minimum-trimmed-read-length=8 \
        --mask-short-adapter-reads=8 \
        --ignore-missing-positions \
        --ignore-missing-controls \
        --ignore-missing-filter \
        --ignore-missing-bcls \
        -R ${IN} \
        --sample-sheet ${IN}/SampleSheet.csv \
        -o ${OUT}/ \
        --no-lane-splitting
```

Sample sheet generation was using 10X website’s [generator](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/bcl2fastq-direct)

*PROBLEM:*

This sample sheet format generates FASTQs that are have these issues:

a) organised by index `SI-GA-*` rather by sample; 

b) we neglected lane number in sample sheet and bcl2fastq removed the `_L00*_` token from file name that is required by _CellRanger_ later.

*SOLUTION:*

Create symbolic links for FASTQs to have names acceptable by _CellRanger_


## Running LongRanger/CellRanger on Spartan

A cluster scheduler template is created for _CellRanger_ AND _LongRanger_ cluster mode job submission (see [`slurm.template`](slurm.template)).



