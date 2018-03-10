# 10X Dumping Ground

Keep 10X notes, reports, results, (small) datasets etc. as you wish.
Basically instead of dumping everything in the
[Google doc](https://docs.google.com/document/d/1EhqPusGRCDKdK5tx5RpEhgwj_LCAi7plb2B62VvbaG4/edit),
feel free to dump it in here.

You can either edit/add files online through GitHub, or clone this
repo locally (`git clone https://github.com/umccr/10x`) and go nuts.
Occasional `git pull`s are recommended. Ask around if you need help ;-)

Just try not to add any 'big' files e.g. > 5Mb. I mean, you _can_, but maybe
you _shouldn't_, since it will take its toll on the repo pulling/pushing
times. I think. Instead you can upload to DropBox or something and
simply link to it on Slack.



## FASTQ Generation - bcl2fastq
10X recommended to generate FASTQ files (at least for single-cell 3', depend on length of reads).

```
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

A cluster scheduler template is created for _CellRanger_ AND _LongRanger_ cluster mode job submission (see [`slurm.template`](scripts/slurm.template)).


