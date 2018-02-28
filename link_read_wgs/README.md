## 10X Whole Genome Sequencing

Test data is downloaded from 10X from [here](https://support.10xgenomics.com/de-novo-assembly/datasets/2.0.0/wfu).

Running _LongRanger_ on test data works using the command in script [`longranger-test.sh`](longranger-test.sh). This runs _LongRanger_ in cluster mode and successfully launched tasks from the master job.

*KNOWN ISSUE:*
This runs _LongRanger_ in `wgs` mode on only ONE flowcell of data (two flowcells in downloaded data). To use data from multiple flowcells (unlikely our use-case), refer to 10X instructions [here](https://support.10xgenomics.com/genome-exome/software/pipelines/latest/advanced/multi-flowcell).

