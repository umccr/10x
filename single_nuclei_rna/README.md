## Single Nuclei RNA-seq Data

Run `180226_A00130_0038_BH3JKFDMXX` contains single nuclei RNA-seq data.

Job submission script for performing *`cellranger count`* is created, but uses array tasks as oppose to _CellRanger_ cluster mode (see [`cellranger-count.sh`](single_nuclei_ran/cellranger-count.sh)).

*NOTE FROM LUCIANO:* Two of the samples, `13_nc_5prime_CUP_1184_B-cell` and `14_nc_5prime_CUP_1184_T-cell` are *VDJ* libraries that requires 2x150 reads whereas 26+98 were sequenced. So they are excluded from CellRanger analysis.

Run time for *`cellranger count`* is approximately 3-4 hours for this data.

*`cellranger aggr`* is run after *`cellranger count`* for all samples are completed. Tasks are separated for 5' and 3' libraries (see [`cellranger-aggr-3p.sh`](single_nuclei_ran/cellranger-aggr-3p.sh) and [`cellranger-aggr-5p.sh`](single_nuclei_ran/cellranger-aggr-5p.sh)).

*PROBLEM:*
After the analyses are completed. `web_summary.html` gives serious warnings about "Low Fraction Reads in Cells".

*SOLUTION:*
Oliver and Lavinia picked up that a pre-mRNA reference needs to be used for single nuclei RNA data to allow reads from introns of pre-splicing transcripts to be counted as "in cell".

Lavinia found the official 10X guide to ["How to make pre-mRNA reference"](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references#premrna) for single nuclei RNA-seq.
