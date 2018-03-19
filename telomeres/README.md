# Brief background on telomeres.

Telomeres are nucleoprotein structures at the ends of all linear
eukaryotic chromosomes. In mammals, telomeres are composed
of the shelterin protein complex and a double-stranded
tract of short tandem repeats of TTAGGG that ends in a singlestranded
3 overhang on the G-strand. Telomeres cap and
protect chromosome ends from eliciting a DNA damage
response and illegitimate recombination events.

As a cell begins to become cancerous, it divides more often and its telomeres are shortened at a higher rate. Cancerous cells escape senescence/death and conversely, they become immortal with the ability to indefinitely replicate, even with shortened telomeres, resulting in uncontrolled tumor growth. It has been demonstrated that this ability to escape the natural cell outcome happens because cancer cells prevent telomeres from becoming critically short by the reactivation of telomerase (phenomenon occurring in 85% of cancers) or by ALT (alternative lengthening of telomeres; 15%).

# Telomeric regions (BED file)

How to find the positions of telomeres? - https://genome.ucsc.edu/FAQ/FAQtracks#tracks20

# Research question(s)

1. Using 10X, is it possible to find the telomere length and/or structure on a chromosome level basis?

Context: Current tools like TelomereHunter et al assess the length on a whole genome basis when using (traditional?) Illumina Truseq sequencing, but we might be able to have a closer view of those by using [linked reads][linked_reads].

# Project organization

There are two directories available from:

`/data/cephfs/punim0010/projects/Valls_10X_Telomeres/*`

`software`: Containing locally edited/patched software `pip install -e .`, like telomerecat. All conda-installed on root environment.
`work`: Data processing and batch scripts, with the following structure (partly resembling bcbio):

```
$ ls -alh
total 4.0K
drwxr-xr-x 1 brainstorm vho 534G Mar 19 15:32 .
drwxr-xr-x 1 brainstorm vho 534G Mar 19 15:29 ..
drwxr-xr-x 1 brainstorm vho 534G Mar 19 15:17 bams             <--- intermediate bams
drwxr-xr-x 1 brainstorm vho 7.1K Mar 19 15:07 bed              <--- intermediate beds
drwxr-xr-x 1 brainstorm vho  13M Mar 19 15:32 in               <--- original input datasets
drwxr-xr-x 1 brainstorm vho 5.5M Mar 19 15:33 final            <--- "final" results
drwxr-xr-x 1 brainstorm vho 3.0K Mar 19 15:21 scripts          <--- slurm scripts used to generate above data
drwxr-xr-x 1 brainstorm vho  14K Mar 19 15:15 slurm_out_err    <--- output/error results from slurm jobs from scripts
```

The input datafiles used **will** be 10X, but we are using COLO829 Truseq to assess current tools:

```
$ ls -alh in/
total 13M
drwxr-xr-x 1 brainstorm vho  13M Mar 19 15:32 .
drwxr-xr-x 1 brainstorm vho 534G Mar 19 15:35 ..
lrwxrwxrwx 1 brainstorm vho  101 Feb 28 11:04 17MHP002Bld-ready.bam -> /data/cephfs/punim0010/data/Results/Avner/2016.249.17.MH.P002/final/17MHP002Bld/17MHP002Bld-ready.bam
-rw-r--r-- 1 brainstorm vho  13M Feb 28 13:38 17MHP002Bld-ready.bam.bai
lrwxrwxrwx 1 brainstorm vho  118 Feb 28 14:24 COLO829BL_TGEN_bwa-ready.bam -> /data/cephfs/punim0010/projects/Saveliev_COLO829_Craig/bcbio_bwa/final/COLO829BL_TGEN_bwa/COLO829BL_TGEN_bwa-ready.bam
lrwxrwxrwx 1 brainstorm vho  122 Feb 28 14:25 COLO829BL_TGEN_bwa-ready.bam.bai -> /data/cephfs/punim0010/projects/Saveliev_COLO829_Craig/bcbio_bwa/final/COLO829BL_TGEN_bwa/COLO829BL_TGEN_bwa-ready.bam.bai
lrwxrwxrwx 1 brainstorm vho  114 Feb 28 14:27 COLO829_TGEN_bwa-ready.bam -> /data/cephfs/punim0010/projects/Saveliev_COLO829_Craig/bcbio_bwa/final/COLO829_TGEN_bwa/COLO829_TGEN_bwa-ready.bam
lrwxrwxrwx 1 brainstorm vho  118 Feb 28 14:27 COLO829_TGEN_bwa-ready.bam.bai -> /data/cephfs/punim0010/projects/Saveliev_COLO829_Craig/bcbio_bwa/final/COLO829_TGEN_bwa/COLO829_TGEN_bwa-ready.bam.bai
lrwxrwxrwx 1 brainstorm vho   83 Mar  8 13:09 GRCh37.fa.fai -> /data/cephfs/punim0010/local/stable/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa.fai
```

TODO: Migrate SLURM scripts to proper Snakemake for reproducibility.

# Plan

At least as Roman imagines it goes at this point... XXX's represent gaps/doubts.

## Calculate number of reads in telomeric regions (absolute, relative to coverage)

XXX: bedtoolsCov/mosdepth/samtools/bbmap/bamutils... is there a tool that just does that and fast?

XXX: Regarding input datasets, we have COLO829, `17MHP002Bld-ready.bam`, any preference on which to attack first?

## Use kmers if the above looks good

XXX: Why?

## Run counts for 10X dataset

XXX: Not sure how to leverage/use the linked reads info, need to read up on tools enabling this analysis

## Segment/mask the 10X dataset per chromosome and run telomerecat and other tools against those segments/bams?

XXX: Not possible since telomerecat needs context, then we should just use a "mask and re-run" approach **per chromosome and telomere region(s)**, computationally too expensive on big BAM files?

XXX: Again, can we leverage some 10X trick to make this easier?


# Telomerecat output reference

Working off [readthedocs](http://telomerecat.readthedocs.io/en/latest/understanding_output.html).


Sample|F1|F2|F4|Psi|Insert_mean|Insert_sd|Read_length|Initial_read_length|F2a|F2a_c|Length
------|--|--|--|---|-----------|---------|-----------|-------------------|---|-----|-------
COLO829BL_TGEN_bwa-ready_except_telomere_regions-output.bam|49532|7350|2962|3.034|346.0|84.396|112|112|4388|4388|2972.7
COLO829_TGEN_bwa-ready_except_telomere_regions-output.bam|39034|7673|2342|2.615|397.0|97.21799999999999|112|112|5331|5331|2476.5
COLO829_TGEN_bwa-ready.bam|38605|8012|2528|2.573|397.0|96.979|112|112|5484|5484|2403.0
COLO829BL_TGEN_bwa-ready.bam|49507|7545|3066|2.9960000000000004|346.0|84.266|112|112|4479|4479|2928.2
sample name | \# of reads which are completely telomeric | \# of reads where one end is completely telomeric, and the other is not, and telomeric end is CCCTAA | as previous column, and end is TTAGGG | Measure of fidelity | Insert size | SD of insert size | | | F2 - F4 

F2-F4 = F2a, the estimated number of reads covering the boundary between telomere and nontelomere.
 
The greater the measure of fidelity, the more we believe the observed measurement of F2a. (Question - so how large can this value be?, how large is this value with a clean sample?)
 
F2a value after undergoing correction (see paper).
 
Length - telomere length as estimated by telomerecat.


[linked_reads]: https://www.10xgenomics.com/linked-reads/
