# Brief background on telomeres.
# XXX: does sphinx support markdown?

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

See README.md

# Plan

At least as Roman imagines it goes at this point... XXX's represent gaps/doubts.

## Calculate number of reads in telomeric regions (absolute, relative to coverage)

1. We will use hg38 to avoid problems with [missing chromosome 17 and other assembly artifacts][chr17_hg19].
2. The UCSC browser assumption that telomeric reads are within 1 to 10000bp is wrong. In reality telomeric hexamers are placed in varying positions after unmapped reads. Those should be determined first, disregarding UCSC's generated BED file.

`bedtools intersect -a data/telomere_regions_GRCh37.bed -b data/EGAZ00001226275_COLO_829_TGEN_IlluminaPipe.bam -c > telomer_counts.txt`

### COLO829BL (regular Illumina run)
#### Multimapped reads

```
$ samtools view final/hg38/COLO829BL-hg38-ready.bam -L bed/hg38_noalt.bed -b > final/hg38/COLO829BL-hg38-ready-without_randoms.bam
$ samtools view -H final/hg38/COLO829BL-hg38-ready-without_randoms.bam | grep -v GL | grep -v KI > COLO829BL-hg38-ready-without_randoms-torehead.sam
$ samtools reheader COLO829BL-hg38-ready-without_randoms-torehead.sam COLO829BL-hg38-ready-without_randoms.bam > COLO829BL-hg38-ready-without_randoms-reheaded.bam
$ samtools view COLO829BL-hg38-ready-without_randoms-reheaded.bam | cut -f 3 | sort | uniq -c | sort -rn
18779 chr5
6604 chr1
3104 chrX
2852 chr20
2654 chr12
2623 chr10
2049 chr4
2008 chr21
1932 chr3
1917 chr15
1699 chr18
1680 chr16
1599 chr22
1505 chr6
1303 chr7
1268 chr13
1206 chr9
 975 chr2
 926 chr17
 892 chr11
 827 chr19
 786 chr14
  77 chr8
```

## Use kmers if the above looks good

I wrote a brute script that assesses the exact hexamer count for each chromosome in GRCh38 (2013) reference assembly, example outputs can be found on `latest_output.txt` on this repo.

## Run counts for 10X dataset

The hypothesis is that using 10X technology we can reduce the uncertainty of multimapped reads (`NH:i:n, n>2`).

So having a base multimapping read count from Truseq will inform how good 10X performs.


## Segment/mask the 10X dataset per chromosome and run telomerecat and other tools against those segments/bams?

A "mask and re-run" approach **per chromosome and telomere region(s)** might be computationally too expensive on big BAM files. On the other hand, if we just focus on the multimapping of reads for QC purposes, for instance, it can run fast.

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
[chr17_hg19]: https://www.biostars.org/p/72730/#72759