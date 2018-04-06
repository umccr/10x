# Exploring novel insertions

### Contents

- [Oncoviruses in DiploidNeverResponder](#oncoviruses-in-diploidneverresponder)
  * [Extracting unmapped reads](#extracting-unmapped-reads)
  * [Exploring taxonomic content](#exploring-taxonomic-content)
    + [Mash](#mash)
    + [BWA-MEM](#bwa-mem)
  * [De-novo assembling HPV18 region](#de-novo-assembling-hpv18-region)
  * [Looking for integration sites](#looking-for-integration-sites)
  * [chr8 integration site](#chr8-integration-site)
- [Going with 10x](#going-with-10x)
  * [Strategies](#strategies)
  * [Datasets](#datasets)
  * [Unmapped read extraction](#unmapped-read-extraction)
  * [GEMs with a high unmapped rate](#gems-with-a-high-unmapped-rate)
  * [Trying de novo assembly of unmapped reads with Supernova](#trying-de-novo-assembly-of-unmapped-reads-with-supernova)
  * [Extracting viral sequences](#extracting-viral-sequences)
    + [Mapping to HPV71](#mapping-to-hpv71)
    + [Trying reads with BX tags and high quality](#trying-reads-with-bx-tags-and-high-quality)
  * [Viruses in COLO829](#viruses-in-colo829)

# Oncoviruses in DiploidNeverResponder

## Extracting unmapped reads

We want to eplore viruses that might be intergrated into the genome. In order to do that, we start with getting reads unmapped to the human genome, or reads with unmapped mate which can support integration breakpoints.

```
sambamba view -f bam -F "unmapped or mate_is_unmapped" -t 20 diploid_tumor-ready.bam | samtools sort diploid_tumor-unmapped_or_mate.unsorted.bam -n -@ 10 > diploid_tumor-unmapped_or_mate.namesorted.bam
```

Then we filter high quality reads:
- all(q >= 10 for q in read.query_qualities) and mean(read.query_qualities) >= 25
- read.query_length >= 125
- paired

```
./filter_bam_highqual.py diploid_tumor-unmapped_or_mate.namesorted.bam
Total:              28,564,510 out of 2,230,924,710
lng_i:              27,894,774
hqual_i:            15,326,411
lng_hqual_i:        14,878,067
paired_i:           28,564,510
paired_lng_hqual_i: 14,878,067
```

Output written to `lng_hqual/diploid_tumor-unmapped_or_mate.namesorted.bam`. Extracting reads:

```
cd lng_hqual
samtools fastq diploid_tumor-unmapped_or_mate.namesorted.bam -1 diploid_tumor.R1.fq -2 diploid_tumor.R2.fq -s diploid_tumor.single.fq
[M::bam2fq_mainloop] discarded 11217916 singletons
[M::bam2fq_mainloop] processed 14878054 reads
```

## Exploring taxonomic content

Querying the target reads against the oncovirus database. First using a fast k-mer approach to give a quick outline, then using exact BWA-MEM alignment to get reads aligned to particular species.

### Mash

Applying a fast [k-mer based approach](http://mash.readthedocs.io/en/latest/tutorials.html?highlight=screen#screening-a-read-set-for-containment-of-refseq-genomes). First, querying against a oncoviral database from GDC:

```
mkdir mash
mash sketch /g/data3/gx8/extras/vlad/bcbio/genomes/Hsapiens/GRCh37/viral/gdc-viral.fa -i -o mash/gdc-viral
mash screen mash/gdc-viral.msh diploid_tumor.R1.fq diploid_tumor.R2.fq | sort -gr > mash/screen_gdc.tab
head mash/screen_gdc.tab

0.99726   944/1000   100  0             HPV18
0.879223  67/1000     49  4.56495e-199  HPV45
0.876649  63/1000     29  1.50084e-185  HPV97
0.835578  23/1000      5  1.51569e-58   HPV59
0.816057  14/1000     91  3.39476e-33   HPV68
0.816057  14/1000      6  3.39476e-33   HPV85
0.816057  14/1000   3712  3.39476e-33   HPV39
0.816057  14/1000      2  3.39476e-33   HPV20
0.806739  11/1000      4  3.07798e-25   HPV70
0.806739  11/1000      4  3.07798e-25   HPV66
```

HPV18 is an obvious candidate. Querying against a full RefSeq just in case:

```
mash screen mash/refseq.genomes.k21s1000.msh diploid_tumor.R1.fq diploid_tumor.R2.fq | sort -gr > mash/screen_refseq.tab
(took 4:14.57 time with 100% CPU total)

0.996652  932/1000  98  0             GCF_000865665.1_ViralProj15506_genomic.fna.gz   NC_001357.1 Human papillomavirus - 18, complete genome
0.97558   595/1000   1  0             GCF_000819615.1_ViralProj14015_genomic.fna.gz   NC_001422.1 Enterobacteria phage phiX174 sensu lato, complete genome
0.932215  229/1000   1  0             GCF_000913595.1_ViralProj222261_genomic.fna.gz  NC_022518.1 Human endogenous retrovirus K113 complete genome
0.916976  162/1000   2  0             GCF_000287295.1_ASM28729v1_genomic.fna.gz       NC_018417.1 Candidatus Carsonella ruddii HT isolate Thao2000, complete genome
0.887203  81/1000    2  3.28623e-247  GCF_000829155.1_ASM82915v1_genomic.fna.gz       NZ_AP013293.1 Candidatus Sulcia muelleri PSPU DNA, complete genome
```

Confirmed HPV18.

### BWA-MEM

First aligning to the GDC oncoviruses fasta:

```
mkdir viral_mapping
bwa mem -t 10 /g/data3/gx8/extras/vlad/bcbio/genomes/Hsapiens/GRCh37/viral/gdc-viral.fa diploid_tumor.R1.fq diploid_tumor.R2.fq | bamsort inputthreads=10 outputthreads=10 inputformat=sam index=1 indexfilename=viral_mapping/diploid_tumor_viral.bam.bai O=viral_mapping/diploid_tumor_viral.bam
[main] Real time: 38.855 sec; CPU: 333.034 sec
[V] wrote 3666269 alignments

samtools idxstats viral_mapping/diploid_tumor_viral.bam | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ($2 != 0) ? $3/$2 : 0 }' | sort -nr -k4,4 | head | cols

name     len  mapped  coverage           normilize by dividing by 1173944366/1000000000 = 1.174 billion reads in genome
HPV18    7857  151264  19.2521           16.3987223169
HCV-2    9711   12828  1.32098
HPV71    8037    2661  0.331094
HCV-1    9646    2862  0.296703
HPV82    7870     242  0.0307497
HPV19    7685     212  0.0275862
HPV21    7779      59  0.00758452
HPV20    7757      37  0.00476989
HPV25    7713      34  0.00440814
```

Aligning to HPV18 only to make sure none of the related virus didn't hijacked any reads:

```
samtools faidx /g/data3/gx8/extras/vlad/bcbio/genomes/Hsapiens/GRCh37/viral/gdc-viral.fa HPV18 > viral_mapping/HPV18.fa
bwa index viral_mapping/HPV18.fa

bwa mem -t 10 viral_mapping/HPV18.fa diploid_tumor.R1.fq diploid_tumor.R2.fq | bamsort inputthreads=10 outputthreads=10 inputformat=sam index=1 indexfilename=viral_mapping/to_HPV18.bam.bai O=viral_mapping/to_HPV18.bam

samtools idxstats viral_mapping/to_HPV18.bam | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ($2 != 0) ? $3/$2 : 0 }'
# HPV18   7857    151271  19.253
```

Getting same 19x coverage, which is a very strong evidence.

## De-novo assembling HPV18 region

Getting reads that map to HPV18 and all their mate that might point us to an integration site:

```
cd /g/data3/gx8/projects/Saveliev_Diploid/lng_hqual/viral_mapping

sambamba view -F "not unmapped or not mate_is_unmapped" -f bam -o to_HPV18.mapped_or_mate.bam to_HPV18.bam
# will fail with segmentation fault, however will generate the output. Alternative option is to `samtools view HPV18` - should do the same.

# Make sure to name-sort reads before feeding into samtools fastq - otherwise the output will be garbage:
samtools sort -n to_HPV18.mapped_or_mate.bam -O bam -o to_HPV18.mapped_or_mate.namesorted.bam

samtools fastq to_HPV18.mapped_or_mate.namesorted.bam -1 to_HPV18.mapped_or_mate.R1.fq -2 to_HPV18.mapped_or_mate.R2.fq -s to_HPV18.mapped_or_mate.single.fq
[M::bam2fq_mainloop] discarded 0 singletons
[M::bam2fq_mainloop] processed 161716 reads
```

Assembling with SPAdes:

```
spades.py --only-assembler -1 to_HPV18.mapped_or_mate.R1.fq -2 to_HPV18.mapped_or_mate.R2.fq -s to_HPV18.mapped_or_mate.single.fq -o spades

# QC the assembly - stats and alignment back to HPV18
quast.py spades/contigs.fasta -R HPV18.fa -o spades/quast --ref-bam to_HPV18.mapped_or_mate.namesorted.bam --no-read-stats --no-sv -1 to_HPV18.mapped_or_mate.R1.fq -2 to_HPV18.mapped_or_mate.R2.fq --debug

# Align contigs to the reference:
minimap2 -a HPV18.fa spades/contigs.fasta | samtools sort > spades/contigs_to_HPV18.bam && samtools index spades/contigs_to_HPV18.bam
```

Reads getting assembled pretty well into 4 long contigs:

![quast](assemble/diploid/img/quast_HPV18.png)

![icarus](assemble/diploid/img/icarus_HPV18.png)

All 4 contigs touch each other and cover the assembly fully uniformly, which is a strong evidence of the virus being present in the dataset with possible breakpoints in integration sites. NODE_2, NODE_4 and NODE_4 are amplified heavily up to 5kX coverage, and NODE_1 has a much smaller coverage. The genome is clearly circular - NODE_2 spans the edges, and the breakpoing supported by discordant pairs:

![igv](assemble/diploid/img/igv_HPV18.png)

## Looking for integration sites

To identify if the virus is integrated into human, we need to review the viral read whose mates map to the human genome. We would expect reads of forward orienation to pile up to the left of the breakpoint, and reads of the opposite orientation to pile up to the right of the breakpoint; a simmetrical picture would be expected in the human chromosome as well. In order to figure that out, we are creating a GRCh37-HPV18 spike reference, and remap reads against this artificial genome:

```
mkdir GRCh37_HPV18
cat /g/data3/gx8/local/development/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa > GRCh37_HPV18/GRCh37_HPV18.fa
cat HPV18.fa >> GRCh37_HPV18/GRCh37_HPV18.fa

# Align
bwa index GRCh37_HPV18/GRCh37_HPV18.fa
bwa mem -t 10 GRCh37_HPV18/GRCh37_HPV18.fa to_HPV18.mapped_or_mate.R1.fq to_HPV18.mapped_or_mate.R2.fq | bamsort inputthreads=10 outputthreads=10 inputformat=sam index=1 indexfilename=GRCh37_HPV18/to_GRCh37_HPV18.bam.bai O=GRCh37_HPV18/to_GRCh37_HPV18.bam

OR EVEN (doesn't require index):
minimap2 -ax sr GRCh37_HPV18/GRCh37_HPV18.fa to_HPV18.mapped_or_mate.R1.fq to_HPV18.mapped_or_mate.R2.fq | bamsort inputformat=sam index=1 indexfilename=to_GRCh37_HPV18_mm2.bam.bai O=to_GRCh37_HPV18_mm2.bam
```

HPV18 region of the spiked reference in IGV, groupped by the chromosome of mate (based on [the color scheme](https://software.broadinstitute.org/software/igv/interpreting_insert_size) - the small green group is from chromosome 3, large purple group is from chromosome 8, and all gray read pairs are entirely mapped to HPV18). Read direction is annotated with arrows.

![igv](assemble/diploid/img/igv_spike_HPV18.png)

There is a quite evident support for human genome integration to chromosome 8 (with purple reads), and a smaller support for the chromosome 3 integration with the green reads. All those colored read piles are consitently one direction and their mates simmetrically cluster in GRCh37.

The leftmost purple reads mates are mapping mapping to 8:128,303,800-128,304,900 

![igv](assemble/diploid/img/igv_spike_integration_site_chr8_leftmost.png)

The second purple pile in the beginning of NODE_3 has its mates mapping 16kb apart - to 8:128,319,000-128,321,00:

![igv](assemble/diploid/img/igv_spike_integration_site_chr8_right_1.png)

The rightmost purple pile mates is mapping nearby the leftmost pile - to 8:128,307,000-128,308,000

![igv](assemble/diploid/img/igv_spike_integration_site_chr8_right_2.png)

This suggest that the deeply covered NODE_2-NODE_3-NODE_4 part of the virus got inserted into chromosome 8.

Also, the green reads map to 3:186,691,636-186,699,490, suggesting possible another integration site. There might be also another less abundant strain of virus, worth exploring other viral sequences to confirm - spiking GRCh37 with the whole GDC database and realigning all human-unmapped reads.

Another interesting obervation comes from blasting the contigs against all human sequences: it reports 18% coverage in a chr8 alternative assembly, and same 18% in hg38:

![blast](assemble/diploid/img/blast_HPV18.png)

It would be interesting to repeate whole experiment with hg38.

## chr8 integration site

The read orientations spanning the breakpoints suggests a quite complex event rather than a simple insertion of the virus. Exploring the full 16kb GRCh37 region covering all of 3 breakpoints:

```
~/bin/sambamba slice 8:128318910-128321289 diploid_tumor-ready.bam > diploid_tumor-chr8_HPV18.bam
```

![igv](assemble/diploid/img/igv_GRCh37_full_region.png)

We can see that 2 of the breakpoint positions show up very clearly, and it's also evident that:

- The NODE_2-NODE_3-NODE_4 region is amplified heavily instead of being deleted, meaning that the virus didn't replace the region with itself resulting in some more complex event;
- The read orientations tell us that the viral amplified region got attached to the rightmost chr8 breakpoint from the right side, and to the left breakpoints from the left side, which rejects the idea of a simple insertion;
- chr8 region between the breakpoints is heavily amplified, at roughly the same coverage as the viral amplified regon.

That suggests that the virus and created a loop by attaching to the leftmost and rightmost breakpoints, and this loop went around many times, heavily amplifying the chr8 16kb region as well as viral region NODE_2-NODE_3-NODE_4. While looping, it also likely occasionally attached to the inner left breakpoint as well. Indeed, things like this are typical for HPV viruses:

![igv](assemble/diploid/img/HPV_loop_from_paper.png)

In our case, the integration events might have led to the following sequence (letters are from the IGV screenshot above): 

A - B - C - (NODE_4-NODE_3-NODE-2 - B - C)*10 - NODE_4-NODE_3 - C - D

The integration site overlaps long non-coding RNA genes CASC21 (Cancer Susceptibility 21, CARLo-2) and CASC8 (Cancer Susceptibility 8, CARLo-1) in their introns. Both genes are associated with cancer, and located in a region 8q24.21 nearby the oncogene MYC.

8q24.21 is well known as HPV integration site hotspot: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4695887/, https://www.ncbi.nlm.nih.gov/pubmed/1649348/: "We have determined the chromosomal localization of integrated HPV type 16 (HPV-16) or HPV-18 genomes in genital cancers... In three cancers, HPV sequences were localized in chromosome band 8q24.1, in which the c-myc gene is mapped... In three of the four cases, the proto-oncogene located near integrated viral sequences was found to be structurally altered and/or overexpressed. These data indicate that HPV genomes are preferentially integrated near myc genes in invasive genital cancers and support the hypothesis that integration plays a part in tumor progression via an activation of cellular oncogenes.", also https://www.nature.com/articles/1207006: "RS–PCR of HPV18-positive tumors revealed a single large cluster at 8q24."

----

## Playground

#### chr3 integration site

```
~/bin/sambamba slice diploid_tumor-ready.bam 3:186,691,636-186,699,490 > diploid_tumor-chr3_HPV18.bam
```

![igv](assemble/diploid/img/igv_GRCh37_chr3_integration_site.png)

Not clear about this one.

#### Aligning to GRCh38

#### Structural variants in this genome location

#### Expreriment with 10x data

# Going with 10x

## Strategies

Next we want to explore if we can more efficiently identify novel insertions using the 10x barcoded data. We explore the following strategies.

1. By Oliver
	1. Extract read group (GEM) tags from unmapped reads
	2. Calculate #/% reads mapped/unmapped for each GEM in that list
	3. Identify outlier GEMs with high % of unmapped reads
	4. Blast/align reads against viral / retrotransposon reference sequences
	5. See if we can identify the integration site for any hits

2. Based on the LongRanger SV caller results (that looks for distant pairs of loci in the genome that share many more barcode than would be expected by chance, then refines the candidate SVs by comparing the layout of reads and barcodes around the event the patterns expected in deletions, inversions, duplications, and translocations to identify the SV type and find the maximum-likelihood breakpoints.)
	1. Find breakpoints with LongRanger SV caller.
	2. List reads that map in proximity with a breakpoint (in both directions), but have an unmapped mate
	3. List BX of those reads
	4. Find all unmapped reads (or reads with unmapped mate) with those barcodes
	5. De novo assemble those reads (with 10x-aware [Supernova](https://support.10xgenomics.com/de-novo-assembly/software/pipelines/latest/using/running) or Athena, or just with SPAdes)

3. Directly assembly all unmapped reads with Supernova and see what happens.

4. By https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-017-0447-8: "We sought to determine whether we could resolve and thus validate the rearranged structure by de novo assembly: 
	1. We extracted all the sequence reads containing SV-specific barcodes from the linked read fastq files 
	2. and then used these subset fastq files as input to the Supernova de novo assembly program to generate contig sequences [29]. This assembler has been recently demonstrated to generate full diploid assemblies. 
	3. We visualized the structures of the resulting contigs by plotting the mapping position of each SV-specific read in the genome versus its mapping position in the contig."

We are going to start exploring with the first strategy.


## Datasets

We are working with 2 10x samples: [NA12878_WGS downloaded from the 10x website](https://support.10xgenomics.com/genome-exome/datasets/2.2.1/NA12878_WGS_v2), and COLO829-10x sequenced in UMCCR. Both datasets were pre-processed with LongRanger to get `NA12878_WGS_v2_phased_possorted_bam.bam` and `COLO829-10x_phased_possorted_bam.bam` barcode-aware GRCh37 BAM files.

We are using the following conda environment for the following processing: bioconda's `sambamba samtools pigz parallel bxtools blast snakemake spades bioawk` + `supernova` and `ema` (installed in extras). We work in the following folder on Raijin:

```
cd /data/cephfs/punim0010/projects/Saveliev_10X/COLO829-10x/insertions
cd /data/cephfs/punim0010/projects/Saveliev_10X/NA12878_WGS/insertions
```

## Unmapped read extraction

Extract unmapped reads or reads with unmapped mate:

```
mkdir unmapped_or_mate_is_unmapped

NA12878:
sambamba view -f bam -F "unmapped or mate_is_unmapped" -t 20 /data/cephfs/punim0010/projects/Saveliev_10X/NA12878-10x/NA12878_WGS_v2_phased_possorted_bam.bam | sambamba sort -n -t 20 -m 100G /dev/stdin --show-progress > unmapped_or_mate_is_unmapped/NA12878_WGS.bam

COLO829:
sambamba view -f bam -F "unmapped or mate_is_unmapped" -t 20 COLO829-10x_phased_possorted_bam.bam | sambamba sort -n -t {threads} -m 100G /dev/stdin --show-progress > unmapped_or_mate_is_unmapped/COLO829-10x.bam
```

Split by BX tag:

```
cd unmapped_or_mate_is_unmapped
mkdir bx_split
bxtools split NA12878_WGS.bam -a bx_split/NA12878_WGS -m 1000 > bx_split/NA12878_WGS_bx_counts.txt
bxtools split COLO829-10x.bam -a bx_split/COLO829-10x -m 1000 > bx_split/COLO829-10x_bx_counts.txt
```

Extract tags with at least 1000 reads:

```
sort -k 2 -r -n NA12878_WGS_bx_counts.txt | py -x "x if int(x.split()[1]) >= 1000 else None" > NA12878_WGS_bx_counts_cnt1000.txt
sort -k 2 -r -n COLO829-10x_bx_counts.txt | py -x "x if int(x.split()[1]) >= 1000 else None" > COLO829-10x_bx_counts_cnt1000.txt
```

Count all reads in full BAM by BX:

```
bxtools stats /data/cephfs/punim0010/projects/Saveliev_10X/NA12878-10x/NA12878_WGS_v2_phased_possorted_bam.bam > /data/cephfs/punim0010/projects/Saveliev_10X/NA12878-10x/qc/NA12878_WGS_v2_phased_possorted_bam.bxstats.txt

bxtools stats /data/cephfs/punim0010/projects/Saveliev_10X/NA12878-10x/NA12878_WGS_v2_phased_possorted_bam.bam > /data/cephfs/punim0010/projects/Saveliev_10X/NA12878-10x/qc/NA12878_WGS_v2_phased_possorted_bam.bxstats.txt

```

Extract BX with many unmapped:

```
grep -f <(cut -f1 NA12878_WGS_bx_counts_cnt1000.txt) /data/cephfs/punim0010/projects/Saveliev_10X/NA12878-10x/qc/NA12878_WGS_v2_phased_possorted_bam.bxstats.txt > NA12878_WGS_full_bx_stats_cnt1000.txt
```

(Optionally, we can split the full BAM by BX tag:)

```
bxtools split NA12878_WGS_v2_phased_possorted_bam.bam -a bx_split/NA12878_WGS -m 1000 > bx_split/NA12878_WGS_bx_counts.txt
```

## GEMs with a high unmapped rate

The following barcodes have at least 1000 unmapped reads or reads with unmapped mate:

```
                     Aln  Unaln | Mate unaln  Unaln  Aln & Mate unaln   IS  MQ
CATGTAGAGTCGTACT-1  4037              94.70% 94.43%                11  237  60
GCCTGTTAGCCCGTGT-1  1342              86.44% 85.77%                 9  348  60
ACTGTAGGTTCGTAAC-1  1647              84.03% 83.79%                 4  347  60
CATGTAGGTAGATGTA-1  1250              84.00% 83.28%                 9   19  19
GCCTGTTGTTCGCGGT-1  1921              73.30% 72.98%                 6  348  60
CATGTAGAGAGTCTCT-1  1828              72.87% 72.10%                14   20  15
GCCTGTTGTCTTACCC-1  1827              72.03% 71.65%                 7  374  60
GCCTGTTAGTCGAGTG-1  1515              69.37% 68.84%                 8  349  60
TCAGTAGGTCATGTAC-1  2070              69.47% 68.41%                22  304  60
GCCTGTTAGCTCCGTG-1  1714              68.14% 67.74%                 7  368  60
ACTGTAGAGGTAAACT-1  1507              67.09% 66.42%                10   20  19
GCCTGTTAGGCCCTCA-1  3324              63.96% 63.45%                17  308  60
CATGTAGGTTGATCAC-1  1649              63.61% 62.77%                14  315  60
GCCTGTTGTTTCGTTT-1  1933              63.68% 62.44%                24  322  60
GCCTGTTAGTCGGCTC-1  2126              61.90% 61.43%                10  339  60
GCCTGTTGTTTGGCTA-1  1921              64.29% 56.95%               141  314  60
TGCGCAGGTACTGTAC-1  2327              56.81% 55.82%                23  334  60
GCCTGTTGTCCGGGTT-1  2053              52.70% 50.02%                55  288  49
GCCTGTTGTACGTTTG-1  2316              50.69% 50.00%                16  286  60
AAGCCCGCATCCCACT-1  2524              41.32% 37.08%               107  264  20
GTGCCTTCAGATCTGT-1  7461              21.61% 20.68%                69  336  60
GTTGCAACAGATCTGT-1  5943              19.60% 15.01%               273  329  60
```

For each tag, we extract mapped reads with unmapped mates:

```python
!cd /data/cephfs/punim0010/projects/Saveliev_10X/COLO829-10x/insertions/unmapped_or_mate_is_unmapped/split_by_bx/

for tag in df.index:
    !sambamba view -F "not unmapped and mapping_quality >= 60" NA12878_WGS.{tag}-1.bam -f bam | samtools sort > NA12878_WGS.{tag}-1.mapped.bam
```

Manually reviewing the reads, we can find one BX with mapped reads spanning consistently more or less one region (possibly a breakpoint?)

```
samtools view NA12878_WGS.GCCTGTTGTTTCGTTT-1.mapped.bam

ST-E00273:259:H7WY3ALXX:1:1217:21207:15479  121 chr2  71844818  60  73M1D11M2I30M4D12M  chr7  128046340  0
ST-E00273:259:H7WY3ALXX:1:1220:8156:43677    73 chr2  71844971  60  47M6I7M2I66M        chr10    654413  0
ST-E00273:259:H7WY3ALXX:2:1118:11373:44152   73 chr2  71845284  60  52M13I2M1D61M       chr21  44035538  0
ST-E00273:259:H7WY3ALXX:1:2121:15006:67814   89 chr2  71845342  60  14S54M7I53M         =      33141434  0
ST-E00273:259:H7WY3ALXX:3:2117:1499:37664    89 chr2  71845402  60  77M2I22M1I17M1D9M   =      71845700  0
ST-E00273:259:H7WY3ALXX:1:1112:18913:57038  169 chr2  71845409  60  7M2I66M76S          *             0  0
ST-E00273:259:H7WY3ALXX:1:2124:12134:1977    73 chr2  71845440  60  79M2D21M1D24M4S     chr20  25048324  0
ST-E00273:259:H7WY3ALXX:2:1106:10602:14828  153 chr2  71845521  60  30S25M2I40M1D13M41S *             0  0
ST-E00273:259:H7WY3ALXX:2:1105:1458:67937    89 chr2  71845525  60  35S93M              chr4   72052386  0
ST-E00273:259:H7WY3ALXX:1:2203:29985:65001   89 chr2  71845588  60  9S78M2D29M2D12M     =      71844822  0
ST-E00273:259:H7WY3ALXX:3:2207:11921:48230  153 chr2  71845590  60  75S54M7I15M         *             0  0
ST-E00273:259:H7WY3ALXX:1:2120:18832:35801   89 chr2  71845617  60  46M2D79M3S          =      71845704  0
ST-E00273:259:H7WY3ALXX:1:2204:30888:68430   73 chr2  71845848  60  6S86M36S            hs37d5 31586471  0
ST-E00273:259:H7WY3ALXX:1:2204:18913:24198   73 chr4  53065963  60  23M4I19M3D38M44S    chr15  79384471  0
ST-E00273:259:H7WY3ALXX:1:1104:28270:7181   121 chr14 32925056  60  128M                chr1   30986160  0
ST-E00273:259:H7WY3ALXX:2:1101:22110:57899   73 chr14 32966430  60  128M                chr13  52303476  0
ST-E00273:259:H7WY3ALXX:2:2103:30452:36311   89 chr14 33028682  60  128M                chr18  71871932  0
ST-E00273:259:H7WY3ALXX:3:2201:23977:47316  105 chr21 28080947  60  128M                chr5   77330570  0
```

Now get all reads for that tag to see if they localized at the same region, then we'll see it in IGV:

```
cd /data/cephfs/punim0010/projects/Saveliev_10X/NA12878-10x
bamtools filter -tag BX:GCCTGTTGTTTCGTTT-1 -in NA12878_WGS_v2_phased_possorted_bam.bam -out NA12878_WGS_v2_phased_possorted_bam-GCCTGTTGTTTCGTTT-1.bam
```

Extract that region to download locally on IGV:

```
sambamba slice NA12878_WGS_v2_phased_possorted_bam.bam chr2:71844818-71847848 -o NA12878_WGS_v2_phased_possorted_chr2_breakpoint.bam
```

Extracting this one to for comparison:

```
sambamba slice NA12878_WGS_v2_phased_possorted_bam.bam chr19:33788840-33791956 -o NA12878_WGS_v2_phased_possorted_CEBPA.bam

```

In the large SVs file downloaded from the 10x website, we can find those SVs called with LOWQ: `/data/cephfs/punim0010/projects/Saveliev_10X/NA12878-10x/longranger/downloaded/NA12878_WGS_v2_large_svs.vcf.gz`

```
chr2   72650821  call_2031    A  <INV>                 5  LOWQ  END=72704553;CIEND=0,1;CIPOS=0,1;SVTYPE=INV;SVLEN=53732;PS=71845775;HAP_ALLELIC_FRAC=0.25;ALLELIC_FRAC=0.0200501253133;PAIRS=1;SPLIT=0;WildCov=.;MolTtl=.;MolTtlNoR=.;MolDel=.
chr2   72790392  call_3830    A  <DUP:TANDEM>          3  LOWQ  END=72858227;CIEND=0,1;CIPOS=0,1;SVTYPE=DUP;SVLEN=67835;PS=71845775;HAP_ALLELIC_FRAC=0.230769230769;ALLELIC_FRAC=0.00617283950617;PAIRS=1;SPLIT=0;WildCov=.;MolTtl=.;MolTtlNoR
chr2   73124465  call_2875    G  <INV>                 4  LOWQ  END=73195439;CIEND=-24,25;CIPOS=-24,25;SVTYPE=INV;SVLEN=70974;PS=71845775;HAP_ALLELIC_FRAC=0.153846153846;ALLELIC_FRAC=0.0128755364807;PAIRS=2;SPLIT=0;WildCov=.;MolTtl=.;MolT
chr2   73877193  call_2032    A  <UNK>                 5  LOWQ  END=73940170;CIEND=-8,8;CIPOS=-6,6;SVTYPE=UNK;SVLEN=62977;PS=71845775;HAP_ALLELIC_FRAC=0.172413793103;ALLELIC_FRAC=0.0178173719376;PAIRS=3;SPLIT=0;WildCov=.;MolTtl=.;MolTtlNo
chr2   76051830  call_840     C  <INV>                 8  LOWQ  END=76110415;CIEND=0,1;CIPOS=0,1;SVTYPE=INV;SVLEN=58585;PS=71845775;HAP_ALLELIC_FRAC=0.228571428571;ALLELIC_FRAC=0.0220264317181;PAIRS=1;SPLIT=0;WildCov=.;MolTtl=.;MolTtlNoR=
chr2   76213102  call_2033    A  <INV>                 5  LOWQ  END=76259489;CIEND=0,1;CIPOS=0,1;SVTYPE=INV;SVLEN=46387;PS=71845775;HAP_ALLELIC_FRAC=0.192307692308;ALLELIC_FRAC=0.0216216216216;PAIRS=1;SPLIT=0;WildCov=.;MolTtl=.;MolTtlNoR=
chr2   78278376  call_3831    T  <INV>                 3  LOWQ  END=78322534;CIEND=0,1;CIPOS=0,1;SVTYPE=INV;SVLEN=44158;PS=71845775;HAP_ALLELIC_FRAC=0.103448275862;ALLELIC_FRAC=0.0113314447592;PAIRS=1;SPLIT=0;WildCov=.;MolTtl=.;MolTtlNoR=
chr2   78415466  call_2876    A  <UNK>                 4  LOWQ  END=78476352;CIEND=-9,9;CIPOS=-9,9;SVTYPE=UNK;SVLEN=60886;PS=71845775;HAP_ALLELIC_FRAC=0.137931034483;ALLELIC_FRAC=0.0176767676768;PAIRS=3;SPLIT=0;WildCov=.;MolTtl=.;MolTtlNo
chr2   78887169  call_3832    G  <INV>                 3  LOWQ  END=78950399;CIEND=0,1;CIPOS=0,1;SVTYPE=INV;SVLEN=63230;PS=71845775;HAP_ALLELIC_FRAC=0.111111111111;ALLELIC_FRAC=0.00973236009732;PAIRS=1;SPLIT=0;WildCov=.;MolTtl=.;MolTtlNoR
chr2   79106744  call_3833    A  <UNK>                 3  LOWQ  END=79159256;CIEND=-30,31;CIPOS=-30,31;SVTYPE=UNK;SVLEN=52512;PS=71845775;HAP_ALLELIC_FRAC=0.2;ALLELIC_FRAC=0.0171568627451;PAIRS=2;SPLIT=0;WildCov=.;MolTtl=.;MolTtlNoR=.;Mol
chr2   80313256  call_1487    T  <INV>                 6  LOWQ  END=80353748;CIEND=-15,15;CIPOS=-15,15;SVTYPE=INV;SVLEN=40492;PS=71845775;HAP_ALLELIC_FRAC=0.25;ALLELIC_FRAC=0.0263157894737;PAIRS=2;SPLIT=0;WildCov=.;MolTtl=.;MolTtlNoR=.;Mo
```

Viewing in IGV: `chr2:71,844,732-71,845,762`

![alt text](img/NA12878_WGS_v2_phased_possorted_chr2_IGV.png)

"Show BX" mode:
![alt text](img/NA12878_WGS_v2_phased_possorted_chr2_IGV_showBX.png)

For comparison, this is how a regular region looks:

![alt text](img/NA12878_CEBPA_igv_snapshot.png)

![alt text](img/NA12878_CEBPA_igv_snapshot_viewBX.png)


## Trying de novo assembly of unmapped reads with Supernova

Extract reads:

```
bedtools bamtofastq -i COLO829-10x.bam -fq COLO829-10x_S1_L001_R1_001.fastq -fq2 COLO829-10x_S1_L001_R2_001.fastq 2>COLO829-10x.unpaired
bedtools bamtofastq -i NA12878_WGS.bam -fq NA12878_WGS_S1_L001_R1_001.fastq -fq2 NA12878_WGS_S1_L001_R2_001.fastq 2>NA12878_WGS.unpaired
```

Try assemble with supernova:

```
cd /g/data3/gx8/projects/Saveliev_10X/NA12878-10x/insertions/unmapped_or_mate_is_unmapped/assemble
qsub run.sh
```

However, it fails with error:

```
[error] The fraction of input reads having valid barcodes is 1.68 pct, whereas the ideal is at least 80 pct.  This condition could have multiple causes including wrong library type, failed library construction and low sequence quality on the barcode bases.  This could have a severe effect on assembly performance, and Supernova has not been tested on data with these properties, so execution will be terminated.

We observe only 36.67 pct of bases on read two with quality scores at least 30. Ideally, we expect at least 75 pct. Data quality issues of this type are difficult to diagnose, but might be caused by defects in sequencing reagents or sequencing instrument condition. Use of low quality reads usually reduces assembly quality.
```

The reason is that extracted reads do no preserve barcodes. We need to extract those reads from the original `fastq` files. To save computational time, we first want to look for particular retroviruses as play with only reads mapping to one virus.

## Extracting viral sequences

Searching viral sequences using 3 different approaches.

- Use Mash (http://mash.readthedocs.io/en/latest/tutorials.html)

```
mash sketch /g/data3/gx8/extras/vlad/bcbio/genomes/Hsapiens/GRCh37/viral/gdc-viral.fa -i -o mash/gdc-viral
mash screen gdc-viral.msh ../NA12878_WGS.merged.fq > screen_gdc.tab
sort -gr screen_gdc.tab | head
0.833811        22/1000 2       2.09944e-33     HPV20
0.828011        19/1000 4       7.40488e-28     HPV6
0.821262        16/1000 3       1.62953e-22     HPV83
0.816057        14/1000 6       4.38911e-19     HPV114
0.816057        14/1000 4       4.38911e-19     HPV19
0.816057        14/1000 3       4.38911e-19     HPV77
0.816057        14/1000 1       4.38911e-19     HPV87
0.813182        13/1000 1       2.05593e-17     HPV29
0.810088        12/1000 49      8.93544e-16     HPV12
0.810088        12/1000 2       8.93544e-16     HPV71
```

Also against RefSeq:

```
mash screen refseq.genomes.k21s1000.msh ../NA12878_WGS.merged.fq > screen_refseq.tab
```

- Kraken against minikraken (http://ccb.jhu.edu/software/kraken/MANUAL.html)

```
kraken --preload --db /g/data3/gx8/extras/kraken/minikraken_20171019_8GB --fastq-input fastq/NA12878_WGS.R1.fq --threads 9 --out kraken/kraken_out --min-hits 2 --quick
```

- Align against the viral database:

```
bwa mem -t 9 /g/data3/gx8/extras/vlad/bcbio/genomes/Hsapiens/GRCh37/viral/gdc-viral.fa NA12878_WGS.R1.fq NA12878_WGS.R2.fq | bamsort inputthreads=9 outputthreads=9 inputformat=sam index=1 indexfilename=viral_mapping/NA12878_viral.bam.bai O=viral_mapping/NA12878_viral.bam
# Real time: 180.500 sec; CPU: 1553.289 sec

samtools idxstats viral_mapping/NA12878_viral.bam | py -x "'\t'.join([x.split()[0], x.split()[1], x.split()[2], str(int(x.split()[2]) / int(x.split()[1]))]) if int(x.split()[1]) > 0 else None" | sort -nr -k4,4 | head | cols

name    mapped coverage
HPV71   87280  10.8597735473
HPV47   45913  5.94266114419
HPV82   20185  2.56480304956
HPV42   10600  1.33889099406
HCV-1   11845  1.22797014306
HCV-2   10985  1.1311914324
HPV51    8678  1.11142418033
HPV21    8375  1.07661653169
HPV105   5904  0.770053475936
HPV20    5500  0.70903699884
```

### Mapping to HPV71

According to the viral mapping, the top hit is HPV71 with 10x coverage. Aligning reads to HPV71 only:

```
samtools faidx /g/data3/gx8/extras/vlad/bcbio/genomes/Hsapiens/GRCh37/viral/gdc-viral.fa HPV71 > viral_mapping/HPV71.fa

bwa index viral_mapping/HPV71.fa

bwa mem -t 9 viral_mapping/HPV71.fa NA12878_WGS.R1.fq NA12878_WGS.R2.fq | bamsort inputthreads=9 outputthreads=9 inputformat=sam index=1 indexfilename=viral_mapping/NA12878_HPV71.bam.bai O=viral_mapping/NA12878_HPV71.bam

samtools idxstats viral_mapping/NA12878_HPV71.bam | py -x "'\t'.join([x.split()[0], x.split()[2], str(int(x.split()[2]) / int(x.split()[1])), str(int(x.split()[3]) / int(x.split()[1]))]) if int(x.split()[1]) > 0 else None" | sort -nr -k3,3 | head | cols
# HPV71  102573  12.7625979843  4.10240139355

samtools stats viral_mapping/NA12878_HPV71.bam > viral_mapping/NA12878_HPV71.stats
SN      bases mapped (cigar):   3815054 # more accurate
SN      bases trimmed:  0
SN      bases duplicated:       0
SN      mismatches:     180824  # from NM fields
SN      error rate:     4.739749e-02    # mismatches / bases mapped (cigar)
SN      average length: 135
SN      maximum length: 151
SN      average quality:        26.4
SN      insert size average:    29.3
SN      insert size standard deviation: 10.4
SN      inward oriented pairs:  14357
SN      outward oriented pairs: 11592
SN      pairs with other orientation:   2819
SN      pairs on different chromosomes: 0
```

### Trying reads with BX tags and high quality

NA12878_WGS read counts:

```
unmapped_or_mate_is_unmapped/NA12878_WGS.bam:  20,876,086 sequences
samtools view -c                               21,750,049  
```

Can we merge all bins from `bxtools split`?

```
cd unmapped_or_mate_is_unmapped
bxtools split NA12878_WGS.bam -a bx_split/NA12878_WGS -m 1000 > bx_split/NA12878_WGS_bx_counts.txt
cd bx_split
parallel -j30 "bedtools bamtofastq -i {} -fq {/.}.R1.fq -fq2 {/.}.R2.fq 2>unpaired/{/.}.unpaired" ::: *-1.bam
cd ..
cat bx_split/*R1.fq > NA12878_WGS.bx_only.R1.fq
cat bx_split/*R2.fq > NA12878_WGS.bx_only.R2.fq
```

Getting total 6_401_804 out of 10_438_043 read pairs - looks quite low.

Aligning to viral fasta:

```
<see commands above>
# Real time: 36.679 sec; CPU: 361.046 sec

name   mapped  coverage         divide by 809965400/1000000000 = 0.8099654 billion reads in genome
HPV71   23552  2.93044668409    3.6179899587
HPV47   10623  1.37496764173
HPV82    7395  0.939644218551
HPV42    4697  0.593280282935
HPV51    4107  0.52599897541
HPV105   1973  0.257336637537
HPV20    1785  0.230114735078
HPV102   1522  0.188412973508
HPV21    1372  0.176372284355
CMV     40204  0.170611849978
```

And to HPV71 only:

```
<see commands above>
# HPV71  30645  3.81298992161
```

Unfortunately, getting a much lower coverage from BX reads (3,8 vs 12,7).

Perhaps we shouldn't be relying on bxtools stats method to get reads with BX. It's weird that summing up values from `bxtools stats` do not match the value above (`6401804 * 2 = 12_803_608` vs `19_613_313` by `cut -f2 NA12878_WGS.bam_bx_stats.txt | awk '{s+=$1} END {print s}'`)

Trying to use pysam to filter reads with BX:

```
mkdir with_bx ; mkdir with_bx_lng_hqual
./filter_bam.py NA12878_WGS.bam

Total:                  21,750,049
bx_i:                   19,613,313
lng_i:                  20,876,724
hqual_i:                 5,639,761
bx_hqual_i:              5,202,715
lng_hqual_i:             4,998,924
bx_lng_hqual_i:          4,632,339
paired_i:               21,750,049
paired_hqual_i:          5,639,761
paired_lng_hqual_i:      4,998,924
paired_bx_hqual_i:       5,202,715
paired_bx_lng_hqual_i:   4,632,339

cd with_bx
samtools fastq -f 1 -F 1536 NA12878_WGS.bam -1 NA12878_WGS.R1.fq -2 NA12878_WGS.R2.fq -s NA12878_WGS.single.fq
NA12878_WGS.R1.fq      - 1,083,568 reads
NA12878_WGS.R2.fq      - 1,083,568 reads
NA12878_WGS.single.fq - 15,865,113 reads

cd with_bx_lng_hqual
samtools fastq -f 1 -F 1536 NA12878_WGS.bam -1 NA12878_WGS.R1.fq -2 NA12878_WGS.R2.fq -s NA12878_WGS.single.fq
NA12878_WGS.R1.fq -         56,907 reads 
NA12878_WGS.R2.fq -         56,907 reads
NA12878_WGS.single.fq -  4,264,669 reads  
```

Trying aligning to viral data again. With BX only:

```
HPV71   3622  0.45066567127        0.5555792877
HPV82   1825  0.231893265565
HPV47   1378  0.178358788506
HPV51    859  0.110015368852
HPV105   706  0.0920829529151
HPV42    684  0.0863963622584
HPV102   604  0.0747709829166
HPV6     308  0.0385192596298
KSHV    5184  0.0375736578507
HPV21    281  0.0361228949736
```

With BX and long and high quality:

```
HPV71   160  0.019907925843        0.0248849073
HPV82   155  0.0196950444727
HPV6    125  0.0156328164082
HPV102   90  0.0111413716266
HPV42    58  0.00732600732601
HPV51    50  0.00640368852459
HPV34    26  0.00336656739609
HPV32    24  0.00301469664615
HPV59    20  0.00253292806484
HPV68    18  0.00230120173869
```

Long and high quality only:

```                                
HPV71   186  0.0231429637925       0.02857156024
HPV82   161  0.020457433291
HPV6    146  0.0182591295648
HPV102  110  0.0136172319881
HPV51    64  0.00819672131148
HPV42    64  0.00808387015284
HPV34    30  0.00388450084164
HPV68    26  0.00332395806699
HPV32    26  0.00326592136666
HPV73    20  0.0025974025974
```

Getting highly mapped reads from the HPV71 bin.

With BX only:

```
HPV71  2529  0.314669652856  
```

With BX and high quality:

```
HPV71  20  0.00248849073037  
```

With BX and long and high quality:

```
HPV71  0         
```

With BX only:

```
HPV71  24  0.00298618887645  
```

## Viruses in COLO829

Merge fastq

```
cd /g/data3/gx8/projects/Saveliev_10X/COLO829-10x/insertions/unmapped_or_mate_is_unmapped/fastq
gunzip -c COLO829-10x_S1_L001_R1_001.fastq.gz >> COLO829-10x_merged.fq
gunzip -c COLO829-10x_S1_L001_R2_001.fastq.gz >> COLO829-10x_merged.fq
  COLO829-10x_merged.fq - 121,428,506 reads
```

Searching viral sequences using 3 different approaches.

- Use Mash
```
mash sketch /g/data3/gx8/extras/vlad/bcbio/genomes/Hsapiens/GRCh37/viral/gdc-viral.fa -i -o mash/gdc-viral
mash screen mash/gdc-viral.msh COLO829-10x_merged.fq | sort -gr > mash/screen_gdc.tab
head mash/screen_gdc.tab
```

Also against RefSeq:

```
mash screen refseq.genomes.k21s1000.msh ../NA12878_WGS.merged.fq | sort -gr > screen_refseq.tab
```

- Align against the viral database:

```
bwa mem -t 9 /g/data3/gx8/extras/vlad/bcbio/genomes/Hsapiens/GRCh37/viral/gdc-viral.fa COLO829.R1.fq COLO829.R2.fq | bamsort inputthreads=9 outputthreads=9 inputformat=sam index=1 indexfilename=viral_mapping/COLO829_viral.bam.bai O=viral_mapping/COLO829_viral.bam

samtools idxstats viral_mapping/COLO829_viral.bam | py -x "'\t'.join([x.split()[0], x.split()[1], x.split()[2], str(int(x.split()[2]) / int(x.split()[1]))]) if int(x.split()[1]) > 0 else None" | sort -nr -k4,4 | head | cols

name   mapped  coverage        divide by 2033578410/1000000000 = 2.03 billion reads in genome
HPV71  268404  33.3960432997   16.451252857
HPV47  170086  22.0147553715
HCV-2  178789  18.4109772423
HCV-1  110412  11.4464026539
HPV82   87949  11.1752223634
HPV42   25999  3.28394593912
HPV20   25253  3.25551115122
HIV-1   28373  3.09040409541
HPV73   20907  2.71519480519
HPV21   20254  2.60367656511
```

- Filtering reads

```
./filter_bam.py COLO829-10x.bam

reads_without_quality:  1,664,954
Total:                128,214,100
bx_i:                 112,376,814
lng_i:                 99,917,997
hqual_i:               83,739,552
lng_hqual_i:           57,912,594
bx_lng_hqual_i:        51,738,038
paired_i:             126,549,146
paired_lng_hqual_i:    57,912,594
paired_bx_lng_hqual_i: 51,738,038

cd with_bx_lng_hqual
samtools fastq -f 1 -F 1536 COLO829-10x.bam -1 COLO829.R1.fq -2 COLO829.R2.fq -s COLO829.single.fq

seq     len  reads  coverage         divide by 2033578410/1000000000 = 2.03 billion reads in genome
HPV71  8037  70862  8.81697150678    4.3356929162
HPV47  7726  30614  3.9624644059
HCV-2  9711  36841  3.7937390588
HPV82  7870  26405  3.35514612452
HCV-1  9646  19666  2.03877254821
HIV-1  9181  12356  1.34582289511
HPV42  7917   9477  1.19704433498
HPV51  7808   7293  0.934042008197
HPV20  7757   6415  0.826994972283
```

- Only HPV71

```
cd /g/data3/gx8/projects/Saveliev_10X/COLO829-10x/insertions/unmapped_or_mate_is_unmapped/with_bx_lng_hqual
samtools faidx /g/data3/gx8/extras/vlad/bcbio/genomes/Hsapiens/GRCh37/viral/gdc-viral.fa HPV71 > viral_mapping/HPV71.fa

bwa index viral_mapping/HPV71.fa

bwa mem -t 10 viral_mapping/HPV71.fa COLO829.R1.fq COLO829.R2.fq | bamsort inputthreads=10 outputthreads=10 inputformat=sam index=1 indexfilename=viral_mapping/COLO829_HPV71.bam.bai O=viral_mapping/COLO829_HPV71.bam

samtools idxstats viral_mapping/COLO829_HPV71.bam | py -x "'\t'.join([x.split()[0], x.split()[2], str(int(x.split()[2]) / int(x.split()[1])), str(int(x.split()[3]) / int(x.split()[1]))]) if int(x.split()[1]) > 0 else None" | sort -nr -k3,3 | head | cols

seq    reads  coverage       divide by 2033578410/1000000000 = 2.03 billion reads in genome
HPV71  91356  11.3669279582  5.5994718966
```