10X Notes
=========

<!-- vim-markdown-toc GFM -->

* [Chromium Technology](#chromium-technology)
    * [Example Calculations](#example-calculations)
* [Long Ranger](#long-ranger)
* [SV Dot Plots in R](#sv-dot-plots-in-r)
* [DNA Molecule Length Distribution](#dna-molecule-length-distribution)
* [Loupe](#loupe)
* [IGV](#igv)

<!-- vim-markdown-toc -->

Chromium Technology
-------------------
* Basic Idea
    - Take a long DNA fragment, chop it into smaller reads and tag each read
      with same molecular barcode. Reads that share the same barcode
      can be grouped as coming from a single long input DNA chunk.
* Load following into microfluidic chip:
    - gel beads coated with millions oligo primers
      `<R1>|<16bp Barcode>|<random 6mer>`
          - 4.4 million different gel beads, each with different
            16bp barcode unique to that bead
            (different colours = different barcodes)
    - HMW (high molecular weight) genomic DNA + enzyme
    - paritioning oil
* Mixture is partitioned into aqueous droplets called GEMs (Gel beads in EMulsion)
* Each GEM contains a gel bead of a different 'colour' representing barcode unique
  to that GEM.
* After the GEMs are formed, a reagent in the GEM dissolves the gel bead, releasing
  its primers into the GEM solution
* In isothermal incubation, primers bind to HMW gDNA and extend to form initial
  portion of library construct
* All resulting amplicons within GEM will contain same unique barcode
* GEMs are then broken: oil is removed, and all different
  amplicons are pooled together
* Because of the limiting amount of input DNA, the chances that two fragments from
  the same genomic region end up in the same GEM is very small (1/7000 -- see below)
* 99.98% unique barcodes at any genomic locus

### Example Calculations
- See Genome Sequencing Chapters 4 and 10 at
  [this link](https://www.10xgenomics.com/10x-university/).
- Input: HMW gDNA with size of 50,000 bases or longer
- 1.25 ng of gDNA into chip (nano = 1 billionth of a gram)
- 1 human genome copy is 3.5 picograms (pico =  1 trillionth of a gram)
- So 1.25ng / 3.5pg = (1.25 * 10^-9 g) / (3.5 * 10^-12 g) =
  1.25 / (3.5 * 10^-3) = 0.35 * 10^3 = 350 genome (haploid) copies
- But 40% capture efficiency: 1.25ng * 0.4 = 0.5ng input DNA / 3.5pg =
  0.5 * 10-9 / 3.5 * 10-12 = 0.14 * 1000 = 150 genome (haploid) copies
- Let's say 1 million GEMs generated, so the 150 genome copies will
  have been partitioned into the 1 million GEMS = 0.00015 (1/7000)
  of genome in each GEM
- 3.2 Gigabases of human genome, so 3,200,000,000 / 7,000 = 500,000 bases
  per GEM
- Since each input DNA molecule is 50,000 bases in length, that gives
  us around 10 HMW gDNA molecules per GEM
- If each input DNA molecule is say 100,000 bases in length, then you'll
  get less (e.g. 5) molecules per GEM
- subset of barcoded template is sequenced (sample more DNA molecules
  instead of more reads per DNA molecule)
- at 30X read coverage, you'll get 35 reads sequenced for each 50kb HMW molecule
- so 35 * 2 PE reads * 150b read length = 10,000b of sequence info per HMW molecule
- so 10kb/50kb = 0.2 read coverage per HMW molecule
- physical coverage: number of gDNA molecules that go into the system and that
  cover a given genome locus
- read coverage: number of sequencing reads generated during Illumina sequencing
  and then mapped to a given genome locus

Long Ranger
-----------
* [10X docs](https://support.10xgenomics.com/genome-exome/software/pipelines/latest/what-is-long-ranger)

* Pipelines:
    - `mkfastq`: basically bcl2fastq
    - `wgs`: take FASTQs and align (Lariat), de-duplicate, filter,
      call and phase SNPs, indels and SVs
    - `targeted`: as above, but for WES. Reports stats on targeted regions in BED
    - `basic`: barcode processing, produces unaligned BAM/FASTQ file. Includes
      read trimming, one-base error correction, whitelist filtering, attachment of
      barcodes to reads. TR/TQ tags preserve trimmed sequences to support lossless
      BAM to FASTQ conversion.
    - `align`: align with BWA

* System requirements:
    - local mode: 16 cores, 96GB RAM, 2TB disk
    - cluster mode: run on SGE/LSF, 8 cores per node, 6GB RAM per core

* BAM tags:
    - BX: barcode sequence after error-correction and whitelist filtering.
    - MI: molecular identifier. Reads with same MI tag are linked to same molecule.
        * If you have 10 molecules in a GEM, then for each GEM you'd expect to see
          one BX tag, and 10 MI tags associated with that BX tag.
    - TR/TQ: sequence and qualities of 7 trimmed bases following barcode at
      the start of R1. Can be used to go back to original R1.
    - PC: Phred-scaled confidence that the read was phased correctly
    - PS: Phase set containing the read
    - HP: Haplotype of molecule that generated the read


SV Dot Plots in R
-----------------
* [10X forum link](https://community.10xgenomics.com/t5/Genome-Exome-Forum/Loupe-Style-Dot-Plots-in-R/m-p/375#M57)
* Check out 'Export to CSV' on SV tab in Loupe
* Some R code:

```r
library(gplots)

# load data from csv, convert to matrix, rotate/transform to Loupe orientation
df <- read.csv("/your/path/to/loupe/output/loupe-sv-barcode-matrix.csv")
m <- data.matrix(df)
rtm <- apply(m, 2, rev)

# define custom colors and color breaks
custcol <- colorRampPalette(c("white", "yellow", "red", "black"))(n = 23)
colbreaks = c(seq(0,2,length=6),
  seq(3,5,length=6),
  seq(6,15,length=6),
  seq(16,100,length=6))

# plot matrix data
heatmap.2(rtm, labRow = FALSE, labCol = FALSE, 
          density.info = "none", key = FALSE, trace = "none",
          col = custcol, breaks = colbreaks, dendrogram = "none",
          Rowv ="NA", Colv="NA", lmat = rbind(c(0, 0), c(0, 1)),
          lwid = c(0.5, 4), lhei = c(.5, 4))
```

DNA Molecule Length Distribution
--------------------------------
* [10X forum link](https://community.10xgenomics.com/t5/Genome-Exome-Forum/Length-distribution-of-input-DNA-molecule/m-p/265#M50)
* Start out with `>200kb` on a pulsed-field gel
* After running Chromium Genome assay: 
    - length-weighted mean is `~90kb`
    - mode/peak is `~55kb`
    - right-skewed

Loupe
-----

* [Application Examples](https://community.10xgenomics.com/t5/Genome-Exome-Forum/Application-examples-using-Loupe-to-visualize-and-interpret/m-p/34#M46)

IGV
---
* [10X blog video](https://community.10xgenomics.com/t5/10x-Blog/10x-pert-Workshop-Structural-Variant-and-Haplotype-Analysis-with/ba-p/563)
* If you want to see small variants and how they link to each other, use IGV
* Two versions (beta 2.4.0 [released in September 2017](https://software.broadinstitute.org/software/igv/ReleaseNotes/2.4.x)):
    - beta: supports linked reads. Right click on read groups, then visualise by
      linked reads/barcode. Thin line represents barcode, with reads aligned to it.
    - old: colour by tag (HP or haplotype tag) and it will show how they are phased
      relatively to each other (blue and pink, vs grey). Can also group by tag.
