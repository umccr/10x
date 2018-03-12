10X Notes
=========

<!-- vim-markdown-toc GFM -->

* [Intro to Linked Reads](#intro-to-linked-reads)
* [SV Dot Plots in R](#sv-dot-plots-in-r)
* [DNA Molecule Length Distribution](#dna-molecule-length-distribution)
* [Loupe](#loupe)
* [IGV](#igv)

<!-- vim-markdown-toc -->

Intro to Linked Reads
---------------------
* [10X blog link](https://community.10xgenomics.com/t5/10x-Blog/A-basic-introduction-to-linked-reads/ba-p/95)
* Take a long DNA fragment, chop it into smaller reads and tag each read
  with same molecular barcode.

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
* Two versions:
    - beta: supports linked reads. Right click on read groups, then visualise by
      linked reads/barcode. Thin line represents barcode, with reads aligned to it.
    - old: colour by tag (HP or haplotype tag) and it will show how they are phased
      relatively to each other (blue and pink, vs grey). Can also group by tag.
