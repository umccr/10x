## Part 1 - check for existing HLA assemblies

PacBio HLA data is available, but not assembled - https://github.com/PacificBiosciences/DevNet/wiki/HLA-Multiplexed-GenDx-Amplicons-HLA-A-,-B,--C,--DQB1,-and--DRB1

Nanopore data is also available, and not assembled - https://figshare.com/articles/Nanopore_reads_and_alignments/1289717 (ref https://f1000research.com/articles/4-17/v2)

Both of these sets of data could be fed into part 2 (optional).

## Part 2 - generate germline HLA consensus

The HLA region is
6   29719561 32883508 
(as defined by the HLA.hg19 data in GWASTools (https://rdrr.io/bioc/GWASTools/man/HLA.html))

Test data downloaded from 10X (https://support.10xgenomics.com/genome-exome/datasets/2.0.0/HCC1954N_WGS)

bcftools creates a consensus by applying VCF variants
https://samtools.github.io/bcftools/howtos/consensus-sequence.html

Input - 10X BAM
      - reference FASTA
samtools mpileup -uf chr6_10X_reference.fasta chr6_10X.bam | bcftools call -c | vcfutils.pl vcf2fq > chr06_10X-cns.fastq

samtools faidx ref.fa 6:29719561-32883508 | bcftools consensus in.vcf.gz -o out.fa

## Part 2 (optional) - generate de novo germline HLA assembly

First suggestion, ARC, see here (https://github.com/ibest/ARC) and here (http://ibest.github.io/ARC/)

Installed ARC locally - Ubuntu 16 Xenia on Windows 10 (yes it is possible, yes it does work), for sandbox purposes.
Instructions below include some dependencies specific for this new installation of Ubuntu.

ARC needs a mapper (bowtie) and an assembler (spades)

```
sudo apt install bowtie2

wget http://cab.spbu.ru/files/release3.11.1/SPAdes-3.11.1.tar.gz
tar -xvzf SPAdes-3.11.1.tar.gz
cd SPAdes-3.11.1

sudo pip install cmake
wget http://www.zlib.net/zlib-1.2.11.tar.gz
tar -xvzf zlib-1.2.11.tar.gz 
cd zlib-1.2.11
./configure
make
sudo apt install libbz2-dev

wget http://www.bzip.org/1.06/bzip2-1.0.6.tar.gz
tar -xvzf bzip2-1.0.6.tar.gz
cd bzip2-1.0.6/

./spades_compile.sh
./spades.py --test

sudo pip install biopython
git clone git://github.com/ibest/ARC.git
sudo python setup.py install
cd ARC/test_data
./runarc
cat log.txt
```

Installed ARC on Spartan

```
Installed SPAdes (http://spades.bioinf.spbau.ru/release3.11.1/SPAdes-3.11.1-Linux.tar.gz)
added SPAdes to PATH
module load Bowtie2
module load Biopthon
cd /Opt/ARC/test_data
./runarc
cat log.txt
```

ARC needs a target/set of targets to map reads against for assembly.  
1. Begin with the HLA region as defined above

No joy, didn't complete.  In hindsight ARC was probably not a smart choice, begins by [aligning the reads](https://ibest.github.io/ARC/#ARCAlgorithm) but has no ability to take advantage of the read cloud attributes of 10X data, so EMA or LongRanger would generate superior alignments.  Supernova would probably be a better option, but the resources were not available.

## Part 3 - call somatic variants against new HLA reference

Using variants called on Raijin (COLO829 10X dilution series run with EMA, variants called with bcbio)
On Raijin, 10X data is here - /g/data/gx8/data/10X/10X_EMA

## Part 4 - compare standard somatic calls vs 10X calls vs new reference calls

Compare COLO829 bwa bcbio variants vs COLO829 EMA bcbio variants
