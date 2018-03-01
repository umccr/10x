## Part 1 - check for existing HLA assemblies

PacBio HLA data is available, but not assembled - https://github.com/PacificBiosciences/DevNet/wiki/HLA-Multiplexed-GenDx-Amplicons-HLA-A-,-B,--C,--DQB1,-and--DRB1
Nanopore data is also available, and not assembled - https://figshare.com/articles/Nanopore_reads_and_alignments/1289717 (ref https://f1000research.com/articles/4-17/v2)

## Part 2 - generate germline HLA consensus

The HLA region is
6   29719561 32883508 
(as defined by the HLA.hg19 data in GWASTools (https://rdrr.io/bioc/GWASTools/man/HLA.html))

https://samtools.github.io/bcftools/howtos/consensus-sequence.html

Input - 10X BAM
      - reference FASTA
samtools mpileup -uf chr6_10X_reference.fasta chr6_10X.bam | bcftools call -c | vcfutils.pl vcf2fq > chr06_10X-cns.fastq

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

## Part 3 - call somatic variants against new HLA reference

## Part 4 - compare standard somatic calls vs 10X calls vs new reference calls
