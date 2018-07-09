Having a mutect error when running on EMA bam files.

```[2018-03-27T06:57Z] spartan-bm034.hpc.unimelb.edu.au: Ignoring SAM validation error: ERROR: Record 1761761, Read name A00130:39:H5GWHDMXX:2:2256:18005:32236, Zero-length read without FZ, CS or CQ tag```

This error is preceded by tons of identical message.

Probably means that zero alignments somehow got extracted into empty reads (was realigning LongRanger's BAMs with BWA), and Mutect doesn't like those, not sure if that's related to the memory issue, again.


Exploring this read:


```
samtools view Colo829_10x_80pc_BWA.bam | grep A00130:39:H5GWHDMXX:1:1235:11026:33739
A00130:39:H5GWHDMXX:1:1235:11026:33739  69      4       49275711        0       *           =       49275711        0       *                                                                                                                                                           *                                                                                                                                                                               AS:i:0  XS:i:0  RG:Z:Colo829_80pc_BWA   MQ:i:0  ms:i:2536       mc:i:49275695   MC:Z:16S52M83S
A00130:39:H5GWHDMXX:1:1235:11026:33739  1161    4       49275711        0       16S52M83S   =       49275711        0       GGAGGCAAAAAAAAAAAAAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAACAAAAAAAAAAAAAAAAAATAAGTATGGATGATTGGTTACGCGATGGCGTGAGTAGGTGGGCGGGAGGAGTGGTGGTGAGGGAGGGTGTCGGGGGGGTTG     :,,,,,,,:,F:,:FF,FFF,,F:,FFFFF,:F,:F:F,FFF,FFFFFFFF,FFFF,FF:::F,,FFF::,F,,F::,,,,,,,,,,,,,:,,F,,F,,,,,,,F,F:FFF,,,:,,,:,,,:,,FF,,,FF,:F,,FF,:,,,:,,,,,F NM:i:1  MD:Z:17A34      AS:i:47 XS:i:46       RG:Z:Colo829_80pc_BWA   SA:Z:21,41612277,-,121S30M,0,0; ms:i:0  mc:i:49275711
A00130:39:H5GWHDMXX:1:1235:11026:33739  441     21      41612277        0       121H30M     =       41612277        0       TTTTTTTTTTTTTTTTTTTTTTTTGCCTCC                                                                                                                              FFFFF,:F,,FFF,FF:,:F,:,,,,,,,:                                                                                                                          NM:i:0  MD:Z:30 AS:i:30 XS:i:30 RG:Z:Colo829_80pc_BWASA:Z:4,49275711,+,16S52M83S,0,1; XA:Z:7,+19341757,30M121S,0;17,-986503,121S30M,0;3,+42546389,30M121S,0;9,+126560999,30M121S,0;

samtools view Colo829_10x_80pc_LongRanger.bam | grep A00130:39:H5GWHDMXX:1:1235:11026:33739
A00130:39:H5GWHDMXX:1:1235:11026:33739  137     chr4    49275711        24      16S52M83S   *       0               0       GGAGGCAAAAAAAAAAAAAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAACAAAAAAAAAAAAAAAAAATAAGTATGGATGATTGGTTACGCGATGGCGTGAGTAGGTGGGCGGGAGGAGTGGTGGTGAGGGAGGGTGTCGGGGGGGTTG     :,,,,,,,:,F:,:FF,FFF,,F:,FFFFF,:F,:F:F,FFF,FFFFFFFF,FFFF,FF:::F,,FFF::,F,,F::,,,,,,,,,,,,,:,,F,,F,,,,,,,F,F:FFF,,,:,,,:,,,:,,FF,,,FF,:F,,FF,:,,,:,,,,,F QT:Z:,:FF,:F:   BC:Z:AGTAGTCT   QX:Z:########   AM:A:0  XM:A:0  RX:Z:NNNNNNNN   AS:f:-65.5      XS:f:-66        XT:i:1  OM:i:5     RG:Z:10X-COLO829:LibraryNotSpecified:1:unknown_fc:0
```


```
00130:39:H5GWHDMXX:1:2131:22046:27242  99      HPV71   94      51      23M104S =       96      23
ACAAAGCAATGCTTGGTTCAGCAATGCTCCATCGAGCAATGCTTGCGATAGCAATGCTCAGTGTAGCAATGCTAGCTTGAGCAATGCTCCACAAAGCACTGCTTGGTTCAGCAATGCTCCATCGAGC
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFF:FFFFFFFFFFFFFF,FFFFFFFFFFFFFF,FFFFFFF,FFFFFFFFFFFFFFF:,FFF
NM:i:1  MD:Z:1T21       MC:Z:94S21M36S  AS:i:21 XS:i:0

A00130:39:H5GWHDMXX:1:2131:22046:27242  125     *       0       0       32S19M76S       chr8    117107344       0
GCTCGATGGAGCATTGCTGAACCAAGCAGTGCTTTGTGGAGCATTGCTCAAGCTAGCATTGCTACACTGAGCATTGCTATCGCAAGCATTGCTCGATGGAGCATTGCTGAACCAAGCATTGCTTTGT
FFF,:FFFFFFFFFFFFFFF,FFFFFFF,FFFFFFFFFFFFFF,FFFFFFFFFFFFFF:FFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
RX:Z:AACCCTCAGCATTGCT   QX:Z:FFFFFFFFFFFFFFFF   TR:Z:CCTCGTT    TQ:Z:FF:FFFF    BC:Z:CAGTACTGQT:Z:F:FFFFFF    XS:f:-138.5     XC:Z:   AC:Z:   AS:f:-138.5     XM:A:0  AM:A:0  XT:i:0  BX:Z:AACCCTCAGCATTGCT-1 RG:Z:10X-COLO829:LibraryNotSpecified:1:unknown_fc:0

A00130:39:H5GWHDMXX:1:2131:22046:27242  189     *       0       0       20M131S *       0       0
CCACCAAGCACTGCTTGGTTCAGCAAGGCTCCATCGAGCACTGCTGGCGAGAGCAATGCTCAGTGTAGCAATGCGAGCTTGAGCAATGCTCCACAAAGCAATGCTTGGTTCAGCAATGCTCCATCGAGCAATGCTATCGCAAGCATTGCTC
F:::,:FFF,,FF:FFFF,F:,F:F,,F::FF:,FFF::F,FFFF,FFFF,F:FFFF:FFF:,,:FFFF:,::F,:F:F,::FFF,FFFFFF:FFF::FFF:FFFFFFF:F:FFFF,FFFFFFFFF,FFFFFFFF:FFFF:FFFF:,FFFF
RX:Z:AACCCTCAGCATTGCT   QX:Z:FFFFFFFFFFFFFFFF   BC:Z:CAGTACTGQT:Z:F:FFFFFF    XS:f:-140       XC:Z:   AC:Z:   AS:f:-138.5     XM:A:0  AM:A:0  XT:i:0  SA:Z:chrX,89395382,+,54H19M78H,19,0;    BX:Z:AACCCTCAGCATTGCT-1 RG:Z:10X-COLO829:LibraryNotSpecified:1:unknown_fc:0
```


```
Read name = A00130:36:H32JNDSXX:2:2422:21097:19742
Read length = 151bp
----------------------
Mapping = Primary @ MAPQ 60
Reference span = NODE_1_length_4594_cov_1709.398716:1-151 (+) = 151bp
Cigar = 151M
Clipping = None
----------------------
Mate is mapped = yes
Mate start = NODE_1_length_4594_cov_1709.398716:608 (-)
Insert size = 759
First in pair
Pair orientation = F1R2
----------------------
MC = 151M
NM = 0
AS = 151
XS = 0
Hidden tags: MD<hr>Location = NODE_1_length_4594_cov_1709.398716:62
Base = T @ QV 37

Alignment start position = NODE_1_length_4594_cov_1709.398716:1
TTGGGGAGGGATGTGGAAAGCATAGACAGAAGGCACTGCAGAAGAGGAGAGTAAATGACTTTAGGTCATTGAACTGGTTTTTGATAGAAATCAAGTTAAGAAAAATAGAAGTCAATAGTACTAAGTTTGAACTTTATCTTTCTAAAAACTG
```

```
Read name = A00130:36:H32JNDSXX:2:1446:9426:16971
Read length = 151bp
----------------------
Mapping = Primary @ MAPQ 60
Reference span = NODE_1_length_4594_cov_1709.398716:78-228 (+) = 151bp
Cigar = 151M
Clipping = None
----------------------
Mate is mapped = yes
Mate start = NODE_1_length_4594_cov_1709.398716:651 (-)
Insert size = 724
First in pair
Pair orientation = F1R2
----------------------
MC = 150M
NM = 0
AS = 151
XS = 0
Hidden tags: MDLocation = NODE_1_length_4594_cov_1709.398716:159
Base = G @ QV 37
```