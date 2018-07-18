## HapCUT2

```
cd /home/563/vs2870/gx8/data/10X/Phasing/test_small
```

### Pipeline

Script `run.sh` will phase a 10x VCF file.

It contains of several HAPCUT2 steps according to https://github.com/vibansal/HapCUT2#10x-genomics-linked-reads
Plus a fgbio tool to convert HAPCUT2 output back into a phased VCF: https://github.com/fulcrumgenomics/fgbio#list-of-tools -> HapCutToVcf

### Testing on chr21

Subsetting a 100pc EMA BAM file

```
ln -s /home/563/vs2870/gx8/data/10X/10X_EMA/bcbio_original/final/2018-07-14_10X-Collaboration_EMA/COLO829_100pc-batch-ensemble-annotated.vcf.gz EMA_COLO829_100pc-ensemble.vcf.gz

# Uncompressed VCF file is required. Tricks like <(gunzip -c ...vcf.gz) or <(bcftools view ...vcf.gz) for some reason won't work
bcftools view -f.,PASS EMA_COLO829_100pc-ensemble.vcf.gz 21 > EMA_COLO829_100pc-ensemble.21.vcf

ln -s /home/563/vs2870/gx8/projects/Saveliev_10X/COLO829-10x/ema_trimmed_runs/ema_COLO829_100pc/COLO829_100pc.bam
samtools view -O BAM COLO829_100pc.bam 21 > COLO829_100pc.21.bam
samtools index COLO829_100pc.21.bam
```

Running

```
bash run.sh results/EMA_ EMA_COLO829_100pc-ensemble.21.vcf EMA_COLO829_100pc.21.bam
```

Does't work on ensemble calls. Trying strelka2:

```
ln -s /home/563/vs2870/gx8/data/10X/10X_EMA/bcbio_original/final/2018-07-14_10X-Collaboration_EMA/COLO829_100pc-batch-strelka2-annotated.vcf.gz EMA_COLO829_100pc-strelka2.vcf.gz
bcftools view -f.,PASS EMA_COLO829_100pc-strelka2.vcf.gz 21 > TruSeq_COLO829_100pc-strelka2.21.vcf

bash run.sh results/EMA_ EMA_COLO829_100pc-strelka2.21.vcf EMA_COLO829_100pc.21.bam
```

Works good.

Now comparing against a non-10x run. Subsetting:

```
ln -s /home/563/vs2870/gx8/data/10X/TruSeq/bcbio_original/final/2018-07-15_10X-Collaboration_WGS-merged/COLO829_100pc-batch-strelka2-annotated.vcf.gz TruSeq_COLO829_100pc-strelka2.vcf.gz
bcftools view -f.,PASS TruSeq_COLO829_100pc-strelka2.vcf.gz 21 > TruSeq_COLO829_100pc-strelka2.21.vcf

ln -s /g/data3/gx8/data/10X/TruSeq/bcbio_original/BAMs/Colo829-ready.bam
samtools view Colo829-ready.bam 21 -O BAM > Colo829-ready.21.bam
samtools index Colo829-ready.21.bam
```

Running sligrly different workflow `run_TruSeq.sh`:

```
bash run_TruSeq.sh results/TruSeq_ TruSeq_COLO829_100pc-strelka2.21.vcf TruSeq_COLO829_100pc.21.bam
```

Doesn't work on TruSeq calls. Maybe need some clean up, like ensemble. Trying on EMA calls:

```
bash run_TruSeq.sh results/No10x_ EMA_COLO829_100pc-strelka2.21.vcf EMA_COLO829_100pc.21.bam
```

All good.

Without using linked reads, the output looks like this, e.g. many discontiguous haplotype blocks, based on paired read evidence. 1221 blocks for one 24M-long chr14 chunck:

```
BLOCK: offset: 1 len: 3 phased: 3 SPAN: 112 fragments 40
1                                                               0  1  chr14  19057082  A                          G                          1/0:51:37:14,37:0.7255:6,8:20,17    0  .  100.00
2                                                               1  0  chr14  19057114  C                          T                          0/1:39:7:31,7:0.1795:13,18:5,2      0  .  100.00
3                                                               1  0  chr14  19057194  G                          A                          0/1:43:14:29,14:0.3256:16,13:9,5    0  .  28.98
********
BLOCK: offset: 4 len: 2 phased: 2 SPAN: 106 fragments 11
4                                                               0  1  chr14  19057714  T                          G                          1/0:46:28:18,28:0.6087:9,9:16,12    0  .  29.90
5                                                               0  1  chr14  19057820  A                          G                          0/1:38:7:29,7:0.1842:13,16:3,4      0  .  29.90
********
BLOCK: offset: 6 len: 10 phased: 8 SPAN: 2110 fragments 53
6                                                               0  1  chr14  19058650  G                          A                          0/1:52:13:39,13:0.25:21,18:4,9      0  .  100.00
7                                                               0  1  chr14  19058697  A                          G                          0/1:38:8:29,8:0.2105:15,14:2,6      0  .  100.00
8                                                               1  0  chr14  19059213  A                          T                          1/0:53:41:11,41:0.7736:3,8:15,26    0  .  100.00
10                                                              0  1  chr14  19059626  C                          T                          0/1:51:13:38,13:0.2549:16,22:5,8    0  .  100.00
11                                                              1  0  chr14  19059977  T                          G                          1/0:32:21:10,21:0.6563:7,3:12,9     0  .  100.00
13                                                              1  0  chr14  19060491  C                          T                          0/1:24:5:19,5:0.2083:5,14:1,4       0  .  100.00
14                                                              1  0  chr14  19060539  A                          T                          0/1:38:10:28,10:0.2632:15,13:7,3    0  .  100.00
15                                                              1  0  chr14  19060760  G                          A                          0/1:30:5:25,5:0.1667:10,15:3,2      0  .  100.00
********
BLOCK: offset: 16 len: 43 phased: 38 SPAN: 5367 fragments 342
16                                                              0  1  chr14  19061333  A                          G                          0/1:44:19:24,19:0.4318:9,15:14,5    0  .  20.42
18                                                              0  1  chr14  19061666  T                          C                          0/1:33:6:27,6:0.1818:14,13:1,5      0  .  100.00
19                                                              0  1  chr14  19061858  A                          G                          0/1:31:4:27,4:0.129:6,21:4,0        0  .  100.00
20                                                              0  1  chr14  19061860  T                          A                          0/1:31:4:27,4:0.129:6,21:4,0        0  .  100.00
21                                                              0  1  chr14  19061985  C                          T                          0/1:30:14:16,14:0.4667:6,10:11,3    0  .  100.00
22                                                              0  1  chr14  19062052  G                          A                          0/1:34:7:26,7:0.2059:18,8:6,1       0  .  100.00
23                                                              1  0  chr14  19062096  T                          A                          0/1:34:11:23,11:0.3235:11,12:7,4    0  .  100.00
24                                                              0  1  chr14  19062363  G                          A                          0/1:31:13:18,13:0.4194:7,11:5,8     0  .  100.00
26                                                              0  1  chr14  19062531  C                          G                          0/1:30:10:20,10:0.3333:10,10:5,5    0  .  100.00
27                                                              0  1  chr14  19062640  A                          G                          0/1:34:7:26,7:0.2059:14,12:5,2      0  .  100.00
28                                                              0  1  chr14  19062945  G                          A                          0/1:25:3:21,3:0.12:7,14:1,2         0  .  100.00
...
```

However, when adding the linked reads information, getting one single haplotype block for that chunk:

```
BLOCK: offset: 1 len: 17194 phased: 11830 SPAN: 4863199 fragments 26060
1                                                                        1  0  chr14  19057082  A                          G                          1/0:51:37:14,37:0.7255:6,8:20,17    0  .  100.00
2                                                                        1  0  chr14  19057114  C                          T                          0/1:39:7:31,7:0.1795:13,18:5,2      0  .  100.00
3                                                                        1  0  chr14  19057194  G                          A                          0/1:43:14:29,14:0.3256:16,13:9,5    0  .  100.00
4                                                                        1  0  chr14  19057714  T                          G                          1/0:46:28:18,28:0.6087:9,9:16,12    0  .  100.00
5                                                                        1  0  chr14  19057820  A                          G                          0/1:38:7:29,7:0.1842:13,16:3,4      0  .  100.00
6                                                                        1  0  chr14  19058650  G                          A                          0/1:52:13:39,13:0.25:21,18:4,9      0  .  100.00
7                                                                        1  0  chr14  19058697  A                          G                          0/1:38:8:29,8:0.2105:15,14:2,6      0  .  100.00
8                                                                        1  0  chr14  19059213  A                          T                          1/0:53:41:11,41:0.7736:3,8:15,26    0  .  100.00
10                                                                       1  0  chr14  19059626  C                          T                          0/1:51:13:38,13:0.2549:16,22:5,8    0  .  100.00
11                                                                       1  0  chr14  19059977  T                          G                          1/0:32:21:10,21:0.6563:7,3:12,9     0  .  100.00
13                                                                       1  0  chr14  19060491  C                          T                          0/1:24:5:19,5:0.2083:5,14:1,4       0  .  100.00
14                                                                       1  0  chr14  19060539  A                          T                          0/1:38:10:28,10:0.2632:15,13:7,3    0  .  100.00
15                                                                       1  0  chr14  19060760  G                          A                          0/1:30:5:25,5:0.1667:10,15:3,2      0  .  100.00
16                                                                       1  0  chr14  19061333  A                          G                          0/1:44:19:24,19:0.4318:9,15:14,5    0  .  100.00
18                                                                       1  0  chr14  19061666  T                          C                          0/1:33:6:27,6:0.1818:14,13:1,5      0  .  100.00
19                                                                       1  0  chr14  19061858  A                          G                          0/1:31:4:27,4:0.129:6,21:4,0        0  .  100.00
20                                                                       1  0  chr14  19061860  T                          A                          0/1:31:4:27,4:0.129:6,21:4,0        0  .  100.00
21                                                                       1  0  chr14  19061985  C                          T                          0/1:30:14:16,14:0.4667:6,10:11,3    0  .  100.00
22                                                                       1  0  chr14  19062052  G                          A                          0/1:34:7:26,7:0.2059:18,8:6,1       0  .  100.00
23                                                                       0  1  chr14  19062096  T                          A                          0/1:34:11:23,11:0.3235:11,12:7,4    0  .  100.00
24                                                                       1  0  chr14  19062363  G                          A                          0/1:31:13:18,13:0.4194:7,11:5,8     0  .  100.00
26                                                                       1  0  chr14  19062531  C                          G                          0/1:30:10:20,10:0.3333:10,10:5,5    0  .  100.00
27                                                                       1  0  chr14  19062640  A                          G                          0/1:34:7:26,7:0.2059:14,12:5,2      0  .  100.00
28                                                                       1  0  chr14  19062945  G                          A                          0/1:25:3:21,3:0.12:7,14:1,2         0  .  100.00
```
