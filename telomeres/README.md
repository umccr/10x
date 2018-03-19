# Brief background on telomeres.

Telomeres are nucleoprotein structures at the ends of all linear
eukaryotic chromosomes. In mammals, telomeres are composed
of the shelterin protein complex and a double-stranded
tract of short tandem repeats of TTAGGG that ends in a singlestranded
3 overhang on the G-strand. Telomeres cap and
protect chromosome ends from eliciting a DNA damage
response and illegitimate recombination events.

As a cell begins to become cancerous, it divides more often and its telomeres are shortened at a higher rate. Cancerous cells escape senescence/death and conversely, they become immortal with the ability to indefinitely replicate, even with shortened telomeres, resulting in uncontrolled tumor growth. It has been demonstrated that this ability to escape the natural cell outcome happens because cancer cells prevent telomeres from becoming critically short by the reactivation of telomerase (phenomenon occurring in 85% of cancers) or by ALT (alternative lengthening of telomeres; 15%).

# Notes

How to find the positions of telomeres? - https://genome.ucsc.edu/FAQ/FAQtracks#tracks20

# Multimapping of telomeric regions

The 10x dataset used (named 10x.bam here) is `/data/cephfs/punim0010/data/External/Reference/NA12878-10x-2018/NA12878_WGS_v2_phased_possorted_bam.bam`:

1. `samtools view -h -L telomere.bed 10x.bam > 10x.telomeres.sam`
2. `samtools view -S -b 10x.telomeres.sam > 10x.telomeres.sam.bam`
3. `samtools bam2fq 10x.telomeres.sam.bam > 10x.telomeres.fq`
4. `bbmap.sh ambiguous=all ref=/data/projects/punim0010/local/share/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa in=10x.telomeres.fq out=10x_telomeres_GRCh37.sam`

Reads Used:9486(1299128 bases)

Mapping:   62.736 seconds.
Reads/sec: 151.20
kBases/sec:20.71

|Read 1 data|pct reads|num reads|pct bases|num bases|
|-----------|---------|---------|---------|---------|
|mapped| 65.4333 |     6207 | 64.4140 |      836820|
|unambiguous| 12.0704 |     1145 | 12.0475 |      156513|
|ambiguous| 53.3629 |     5062 | 52.3664 |      680307|
|low-Q discards|  0.0000 |        0 |  0.0000 |           0|
|perfect best site|  1.8659 |      177 |  1.4290 |       18564|
|semiperfect site|  8.9922 |      853 |  8.1521 |      105906|
|Match Rate|      NA |       NA | 55.4591 |      751906|
|Error Rate| 59.2866 |     5385 | 42.7044 |      578980|
|Sub Rate| 58.5930 |     5322 |  4.1320 |       56021|
|Del Rate| 19.8723 |     1805 | 38.2778 |      518964|
|Ins Rate | 15.0611 |     1368 |  0.2947 |        3995|
|N Rate| 21.8100 |     1981 |  1.8364 |       24898|


5. `grep -oP 'NH:i:\K.*' 10x_telomeres.sam | sort | uniq -c > grep.results`

1190 1

6654 2

1329 3

1280 4

4635 5

# BBMAP's `pileup.sh` coverage

```
$ head COLO829_TGEN_bwa-ready_except_telomere_regions-output.bbmap_coverage
#ID     Avg_fold        Length  Ref_GC  Covered_percent Covered_bases   Plus_reads      Minus_reads     Read_GC Median_fold     Std_Dev
1       78.3707 249250621       0.0000  90.1977 224818427       88049200        87993233        0.4090  65      160.34
2       86.6427 243199373       0.0000  97.8991 238089979       95298296        95600531        0.3980  84      343.57
3       104.8785        198022430       0.0000  98.3519 194758850       93407852        93372869        0.3910  109     63.17
4       100.0256        191154276       0.0000  98.0499 187426510       87635386        87647137        0.3777  102     380.29
5       55.8516 180915260       0.0000  98.1176 177509692       45512490        45505018        0.3902  56      35.78
6       70.0852 171115067       0.0000  97.7692 167297822       53974270        53982658        0.3934  64      275.87
7       114.2343        159138663       0.0000  97.5917 155306175       81960373        81974789        0.4005  112     226.43
8       89.4962 146364022       0.0000  97.5911 142838265       60115957        60130863        0.3826  84      293.91
9       74.3349 141213431       0.0000  85.0412 120089528       47342464        47359013        0.4107  82      199.85


$ head COLO829BL_TGEN_bwa-ready_except_telomere_regions-output.bbmap_coverage
#ID     Avg_fold        Length  Ref_GC  Covered_percent Covered_bases   Plus_reads      Minus_reads     Read_GC Median_fold     Std_Dev
1       80.8254 249250621       0.0000  90.3430 225180539       91288844        91215087        0.4102  82      163.86
2       85.8845 243199373       0.0000  97.8725 238025419       94910965        95200781        0.3968  83      335.98
3       82.4909 198022430       0.0000  98.3520 194758966       73909292        73869763        0.3912  83      48.80
4       85.6878 191154276       0.0000  98.0497 187426168       75340085        75360032        0.3790  82      323.10
5       82.4009 180915260       0.0000  98.1990 177657056       67466386        67447485        0.3892  83      37.43
6       83.2797 171115067       0.0000  97.7710 167300907       64466919        64473102        0.3901  83      273.40
7       82.4560 159138663       0.0000  97.6000 155319315       59537310        59544839        0.3999  82      209.52
8       87.9029 146364022       0.0000  97.5910 142838157       59420561        59416795        0.3812  82      288.08
9       71.8934 141213431       0.0000  85.0513 120103828       46017017        46053224        0.4058  79      201.66
```

# Telomerecat output

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
