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

# Telomerecat output

Working off [readthedocs](http://telomerecat.readthedocs.io/en/latest/understanding_output.html).


Sample|F1|F2|F4|Psi|Insert_mean|Insert_sd|Read_length|Initial_read_length|F2a|F2a_c|Length
------|--|--|--|---|-----------|---------|-----------|-------------------|---|-----|-------
COLO829_TGEN_bwa-ready.bam|38605|8012|2528|2.573|397.0|96.979|112|112|5484|5484|2403.0
COLO829BL_TGEN_bwa-ready.bam|49507|7545|3066|2.9960000000000004|346.0|84.266|112|112|4479|4479|2928.2
sample name | \# of reads which are completely telomeric | \# of reads where one end is completely telomeric, and the other is not, and telomeric end is CCCTAA | as previous column, and end is TTAGGG | Measure of fidelity | Insert size | SD of insert size | | | F2 - F4 

F2-F4 = F2a, the estimated number of reads covering the boundary between telomere and nontelomere.
 
The greater the measure of fidelity, the more we believe the observed measurement of F2a. (Question - so how large can this value be?, how large is this value with a clean sample?)
 
F2a value after undergoing correction (see paper).
 
Length - telomere length as estimated by telomerecat.
