# Brief background on telomeres.

Telomeres are nucleoprotein structures at the ends of all linear
eukaryotic chromosomes. In mammals, telomeres are composed
of the shelterin protein complex and a double-stranded
tract of short tandem repeats of TTAGGG that ends in a singlestranded
3 overhang on the G-strand. Telomeres cap and
protect chromosome ends from eliciting a DNA damage
response and illegitimate recombination events.

As a cell begins to become cancerous, it divides more often and its telomeres are shortened at a higher rate. Cancerous cells escape senescence/death and conversely, they become immortal with the ability to indefinitely replicate, even with shortened telomeres, resulting in uncontrolled tumor growth. It has been demonstrated that this ability to escape the natural cell outcome happens because cancer cells prevent telomeres from becoming critically short by the reactivation of telomerase (phenomenon occurring in 85% of cancers) or by ALT (alternative lengthening of telomeres; 15%). 

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
