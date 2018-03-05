# Other stuff
# Telomerecat output

Working off [readthedocs](http://telomerecat.readthedocs.io/en/latest/understanding_output.html).


Sample|F1|F2|F4|Psi|Insert_mean|Insert_sd|Read_length|Initial_read_length|F2a|F2a_c|Length
------------------------------------------------------------------------------------------
COLO829_TGEN_bwa-ready.bam|38605|8012|2528|2.573|397.0|96.979|112|112|5484|5484|2403.0
 | \# of reads which are completely telomeric | \# of reads where one end is completely telomeric, and the other is not, and telomeric end is CCCTAA | as previous column, and end is TTAGGG | Measure of fidelity | Insert size | SD of insert size | F2 - F4
 

 F2-F4 = F2a, the estimated number of reads covering the boundary between telomere and nontelomere.
 
  The greater the measure of fidelity, the more we believe the observed measurement of F2a.
 
 F2a value after undergoing correction (see paper).
 
 Length - telomere length as estimated by telomerecat.
