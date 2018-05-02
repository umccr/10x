all: igv bedfile 
igv:
	bedtools igv -slop -400 -path ${PWD}/igv -i GRC38h_real_coords.bed > telomere_igv_screenshots.batch
bedfile:
	tail -48 latest_output.txt > GRC38h_real_coords.bed
