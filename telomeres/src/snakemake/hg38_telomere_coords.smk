rule telomere_coords:
    output: "GRC38h_telomeric_regions.bed"
    shell:
        "src/10x_telomeres/telomere_coords.py data/external/hg38.fa.gz"

# XXX: Tweak telomere_coords.py accordingly to output a well formed bedfile