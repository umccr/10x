rule telomere_coords:
    output: "data/processed/telomere_coords.bed"
    params:
        genome_path = "data/external/hg38.fa.gz"
    shell:
        "src/10x_telomeres/telomere_coords.py {params.genome_path} > logs/latest_output.txt"