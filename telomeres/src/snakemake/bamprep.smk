# XXX: Pending pullrequests:
# https://bitbucket.org/brainstorm/snakemake-wrappers/commits/d7936a42db7139f3ec41fc25100abf9d66b0e2a8
# https://github.com/bioconda/bioconda-recipes/pull/8787
#rule sambamba_slice:
#    input:
#        "../../bams/hg38/COLO829BL-ready.bam"
#    output:
#        "../../bams/hg38/COLO829BL-ready-hg38_telomeric_regions.bam"
#    params:
#        "-L ../../bed/hg38_noalt.bed"
#    wrapper:
#       "0.23.1/bio/sambamba/slice"

rule count_multimap:
    input: data/processed/{sample}.bam
    shell:
        samtools view data/processed/{sample}.bam | cut -f 3 | sort | uniq -c | sort -rn