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

rule samtools_view:
    input:
        "../../bams/hg38/COLO829BL-ready.bam"
    output:
        "../../bams/hg38/COLO829BL-ready-telomeres.bam"
    params:
        "-b",
        "-L ../../bed/hg38_noalt.bed"
#        samtools view {sample} -L "../../bed/hg38_noalt.bed" -b > "final/hg38/COLO829BL-hg38-ready-without_randoms.bam"
#        samtools view -H final/hg38/COLO829BL-hg38-ready-without_randoms.bam | grep -v GL | grep -v KI > COLO829BL-hg38-ready-without_randoms-torehead.sam
#        samtools reheader COLO829BL-hg38-ready-without_randoms-torehead.sam COLO829BL-hg38-ready-without_randoms.bam > COLO829BL-hg38-ready-without_randoms-reheaded.bam
#        samtools view COLO829BL-hg38-ready-without_randoms-reheaded.bam | cut -f 3 | sort | uniq -c | sort -rn
    wrapper:
       "0.23.1/bio/samtools/view"