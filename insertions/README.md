See [notebook](10x_novel_insertions.ipynb)

Misc commands:

```
bedtools bamtofastq -i COLO829-10x.sorted.bam -fq COLO829-10x.R1.fq -fq2 COLO829-10x.R2.fq 2>COLO829-10x.unmpaired

bxtools split NA12878_WGS.sorted.bam -a NA12878_WGS -m 100 | sort -k 2 -r -n > NA12878_WGS_bx_counts.txt

sambamba view -f bam -F "unmapped" -t 10 NA12878_WGS.sorted.bam > ../unmapped/NA12878_WGS.sorted.bam

bxtools split COLO829-10x_phased_possorted_bam.bam -a bx_split/COLO829-10x -m 1000 > bx_split/COLO829-10x_bx_counts.txt
```

In `/data/cephfs/punim0010/projects/Saveliev_10X/NA12878-10x/insertions/unmapped_or_mate_is_unmapped/NA12878_WGS_bx_counts.txt`: 
Barcodes: 1,460,274
Reads: 19,613,313

EMA count command on unmapped_or_mate_unmapped reads:
R1: :: Reads with OK barcode: 12,195 out of 5,232,233
R2: :: Reads with OK barcode: 23,052 out of 5,232,233


TODO: check why barcodes are bad? Howe to extract on good barcodes for denovo assembly?