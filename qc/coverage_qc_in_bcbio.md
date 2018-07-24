### Running coverage QC in bcbio to compare 10x aligners

1. `cd /data/cephfs/punim0010/projects/Saveliev_10X/NA12878-10x/bcbio_grch37`

2. Put [QC config](bcbio_na12878_spartan_qc.yaml) into `config` (for extranal BAMs, set `bam_clean: fixrg` to make that `RG` equals the `description` you use, and turn off `mark_duplicates`)

3. Run bcbio

3. For external samples, `goleft indexcov` and `mosdepth`'s dist [wont run](https://github.com/bcbio/bcbio-nextgen/issues/2343). For that reason, generating those manually:

```
cd /data/cephfs/punim0010/projects/Saveliev_10X/NA12878-10x/bcbio_grch37/final_qc

goleft indexcov --directory ema_final/qc/coverage/ema_final --sex X,Y -- ../../ema_final.bam && mv ema_final/qc/coverage/ema_final ema_final/qc/coverage/indexcov

goleft indexcov --directory longranger_decoy_pos_sorted_bam/qc/coverage/longranger_decoy_pos_sorted_bam --sex X,Y -- /data/cephfs/punim0010/projects/Saveliev_10X/NA12878-10x/bcbio_grch37/bams/longranger_decoy_pos_sorted.bam && mv longranger_decoy_pos_sorted_bam/qc/coverage/longranger_decoy_pos_sorted_bam longranger_decoy_pos_sorted_bam/qc/coverage/indexcov

goleft indexcov --directory longranger_pos_sorted/qc/coverage/longranger_pos_sorted --sex X,Y -- /data/cephfs/punim0010/projects/Saveliev_10X/NA12878-10x/bcbio_grch37/bams/longranger_pos_sorted.bam && mv longranger_pos_sorted/qc/coverage/longranger_pos_sorted longranger_pos_sorted/qc/coverage/indexcov

mosdepth -t 30 -F 1804  --no-per-base --by /data/cephfs/punim0010/projects/Saveliev_10X/NA12878-10x/bcbio_grch37/work/bedprep/cov-2018-02-16_PMCC_Panel_Regions-merged.bed  longranger_decoy_pos_sorted_bam/qc/coverage/longranger_pos_sorted-coverage /data/cephfs/punim0010/projects/Saveliev_10X/NA12878-10x/bcbio_grch37/work/bamclean/longranger_pos_sorted/longranger_pos_sorted-fixrg.bam -T 1,5,10,20,50,100,250,500,1000,5000,10000,50000

mosdepth -t 30 -F 1804  --no-per-base --by /data/cephfs/punim0010/projects/Saveliev_10X/NA12878-10x/bcbio_grch37/work/bedprep/cov-2018-02-16_PMCC_Panel_Regions-merged.bed  longranger_decoy_pos_sorted_bam/qc/coverage/longranger_decoy_pos_sorted_bam-coverage /data/cephfs/punim0010/projects/Saveliev_10X/NA12878-10x/bcbio_grch37/work/bamclean/longranger_decoy_pos_sorted_bam/longranger_decoy_pos_sorted-fixrg.bam -T 1,5,10,20,50,100,250,500,1000,5000,10000,50000

mosdepth -t 30 -F 1804  --no-per-base --by /data/cephfs/punim0010/projects/Saveliev_10X/NA12878-10x/bcbio_grch37_ema/work/bedprep/cov-2018-02-16_PMCC_Panel_Regions-merged.bed  ema_final/qc/coverage/ema_final-coverage ../../ema_final.bam -T 1,5,10,20,50,100,250,500,1000,5000,10000,50000

# Convert .bed.gz to .tsv because otherwise MultiQC won't handle those:

gunzip -c longranger_decoy_pos_sorted_bam/qc/coverage/indexcov/longranger_decoy_pos_sorted_bam-indexcov.bed.gz > longranger_decoy_pos_sorted_bam/qc/coverage/indexcov/longranger_deco
y_pos_sorted_bam-indexcov.tsv

gunzip -c longranger_pos_sorted/qc/coverage/indexcov/longranger_pos_sorted-indexcov.bed.gz > longranger_pos_sorted/qc/coverage/indexcov/longranger_pos_sorted-indexcov.tsv

gunzip -c NA12878_10x_EMA/qc/coverage/indexcov/NA12878_10x_EMA-indexcov.bed.gz > NA12878_10x_EMA/qc/coverage/indexcov/NA12878_10x_EMA-indexcov.tsv
```

3. Copy `cp 2018-03-20_project/multiqc/list_files_final.txt .` info your final directory and add newly generated files into it:

```
ema_final/qc/coverage/indexcov/ema_final-indexcov.ped
ema_final/qc/coverage/indexcov/ema_final-indexcov.roc
ema_final/qc/coverage/indexcov/ema_final-indexcov.tsv
ema_final/qc/coverage/ema_final-coverage.mosdepth.dist.txt
longranger_pos_sorted/qc/coverage/indexcov/longranger_pos_sorted-indexcov.ped
longranger_pos_sorted/qc/coverage/indexcov/longranger_pos_sorted-indexcov.roc
longranger_pos_sorted/qc/coverage/indexcov/longranger_pos_sorted-indexcov.tsv
longranger_pos_sorted/qc/coverage/longranger_pos_sorted-coverage.mosdepth.dist.txt
longranger_decoy_pos_sorted_bam/qc/coverage/indexcov/longranger_decoy_pos_sorted_bam-indexcov.ped
longranger_decoy_pos_sorted_bam/qc/coverage/indexcov/longranger_decoy_pos_sorted_bam-indexcov.roc
longranger_decoy_pos_sorted_bam/qc/coverage/indexcov/longranger_decoy_pos_sorted_bam-indexcov.tsv
longranger_decoy_pos_sorted_bam/qc/coverage/longranger_decoy_pos_sorted_bam-coverage.mosdepth.dist.txt
```

4. Copy `cp 2018-03-20_project/multiqc/multiqc_config.yaml .` and add fix the `table_columns_visible` section with the following:

```
bcftools:
  write_separate_table: true
module_order:
- bcbio
- samtools
- goleft_indexcov
- peddy
- bcftools
- picard
- qualimap
- snpeff
- fastqc
- preseq
table_columns_visible:
  FastQC:
    percent_gc: false
  bcbio:
    Average_insert_size: true
    Usable_pct: false
    Ontarget_pct: false
  Peddy:
    error_sex_check: false
    family_id: false
    sex_het_ratio: false
    percent_gc: false
  QualiMap:
    1_x_pc: true
    5_x_pc: true
    10_x_pc: true
    30_x_pc: true
    50_x_pc: true
  "Samtools Stats":
    error_rate: false
```

5. Run MultiQC (load bcbio if needed):

```
multiqc -f -l list_files_final.txt -c multiqc_config.yaml -o multiqc 
```










