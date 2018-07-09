Exploring tricky regions

Looking at tricky variants from the EMA article https://github.com/arshajii/ema-paper-data/blob/master/experiments.ipynb

```
C4A 6:31965242
AMY1 1:104197843
CYP2D7 22:42537120
```

Making a bed file for the regions

```tricky_3.bed
1	104090000	104310000	AMY1/2
6	31940000	32010000	C4A/B
22	42521302	42549900	CY2D6/7
```

And slicing BAMs

```
cd /g/data3/gx8/projects/Saveliev_10X/COLO829-10x/ema_trimmed_runs/sliced_bams

~/bin/sambamba slice -L tricky_3.bed ../ema_COLO829_100pc/COLO829_100pc.bam -o COLO829_100pc.SLICED_3.bam
~/bin/sambamba slice -L tricky_3.bed ../ema_COLO829_20pc/COLO829_20pc.bam -o COLO829_20pc.SLICED_3.bam
~/bin/sambamba slice -L tricky_3.bed ../ema_COLO829_40pc/COLO829_40pc.bam -o COLO829_40pc.SLICED_3.bam
~/bin/sambamba slice -L tricky_3.bed ../ema_COLO829BL/COLO829BL.bam -o COLO829BL.SLICED_3.bam
```

Also slicing older EMA runs, Longranger, and BWA

```
~/bin/sambamba slice -L tricky_3.bed /g/data3/gx8/projects/Saveliev_10X/COLO829-10x/ema_older_runs/bams/Colo829_80pc_EMA.bam -o COLO829_80pc.OLD_EMA.SLICED3.bam
#~/bin/sambamba slice -L tricky_3.bed /g/data3/gx8/projects/Saveliev_10X/COLO829-10x/ema_older_runs/bams/Colo829Bl_10x_EMA.bam -o COLO829BL.OLD_EMA.SLICED3.bam

~/bin/sambamba slice -L tricky_3.bed /g/data3/gx8/projects/Saveliev_10X/COLO829-10x/longranger/bams/COLO829_80pc_LR.bam -o COLO829_80pc.LR.SLICED3.bam
#~/bin/sambamba slice -L tricky_3.bed /g/data3/gx8/projects/Saveliev_10X/COLO829-10x/longranger/bams/COLO829BL_LR.bam -o COLO829BL.LR.SLICED3.bam

```