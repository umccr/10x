#!/usr/bin/env python
# Extracts manually curated (eyeball) GRC38h coordinates from IGV screenshot filenames
import re
from pathlib import Path

dataset_root=Path("data/processed/igv_telomere_curated")

for image_path in dataset_root.iterdir():
    fname = image_path.parts[-1]
    if not "notelomeresseen" in fname:
        fname_mdata = fname.split("_")
        
        chrom = fname_mdata[0]
        direction = fname_mdata[1]
        coords = fname_mdata[2].replace(".png", "").split("-")

        print("{}\t{}\t{}".format(chrom, coords[0], coords[1]))
        


