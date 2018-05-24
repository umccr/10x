import sys
import pysam
from collections import defaultdict

vcf_path = sys.argv[1]
bam_path = sys.argv[2]

vcf = pysam.VariantFile(vcf_path)
bam = pysam.AlignmentFile(bam_path, "rb" )

vcf.header.filters.add('BX_support', None, None, 'The alternate allele supported by less than 2 GEMs')
vcf.header.info.add('BX_support', 1, 'Integer', 'The number of GEMs supporting the alternate allele')
vcf.header.info.add('BX_support_list', 1, 'String', 'The list GEMs supporting the alternate allele')

sys.stdout.write(str(vcf.header))

variants_total = 0
variants_filtered = 0

for rec in vcf:
    variants_total += 1

    pileupcolumns = bam.pileup(rec.chrom, rec.pos - 1, rec.pos)
    pileupcolumn = next((pc for pc in pileupcolumns if pc.reference_pos == rec.pos - 1), None)
    if pileupcolumn is None:
        sys.stderr.write(f'No pileup column fpr {rec.chrom}:{rec.pos}\n')
        continue

    barcodes_by_alt = defaultdict(set)

    for pileupread in pileupcolumn.pileups:
        aln = pileupread.alignment
        if aln.is_duplicate or aln.is_qcfail:
            continue

        if aln.has_tag('BX'):
            bx = aln.get_tag('BX')
        else:
            bx = 'None'

        if pileupread.is_del or pileupread.is_refskip:
            base = None  # TODO: better process deletions
        else:
            base = aln.query_sequence[pileupread.query_position]
            barcodes_by_alt[base].add(bx)

    for alt in rec.alts:
        bx_support_list = barcodes_by_alt.get(alt, [])
        rec.info['BX_support_list'] = ','.join(bx_support_list)
        rec.info['BX_support'] = len(bx_support_list)
        if len(bx_support_list) < 2:
            rec.filter.add('BX_support')
            variants_filtered += 1

    sys.stdout.write(str(rec))

    if variants_total % 1000 == 0:
        sys.stderr.write(f'Written {variants_total} variants, filtered {variants_filtered} variants\n')

sys.stderr.write(f'Written {variants_total} variants, filtered {variants_filtered} variants\n')

bam.close()
vcf.close()


# 1:
# NeverResponder filtering running - count
# NeverResponder LR variant calling running /data/cephfs/punim0010/projects/Saveliev_10X/DiploidNeverResponder/diploid_lariat/work
# Review LR calls vs EMA calls vs older normal calls

# WRONG!!!:
# * All *      non-filt     filt
# EMA only:       32181    15914
# non-10x only:    6208     6321
# shared:          5425     5370

# * SNP *      non-filt     filt
# EMA only:        9116     8613 
# non-10x only:    4319     4326 
# shared:          5212     5205 


# 2:
# Running Strelka2 for EMA and LR COLO829: /data/cephfs/punim0010/projects/Saveliev_10X/COLO829-10x/bcbio_vc/work
# After finishing, compare older our COLO829 run vs truth set VCF vs 10x strelka
# Review variants + BAM files from COLO829 10x vs regular
