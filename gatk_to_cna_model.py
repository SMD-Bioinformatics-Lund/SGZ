#!/usr/bin/env python3
"""
gatk_to_cna_model.py  —  Convert GATK4 ModelSegments output to fmiSGZ CNA model format.

Reads:  *.modelFinal.seg  (output of GATK4 ModelSegments run with a PON)
Writes: *.cna_calls.txt   (input to fmiSGZ_v2.py)

Usage
-----
    # With purity known (best accuracy):
    python3 gatk_to_cna_model.py tumor.modelFinal.seg output.cna_calls.txt --purity 0.72

    # With purity unknown (default 0.5, approximate):
    python3 gatk_to_cna_model.py tumor.modelFinal.seg output.cna_calls.txt

    # Diploid tumor, purity 0.8, ploidy 2 (default):
    python3 gatk_to_cna_model.py tumor.modelFinal.seg output.cna_calls.txt --purity 0.8 --ploidy 2

    # Aneuploid tumor, ploidy estimated at 3:
    python3 gatk_to_cna_model.py tumor.modelFinal.seg output.cna_calls.txt --purity 0.65 --ploidy 3

GATK4 ModelSegments output (*.modelFinal.seg)
---------------------------------------------
SAM-style @ header lines (skipped), then a tab-separated data table:

  CONTIG                            chromosome
  START                             1-based start  (already 1-based in GATK4)
  END                               1-based end
  NUM_POINTS_COPY_RATIO             number of CR targets in segment
  NUM_POINTS_ALLELE_FRACTION        number of het SNPs in segment
  LOG2_COPY_RATIO_POSTERIOR_50      median log2 copy ratio  (segLR)
  LOG2_COPY_RATIO_POSTERIOR_10      10th-percentile LR
  LOG2_COPY_RATIO_POSTERIOR_90      90th-percentile LR
  MINOR_ALLELE_FRACTION_POSTERIOR_50 median minor-allele fraction (segMAF)
  MINOR_ALLELE_FRACTION_POSTERIOR_10 10th-percentile MAF
  MINOR_ALLELE_FRACTION_POSTERIOR_90 90th-percentile MAF

When NUM_POINTS_ALLELE_FRACTION == 0, the MAF columns are NaN — these
segments have no informative heterozygous SNPs (e.g. centromeres, very short
segments, or homozygous deletions).

Derived quantities
------------------
GATK4 does NOT call integer copy numbers directly.  This script derives them
from the log2 copy ratio using the mixture model:

    Expected log2 CR = log2((p * CN  +  2 * (1-p)) / baseLevel)

    where  p         = tumour purity
           CN        = integer total copy number in tumour cells
           baseLevel = p * ploidy + 2 * (1-p)   [linear diploid reference]

Solving for CN:
    CN = round( (2^segLR * baseLevel - 2*(1-p)) / p )

Minor-allele copy number is derived from segMAF:
    segMAF = (p * M  +  (1-p)) / (p * CN + 2*(1-p))
    M      = round( (segMAF * (p*CN + 2*(1-p)) - (1-p)) / p )

Model-predicted MAF (mafPred) is then back-computed from the called M and CN:
    mafPred = (p * M  +  (1-p)) / (p * CN + 2*(1-p))

Column mapping (GATK4 → fmiSGZ)
---------------------------------
  CONTIG                              → CHR
  START                               → segStart   (no conversion needed)
  END                                 → segEnd
  LOG2_COPY_RATIO_POSTERIOR_50        → segLR
  MINOR_ALLELE_FRACTION_POSTERIOR_50  → segMAF     (NA when no SNPs)
  derived from segLR + purity         → CN
  derived from segMAF + CN + purity   → numMAtumorPred
  back-computed from CN + M + purity  → mafPred    (NA when MAF absent)
  user-supplied                       → purity
  p * ploidy + 2 * (1-p)             → baseLevel
"""

import argparse
import csv
import math
import sys


FMISGZ_HEADER = [
    'CHR', 'segStart', 'segEnd',
    'mafPred', 'CN', 'segLR', 'segMAF',
    'numMAtumorPred', 'purity', 'baseLevel',
]


def _fmt(value, digits=6) -> str:
    """Format a float to a fixed number of significant digits, or return 'NA'."""
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return 'NA'
    return f'{value:.{digits}g}'


def _derive_cn(seg_lr: float, purity: float, base_level: float) -> int:
    """Integer total copy number from observed log2 ratio.

    Expected LR = log2((p*CN + 2*(1-p)) / base_level)
    → CN = (2^LR * base_level - 2*(1-p)) / p
    """
    linear = (2 ** seg_lr) * base_level
    cn_raw = (linear - 2 * (1 - purity)) / purity
    return max(0, round(cn_raw))


def _derive_minor_cn(seg_maf: float, cn: int, purity: float) -> int:
    """Integer minor-allele copy number from observed BAF/MAF.

    segMAF = (p*M + (1-p)) / (p*CN + 2*(1-p))
    → M = (segMAF * (p*CN + 2*(1-p)) - (1-p)) / p

    M is clamped to [0, CN//2] (minor allele ≤ major allele).
    """
    denom = purity * cn + 2 * (1 - purity)
    m_raw = (seg_maf * denom - (1 - purity)) / purity
    return max(0, min(cn // 2, round(m_raw)))


def _maf_pred(cn: int, minor_cn: int, purity: float) -> float:
    """Back-compute the model-predicted MAF from integer CN and minor CN."""
    denom = purity * cn + 2 * (1 - purity)
    if denom == 0:
        return float('nan')
    return (purity * minor_cn + (1 - purity)) / denom


def convert(input_seg: str, output_cna: str, purity: float, ploidy: float) -> None:
    # Diploid linear reference: value such that log2((p*ploidy + 2*(1-p)) / base_level) = 0
    base_level = purity * ploidy + 2 * (1 - purity)
    purity_str    = f'{purity:.4f}'
    base_level_str = f'{base_level:.6g}'

    n_written   = 0
    n_no_maf    = 0   # segments with no heterozygous SNPs
    n_del       = 0   # homozygous deletions (CN=0)

    with open(input_seg) as fin, open(output_cna, 'w', newline='') as fout:
        writer = csv.DictWriter(fout, fieldnames=FMISGZ_HEADER, delimiter='\t',
                                lineterminator='\n')
        writer.writeheader()

        reader = None
        for line in fin:
            stripped = line.rstrip('\n')

            # Skip SAM-style @ header lines and blank lines.
            if stripped.startswith('@') or not stripped:
                continue

            # First non-@ line is the column header.
            if reader is None:
                cols   = stripped.split('\t')
                reader = True    # sentinel — header consumed
                continue

            # Data rows.
            parts = stripped.split('\t')
            row   = dict(zip(cols, parts))

            chrom     = row['CONTIG']
            seg_start = int(row['START'])     # already 1-based
            seg_end   = int(row['END'])

            # Log2 copy ratio.
            lr_str = row.get('LOG2_COPY_RATIO_POSTERIOR_50', '')
            if not lr_str or lr_str in ('NA', 'NaN', 'nan', '.'):
                seg_lr = None
            else:
                seg_lr = float(lr_str)

            # Minor allele fraction (absent when no informative SNPs).
            maf_str = row.get('MINOR_ALLELE_FRACTION_POSTERIOR_50', '')
            n_snps  = int(row.get('NUM_POINTS_ALLELE_FRACTION', '0') or '0')
            if not maf_str or maf_str in ('NA', 'NaN', 'nan', '.') or n_snps == 0:
                seg_maf = None
                n_no_maf += 1
            else:
                seg_maf = float(maf_str)

            # Derive integer CN from log2 ratio.
            if seg_lr is not None:
                cn = _derive_cn(seg_lr, purity, base_level)
            else:
                cn = None

            if cn == 0:
                n_del += 1

            # Derive minor-allele CN from MAF (only meaningful when MAF present
            # and CN > 0).
            if seg_maf is not None and cn is not None and cn > 0:
                minor_cn  = _derive_minor_cn(seg_maf, cn, purity)
                maf_pred  = _maf_pred(cn, minor_cn, purity)
            else:
                minor_cn  = None
                maf_pred  = None

            writer.writerow({
                'CHR'            : chrom,
                'segStart'       : seg_start,
                'segEnd'         : seg_end,
                'mafPred'        : _fmt(maf_pred),
                'CN'             : cn if cn is not None else 'NA',
                'segLR'          : _fmt(seg_lr),
                'segMAF'         : _fmt(seg_maf),
                'numMAtumorPred' : minor_cn if minor_cn is not None else 'NA',
                'purity'         : purity_str,
                'baseLevel'      : base_level_str,
            })
            n_written += 1

    # Summary
    print(f'Written : {output_cna}')
    print(f'Segments: {n_written}  '
          f'(no-MAF: {n_no_maf},  homoDel CN=0: {n_del})')
    print(f'Purity  : {purity_str}   Ploidy: {ploidy}   baseLevel: {base_level_str}')
    if purity == 0.5:
        print('WARNING : purity defaulted to 0.5 — supply --purity for accurate CN calls.')
    if ploidy != 2:
        print(f'NOTE    : non-diploid ploidy {ploidy} used for baseLevel; '
              'verify with genome-wide LR median.')


def main() -> None:
    parser = argparse.ArgumentParser(
        description='Convert GATK4 ModelSegments *.modelFinal.seg to fmiSGZ CNA model format.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument('input_seg',  help='GATK4 ModelSegments output (*.modelFinal.seg)')
    parser.add_argument('output_cna', help='Output CNA model file for fmiSGZ_v2.py')
    parser.add_argument(
        '--purity', type=float, default=0.5,
        help='Tumour purity [0,1] (default: 0.5). '
             'Estimate from ABSOLUTE, PureCN, or pathology report.',
    )
    parser.add_argument(
        '--ploidy', type=float, default=2.0,
        help='Tumour ploidy (default: 2.0 = diploid). '
             'Use genome-wide log2 ratio median to check: '
             'if median LR >> 0 the tumour may be polyploid.',
    )
    args = parser.parse_args()

    if not (0.0 < args.purity <= 1.0):
        sys.exit('ERROR: --purity must be in (0, 1]')
    if args.ploidy <= 0:
        sys.exit('ERROR: --ploidy must be positive')

    convert(args.input_seg, args.output_cna, args.purity, args.ploidy)


if __name__ == '__main__':
    main()
