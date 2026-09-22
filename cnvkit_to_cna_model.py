#!/usr/bin/env python3
"""
cnvkit_to_cna_model.py  —  Convert a CNVkit .cns file to fmiSGZ CNA model format.

Usage
-----
    python3 cnvkit_to_cna_model.py input.cns output.cna_calls.txt [--purity 0.5]

Column mapping (CNVkit → fmiSGZ)
---------------------------------
    chromosome  → CHR           (kept as-is; chr prefix optional)
    start + 1   → segStart      (BED 0-based → 1-based)
    end         → segEnd        (BED exclusive end = 1-based inclusive)
    log2        → segLR         (observed segment log2-ratio)
    baf         → segMAF        (observed BAF; empty cells become NA)
    cn          → CN            (total integer copy number)
    cn1         → numMAtumorPred (minor-allele copy number from CNVkit call)
    NA          → mafPred       (model-predicted MAF; requires purity — set NA)
    <purity>    → purity        (from --purity flag; default 0.5)
    2.0         → baseLevel     (diploid linear baseline: log2(2/2) = 0)

Notes
-----
* CNVkit uses a BED-style 0-based half-open coordinate system (start is
  0-based, end is exclusive).  fmiSGZ expects 1-based inclusive coordinates,
  so segStart = start + 1 and segEnd = end.

* baseLevel = 2.0 is the linear diploid reference level when CNVkit log2
  ratios are already centred at 0 for diploid segments.  fmiSGZ uses it in
  cn2lr_bl(p, bl, cn) = log2((p*cn + 2*(1-p)) / bl), so bl=2 yields 0 for
  a diploid segment (cn=2) regardless of purity, which is consistent with
  CNVkit's own normalisation.

* mafPred (model-predicted minor-allele frequency) requires purity, CN, and
  minor-allele CN.  When purity is unknown it is set to NA so fmiSGZ will
  report 'ambiguous_CNA_model' for variants where model error checks are
  triggered.  CN and segLR are still used for zygosity and basic call rules.
"""

import argparse
import csv
import sys


FMISGZ_HEADER = [
    'CHR', 'segStart', 'segEnd',
    'mafPred', 'CN', 'segLR', 'segMAF',
    'numMAtumorPred', 'purity', 'baseLevel',
]

# Diploid linear baseline: log2((p*2 + 2*(1-p)) / bl) = 0  →  bl = 2.0
DIPLOID_BASE_LEVEL = '2.0'


def _safe(value: str, fallback: str = 'NA') -> str:
    """Return value if non-empty/non-null, else fallback."""
    v = value.strip()
    if v in ('', '.', 'NA', 'NaN', 'nan', 'None'):
        return fallback
    try:
        float(v)
        return v
    except ValueError:
        return fallback


def _parse_space_cns(data_lines: list) -> list:
    """Parse space-delimited CNVkit .cns rows where some cells may be empty.

    The layout after splitting on whitespace is:
        tokens[0]   = chromosome
        tokens[1]   = start
        tokens[2]   = end
        tokens[3]   = gene          (comma-separated, no internal spaces)
        tokens[-3]  = depth
        tokens[-2]  = probes
        tokens[-1]  = weight
        tokens[4:-3] = middle fields: log2 [baf] ci_hi ci_lo cn [cn1] [cn2]

    The number of middle tokens indicates which optional fields are present:
        4 tokens → log2, ci_hi, ci_lo, cn             (no baf / cn1 / cn2)
        5 tokens → log2, baf, ci_hi, ci_lo, cn         (no cn1 / cn2)
        6 tokens → log2, baf, ci_hi, ci_lo, cn, cn1    (no cn2)
        7 tokens → log2, baf, ci_hi, ci_lo, cn, cn1, cn2
    """
    rows = []
    for line in data_lines:
        t = line.split()
        chrom = t[0]
        start = t[1]
        end   = t[2]
        gene  = t[3] if len(t) > 3 else ''
        # depth, probes, weight are always the last 3
        middle = t[4:-3] if len(t) > 7 else t[4:]

        n = len(middle)
        log2  = middle[0]      if n >= 1 else ''
        baf   = middle[1]      if n >= 5 else ''   # present only when ≥5 middle fields
        ci_hi = middle[2]      if n >= 5 else (middle[1] if n >= 2 else '')
        ci_lo = middle[3]      if n >= 5 else (middle[2] if n >= 3 else '')
        cn    = middle[4]      if n >= 5 else (middle[3] if n >= 4 else '')
        cn1   = middle[5]      if n >= 6 else ''
        cn2   = middle[6]      if n >= 7 else ''

        rows.append({
            'chromosome': chrom, 'start': start, 'end': end, 'gene': gene,
            'log2': log2, 'baf': baf, 'ci_hi': ci_hi, 'ci_lo': ci_lo,
            'cn': cn, 'cn1': cn1, 'cn2': cn2,
        })
    return rows


def convert(input_cns: str, output_cna: str, purity: float) -> None:
    purity_str = f'{purity:.4f}'

    with open(input_cns) as fin, \
         open(output_cna, 'w', newline='') as fout:

        raw_lines   = [l.rstrip('\n') for l in fin if l.strip()]
        header_line = raw_lines[0]
        data_lines  = raw_lines[1:]

        # Detect delimiter: tab-separated (standard CNVkit) or space-padded.
        if '\t' in header_line:
            import csv as _csv
            reader = _csv.DictReader([header_line] + data_lines, delimiter='\t')
            rows = list(reader)
        else:
            # Space-padded fixed-width: use column positions from the header.
            rows = _parse_space_cns(data_lines)

        writer = csv.DictWriter(fout, fieldnames=FMISGZ_HEADER, delimiter='\t',
                                lineterminator='\n')
        writer.writeheader()

        for row in rows:
            chrom    = row['chromosome'].strip()
            start_0  = int(row['start'])          # 0-based
            end_excl = int(row['end'])             # exclusive

            seg_start = start_0 + 1               # → 1-based inclusive
            seg_end   = end_excl                  # already 1-based inclusive

            seg_lr  = _safe(row.get('log2', ''))
            seg_maf = _safe(row.get('baf', ''))
            cn      = _safe(row.get('cn', ''))
            cn1     = _safe(row.get('cn1', ''))   # minor-allele CN

            writer.writerow({
                'CHR'            : chrom,
                'segStart'       : seg_start,
                'segEnd'         : seg_end,
                'mafPred'        : 'NA',           # cannot compute without purity
                'CN'             : cn,
                'segLR'          : seg_lr,
                'segMAF'         : seg_maf,
                'numMAtumorPred' : cn1,
                'purity'         : purity_str,
                'baseLevel'      : DIPLOID_BASE_LEVEL,
            })

    print(f'Written: {output_cna}  (purity={purity_str}, baseLevel={DIPLOID_BASE_LEVEL})')
    if purity == 0.5:
        print('WARNING: purity defaulted to 0.5.  Somatic/germline calls will be '
              'approximate.  Supply a purity estimate with --purity for better accuracy.')
    print('NOTE: mafPred=NA — variants will fall back to CN+logOR rules only '
          '(no MAF model error filter).')


def main() -> None:
    parser = argparse.ArgumentParser(
        description='Convert CNVkit .cns output to fmiSGZ CNA model format.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument('input_cns',  help='CNVkit called segments file (.cns)')
    parser.add_argument('output_cna', help='Output CNA model file for fmiSGZ_v2.py')
    parser.add_argument(
        '--purity', type=float, default=0.5,
        help='Tumour purity estimate [0,1] (default: 0.5). '
             'Use a value from ABSOLUTE, PureCN, or pathology estimate.',
    )
    args = parser.parse_args()

    if not (0.0 < args.purity <= 1.0):
        sys.exit('ERROR: --purity must be in (0, 1]')

    convert(args.input_cns, args.output_cna, args.purity)


if __name__ == '__main__':
    main()
