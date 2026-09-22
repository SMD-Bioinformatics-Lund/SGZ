#!/usr/bin/env python3.13
"""
Parse a VEP-annotated VCF file and output a mutation table compatible with
the mut_aggr.full.txt format used by fmiSGZ_v2 and basicSGZ pipelines.

Mutation annotation uses the MANE SELECT transcript (falls back to CANONICAL).
Frequency = VAF (FORMAT field), Depth = VD (ALT allele count, FORMAT field).
"""

import argparse
import os
import re
import sys
from urllib.parse import unquote

import pysam

# Suppress pysam/htslib contig-not-in-header warnings
pysam.set_verbosity(0)


# ---------------------------------------------------------------------------
# Amino acid 3-letter → 1-letter conversion
# ---------------------------------------------------------------------------
AA3TO1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*", "Sec": "U", "Pyl": "O", "Xaa": "X",
}

# Pre-compile regex that matches any 3-letter AA code
_AA3_RE = re.compile("|".join(AA3TO1.keys()))


def convert_hgvsp(raw: str) -> str:
    """
    Convert a raw HGVSp string (e.g. 'ENSP00000xxx.y:p.Asp1228Glu') to
    1-letter amino-acid notation (e.g. 'p.D1228E').
    URL-encoded characters (%3D, %2A …) are decoded first.
    """
    if not raw:
        return ""
    # Strip protein-ID prefix
    if ":" in raw:
        raw = raw.split(":", 1)[1]
    # URL-decode (%3D → =, %2A → *, etc.)
    raw = unquote(raw)
    # Replace 3-letter codes with 1-letter
    return _AA3_RE.sub(lambda m: AA3TO1[m.group()], raw)


# ---------------------------------------------------------------------------
# Consequence → effect label
# ---------------------------------------------------------------------------
def consequence_to_effect(consequence: str) -> str:
    """Map a VEP Consequence string to the short effect label."""
    terms = set(consequence.split("&"))
    if "frameshift_variant" in terms:
        return "frameshift"
    if terms & {"stop_gained", "stop_lost", "start_lost"}:
        return "nonsense"
    if terms & {
        "splice_donor_variant", "splice_acceptor_variant",
        "splice_region_variant", "splice_donor_5th_base_variant",
        "splice_polypyrimidine_tract_variant",
    }:
        return "splice"
    if "missense_variant" in terms:
        return "missense"
    if terms & {"inframe_insertion", "inframe_deletion"}:
        return "inframe_indel"
    if "synonymous_variant" in terms:
        return "synonymous"
    return consequence.split("&")[0]  # fallback: first term


# ---------------------------------------------------------------------------
# Mutation string builder
# ---------------------------------------------------------------------------
def build_mutation(symbol: str, nm_tx: str, hgvsc: str,
                   hgvsp_raw: str, consequence: str) -> str:
    """
    Build the mutation string in the form:
      GENE:NM_xxxxx:c.XXX_p.YYY         (missense / synonymous / nonsense)
      GENE:NM_xxxxx:c.XXX:splice         (splice, no protein change)
      GENE:NM_xxxxx:c.XXX_p.YYYfs*N:frameshift  (frameshift)
    """
    # Strip transcript/protein prefixes from HGVS fields
    if ":" in hgvsc:
        hgvsc = hgvsc.split(":", 1)[1]

    hgvsp = convert_hgvsp(hgvsp_raw)

    effect = consequence_to_effect(consequence)

    if effect == "splice" and not hgvsp:
        return f"{symbol}:{nm_tx}:{hgvsc}:splice"

    if hgvsp:
        base = f"{symbol}:{nm_tx}:{hgvsc}_{hgvsp}"
        if effect == "frameshift":
            return base + ":frameshift"
        return base

    # Fallback: no protein change available
    return f"{symbol}:{nm_tx}:{hgvsc}"


# ---------------------------------------------------------------------------
# MANE SELECT / CANONICAL transcript selector
# ---------------------------------------------------------------------------
def select_transcript(csq_list: list[dict]) -> dict | None:
    """
    Return the CSQ annotation dict for the MANE SELECT transcript.
    Falls back to CANONICAL (YES) if no MANE SELECT is found.
    Only considers protein-coding transcripts.
    """
    mane = None
    canonical = None
    for ann in csq_list:
        if ann.get("BIOTYPE") != "protein_coding":
            continue
        if ann.get("MANE_SELECT"):
            mane = ann
            break
        if ann.get("CANONICAL") == "YES" and canonical is None:
            canonical = ann
    return mane or canonical


# ---------------------------------------------------------------------------
# Parse CSQ INFO field
# ---------------------------------------------------------------------------
def parse_csq(record: pysam.VariantRecord, csq_fields: list[str]) -> list[dict]:
    """Return list of per-transcript annotation dicts from the CSQ INFO tag."""
    raw = record.info.get("CSQ")
    if not raw:
        return []
    result = []
    for entry in raw:
        vals = entry.split("|")
        # Pad to expected length in case trailing empty fields are stripped
        vals += [""] * (len(csq_fields) - len(vals))
        result.append(dict(zip(csq_fields, vals)))
    return result


# ---------------------------------------------------------------------------
# NM transcript ID (strip version suffix)
# ---------------------------------------------------------------------------
def nm_id(mane_select: str) -> str:
    """
    Extract the NM_ ID without version from a MANE SELECT value.
    E.g. 'NM_004958.3' → 'NM_004958'
    If the value contains a space-separated ENST, take only the NM_ part.
    """
    # MANE_SELECT may be like 'NM_004958.3' or 'NM_004958.3 ENST00000...'
    nm = mane_select.split()[0]
    return nm.split(".")[0]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def parse_vcf(vcf_path: str, sample_name: str, out_fh,
              include_synonymous: bool = False,
              min_vaf: float | None = None):
    vcf = pysam.VariantFile(vcf_path)

    # Extract CSQ field names from the header description
    csq_desc = vcf.header.info["CSQ"].description
    # Format: "... Format: Field1|Field2|..."
    csq_fields = csq_desc.split("Format: ", 1)[1].split("|")

    header = ["#sample", "mutation", "frequency", "depth", "pos",
              "status", "strand", "effect"]
    out_fh.write("\t".join(header) + "\n")

    for record in vcf.fetch():
        chrom = record.chrom
        if not chrom.startswith("chr"):
            chrom = "chr" + chrom
        pos = record.pos  # pysam returns 1-based POS

        # Parse format fields from the single sample
        sample_data = record.samples[list(record.samples)[0]]
        vaf = sample_data.get("VAF")
        vd = sample_data.get("DP")

        if vaf is None or vd is None:
            continue

        # --------------- VAF FILTER ----------------------------------------
        # Toggle with --min-vaf / --no-vaf-filter at the command line.
        # Set min_vaf=None to disable entirely.
        if min_vaf is not None and float(vaf) < min_vaf:
            continue
        # -------------------------------------------------------------------

        # Parse CSQ and pick the MANE SELECT / CANONICAL transcript
        csq_list = parse_csq(record, csq_fields)
        ann = select_transcript(csq_list)
        if ann is None:
            continue

        # Skip non-coding / synonymous consequences
        consequence = ann.get("Consequence", "")
        effect = consequence_to_effect(consequence)
        allowed = {"missense", "nonsense", "splice", "frameshift", "inframe_indel"}
        if include_synonymous:
            allowed.add("synonymous")
        if effect not in allowed:
            continue

        symbol = ann.get("SYMBOL", "")
        mane_select = ann.get("MANE_SELECT", "")
        feature = ann.get("Feature", "")
        hgvsc = ann.get("HGVSc", "")
        hgvsp_raw = ann.get("HGVSp", "")
        strand_val = ann.get("STRAND", "")

        # Determine transcript ID (prefer MANE SELECT NM_, else Feature)
        if mane_select and mane_select.startswith("NM_"):
            tx = nm_id(mane_select)
        elif feature and feature.startswith("NM_"):
            tx = nm_id(feature)
        else:
            # No NM transcript available – skip
            continue

        if not symbol or not hgvsc:
            continue

        mutation = build_mutation(symbol, tx, hgvsc, hgvsp_raw, consequence)

        strand = "+" if strand_val == "1" else ("-" if strand_val == "-1" else ".")
        freq = round(float(vaf), 4)
        depth = int(vd)   # VD = ALT allele observation count (as requested)

        row = [
            sample_name,
            mutation,
            str(freq),
            str(depth),
            f"{chrom}:{pos}",
            ann.get("IMPACT", "unknown").lower() or "unknown",
            strand,
            effect,
        ]
        out_fh.write("\t".join(row) + "\n")

    vcf.close()


def main():
    parser = argparse.ArgumentParser(
        description="Convert VEP-annotated VCF to mutation table (mut_aggr.full.txt format)"
    )
    parser.add_argument("vcf", help="Input VCF file (VEP-annotated)")
    parser.add_argument(
        "-s", "--sample",
        help="Sample name to use in output (default: VCF filename stem)",
    )
    parser.add_argument(
        "-o", "--output",
        help="Output file path (default: stdout)",
    )
    parser.add_argument(
        "--synonymous", action="store_true",
        help="Also include synonymous variants (excluded by default)",
    )

    vaf_group = parser.add_mutually_exclusive_group()
    vaf_group.add_argument(
        "--min-vaf", type=float, default=0.5, metavar="FLOAT",
        help="Filter out variants with VAF < FLOAT (default: 0.5)",
    )
    vaf_group.add_argument(
        "--no-vaf-filter", action="store_true",
        help="Disable the VAF filter entirely",
    )

    args = parser.parse_args()

    sample_name = args.sample or os.path.basename(args.vcf).split(".")[0]

    min_vaf = None if args.no_vaf_filter else args.min_vaf

    out_fh = open(args.output, "w") if args.output else sys.stdout
    try:
        parse_vcf(args.vcf, sample_name, out_fh,
                  include_synonymous=args.synonymous,
                  min_vaf=min_vaf)
    finally:
        if args.output:
            out_fh.close()


if __name__ == "__main__":
    main()
