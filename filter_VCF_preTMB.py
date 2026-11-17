#!/usr/local/bin/python

"""
Description: this script, is based on post_annotation_filtering.py from the SPP, and it aims at removing variants from the
VEP annotated final.filtered.vcf (vcf uploaded to Coyote, output from SNV_ANNOTATE), 
that has a gnomAD above --max_freq (default 0.001) or a user defined value. 
It will keep variants that have no gnomAD annotation.
Variants are excluded if their gnomAD allele frequency is greater
than or equal to --max_freq (default 0.001, or 0.1%), or if their
FILTER field matches any of the patterns provided via --filters.

Variants without a usable gnomAD frequency annotation are kept,
provided they do not match any of the specified FILTER patterns.

The script checks gnomAD exome (gnomAD_AF) and genome (gnomADg_AF)
frequencies. If both are available, the highest frequency is used.
Multiple values separated by '&' are handled by taking the maximum.

Usage: python filter_gnomad.py \
    -f input.vcf \
    -o filtered.vcf \
    --max_freq 0.001 \
    --filters "LowQual,VQSR*,StrandBias"
Arguments:
--vcf, -f:
Path to the input VCF file (required).
--output, -o:
    Path to the output VCF file (required).

--max_freq:
    Maximum allowed gnomAD allele frequency. Variants with a
    frequency greater than or equal to this value are excluded.
    Optional; defaults to 0.001 (0.1%). Accepts a float between
    0 and 1.

--filters:
    Optional, comma-separated list of VCF FILTER strings or patterns
    to exclude. Use of wildcards possible, e.g. "LowQual,VQSR*".
    If omitted, no FILTER-string exclusions are applied.

"""

import argparse
import fnmatch
import re
import sys


def max_gnomad(value):
    """Return the maximum numeric value from a '&'-separated annotation."""
    if not value or value == ".":
        return None

    try:
        values = [
            float(x)
            for x in value.split("&")
            if x and x != "."
        ]
        return max(values) if values else None
    except (ValueError, TypeError):
        return None


def get_csq_fields(header_lines):
    """Extract VEP CSQ field names from the VCF header."""
    for line in header_lines:
        if line.startswith("##INFO=<ID=CSQ,"):
            match = re.search(r'Format:\s*([^"]+)', line)

            if match:
                return match.group(1).rstrip(">").split("|")

    return []


def get_annotations(info, csq_fields):
    """Parse VEP CSQ annotations from the INFO column."""
    csq_value = None

    for item in info.split(";"):
        if item.startswith("CSQ="):
            csq_value = item[4:]
            break

    if not csq_value or not csq_fields:
        return []

    annotations = []

    for entry in csq_value.split(","):
        values = entry.split("|")

        annotation = {
            field: values[i] if i < len(values) else ""
            for i, field in enumerate(csq_fields)
        }

        annotations.append(annotation)

    return annotations


def check_gnomad_vaf(annotations, af_cutoff):
    """Return False if any available gnomAD AF reaches the cutoff."""
    for annotation in annotations:
        frequencies = [
            max_gnomad(annotation.get("gnomAD_AF")),
            max_gnomad(annotation.get("gnomADg_AF")),
        ]

        frequencies = [
            af for af in frequencies if af is not None
        ]

        if frequencies and max(frequencies) >= af_cutoff:
            return False

    return True


def filter_flags(filter_column, filters):
    """Return False if a FILTER value matches an exclusion pattern."""
    if not filters:
        return True

    filter_values = filter_column.split(";")

    return not any(
        fnmatch.fnmatch(value, pattern)
        for value in filter_values
        for pattern in filters
    )


def read_vcf(infile, outfile, filters, af_cutoff):
    """Filter a VEP-annotated VCF using text-based I/O."""

    total = 0
    excluded_filter = 0
    excluded_gnomad = 0
    kept = 0

    with open(infile, "r") as source, open(outfile, "w") as output:

        # Read and preserve the complete original header.
        header_lines = []

        for line in source:
            if line.startswith("#"):
                header_lines.append(line)
            else:
                break

        output.writelines(header_lines)

        csq_fields = get_csq_fields(header_lines)

        if not csq_fields:
            print(
                "WARNING: Could not identify VEP CSQ fields in the header. "
                "Variants will be retained based on FILTER only.",
                file=sys.stderr,
            )

        elif not any(
            field in csq_fields
            for field in ("gnomAD_AF", "gnomADg_AF")
        ):
            print(
                "WARNING: Neither gnomAD_AF nor gnomADg_AF was found "
                "in the CSQ header. gnomAD filtering cannot be applied.",
                file=sys.stderr,
            )

        # Process variant records as original text lines.
        for line in source:
            if not line.strip():
                continue

            total += 1
            fields = line.rstrip("\r\n").split("\t")

            if len(fields) < 8:
                print(
                    f"WARNING: Skipping malformed record: {line[:100]}",
                    file=sys.stderr,
                )
                continue

            # Apply FILTER exclusions first.
            if not filter_flags(fields[6], filters):
                excluded_filter += 1
                continue

            # Parse and apply gnomAD filtering.
            annotations = get_annotations(fields[7], csq_fields)

            if not check_gnomad_vaf(annotations, af_cutoff):
                excluded_gnomad += 1
                continue

            # Preserve the original variant record unchanged.
            output.write(line)
            kept += 1

    print(f"Total input variants:       {total}", file=sys.stderr)
    print(f"Excluded by FILTER:         {excluded_filter}", file=sys.stderr)
    print(f"Excluded by gnomAD:         {excluded_gnomad}", file=sys.stderr)
    print(f"Retained variants:           {kept}", file=sys.stderr)
    print(f"Output written to:           {outfile}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Filter a VEP-annotated VCF by gnomAD frequency "
            "and FILTER strings."
        )
    )

    parser.add_argument(
        "--vcf", "-f",
        required=True,
        help="Input VCF file",
    )

    parser.add_argument(
        "--output", "-o",
        required=True,
        help="Output VCF file",
    )

    parser.add_argument(
        "--max_freq",
        type=float,
        default=0.001,
        help="Maximum gnomAD frequency (default: 0.001 = 0.1%%)",
    )

    parser.add_argument(
        "--filters",
        default="",
        help="Comma-separated FILTER patterns to exclude, e.g. FAIL*,*PON*,Q10",
    )

    args = parser.parse_args()

    if not 0 <= args.max_freq <= 1:
        parser.error("--max_freq must be between 0 and 1")

    filters = [
        pattern.strip()
        for pattern in args.filters.split(",")
        if pattern.strip()
    ]

    read_vcf(
        args.vcf,
        args.output,
        filters,
        args.max_freq,
    )


if __name__ == "__main__":
    main()