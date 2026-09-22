#!/usr/bin/env python3
"""
basicSGZ_v2 - Basic Somatic/Germline Zygosity (SGZ) classifier.

This module implements a simplified SGZ method for classifying short variants
(SNVs and small indels) detected by targeted next-generation sequencing as
somatic or germline, without requiring a copy-number (CNA) model.

Algorithm overview
------------------
For each short variant the algorithm:

1. Excludes chrX variants (sex chromosomes are not modelled).
2. Classifies high-frequency variants (AF > 0.95) as homozygous.
   If tumour purity is available and below *pthresh*, they are called germline.
3. For remaining variants, performs a two-sided binomial exact test against
   H0: AF = 0.5.  Variants with p < alpha are called somatic (allelic
   imbalance inconsistent with a diploid germline heterozygous state);
   the remainder are called germline heterozygous.
4. Optionally computes Tumor Mutation Burden (TMB) as the count of somatic
   variants divided by the effective panel size in megabases.

Output
------
One tab-separated output file per sample:

* ``*.basic.sgz.txt`` – per-variant somatic/germline calls with zygosity.
* ``*.tmb.txt``        – TMB summary (written only when ``--panel-size`` is given).

Usage
-----
::

    basicSGZ_v2.py [options] aggregated_mutations_file
    basicSGZ_v2.py [-h | --help]
    basicSGZ_v2.py [--version]

Dependencies
------------
* Python >= 3.8
* numpy
* pandas
* scipy >= 1.7  (uses ``scipy.stats.binomtest``, replacing the deprecated
                 ``scipy.stats.binom_test``)
"""

import argparse
import logging
import os
import sys
from typing import List, Optional, Union

import numpy as np
import pandas as pd
from scipy.stats import binomtest

__author__  = ''
__version__ = '2.0.0'
__doc__     = 'Basic SGZ method to evaluate variant origin and zygosity.'

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Purity helper
# ---------------------------------------------------------------------------

def read_pathology_purity(pathology_purity_file: Optional[str]) -> Union[float, str]:
    """Read tumour purity from a single-value text file.

    Parameters
    ----------
    pathology_purity_file : str or None
        Path to a plain-text file whose first line contains the percent
        tumour nuclei (e.g. ``"75"`` for 75 %).  Pass ``None`` to skip.

    Returns
    -------
    float
        Purity as a fraction in [0, 1] (raw value divided by 100).
    str
        ``"NA"`` when the file is absent, unreadable, or non-numeric.
    """
    if not pathology_purity_file or not os.path.isfile(pathology_purity_file):
        return "NA"

    with open(pathology_purity_file) as f:
        value = f.readline().strip()

    try:
        return float(value) / 100
    except ValueError:
        logger.warning("Cannot parse purity value %r from %s – using NA",
                       value, pathology_purity_file)
        return "NA"


# ---------------------------------------------------------------------------
# Statistical helpers
# ---------------------------------------------------------------------------

def compute_pvalues(freq: pd.Series, depth: pd.Series) -> np.ndarray:
    """Run vectorised two-sided binomial exact tests against H0: AF = 0.5.

    For each variant the observed alt-allele count is estimated as
    ``round(freq * depth)`` and tested against a Binomial(n=depth, p=0.5)
    distribution.  A significant result (small p-value) indicates that the
    observed frequency is inconsistent with a diploid heterozygous germline
    state and therefore supports a somatic call.

    Parameters
    ----------
    freq : pd.Series
        Observed allele frequencies (values clipped to [0, 1]).
    depth : pd.Series
        Total sequencing depths.  Values are cast to ``int`` before use.

    Returns
    -------
    np.ndarray
        Array of two-sided p-values, one per variant.
    """
    depth_int = depth.astype(int)
    successes = np.round(freq * depth_int).astype(int)

    pvals = [
        binomtest(int(k), int(n), p=0.5).pvalue
        for k, n in zip(successes, depth_int)
    ]

    return np.array(pvals)


# ---------------------------------------------------------------------------
# TMB
# ---------------------------------------------------------------------------

def compute_tmb(result_df: pd.DataFrame, panel_size_mb: float) -> float:
    """Compute Tumor Mutation Burden (TMB) from classified variants.

    TMB is defined as the number of somatic variants divided by the
    effective genomic panel size in megabases (mut/Mb).  Only variants
    whose ``germline/somatic`` column equals ``"somatic"`` are counted.

    Parameters
    ----------
    result_df : pd.DataFrame
        Output of :func:`run_sgz`.  Must contain a ``"germline/somatic"``
        column with values ``"somatic"``, ``"germline"``, or ``"ambiguous"``.
    panel_size_mb : float
        Effective size of the sequenced/analysed genomic region in Mb
        (e.g. 1.0 for a 1 Mb targeted panel, ~30 for whole-exome sequencing).

    Returns
    -------
    float
        TMB in mutations per megabase (mut/Mb).

    Raises
    ------
    ValueError
        If *panel_size_mb* is not positive.
    """
    if panel_size_mb <= 0:
        raise ValueError(f"panel_size_mb must be positive, got {panel_size_mb}")

    n_somatic = (result_df["germline/somatic"] == "somatic").sum()
    return float(n_somatic) / panel_size_mb


# ---------------------------------------------------------------------------
# Core SGZ algorithm
# ---------------------------------------------------------------------------

def run_sgz(
    df: pd.DataFrame,
    pathology_purity: Union[float, str],
    alpha: float,
    pthresh: float,
) -> pd.DataFrame:
    """Classify each variant as somatic or germline using the basic SGZ model.

    The function applies a two-step decision rule:

    1. **High-frequency rule** (AF > 0.95): variant is called *homozygous*.
       If tumour purity is known and below *pthresh*, it is additionally
       labelled *germline*.
    2. **Binomial test** (AF <= 0.95): a two-sided exact binomial test is
       performed against H0: AF = 0.5.  Significant variants (p < *alpha*)
       are called *somatic*; the rest are called *germline het*.

    chrX variants are silently excluded prior to classification.

    Parameters
    ----------
    df : pd.DataFrame
        Aggregated mutations table.  Required columns: ``mutation``,
        ``pos`` (``"chrN:position"`` format), ``depth``, ``frequency``.
    pathology_purity : float or "NA"
        Tumour purity as a fraction in [0, 1], or the string ``"NA"``
        when purity is unavailable.
    alpha : float
        Significance threshold for the two-sided binomial test (default 0.05).
    pthresh : float
        Purity threshold below which high-frequency variants are called
        germline (default 0.9).  If purity >= pthresh the call stays
        ``"ambiguous"`` because extreme purity makes germline and somatic
        homozygous indistinguishable without CNA data.

    Returns
    -------
    pd.DataFrame
        Result table with columns: ``mutation``, ``pos``, ``depth``,
        ``frequency``, ``zygosity``, ``germline/somatic``.
    """
    df = df.copy()

    # Remove chrX variants – sex chromosomes require allele-dosage correction
    # that is outside the scope of the basic (no-CNA) model.
    df["chr"] = df["pos"].str.split(":").str[0].str.replace("chr", "", regex=False)
    df = df[df["chr"] != "X"].copy()

    freq  = df["frequency"].astype(float).clip(upper=1.0)
    depth = df["depth"].astype(float)

    df["frequency"] = freq
    df["depth"]     = depth

    # Initialise output columns.
    df["zygosity"]         = "ambiguous"
    df["germline/somatic"] = "ambiguous"

    # ------------------------------------------------------------------
    # CASE 1: AF > 0.95 — likely homozygous
    # ------------------------------------------------------------------
    # An AF above 0.95 is consistent with a homozygous state.  When purity
    # is available and low (< pthresh) the read mixture is mostly normal
    # diploid, so a near-100 % AF must come from a germline homozygous SNP.
    high_freq = freq > 0.95
    df.loc[high_freq, "zygosity"] = "homozygous"

    if pathology_purity != "NA":
        germline_mask = high_freq & (pathology_purity < pthresh)
        df.loc[germline_mask, "germline/somatic"] = "germline"

    # ------------------------------------------------------------------
    # CASE 2: AF <= 0.95 — test against diploid heterozygous expectation
    # ------------------------------------------------------------------
    # Under a diploid germline heterozygous model the expected AF is 0.5.
    # Variants with AF significantly different from 0.5 (p < alpha) are
    # inconsistent with that model and are called somatic.
    low_freq     = ~high_freq
    pvals        = compute_pvalues(freq[low_freq], depth[low_freq])
    somatic_mask = pvals < alpha
    low_freq_idx = df.index[low_freq]

    df.loc[low_freq_idx[somatic_mask],  "germline/somatic"] = "somatic"
    df.loc[low_freq_idx[~somatic_mask], "germline/somatic"] = "germline"
    df.loc[low_freq_idx[~somatic_mask], "zygosity"]         = "het"

    # ------------------------------------------------------------------
    # Format and select output columns
    # ------------------------------------------------------------------
    out_df = df[["mutation", "pos", "depth", "frequency",
                 "germline/somatic"]].copy()
    out_df["depth"]     = out_df["depth"].astype(int)
    out_df["frequency"] = out_df["frequency"].round(2)

    return out_df


# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    """Build the CLI argument parser and parse ``sys.argv``.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "aggregated_mutations_file",
        help="Full aggregated mutations text file (tab-separated)",
    )
    parser.add_argument(
        "-f",
        dest="pathology_purity_file",
        default=None,
        metavar="PURITY_FILE",
        help="File containing percent tumor nuclei (single numeric value)",
    )
    parser.add_argument(
        "-a",
        dest="alpha",
        type=float,
        default=0.05,
        metavar="ALPHA",
        help="Significance threshold for binomial test (default: 0.05)",
    )
    parser.add_argument(
        "-p",
        dest="pthresh",
        type=float,
        default=0.9,
        metavar="PTHRESH",
        help="Purity threshold for high-frequency germline call (default: 0.9)",
    )
    parser.add_argument(
        "-o",
        dest="output_header",
        default="",
        metavar="OUT_PREFIX",
        help="Output file prefix (default: input filename stem)",
    )
    parser.add_argument(
        "--panel-size",
        dest="panel_size_mb",
        type=float,
        default=None,
        metavar="MB",
        help=(
            "Effective panel size in megabases for TMB calculation.  "
            "If omitted, TMB is not computed."
        ),
    )
    parser.add_argument(
        "--version",
        action="version",
        version=__version__,
    )

    return parser.parse_args()


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def main() -> int:
    """Run the basic SGZ pipeline end-to-end.

    Reads the aggregated mutations file, runs :func:`run_sgz`, writes the
    per-variant classification file (``*.basic.sgz.txt``), and optionally
    computes TMB and writes ``*.tmb.txt``.

    Returns
    -------
    int
        Exit code (0 = success).
    """
    args = parse_args()

    fname      = args.aggregated_mutations_file
    stem, _    = os.path.splitext(fname)
    out_header = args.output_header if args.output_header else stem
    outfile    = f"{out_header}.basic.sgz.txt"

    pathology_purity = read_pathology_purity(args.pathology_purity_file)
    logger.info("Tumour purity: %s", pathology_purity)

    df = pd.read_csv(fname, sep="\t")
    logger.info("Loaded %d variants from %s", len(df), fname)

    result_df = run_sgz(df, pathology_purity, args.alpha, args.pthresh)

    n_somatic   = (result_df["germline/somatic"] == "somatic").sum()
    n_germline  = (result_df["germline/somatic"] == "germline").sum()
    n_ambiguous = (result_df["germline/somatic"] == "ambiguous").sum()
    logger.info(
        "Classification: %d somatic, %d germline, %d ambiguous",
        n_somatic, n_germline, n_ambiguous,
    )

    result_df.to_csv(outfile, sep="\t", index=False, float_format='%.2f')
    logger.info("SGZ calls written to %s", outfile)

    # ------------------------------------------------------------------
    # Tumor Mutation Burden (TMB)
    # ------------------------------------------------------------------
    if args.panel_size_mb is not None:
        tmb = compute_tmb(result_df, args.panel_size_mb)
        tmb_file = f"{out_header}.tmb.txt"
        with open(tmb_file, "w") as fout:
            fout.write("sample\tn_somatic\tpanel_size_mb\tTMB\n")
            fout.write(
                f"{out_header}\t{n_somatic}\t{args.panel_size_mb}\t{tmb:.4f}\n"
            )
        logger.info(
            "TMB = %.4f mut/Mb (%d somatic / %.2f Mb panel)",
            tmb, n_somatic, args.panel_size_mb,
        )
        logger.info("TMB written to %s", tmb_file)

    return 0


# ---------------------------------------------------------------------------
# Script entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    logging.basicConfig(
        format="%(levelname)s: [%(asctime)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S %Z",
    )
    logger.setLevel(logging.INFO)
    sys.exit(main())
