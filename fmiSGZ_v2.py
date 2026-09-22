#!/usr/bin/env python3
"""
fmiSGZ_v2 - Foundation Medicine Somatic/Germline Zygosity (SGZ) classifier.

This module implements the FMI SGZ method for evaluating copy-number-aware
somatic vs. germline origin and allelic zygosity of short variants detected
by targeted next-generation sequencing.

Algorithm overview
------------------
For each short variant the algorithm:

1. Locates the copy-number (CN) segment that overlaps the variant locus.
2. Uses the segment's purity estimate (p), total copy number (C), and minor
   allele copy number (M) to compute the expected allele frequencies (EAF)
   under four hypotheses:

   * G1 – germline on the major allele:  AF = (p*M + 1*(1-p)) / (p*C + 2*(1-p))
   * G2 – germline on the minor allele:  AF = (p*(C-M) + 1*(1-p)) / (p*C + 2*(1-p))
   * S1 – somatic on the major allele:   AF = p*M / (p*C + 2*(1-p))
   * S2 – somatic on the minor allele:   AF = p*(C-M) / (p*C + 2*(1-p))

3. Tests each hypothesis with a binomial exact test against the observed
   allele count.
4. Derives a log-odds score  log10(P_germline / P_somatic) and applies
   a series of priority rules (low-frequency subclonal, extreme purity,
   CYP2D6/HLA exclusions, …) to produce a call string and zygosity label.
5. Estimates an *allele burden* (quantitative clonality) and its 95 % CI
   via 1 000-sample bootstrap resampling.

Output
------
Two tab-separated output files are produced for each sample:

* ``*.fmi.sgz.full.txt`` – all intermediate statistics.
* ``*.fmi.sgz.txt``      – compact summary (mutation, call, zygosity).

Usage
-----
::

    fmiSGZ_v2.py [options] aggregated_mutations_file cna_model_file
    fmiSGZ_v2.py [-h | --help]
    fmiSGZ_v2.py [--version]

Dependencies
------------
* Python >= 3.8
* numpy
* scipy >= 1.7   (uses ``scipy.stats.binomtest``, replacing the deprecated
                  ``scipy.stats.binom_test``)
"""

import os
import sys
import csv
import logging
import argparse

import numpy as np
from scipy.stats import binom, binomtest
from scipy.special import logsumexp

__author__  = 'James Sun'
__contact__ = 'jsun@foundationmedicine.com'
__version__ = '2.0.0'
__doc__     = 'FMI SGZ method to evaluate CNA origin and zygosity.'

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# CNA model helpers
# ---------------------------------------------------------------------------

def cn2lr_bl(p: float, bl: float, cn: float) -> float:
    """Convert copy number to log-ratio relative to a baseline.

    Computes the expected log2 ratio of a segment given tumor purity,
    the sample-level baseline (diploid) log-ratio, and the integer total
    copy number of the segment.

    Parameters
    ----------
    p : float
        Tumor purity (fraction of tumor cells), in [0, 1].
    bl : float
        Diploid baseline log-ratio for the sample (accounts for overall
        ploidy and stromal contamination).
    cn : float
        Integer total copy number of the segment in tumor cells.

    Returns
    -------
    float
        Expected log2 ratio: log2((p*cn + 2*(1-p)) / bl).
    """
    return np.log2((p * cn + 2 * (1 - p)) / bl)


# ---------------------------------------------------------------------------
# Core SGZ algorithm
# ---------------------------------------------------------------------------

def core_SGZ(data_CNA: list, short_variants: list) -> list:
    """Classify each short variant as somatic or germline using the SGZ model.

    For every variant in *short_variants* the function:

    * skips sex-chromosome variants (chrX / chrY),
    * finds the overlapping CN segment from *data_CNA*,
    * computes expected allele frequencies under germline (G1, G2) and
      somatic (S1, S2) hypotheses,
    * performs binomial exact tests for each hypothesis,
    * applies a decision tree to produce a call, zygosity, clonality, and
      in-tumor labels,
    * estimates allele burden and its bootstrap 95 % CI.

    Parameters
    ----------
    data_CNA : list of Segment
        Copy-number segments loaded from the CNA model file.  Each element
        must expose the attributes ``chr_``, ``segStart``, ``segEnd``,
        ``purity``, ``numMAtumorPred``, ``CN``, ``baseLevel``, ``mafPred``,
        ``segMAF``, and ``segLR``.
    short_variants : list of SV
        Short variants to classify.  Each element must expose ``chr_``,
        ``position``, ``frequency``, ``depth``, and ``mutation``.

    Returns
    -------
    list of dict
        One result dictionary per variant (sex-chromosome variants are
        silently skipped).  Each dictionary contains:

        ============  ======================================================
        Key           Description
        ============  ======================================================
        mutation      Variant identifier string
        position      Genomic position "chrX:pos"
        depth         Sequencing depth (string)
        AF_obs        Observed allele frequency (string, 2 d.p.)
        AF_E[G]       [AF_G1, AF_G2] – expected germline AFs
        AF_E[S]       [AF_S1, AF_S2] – expected somatic AFs
        purity        Tumor purity estimate
        CN            Total copy number
        M             Minor (or best-fit) allele copy number
        errMAF        |model_MAF - data_MAF|
        errLR         |model_LR - data_LR|
        P(G)          max binomial p-value under germline hypotheses
        P(S)          max binomial p-value under somatic hypotheses
        log(G/S)      log10(P_G / P_S) log-odds score
        call          SGZ classification string
        zygosity      Zygosity label (het / homozygous / homoDel / …)
        clonality     'clonal' | 'subclonal' | 'NA'
        in_tumor      'in tumor' | 'not in tumor' | 'NA'
        burden        Allele burden (float)
        burdenCI      [lower, upper] 95 % bootstrap CI
        ============  ======================================================
    """
    ALPHA = 0.01

    results = []

    for sv in short_variants:
        # Sex chromosomes are excluded from the SGZ model.
        if sv.chr_.upper() in ('X', 'Y'):
            continue

        # Clamp observed frequency to [0, 1].
        if sv.frequency > 1:
            sv.frequency = 1.0

        p = M = C = np.nan

        seg_found = False
        clonality = 'NA'
        in_tumor  = 'NA'

        for seg in data_CNA:
            if seg.chr_ == sv.chr_ and seg.segStart <= sv.position <= seg.segEnd:
                seg_found = True

                p  = seg.purity
                M  = seg.numMAtumorPred
                C  = seg.CN
                bL = seg.baseLevel

                model_MAF = seg.mafPred
                data_MAF  = seg.segMAF
                error_MAF = abs(model_MAF - data_MAF)

                model_LR  = cn2lr_bl(p, bL, C)
                data_LR   = seg.segLR
                error_LR  = abs(model_LR - data_LR)

                # Ensure M refers to the major allele copy number.
                if M == 0:
                    M = C - M

                P_G1 = P_G2 = P_S1 = P_S2 = 0
                AF_G1 = AF_G2 = AF_S1 = AF_S2 = np.nan

                # Skip binomial tests if key CNA values are missing.
                if np.isnan(p) or np.isnan(C) or np.isnan(M):
                    SG_status = 'ambiguous_CNA_model'
                    break

                # Expected AFs under G1 and S1 (major allele hypotheses).
                AF_G1 = (p * M         + 1 * (1 - p)) / (p * C + 2 * (1 - p))
                AF_S1 = (p * M         + 0 * (1 - p)) / (p * C + 2 * (1 - p))

                obs_count = int(sv.depth * sv.frequency + 0.5)
                P_G1 = binomtest(obs_count, int(sv.depth), AF_G1).pvalue
                P_S1 = binomtest(obs_count, int(sv.depth), AF_S1).pvalue

                # Minor allele hypotheses (G2 / S2) only if M != C-M.
                if M != C - M:
                    AF_G2 = (p * (C - M) + 1 * (1 - p)) / (p * C + 2 * (1 - p))
                    P_G2  = binomtest(obs_count, int(sv.depth), AF_G2).pvalue

                    if C - M != 0:
                        AF_S2 = (p * (C - M) + 0 * (1 - p)) / (p * C + 2 * (1 - p))
                        P_S2  = binomtest(obs_count, int(sv.depth), AF_S2).pvalue
                    else:
                        AF_S2 = P_S2 = np.nan
                else:
                    AF_G2 = AF_S2 = P_G2 = P_S2 = np.nan

                # Best p-values across the two allele hypotheses.
                max_prob_germline = np.nanmax([P_G1, P_G2])
                max_prob_somatic  = np.nanmax([P_S1, P_S2])

                # Log-odds: positive → germline favoured.
                if max_prob_somatic == 0 and max_prob_germline == 0:
                    logodds = np.nan
                elif max_prob_somatic == 0:
                    logodds = np.inf
                elif max_prob_germline == 0:
                    # p_germline underflowed to 0; compute log10(P_G) in log-space
                    # to avoid -inf.  Use logsumexp over logpmf values in the tail.
                    _log_pG_list = []
                    for _af, _p in [(AF_G1, P_G1), (AF_G2, P_G2)]:
                        if np.isnan(_af):
                            continue
                        _k = obs_count
                        _n = int(sv.depth)
                        _expected = _n * _af
                        if _k <= _expected:
                            _logpmfs = np.array([binom.logpmf(j, _n, _af)
                                                 for j in range(_k + 1)])
                        else:
                            _logpmfs = np.array([binom.logpmf(j, _n, _af)
                                                 for j in range(_k, _n + 1)])
                        _log_pG_list.append(logsumexp(_logpmfs))
                    log10_pG = (np.log(2) + max(_log_pG_list)) / np.log(10)
                    logodds = log10_pG - np.log10(max_prob_somatic)
                else:
                    logodds = np.log10(max_prob_germline) - np.log10(max_prob_somatic)

                # ----------------------------------------------------------------
                # Decision tree
                # ----------------------------------------------------------------
                if 'CYP2D6' in sv.mutation:
                    SG_status = 'nocall_CYP2D6'

                elif 'HLA' in sv.mutation:
                    SG_status = 'nocall_HLA'

                elif sv.frequency < 0.05 and 0.2 < p < 0.9 and sv.status in ('high', 'moderate'):
                    SG_status = 'subclonal somatic'
                    clonality = 'subclonal'
                    in_tumor  = 'in tumor'

                elif np.isnan(C) or np.isnan(M):
                    SG_status = 'ambiguous_CNA_model'

                elif p > 0.95:
                    SG_status = 'nocall_purity>95%'

                elif sv.frequency > 0.95:
                    SG_status = 'germline'
                    clonality = 'clonal'
                    in_tumor  = 'in tumor'

                elif error_MAF > 0.06 or error_LR > 0.4:
                    SG_status = 'ambiguous_CNA_model'

                elif max_prob_germline > ALPHA and max_prob_somatic < ALPHA:
                    SG_status = 'probable germline' if logodds < 2 else 'germline'
                    clonality = 'clonal'

                    if np.nanargmax([P_G1, P_G2]) == 0:
                        in_tumor = 'in tumor'
                    else:
                        M = C - M
                        if M == 0:
                            in_tumor = 'not in tumor'

                elif max_prob_germline < ALPHA and max_prob_somatic > ALPHA:
                    SG_status = 'probable somatic' if logodds > -2 else 'somatic'
                    clonality = 'clonal'
                    in_tumor  = 'in tumor'

                    if np.nanargmax([P_S1, P_S2]) == 1:
                        M = C - M

                elif max_prob_germline > ALPHA and max_prob_somatic > ALPHA:
                    SG_status = 'ambiguous_both_G_and_S'

                elif max_prob_germline < ALPHA and max_prob_somatic < ALPHA:
                    min_soma_EAF = np.nanmin([AF_S1, AF_S2])
                    min_germ_EAF = np.nanmin([AF_G1, AF_G2])

                    if (p >= 0.3 and sv.frequency < 0.25
                            and sv.frequency < min_soma_EAF / 1.5
                            and min_soma_EAF <= min_germ_EAF):
                        SG_status = 'subclonal somatic'
                        clonality = 'subclonal'
                        in_tumor  = 'in tumor'

                    elif (p >= 0.3 and sv.frequency < 0.25
                            and sv.frequency < min_germ_EAF / 2.0
                            and min_germ_EAF < min_soma_EAF):
                        SG_status = 'subclonal somatic'
                        clonality = 'subclonal'
                        in_tumor  = 'in tumor'

                    elif logodds < -5 and max_prob_somatic > 1e-10:
                        SG_status = 'somatic'
                        clonality = 'clonal'
                        in_tumor  = 'in tumor'

                        if np.nanargmax([P_S1, P_S2]) == 1:
                            M = C - M

                    elif logodds > 5 and max_prob_germline > 1e-4:
                        SG_status = 'germline'
                        clonality = 'clonal'

                        if np.nanargmax([P_G1, P_G2]) == 0:
                            in_tumor = 'in tumor'
                        else:
                            M = C - M
                            if M == 0:
                                in_tumor = 'not in tumor'

                    else:
                        SG_status = 'ambiguous_neither_G_nor_S'

                else:
                    SG_status = 'unknown'

                break  # Stop after the first matching segment.

        # --------------------------------------------------------------------
        # Fallback when no segment overlaps the variant.
        # --------------------------------------------------------------------
        if not seg_found:
            AF_G1 = AF_G2 = AF_S1 = AF_S2 = np.nan
            p = C = M = error_LR = error_MAF = np.nan
            max_prob_germline = max_prob_somatic = logodds = np.nan
            SG_status = 'nocall_segmentMissing'

        # --------------------------------------------------------------------
        # Zygosity assignment
        # --------------------------------------------------------------------
        zygosity = 'NA'
        if np.isnan(C) or np.isnan(M) or p < 0.19:
            zygosity = 'NA'
            if 'germline' in SG_status:
                in_tumor = 'NA'
        elif C == 0 and M == 0: zygosity = 'homoDel'
        elif C == M and M == 1: zygosity = 'homozygous'
        elif C == M and M >= 2: zygosity = 'homozygous'
        elif C >= 1 and M == 0: zygosity = 'not in tumor'
        elif C != M and M != 0: zygosity = 'het'

        # --------------------------------------------------------------------
        # Allele burden (quantitative clonality) with 95 % bootstrap CI
        # --------------------------------------------------------------------
        allele_burden = np.nan
        burden_CI     = [np.nan, np.nan]

        if 'germline' in SG_status or 'somatic' in SG_status:
            if 'germline' in SG_status:
                EAF = (p * M + 1 - p) / (p * C + 2 * (1 - p))
            else:
                EAF = p * M / (p * C + 2 * (1 - p))

            allele_burden = sv.frequency / EAF

            try:
                # Bootstrap 1 000 resamples to obtain the 95 % CI.
                EAF_sim   = binom.rvs(int(sv.depth), EAF, size=1000) / sv.depth
                ab_sim    = sv.frequency / EAF_sim
                burden_CI = np.percentile(ab_sim, [2.5, 97.5])
            except Exception as exc:
                logging.warning('Problem computing EAF_sim for depth=%s EAF=%s: %s',
                                sv.depth, EAF, exc)

        result = {
            'mutation' : sv.mutation,
            'position' : 'chr%s:%d' % (sv.chr_, sv.position),
            'depth'    : '%d' % sv.depth,
            'AF_obs'   : '%0.2f' % sv.frequency,
            'status'   : sv.status,
            'AF_E[G]'  : [AF_G1, AF_G2],
            'AF_E[S]'  : [AF_S1, AF_S2],
            'purity'   : p,
            'CN'       : C,
            'M'        : M,
            'errMAF'   : error_MAF,
            'errLR'    : error_LR,
            'P(G)'     : max_prob_germline,
            'P(S)'     : max_prob_somatic,
            'log(G/S)' : logodds,
            'call'     : SG_status,
            'zygosity' : zygosity,
            'clonality': clonality,
            'in_tumor' : in_tumor,
            'burden'   : allele_burden,
            'burdenCI' : burden_CI,
        }

        results.append(result)

    return results


# ---------------------------------------------------------------------------
# TMB
# ---------------------------------------------------------------------------

def compute_tmb(results: list, panel_size_mb: float) -> float:
    """Compute Tumor Mutation Burden (TMB) from classified variants.

    TMB is defined as the number of somatic variants divided by the
    effective genomic panel size in megabases (mut/Mb).  Variants whose
    ``call`` field contains the substring ``"somatic"`` are counted, which
    covers the categories ``"somatic"``, ``"probable somatic"``, and
    ``"subclonal somatic"``.

    Parameters
    ----------
    results : list of dict
        Output of :func:`core_SGZ`.  Each dict must contain a ``"call"``
        key with the SGZ classification string.
    panel_size_mb : float
        Effective size of the sequenced/analysed genomic region in Mb
        (e.g. 0.8 for the Foundation Medicine CDx panel, ~30 for WES).

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

    n_somatic = sum(1 for sv in results if "somatic" in sv.get("call", ""))
    return float(n_somatic) / panel_size_mb


# ---------------------------------------------------------------------------
# File readers
# ---------------------------------------------------------------------------

def read_cna_model_file(fname: str) -> list:
    """Read a tab-separated CNA model file into a list of Segment objects.

    The file must have a header row with at least the following columns:
    ``CHR``, ``segStart``, ``segEnd``, ``mafPred``, ``CN``, ``segLR``,
    ``segMAF``, ``numMAtumorPred``, ``purity``, ``baseLevel``.

    Parameters
    ----------
    fname : str
        Path to the CNA model text file.

    Returns
    -------
    list of Segment
        One :class:`Segment` object per data row.

    Raises
    ------
    FileNotFoundError
        If *fname* does not exist.
    KeyError
        If a required column is absent from the file header.
    """
    data_CNA = []
    with open(fname, 'r', newline='') as fin:
        reader = csv.DictReader(fin, dialect='excel-tab')
        for line in reader:
            seg = Segment(
                line['CHR'], line['segStart'], line['segEnd'],
                line['mafPred'], line['CN'], line['segLR'],
                line['segMAF'], line['numMAtumorPred'],
                line['purity'], line['baseLevel'],
            )
            data_CNA.append(seg)
    return data_CNA


def read_mut_aggr_full(fname: str) -> list:
    """Read a tab-separated aggregated mutations file into a list of SV objects.

    The file must have a header row containing at least the columns
    ``mutation``, ``frequency``, ``depth``, and ``pos``.

    Parameters
    ----------
    fname : str
        Path to the aggregated mutations text file.

    Returns
    -------
    list of SV
        One :class:`SV` object per data row.

    Raises
    ------
    FileNotFoundError
        If *fname* does not exist.
    KeyError
        If a required column is absent from the file header.
    """
    short_variants = []
    with open(fname, 'r', newline='') as fin:
        reader = csv.DictReader(fin, dialect='excel-tab')
        for line in reader:
            status = line.get('status', '')
            sv = SV(line['mutation'], line['frequency'], line['depth'], line['pos'], status)
            short_variants.append(sv)
    return short_variants


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

class Segment:
    """A copy-number segment from the CNA model.

    Attributes
    ----------
    chr_ : str
        Chromosome name with the 'chr' prefix stripped (e.g. ``'1'``, ``'X'``).
    segStart : int
        Genomic start coordinate of the segment (1-based, inclusive).
    segEnd : int
        Genomic end coordinate of the segment (1-based, inclusive).
    mafPred : float or nan
        Model-predicted minor-allele frequency for the segment.
    CN : int or nan
        Total integer copy number in tumour cells.
    segLR : float or nan
        Observed segment log-ratio.
    segMAF : float or nan
        Observed segment minor-allele frequency.
    numMAtumorPred : float or nan
        Predicted number of minor-allele copies in tumour cells.
    purity : float
        Tumour purity estimate (fraction of tumour cells).
    baseLevel : float
        Diploid baseline log-ratio for the sample.
    """

    def __init__(self, chr_: str, segStart: str, segEnd: str,
                 mafPred: str, CN: str, segLR: str, segMAF: str,
                 numMAtumorPred: str, purity: str, baseLevel: str) -> None:
        super().__init__()
        self.chr_           = chr_.strip('chr')
        self.segStart       = int(segStart)
        self.segEnd         = int(segEnd)
        self.mafPred        = float(mafPred)        if mafPred        != 'NA' else np.nan
        self.CN             = int(float(CN))        if CN             != 'NA' else np.nan
        self.segLR          = float(segLR)          if segLR          != 'NA' else np.nan
        self.segMAF         = float(segMAF)         if segMAF         != 'NA' else np.nan
        self.numMAtumorPred = float(numMAtumorPred) if numMAtumorPred != 'NA' else np.nan
        self.purity         = float(purity)         if purity         != 'NA' else np.nan
        self.baseLevel      = float(baseLevel)


class SV:
    """A short variant (SNV or small indel) from the aggregated mutations file.

    Attributes
    ----------
    mutation : str
        Variant identifier (e.g. HGVS notation or gene name).
    frequency : float
        Observed variant allele frequency (VAF), in [0, 1].
    depth : float
        Total sequencing depth at the variant locus.
    chr_ : str
        Chromosome name with the 'chr' prefix stripped.
    position : int
        Genomic position of the variant.
    """

    def __init__(self, mutation: str, frequency: str, depth: str, position: str,
                 status: str = '') -> None:
        super().__init__()
        self.mutation  = mutation
        self.frequency = float(frequency)
        self.depth     = float(depth)
        self.status    = status.lower()

        parts          = position.split(':')
        self.chr_      = parts[0].strip('chr')
        self.position  = int(parts[1])


# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------

def _arg_parser() -> argparse.ArgumentParser:
    """Build and return the command-line argument parser.

    Returns
    -------
    argparse.ArgumentParser
        Configured parser ready for ``parse_args()``.
    """
    script_version     = globals().get('__version__')
    script_description = globals().get('__doc__')
    script_usage = (
        '%(prog)s [options] aggregated_mutations_file cna_model_file\n'
        '       %(prog)s [-h|--help]\n'
        '       %(prog)s [--version]'
    )

    parser = argparse.ArgumentParser(
        usage=script_usage,
        description=script_description,
        add_help=False,
    )

    g = parser.add_argument_group('Program Help')
    g.add_argument('-h', '--help', action='help',
                   help='show this help message and exit')
    g.add_argument('--version', action='version', version=script_version)

    g = parser.add_argument_group('Required Arguments')
    g.add_argument('aggregated_mutations_file', action='store', type=str,
                   help='Full aggregated mutations text file')
    g.add_argument('cna_model_file', action='store', type=str,
                   help='CNA model calls')

    g = parser.add_argument_group('Optional Arguments')
    g.add_argument('-o', dest='output_header', action='store', type=str,
                   default='',
                   help='Output file folder and prefix (default: input filename stem)')
    g.add_argument('--panel-size', dest='panel_size_mb', action='store',
                   type=float, default=None, metavar='MB',
                   help=(
                       'Effective panel size in megabases for TMB calculation.  '
                       'If omitted, TMB is not computed.  '
                       'Example: 0.8 for Foundation Medicine CDx, ~30 for WES.'
                   ))

    return parser


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def main(args: argparse.Namespace) -> int:
    """Run the SGZ pipeline end-to-end.

    Reads the CNA model and mutations files, runs :func:`core_SGZ`, and
    writes two output files (full statistics and compact summary).

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments with attributes:

        * ``aggregated_mutations_file`` – path to the mutations file.
        * ``cna_model_file``           – path to the CNA model file.
        * ``output_header``            – optional output path prefix.
        * ``panel_size_mb``            – effective panel size in Mb for TMB
                                         (``None`` disables TMB output).

    Returns
    -------
    int
        Exit code (0 = success).
    """
    np.random.seed(10)

    fname_muts    = args.aggregated_mutations_file
    fname_CNmodel = args.cna_model_file
    out_header    = args.output_header

    if not out_header:
        stem, _ = os.path.splitext(fname_muts)
        out_header = stem

    sgz_full = out_header + '.fmi.sgz.full.txt'
    sgz_out  = out_header + '.fmi.sgz.txt'

    # If the CNA model file is empty no segment-level data are available.
    # Produce empty output files and exit cleanly.
    if not os.path.getsize(fname_CNmodel):
        logger.warning('Empty cna_calls.txt file.  Not running SGZ code.')
        for filename in (sgz_full, sgz_out):
            open(filename, 'w').close()
        return 0

    # 1) Read CN model file
    data_CNA = read_cna_model_file(fname_CNmodel)
    # 2) Read mutations file
    short_variants = read_mut_aggr_full(fname_muts)
    # 3) Run SGZ
    y_obj = core_SGZ(data_CNA, short_variants)

    # 4) Write output files
    full_header = [
        'mutation', 'pos', 'depth', 'frequency', 'status',
        'afG1', 'afS1', 'afG2', 'afS2',
        'p', 'C', 'M', 'logOR_G',
        'clonality', 'clonality_CI_low', 'clonality_CI_high',
        'germline/somatic', 'zygosity',
    ]
    summary_header = ['mutation', 'pos', 'depth', 'frequency', 'status', 'C', 'germline/somatic', 'zygosity']

    with open(sgz_full, 'w', newline='') as fout, \
         open(sgz_out,  'w', newline='') as fout2:

        fout.write('\t'.join(full_header)    + '\n')
        fout2.write('\t'.join(summary_header) + '\n')

        for sv in y_obj:
            full_row = [
                sv['mutation'],
                sv['position'],
                sv['depth'],
                sv['AF_obs'],
                sv['status'],
                '%.2f' % sv['AF_E[G]'][0],
                '%.2f' % sv['AF_E[S]'][0],
                '%.2f' % sv['AF_E[G]'][1],
                '%.2f' % sv['AF_E[S]'][1],
                '%.3f' % sv['purity'],
                '%g'   % sv['CN'],
                '%g'   % sv['M'],
                '%.1f' % sv['log(G/S)'],
                '%.2f' % sv['burden'],
                '%.2f' % sv['burdenCI'][0],
                '%.2f' % sv['burdenCI'][1],
                sv['call'],
                sv['zygosity'],
            ]
            fout.write('\t'.join(full_row) + '\n')

            summary_row = [
                sv['mutation'],
                sv['position'],
                sv['depth'],
                sv['AF_obs'],
                sv['status'],
                '%g' % sv['CN'],
                sv['call'],
                sv['zygosity'],
            ]
            fout2.write('\t'.join(summary_row) + '\n')

    logger.info("SGZ calls written to %s and %s", sgz_out, sgz_full)

    # ------------------------------------------------------------------
    # Tumor Mutation Burden (TMB)
    # ------------------------------------------------------------------
    if getattr(args, 'panel_size_mb', None) is not None:
        tmb = compute_tmb(y_obj, args.panel_size_mb)
        n_somatic = sum(1 for sv in y_obj if 'somatic' in sv.get('call', ''))
        tmb_file  = out_header + '.tmb.txt'
        with open(tmb_file, 'w', newline='') as ftmb:
            ftmb.write('sample\tn_somatic\tpanel_size_mb\tTMB\n')
            ftmb.write('%s\t%d\t%.2f\t%.4f\n' % (
                out_header, n_somatic, args.panel_size_mb, tmb,
            ))
        logger.info(
            'TMB = %.4f mut/Mb (%d somatic / %.2f Mb panel)',
            tmb, n_somatic, args.panel_size_mb,
        )
        logger.info('TMB written to %s', tmb_file)

    return 0


# ---------------------------------------------------------------------------
# Script entry point
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    logging.basicConfig(
        format='%(levelname)s: [%(asctime)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S %Z',
    )
    logger.setLevel(logging.INFO)

    args = _arg_parser().parse_args(sys.argv[1:])

    if getattr(args, 'debug', False):
        logger.setLevel(logging.DEBUG)

    sys.exit(main(args))
