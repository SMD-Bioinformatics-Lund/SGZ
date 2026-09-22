# SGZ — Somatic/Germline Zygosity Classifiers

**Version:** 2.0.2

This repository contains two complementary tools for classifying short variants
(SNVs and small indels) detected by targeted next-generation sequencing as
**somatic** or **germline**, and for computing **Tumor Mutation Burden (TMB)**.

| Script | Model | Requires CNA data | TMB |
|--------|-------|:-----------------:|:---:|
| `basicSGZ_v2.py` | Binomial test against AF = 0.5 | No | Yes |
| `fmiSGZ_v2.py`  | Copy-number-aware FMI SGZ model | Yes | Yes |

Use `basicSGZ_v2.py` when no CNA model is available (quick screen).  Use
`fmiSGZ_v2.py` when a CNA model is available for maximum accuracy.

Two helper scripts convert upstream CNA caller output into the CNA model format
required by `fmiSGZ_v2.py`:

| Script | Source tool | Panel-of-normals support |
|--------|-------------|:------------------------:|
| `cnvkit_to_cna_model.py` | CNVkit (`.cns`) | Yes (via CNVkit PoN) |
| `gatk_to_cna_model.py`   | GATK4 ModelSegments (`.modelFinal.seg`) | Yes (via GATK4 PON) |

A VCF converter generates the mutations input file directly from a
VEP-annotated VCF:

| Script | Input | Output |
|--------|-------|--------|
| `vcf_to_mutations.py` | VEP-annotated VCF (VarDict / any caller) | `mut_aggr.full.txt` |

---

## Table of Contents

1. [Requirements](#requirements)
2. [Quick Start](#quick-start)
3. [Preparing Mutation Input from a VCF — `vcf_to_mutations.py`](#preparing-mutation-input-from-a-vcf--vcf_to_mutationspy)
4. [Generating a CNA Model File](#generating-a-cna-model-file)
   - [From CNVkit](#from-cnvkit-cnvkit_to_cna_modelpy)
   - [From GATK4 ModelSegments](#from-gatk4-modelsegments-gatk_to_cna_modelpy)
5. [basicSGZ\_v2 — Basic Classifier](#basicsgz_v2--basic-classifier)
6. [fmiSGZ\_v2 — FMI CNA-Aware Classifier](#fmisgz_v2--fmi-cna-aware-classifier)
7. [Tumor Mutation Burden (TMB)](#tumor-mutation-burden-tmb)
8. [Input File Formats](#input-file-formats)
9. [Output File Formats](#output-file-formats)
10. [Call Strings Reference](#call-strings-reference)
11. [Module API](#module-api)
12. [Regression Tests](#regression-tests)
13. [Changelog](#changelog)

---

## Requirements

| Package | Version | Used by |
|---------|---------|---------|
| Python  | ≥ 3.8 (≥ 3.13 for `vcf_to_mutations.py`) | all scripts |
| numpy   | any recent | `basicSGZ_v2`, `fmiSGZ_v2` |
| pandas  | any recent | `basicSGZ_v2` only |
| scipy   | ≥ 1.7 (requires `scipy.stats.binomtest`) | `basicSGZ_v2`, `fmiSGZ_v2` |
| pysam   | ≥ 0.22 | `vcf_to_mutations.py` only |

```bash
pip install numpy pandas "scipy>=1.7" pysam
```

The CNA converter scripts (`cnvkit_to_cna_model.py`, `gatk_to_cna_model.py`)
use only the Python standard library — no additional packages required.

### Which `python3` runs this?

If your system has more than one Python 3 installed (common on shared
clusters and Ubuntu/Debian systems, where the default `python3` and
`python3.13` can be entirely separate installations with separate
`site-packages`), `pip install` and the script invocation must use the
**same** interpreter, or you'll hit `ModuleNotFoundError: No module named
'numpy'` even though the packages are "installed."  Before running anything,
confirm:

```bash
python3 -c "import numpy, pandas, scipy, pysam; print('OK')"
```

If that fails, either install into that interpreter explicitly
(`python3 -m pip install ...`) or invoke the scripts with the interpreter
that already has the packages (e.g. `python3.13 basicSGZ_v2.py ...`) — the
examples in this README use `python3`/`python3.13` as a shorthand; substitute
whichever interpreter on your system actually satisfies the requirements
above.  A dedicated virtual environment avoids the ambiguity entirely:

```bash
python3.13 -m venv .venv
source .venv/bin/activate
pip install numpy pandas "scipy>=1.7" pysam
```

---

## Quick Start

```bash
# Step 0 — Generate mutations file from a VEP-annotated VCF
python3.13 vcf_to_mutations.py sample.final.filtered.vcf \
    -s SAMPLE_NAME \
    -o sample.mut_aggr.full.txt

# Basic classifier (no CNA data needed)
python3 basicSGZ_v2.py sample.mut_aggr.full.txt

# Basic classifier + purity + TMB on a 1 Mb panel
python3 basicSGZ_v2.py sample.mut_aggr.full.txt \
    -f sample.purity.txt \
    --panel-size 1.0 \
    -o results/sample

# FMI CNA-aware classifier (CNA model file already prepared)
python3 fmiSGZ_v2.py sample.mut_aggr.full.txt sample.cna_calls.txt

# FMI classifier + TMB on MSK-IMPACT 468-gene panel (~1.2 Mb)
python3 fmiSGZ_v2.py sample.mut_aggr.full.txt sample.cna_calls.txt \
    --panel-size 1.2 \
    -o results/sample

# Full pipeline — CNVkit → CNA model → fmiSGZ
python3 cnvkit_to_cna_model.py sample.cns sample.cna_calls.txt --purity 0.72
python3 fmiSGZ_v2.py sample.mut_aggr.full.txt sample.cna_calls.txt \
    --panel-size 1.2 -o results/sample

# Full pipeline — GATK4 → CNA model → fmiSGZ
python3 gatk_to_cna_model.py sample.modelFinal.seg sample.cna_calls.txt \
    --purity 0.72 --ploidy 2
python3 fmiSGZ_v2.py sample.mut_aggr.full.txt sample.cna_calls.txt \
    --panel-size 1.2 -o results/sample
```

---

## Preparing Mutation Input from a VCF — `vcf_to_mutations.py`

### Overview

`vcf_to_mutations.py` converts a **VEP-annotated VCF** file (produced by any
variant caller — VarDict, FreeBayes, GATK Mutect2, etc. — annotated with
Ensembl VEP) into the `mut_aggr.full.txt` tab-separated format consumed by
`basicSGZ_v2.py` and `fmiSGZ_v2.py`.

For each variant the script selects the **MANE SELECT** transcript and
extracts:

- Mutation identifier string (gene, NM accession, HGVSc, HGVSp in 1-letter
  amino-acid notation)
- Genomic position (`chrN:POS`)
- Allele frequency (`VAF` FORMAT field)
- Depth (`VD` FORMAT field — ALT allele observation count)
- Strand and effect category

It requires **Python ≥ 3.13** and **pysam**.

### Transcript selection priority

1. **MANE SELECT** transcript (`MANE_SELECT` CSQ field, `NM_*` accession,
   `protein_coding` biotype only)
2. **CANONICAL** (`CANONICAL == YES`, `protein_coding` only)
3. Variant is skipped if neither is available with an `NM_*` accession

### Mutation string format

| Variant class | Format | Example |
|---------------|--------|---------|
| Missense / nonsense / inframe | `GENE:NM_xxxxx:c.HGVSc_p.HGVSp` | `BRAF:NM_004333:c.1799T>A_p.V600E` |
| Splice (no protein change) | `GENE:NM_xxxxx:c.HGVSc:splice` | `FGFR1:NM_023110:c.1552+2T>C:splice` |
| Frameshift | `GENE:NM_xxxxx:c.HGVSc_p.HGVSp:frameshift` | `TP53:NM_000546:c.817del_p.L273Cfs*3:frameshift` |

HGVSp is always converted from 3-letter to 1-letter amino-acid codes
(`Asp1228Glu` → `D1228E`). URL-encoded characters (`%3D` → `=`, `%2A` → `*`)
are decoded automatically.

### Effect categories

| Output label | VEP consequence terms |
|---|---|
| `missense` | `missense_variant` |
| `nonsense` | `stop_gained`, `stop_lost`, `start_lost` |
| `splice` | `splice_donor_variant`, `splice_acceptor_variant`, `splice_region_variant`, `splice_donor_5th_base_variant`, `splice_polypyrimidine_tract_variant` |
| `frameshift` | `frameshift_variant` |
| `inframe_indel` | `inframe_insertion`, `inframe_deletion` |

Synonymous and non-coding variants are excluded by default. Use `--synonymous`
to include synonymous variants.

### Usage

```
vcf_to_mutations.py [options] vcf_file
vcf_to_mutations.py [-h | --help]
```

### Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `vcf` | *(required)* | VEP-annotated VCF file (uncompressed or bgzipped) |
| `-s / --sample NAME` | VCF filename stem | Sample name written in the `#sample` column |
| `-o / --output FILE` | stdout | Output file path |
| `--synonymous` | off | Include synonymous variants (excluded by default) |
| `--min-vaf FLOAT` | `0.5` | Discard variants with VAF below this threshold. Mutually exclusive with `--no-vaf-filter`. |
| `--no-vaf-filter` | off | Disable the VAF filter entirely (keep all VAF values). Mutually exclusive with `--min-vaf`. |

### Examples

**Minimal — write to stdout:**

```bash
python3.13 vcf_to_mutations.py sample.final.filtered.vcf -s MySample
```

**Full run with explicit output:**

```bash
python3.13 vcf_to_mutations.py \
    sample.final.filtered.vcf \
    -s 25PH01170-0101-DNA \
    -o sample.mut_aggr.full.txt
```

**Including synonymous variants:**

```bash
python3.13 vcf_to_mutations.py sample.final.filtered.vcf \
    -s MySample \
    --synonymous \
    -o sample.mut_aggr.full.synonymous.txt
```

**VAF filter — keep only variants with VAF ≥ 0.5 (default):**

```bash
python3.13 vcf_to_mutations.py sample.final.filtered.vcf \
    -s MySample \
    -o sample.mut_aggr.full.txt
# Equivalent (explicit):
python3.13 vcf_to_mutations.py sample.final.filtered.vcf \
    -s MySample \
    --min-vaf 0.5 \
    -o sample.mut_aggr.full.txt
```

**VAF filter — custom threshold (e.g. keep VAF ≥ 0.05):**

```bash
python3.13 vcf_to_mutations.py sample.final.filtered.vcf \
    -s MySample \
    --min-vaf 0.05 \
    -o sample.mut_aggr.full.txt
```

**Disable the VAF filter entirely:**

```bash
python3.13 vcf_to_mutations.py sample.final.filtered.vcf \
    -s MySample \
    --no-vaf-filter \
    -o sample.mut_aggr.full.txt
```

**Full pipeline — VCF → mutations → fmiSGZ with CNA model:**

```bash
# 1. Convert VCF to mutation table
python3.13 vcf_to_mutations.py sample.final.filtered.vcf \
    -s MySample -o MySample.mut_aggr.full.txt

# 2. Convert CNVkit output to CNA model
python3 cnvkit_to_cna_model.py sample.called.cns sample.cna_calls.txt \
    --purity 0.72

# 3. Run FMI SGZ classifier
python3 fmiSGZ_v2.py MySample.mut_aggr.full.txt sample.cna_calls.txt \
    --panel-size 1.2 -o results/MySample
```

### Output format

The output is a tab-separated file identical to the `mut_aggr.full.txt` files
in `data/samples/`:

| Column | Description |
|--------|-------------|
| `#sample` | Sample name (from `-s` or VCF filename stem) |
| `mutation` | Mutation identifier string (see [Mutation string format](#mutation-string-format) above) |
| `frequency` | VAF — ALT allele fraction from the `VAF` FORMAT field |
| `depth` | VD — ALT allele observation count from the `VD` FORMAT field |
| `pos` | Genomic position in `chrN:POS` format (prefixes `chr` if absent) |
| `status` | Always `unknown` (no ClinVar/hotspot classification performed) |
| `strand` | `+` or `-` from the VEP `STRAND` CSQ sub-field |
| `effect` | Short effect label (see [Effect categories](#effect-categories) above) |

**Example output:**

```
#sample              mutation                              frequency  depth  pos               status   strand  effect
25PH01170-0101-DNA   TNFRSF14:NM_003820:c.316C>T_p.R106C  0.0103     2      chr1:2559834      unknown  +       missense
25PH01170-0101-DNA   MTOR:NM_004958:c.6811-5T>A:splice     0.0161     2      chr1:11121373     unknown  -       splice
25PH01170-0101-DNA   ARID1A:NM_006015:c.3580G>T_p.G1194*   0.0105     2      chr1:26772852     unknown  +       nonsense
25PH01170-0101-DNA   SPEN:NM_015001:c.1639del_p.V547Cfs*5:frameshift  0.0116  2  chr1:15920871  unknown  +  frameshift
```

> **Depth column note:** `depth` contains the **VD** field (ALT allele
> observation count), as requested.  The existing reference sample files in
> `data/samples/` use the **DP** field (total depth) in this column.  To
> switch to total depth, replace `sample_data.get("VD")` with
> `sample_data.get("DP")` on line 189 of `vcf_to_mutations.py`.

### VCF requirements

| Requirement | Detail |
|-------------|--------|
| VEP annotation | CSQ INFO field must be present with `Format:` description in the header |
| FORMAT fields | `VAF` (float) and `VD` (integer) must be present per-sample |
| Single-sample VCF | Only the first sample column is read |
| Chromosome style | Both `chr1` and `1` are accepted; `chr` prefix is added if missing |

---

## Generating a CNA Model File

`fmiSGZ_v2.py` requires a CNA model file with per-segment copy number, purity,
and allele fraction information.  Two converter scripts are provided to produce
this file from common upstream callers.

### From CNVkit — `cnvkit_to_cna_model.py`

**Input:** `*.cns` file from `cnvkit.py call` (copy-number segments with integer
CN and allele fraction calls).

**Typical CNVkit pipeline (MSK-IMPACT example):**

```bash
# Build panel-of-normals reference
cnvkit.py batch --normal normal1.bam normal2.bam ... \
    --targets MSK_IMPACT_targets.bed \
    --fasta hg19.fa \
    --output-reference pon_reference.cnn

# Call CNVs on the tumour sample
cnvkit.py batch tumor.bam \
    --reference pon_reference.cnn \
    --output-dir cnvkit_out/

# Optional: call integer copy numbers
cnvkit.py call cnvkit_out/tumor.cns -o tumor.called.cns
```

**Convert to CNA model:**

```bash
# Purity known (recommended)
python3 cnvkit_to_cna_model.py tumor.called.cns sample.cna_calls.txt \
    --purity 0.72

# Purity unknown (defaults to 0.5 — approximate)
python3 cnvkit_to_cna_model.py tumor.called.cns sample.cna_calls.txt
```

**Arguments:**

| Argument | Default | Description |
|----------|---------|-------------|
| `input_cns` | *(required)* | CNVkit `.cns` called-segments file |
| `output_cna` | *(required)* | Output CNA model file for `fmiSGZ_v2.py` |
| `--purity` | `0.5` | Tumour purity in (0, 1] |

**Column mapping (CNVkit → CNA model):**

| CNVkit column | CNA model column | Notes |
|---------------|-----------------|-------|
| `chromosome` | `CHR` | as-is |
| `start + 1` | `segStart` | CNVkit is 0-based BED; convert to 1-based |
| `end` | `segEnd` | BED exclusive end = 1-based inclusive |
| `log2` | `segLR` | observed log2-ratio |
| `baf` | `segMAF` | minor allele fraction; empty → `NA` |
| `cn` | `CN` | total integer copy number |
| `cn1` | `numMAtumorPred` | minor-allele CN; empty → `NA` |
| `NA` | `mafPred` | requires purity × CN calculation; always `NA` |
| `--purity` | `purity` | user-supplied |
| `2.0` | `baseLevel` | diploid linear baseline (log2=0 → linear=2) |

> **Note:** `mafPred` is left as `NA` because CNVkit does not model tumour
> purity internally.  Variants in segments with `mafPred=NA` fall through to
> the CN and log-odds decision rules in `fmiSGZ`.

---

### From GATK4 ModelSegments — `gatk_to_cna_model.py`

**Input:** `*.modelFinal.seg` file from the GATK4 somatic CNV pipeline run
with a Panel of Normals (PON).

**Typical GATK4 CNV pipeline:**

```bash
# 1. Collect read counts
gatk CollectReadCounts -I tumor.bam -L targets.interval_list \
    --interval-merging-rule OVERLAPPING_ONLY -O tumor.counts.hdf5

# 2. Denoise against PON
gatk DenoiseReadCounts -I tumor.counts.hdf5 \
    --count-panel-of-normals pon.hdf5 \
    --standardized-copy-ratios tumor.standardizedCR.tsv \
    --denoised-copy-ratios tumor.denoisedCR.tsv

# 3. Collect allelic counts at common SNP sites
gatk CollectAllelicCounts -I tumor.bam -L snp_sites.interval_list \
    -R hg38.fa -O tumor.allelicCounts.tsv

# 4. Model segments
gatk ModelSegments \
    --denoised-copy-ratios tumor.denoisedCR.tsv \
    --allelic-counts tumor.allelicCounts.tsv \
    --normal-allelic-counts normal.allelicCounts.tsv \
    --output . --output-prefix tumor
# → produces tumor.modelFinal.seg
```

**Convert to CNA model:**

```bash
# Diploid tumour, purity known
python3 gatk_to_cna_model.py tumor.modelFinal.seg sample.cna_calls.txt \
    --purity 0.72 --ploidy 2

# Aneuploid tumour (ploidy estimated from genome-wide log2 median)
python3 gatk_to_cna_model.py tumor.modelFinal.seg sample.cna_calls.txt \
    --purity 0.65 --ploidy 3

# Purity unknown (defaults to 0.5 — approximate)
python3 gatk_to_cna_model.py tumor.modelFinal.seg sample.cna_calls.txt
```

**Arguments:**

| Argument | Default | Description |
|----------|---------|-------------|
| `input_seg` | *(required)* | GATK4 `*.modelFinal.seg` file |
| `output_cna` | *(required)* | Output CNA model file for `fmiSGZ_v2.py` |
| `--purity` | `0.5` | Tumour purity in (0, 1] |
| `--ploidy` | `2.0` | Tumour ploidy (default diploid) |

**Column mapping (GATK4 → CNA model):**

| GATK4 column | CNA model column | Notes |
|--------------|-----------------|-------|
| `CONTIG` | `CHR` | as-is |
| `START` | `segStart` | already 1-based (no conversion needed) |
| `END` | `segEnd` | as-is |
| `LOG2_COPY_RATIO_POSTERIOR_50` | `segLR` | median posterior log2-ratio |
| `MINOR_ALLELE_FRACTION_POSTERIOR_50` | `segMAF` | `NA` when no het SNPs in segment |
| derived | `CN` | `round((2^segLR × baseLevel − 2(1−p)) / p)` |
| derived | `numMAtumorPred` | `round((segMAF × denom − (1−p)) / p)` |
| back-computed | `mafPred` | `(p×M + (1−p)) / (p×CN + 2(1−p))`; `NA` when MAF absent |
| `--purity` | `purity` | user-supplied |
| `p×ploidy + 2(1−p)` | `baseLevel` | accounts for tumour ploidy |

> **Key difference from CNVkit:** GATK4 does not output integer copy numbers
> directly.  `gatk_to_cna_model.py` derives CN from the log2-ratio using the
> tumour purity + ploidy mixture model.  Accurate purity is therefore more
> important for the GATK4 converter than for the CNVkit converter.

**Estimating ploidy from the GATK4 output:**

```python
# Rough ploidy check — run before converting
import statistics, csv

with open("tumor.modelFinal.seg") as f:
    lrs = []
    for line in f:
        if line.startswith('@') or line.startswith('CONTIG'):
            continue
        val = line.split('\t')[5]   # LOG2_COPY_RATIO_POSTERIOR_50
        if val not in ('NaN', 'NA', ''):
            lrs.append(float(val))

median_lr = statistics.median(lrs)
approx_ploidy = 2 ** (median_lr + 1)
print(f"Median LR = {median_lr:.3f}  →  approximate ploidy = {approx_ploidy:.1f}")
# median ≈ 0.0  → ploidy 2 (diploid)
# median ≈ 0.58 → ploidy 3
# median ≈ 1.0  → ploidy 4
```

**Converter comparison:**

| Feature | `cnvkit_to_cna_model.py` | `gatk_to_cna_model.py` |
|---------|--------------------------|------------------------|
| Coordinate system | 0-based BED → **+1** for `segStart` | 1-based → no conversion |
| Integer CN | read from `cn` column | **derived** from log2 + purity |
| Minor allele CN | read from `cn1` column | **derived** from MAF + CN + purity |
| `mafPred` | always `NA` | back-computed when MAF is available |
| `baseLevel` | fixed `2.0` (diploid) | `p × ploidy + 2(1−p)` |
| `--ploidy` flag | not needed | required for aneuploid tumours |
| No-MAF segments | `baf` column empty | `NUM_POINTS_ALLELE_FRACTION == 0` |

---

## basicSGZ\_v2 — Basic Classifier

### Overview

`basicSGZ_v2.py` classifies variants without a copy-number model.  It is
useful when CNA data are unavailable or as a fast first-pass filter.

### Algorithm

```
For each autosomal variant (chrX excluded):

  IF AF > 0.95
    → zygosity = "homozygous"
    → if purity is known AND purity < pthresh: call = "germline"
    → otherwise:                               call = "ambiguous"

  ELSE (AF ≤ 0.95)
    → two-sided binomial test: H₀: AF = 0.5
    → if p < alpha:  call = "somatic"
    → if p ≥ alpha:  call = "germline", zygosity = "het"
```

The binomial test is the key: a variant with AF significantly different from
0.5 is inconsistent with a diploid heterozygous germline state.

### Usage

```
basicSGZ_v2.py [options] aggregated_mutations_file
basicSGZ_v2.py [-h | --help]
basicSGZ_v2.py [--version]
```

### Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `aggregated_mutations_file` | *(required)* | Tab-separated mutations file |
| `-f PURITY_FILE` | `None` | File with percent tumour nuclei (single number) |
| `-a ALPHA` | `0.05` | Significance threshold for binomial test |
| `-p PTHRESH` | `0.9` | Purity threshold for high-AF germline call |
| `-o OUT_PREFIX` | filename stem | Output file prefix |
| `--panel-size MB` | `None` | Panel size in Mb; enables TMB output |
| `--version` | | Print version and exit |

### Examples

**Minimal run — no purity, no TMB:**

```bash
python3 basicSGZ_v2.py sample.mut_aggr.full.txt
# Output: sample.mut_aggr.basic.sgz.txt
```

**With purity file:**

```bash
echo "72" > sample.purity.txt   # 72 % tumour nuclei
python3 basicSGZ_v2.py sample.mut_aggr.full.txt \
    -f sample.purity.txt
```

**With stricter binomial threshold and custom output prefix:**

```bash
python3 basicSGZ_v2.py sample.mut_aggr.full.txt \
    -a 0.01 \
    -o results/sample11
# Output: results/sample11.basic.sgz.txt
```

**Full run — purity + TMB on a 1.1 Mb panel:**

```bash
python3 basicSGZ_v2.py sample.mut_aggr.full.txt \
    -f sample.purity.txt \
    --panel-size 1.1 \
    -o results/sample11
# Outputs:
#   results/sample11.basic.sgz.txt
#   results/sample11.tmb.txt
```

---

## fmiSGZ\_v2 — FMI CNA-Aware Classifier

**Author:** James Sun (`jsun@foundationmedicine.com`)

### Overview

`fmiSGZ_v2.py` implements the Foundation Medicine SGZ method.  By using the
local copy-number (CN) model it computes the exact expected allele frequency
under four germline/somatic hypotheses, produces a log-odds score, and applies
a priority-ordered decision tree.  It also estimates allele burden and its
95 % bootstrap confidence interval.

### Algorithm

#### 1. Segment lookup

Each variant is matched to the overlapping CN segment.  The segment supplies:

| Symbol | Meaning |
|--------|---------|
| `p` | Tumour purity (fraction of tumour cells) |
| `C` | Total integer copy number in tumour cells |
| `M` | Minor-allele copy number in tumour cells |
| `bL` | Diploid baseline log-ratio for the sample |

Sex-chromosome variants (chrX, chrY) are silently excluded.

#### 2. Expected allele frequencies

| Label | Hypothesis | Expected AF |
|-------|-----------|-------------|
| G1 | Germline on the **major** allele | `(p·M + 1·(1−p)) / (p·C + 2·(1−p))` |
| G2 | Germline on the **minor** allele | `(p·(C−M) + 1·(1−p)) / (p·C + 2·(1−p))` |
| S1 | Somatic on the **major** allele  | `p·M / (p·C + 2·(1−p))` |
| S2 | Somatic on the **minor** allele  | `p·(C−M) / (p·C + 2·(1−p))` |

G2 and S2 are only computed when `M ≠ C − M`.

#### 3. Binomial exact tests

```
max_prob_germline = max(P_G1, P_G2)
max_prob_somatic  = max(P_S1, P_S2)
log_odds          = log10(P_germline / P_somatic)
```

Positive log-odds favour germline; negative log-odds favour somatic.
The internal significance threshold is α = 0.01.

#### 4. Decision tree (priority order)

| Condition | Call |
|-----------|------|
| Gene contains `CYP2D6` | `nocall_CYP2D6` |
| Gene contains `HLA` | `nocall_HLA` |
| VAF < 0.05 and 0.2 < purity < 0.9 **and status ∈ {`high`, `moderate`}** | `subclonal somatic` |
| CN or M is NaN | `ambiguous_CNA_model` |
| purity > 0.95 | `nocall_purity>95%` |
| VAF > 0.95 | `germline` |
| MAF model error > 0.06 or LR model error > 0.4 | `ambiguous_CNA_model` |
| P(G) > α **and** P(S) < α | `germline` / `probable germline` |
| P(G) < α **and** P(S) > α | `somatic`  / `probable somatic`  |
| Both > α | `ambiguous_both_G_and_S` |
| Both < α | subclonal rescue → extreme log-odds rescue → `ambiguous_neither_G_nor_S` |
| No matching segment | `nocall_segmentMissing` |

#### 5. Zygosity

| Condition | Zygosity |
|-----------|----------|
| C = 0, M = 0 | `homoDel` |
| C = M = 1 | `homozygous` |
| C = M ≥ 2 | `homozygous` |
| C ≥ 1, M = 0 | `not in tumor` |
| C ≠ M, M ≠ 0 | `het` |
| purity < 0.19 or NaN | `NA` |

#### 6. Allele burden

For definitively called variants, allele burden = `observed_VAF / EAF`.
The 95 % CI is derived from 1 000 binomial resamples of the observed depth.

### Usage

```
fmiSGZ_v2.py [options] aggregated_mutations_file cna_model_file
fmiSGZ_v2.py [-h | --help]
fmiSGZ_v2.py [--version]
```

### Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `aggregated_mutations_file` | *(required)* | Tab-separated mutations file |
| `cna_model_file` | *(required)* | Tab-separated CNA model calls file |
| `-o OUT_PREFIX` | filename stem | Output path prefix |
| `--panel-size MB` | `None` | Panel size in Mb; enables TMB output |

### Examples

**Standard run:**

```bash
python3 fmiSGZ_v2.py \
    sample11.mut_aggr.full.txt \
    sample11.cna_calls.txt \
    -o results/sample11
# Outputs:
#   results/sample11.fmi.sgz.txt        (compact summary)
#   results/sample11.fmi.sgz.full.txt   (full statistics)
```

**With TMB — Foundation Medicine CDx panel (0.8 Mb):**

```bash
python3 fmiSGZ_v2.py \
    sample11.mut_aggr.full.txt \
    sample11.cna_calls.txt \
    --panel-size 0.8 \
    -o results/sample11
# Additional output:
#   results/sample11.tmb.txt
```

**With TMB — whole-exome sequencing (~30 Mb coding):**

```bash
python3 fmiSGZ_v2.py \
    sample11.mut_aggr.full.txt \
    sample11.cna_calls.txt \
    --panel-size 30.0 \
    -o results/sample11_wes
```

**Empty CNA file (graceful fallback):**

```bash
touch empty.cna_calls.txt
python3 fmiSGZ_v2.py sample11.mut_aggr.full.txt empty.cna_calls.txt
# Writes empty output files and exits 0
```

---

## Tumor Mutation Burden (TMB)

TMB is a genomic biomarker widely used to predict response to immune checkpoint
inhibitors.  It is defined as the total number of somatic mutations per
megabase of the analysed genome.

### Formula

```
TMB (mut/Mb) = number of somatic variants / panel size (Mb)
```

### Which variants count as "somatic"?

| Script | Variants counted |
|--------|-----------------|
| `basicSGZ_v2.py` | `germline/somatic == "somatic"` |
| `fmiSGZ_v2.py`  | `call` contains `"somatic"` → includes `"somatic"`, `"probable somatic"`, `"subclonal somatic"` |

### Choosing a panel size

| Assay | Typical panel size |
|-------|--------------------|
| Foundation Medicine CDx | ~0.8 Mb |
| MSK-IMPACT (468 gene) | ~1.2 Mb |
| TruSight Oncology 500 | ~1.9 Mb |
| Whole-exome sequencing | ~30 Mb |
| Whole-genome sequencing | ~2 800 Mb |

> **Note:** Always use the *effective* (callable) panel size reported for your
> assay, not the nominal target size.  The effective size accounts for
> regions with insufficient coverage.

### TMB output file (`*.tmb.txt`)

Tab-separated, one data row per sample:

| Column | Description |
|--------|-------------|
| `sample` | Output prefix used for this run |
| `n_somatic` | Raw count of somatic variants |
| `panel_size_mb` | Panel size provided via `--panel-size` |
| `TMB` | Mutations per megabase (4 decimal places) |

**Example:**

```
sample            n_somatic  panel_size_mb  TMB
results/sample11  14         0.80           17.5000
```

### TMB clinical thresholds

> These are indicative reference values only; clinical cut-offs are
> assay-specific and indication-specific.

| Interpretation | Common threshold |
|----------------|-----------------|
| TMB-High | ≥ 10 mut/Mb |
| TMB-Low  | < 10 mut/Mb |

---

## Input File Formats

Both tools share the same mutations file format.

### Aggregated mutations file (both tools)

Tab-separated with a header row.

| Column | Type | Description |
|--------|------|-------------|
| `mutation` | string | Variant identifier (e.g. `BRAF p.V600E`) |
| `frequency` | float | Observed variant allele frequency (VAF) in [0, 1] |
| `depth` | int | Total sequencing depth at the locus |
| `pos` | string | Genomic position, format `chr7:140453136` |
| `status` | string | *(optional)* Variant-call confidence: `high`, `moderate`, or `low`. Used by `fmiSGZ_v2` to gate the `subclonal somatic` shortcut. Rows without this column are treated as `status = ''` and fall through to the full binomial test. |

**Example:**

```
mutation          frequency  depth  pos               status
BRAF p.V600E      0.42       312    chr7:140453136    high
TP53 p.R175H      0.51       280    chr17:7675088     high
KRAS p.G12D       0.08       450    chr12:25398284    moderate
EGFR p.L858R      0.31       195    chrX:55259515     low
```

> chrX variants are silently excluded by both tools.

### Purity file (`basicSGZ` only, `-f`)

A plain-text file whose first line is the percent tumour nuclei (integer or
float, **no** fraction — enter `72` for 72 %, not `0.72`).

```
72
```

### CNA model file (`fmiSGZ` only)

Tab-separated with a header row.  Produced by `cnvkit_to_cna_model.py` or
`gatk_to_cna_model.py`.

| Column | Type | Description |
|--------|------|-------------|
| `CHR` | string | Chromosome (with or without `chr` prefix) |
| `segStart` | int | Segment start (1-based, inclusive) |
| `segEnd` | int | Segment end (1-based, inclusive) |
| `mafPred` | float\|`NA` | Model-predicted minor-allele frequency |
| `CN` | int\|`NA` | Total copy number in tumour cells |
| `segLR` | float\|`NA` | Observed segment log2-ratio |
| `segMAF` | float\|`NA` | Observed segment minor-allele frequency |
| `numMAtumorPred` | float\|`NA` | Predicted minor-allele copy number |
| `purity` | float | Tumour purity estimate |
| `baseLevel` | float | Diploid baseline log-ratio (linear scale) |

**Example:**

```
CHR  segStart   segEnd     mafPred  CN  segLR     segMAF   numMAtumorPred  purity  baseLevel
1    150501     121571999  NA       2   0.012     0.341    1               0.7200  2.0
1    121572000  125315459  NA       1   -0.821    NA       NA              0.7200  2.0
chr7 55086725   55324313   NA       5   1.142     NA       NA              0.7200  2.0
chr17 7674220   7676592    0.3676   3   0.523     0.281    1               0.7200  2.0
```

> Segments where `mafPred`, `CN`, or `numMAtumorPred` are `NA` will receive
> an `ambiguous_CNA_model` call from `fmiSGZ`.  This occurs for segments with
> no informative heterozygous SNPs (e.g. centromeres, homozygous deletions).

An empty CNA file is handled gracefully — empty output files are written and
the script exits with code 0.

---

## Output File Formats

### basicSGZ: `*.basic.sgz.txt`

Tab-separated, one row per classified autosomal variant.

| Column | Description |
|--------|-------------|
| `mutation` | Variant identifier |
| `pos` | Genomic position |
| `depth` | Sequencing depth |
| `frequency` | Observed VAF (2 d.p.) |
| `germline/somatic` | `somatic` \| `germline` \| `ambiguous` |

### fmiSGZ: `*.fmi.sgz.txt` (compact summary)

| Column | Description |
|--------|-------------|
| `mutation` | Variant identifier |
| `pos` | Genomic position |
| `depth` | Sequencing depth |
| `frequency` | Observed VAF |
| `status` | Variant-call confidence from the input file (`high` / `moderate` / `low` / empty) |
| `C` | Total copy number |
| `germline/somatic` | SGZ call string |
| `zygosity` | Zygosity label |

### fmiSGZ: `*.fmi.sgz.full.txt` (full statistics)

All columns of the compact summary plus:

| Column | Description |
|--------|-------------|
| `afG1` | Expected AF under G1 (germline, major allele) |
| `afS1` | Expected AF under S1 (somatic, major allele) |
| `afG2` | Expected AF under G2 (germline, minor allele) |
| `afS2` | Expected AF under S2 (somatic, minor allele) |
| `p` | Tumour purity |
| `M` | Best-fit allele copy number |
| `logOR_G` | log₁₀(P_germline / P_somatic) |
| `clonality` | Allele burden estimate (`clonal` \| `subclonal` \| `NA`) |
| `clonality_CI_low` | Lower bound of 95 % allele-burden CI |
| `clonality_CI_high` | Upper bound of 95 % allele-burden CI |

### Both tools: `*.tmb.txt`

Written only when `--panel-size` is supplied.

| Column | Description |
|--------|-------------|
| `sample` | Output prefix |
| `n_somatic` | Raw somatic variant count |
| `panel_size_mb` | Panel size in Mb |
| `TMB` | Mutations per megabase |

---

## Call Strings Reference

| Call string | Tool | Meaning |
|-------------|------|---------|
| `germline` | both | High-confidence germline variant |
| `probable germline` | fmiSGZ | Germline favoured but log-odds < 2 |
| `somatic` | both | High-confidence somatic variant |
| `probable somatic` | fmiSGZ | Somatic favoured but log-odds > −2 |
| `subclonal somatic` | fmiSGZ | Low-VAF variant (< 5 %) consistent with a subclonal tumour mutation; only assigned when input `status` is `high` or `moderate` |
| `ambiguous` | basicSGZ | Cannot classify without purity or CNA data |
| `ambiguous_CNA_model` | fmiSGZ | CN model fit too poor to classify |
| `ambiguous_both_G_and_S` | fmiSGZ | Both germline and somatic hypotheses pass α |
| `ambiguous_neither_G_nor_S` | fmiSGZ | Neither hypothesis passes α; log-odds inconclusive |
| `nocall_CYP2D6` | fmiSGZ | Variant in CYP2D6 — excluded (complex locus) |
| `nocall_HLA` | fmiSGZ | Variant in HLA region — excluded (complex locus) |
| `nocall_purity>95%` | fmiSGZ | Sample purity > 95 % — origins indistinguishable |
| `nocall_segmentMissing` | fmiSGZ | No overlapping CN segment found |

---

## Module API

Both scripts can be imported as Python modules.

### basicSGZ\_v2

```python
from basicSGZ_v2 import (
    read_pathology_purity,
    compute_pvalues,
    run_sgz,
    compute_tmb,
)
import pandas as pd

df               = pd.read_csv("sample.mut_aggr.full.txt", sep="\t")
purity           = read_pathology_purity("sample.purity.txt")   # float or "NA"
result_df        = run_sgz(df, purity, alpha=0.05, pthresh=0.9)
tmb              = compute_tmb(result_df, panel_size_mb=1.0)

print(result_df.head())
print(f"TMB = {tmb:.2f} mut/Mb")
```

### fmiSGZ\_v2

```python
from fmiSGZ_v2 import (
    read_cna_model_file,
    read_mut_aggr_full,
    core_SGZ,
    compute_tmb,
)

data_CNA       = read_cna_model_file("sample.cna_calls.txt")
short_variants = read_mut_aggr_full("sample.mut_aggr.full.txt")
results        = core_SGZ(data_CNA, short_variants)
tmb            = compute_tmb(results, panel_size_mb=0.8)

for r in results:
    print(r["mutation"], r["call"], r["zygosity"])

print(f"TMB = {tmb:.2f} mut/Mb")
```

### Key functions

| Function | Module | Description |
|----------|--------|-------------|
| `run_sgz(df, purity, alpha, pthresh)` | basicSGZ | Classify variants; return DataFrame |
| `compute_pvalues(freq, depth)` | basicSGZ | Vectorised binomial p-values |
| `read_pathology_purity(file)` | basicSGZ | Read purity file → float or `"NA"` |
| `core_SGZ(data_CNA, short_variants)` | fmiSGZ | CNA-aware classifier; return list of dicts |
| `read_cna_model_file(fname)` | fmiSGZ | Parse CNA model TSV → list of `Segment` |
| `read_mut_aggr_full(fname)` | fmiSGZ | Parse mutations TSV → list of `SV` |
| `cn2lr_bl(p, bl, cn)` | fmiSGZ | Convert copy number to expected log₂ ratio |
| `compute_tmb(data, panel_size_mb)` | both | Compute TMB (mut/Mb) |

---

## Regression Tests

The `data/` directory contains a self-contained regression suite for validating
both classifiers against known-good outputs.

### Test data layout

```
data/
├── samples/                          # Input files (4 samples)
│   ├── sample1.mut_aggr.full.txt     # Aggregated mutations
│   ├── sample1.cna_calls.txt         # CNA model
│   ├── sample1.pathology_purity.txt  # Purity (may contain "NA")
│   ├── sample2.*  …
│   ├── sample3.*  …
│   └── sample4.*  …
└── expected_samples_outcome/         # Expected outputs (ground truth)
    ├── sample1.basic.sgz.txt
    ├── sample1.fmi.sgz.txt
    ├── sample1.fmi.sgz.full.txt
    ├── sample2.*  …
    ├── sample3.*  …
    └── sample4.*  …
```

### Sample characteristics

| Sample | Variants | Purity | CNA segments | Notable features |
|--------|:--------:|--------|:------------:|------------------|
| sample1 | 12 | NA | 78 | Mixed somatic / germline; purity unavailable |
| sample2 | 13 | 0.529 | — | Half-integer VAF × depth edge case (FBXW7) |
| sample3 | 18 | NA | — | Near-homogeneous somatic tumour; one germline (TP53) |
| sample4 | 57 | 0.500 | — | Large panel; VAF underflow edge case (NUP98 VAF=0.02) |

### Running the tests

Run both scripts for all four samples and diff against the expected output:

```bash
# Create a temporary output directory
mkdir -p /tmp/sgz_test

# basicSGZ — all four samples
for i in 1 2 3 4; do
  python3 basicSGZ_v2.py \
    -f data/samples/sample${i}.pathology_purity.txt \
    -o /tmp/sgz_test/sample${i} \
    data/samples/sample${i}.mut_aggr.full.txt
done

# fmiSGZ — all four samples
for i in 1 2 3 4; do
  python3 fmiSGZ_v2.py \
    -o /tmp/sgz_test/sample${i} \
    data/samples/sample${i}.mut_aggr.full.txt \
    data/samples/sample${i}.cna_calls.txt
done

# Compare every output against expected
all_pass=true
for i in 1 2 3 4; do
  for ext in basic.sgz.txt fmi.sgz.txt fmi.sgz.full.txt; do
    got=("/tmp/sgz_test/sample${i}.${ext}")
    exp="data/expected_samples_outcome/sample${i}.${ext}"
    if diff -q "$got" "$exp" > /dev/null; then
      echo "PASS  sample${i}.${ext}"
    else
      echo "FAIL  sample${i}.${ext}"
      diff "$got" "$exp"
      all_pass=false
    fi
  done
done
$all_pass && echo "All tests passed."
```

**Expected output (all 12 files):**

```
PASS  sample1.basic.sgz.txt
PASS  sample1.fmi.sgz.txt
PASS  sample1.fmi.sgz.full.txt
PASS  sample2.basic.sgz.txt
PASS  sample2.fmi.sgz.txt
PASS  sample2.fmi.sgz.full.txt
PASS  sample3.basic.sgz.txt
PASS  sample3.fmi.sgz.txt
PASS  sample3.fmi.sgz.full.txt
PASS  sample4.basic.sgz.txt
PASS  sample4.fmi.sgz.txt
PASS  sample4.fmi.sgz.full.txt
All tests passed.
```

### Bugs identified and fixed during test development

Running the regression suite against v2.0.0 revealed three bugs, all fixed in
v2.0.1 (see [Changelog](#changelog)):

| # | Script | Bug | Symptom | Fix |
|---|--------|-----|---------|-----|
| 1 | `basicSGZ_v2.py` | `zygosity` column included in output | Output had 6 columns; expected 5 (`mutation`, `pos`, `depth`, `frequency`, `germline/somatic`) | Removed `zygosity` from the `to_csv` column list |
| 2 | `basicSGZ_v2.py` | Float formatting of `frequency` | Values such as `0.6`, `0.1` were not zero-padded to two decimal places (`0.60`, `0.10`) | Added `float_format='%.2f'` to `to_csv` |
| 3 | `fmiSGZ_v2.py` | Python 3 banker's rounding for `obs_count` | `round(0.55 × 590) = 324` (Python 3) vs `325` (round-half-up), producing `logOR_G = −29.1` instead of `−28.7` for FBXW7 in sample2 | Changed to `int(depth × vaf + 0.5)` (round-half-up) |
| 4 | `fmiSGZ_v2.py` | `logOR_G = −inf` when germline p-value underflows | For very low VAF variants (e.g. NUP98, VAF = 0.02, depth = 1199), `binomtest` returns `p = 0.0` under the germline hypothesis, making `log10(0) = −∞` | Replaced the `−inf` special case with a `logsumexp` computation over individual `binom.logpmf` values in log-space |
| 5 | `fmiSGZ_v2.py` | `ValueError: p (nan) must be in range [0,1]` | When a CNA segment has `CN`, `numMAtumorPred`, or `purity` = `NA`, the computed expected AF is NaN and `binomtest` raises a `ValueError` | Added a NaN guard before the binomial tests; `purity = 'NA'` now stored as `np.nan` like other nullable segment fields |
| 6 | `fmiSGZ_v2.py` | Over-counting of somatic variants due to low-quality calls | With purity ~0.70, every variant with VAF < 5 % was auto-classified `subclonal somatic` regardless of variant-call quality, inflating TMB by including noise calls | The `subclonal somatic` shortcut now requires input `status ∈ {high, moderate}`; `low`-quality variants fall through to full binomial testing |

> **Maintenance note (v2.0.2):** the `status` column was added to
> `fmiSGZ_v2.py`'s output files (see [Changelog](#changelog)), but the
> checked-in `data/expected_samples_outcome/*.fmi.sgz*.txt` fixtures were not
> regenerated at the same time, so the regression suite briefly reported
> false `FAIL`s on an otherwise-correct classifier — the call strings,
> zygosity, log-odds, and confidence intervals were unchanged; only the
> fixtures were stale. **Whenever a script's output schema changes,
> regenerate the fixtures in the same commit** (diff the new output against
> the old fixture field-by-field, not just byte-for-byte, to confirm the
> classification logic itself didn't move).

---

## Changelog

See [CHANGELOG.md](CHANGELOG.md) for the full version history.
