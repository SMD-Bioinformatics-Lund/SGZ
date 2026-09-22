# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [2.0.2] - 2026-09-22

### Added

#### `fmiSGZ_v2.py`

- **`status` column in output** — `*.fmi.sgz.txt` (compact summary) and
  `*.fmi.sgz.full.txt` (full statistics) now include the input file's
  `status` column (variant-call confidence: `high` / `moderate` / `low` /
  empty), positioned after `frequency`. Previously this field was read and
  used internally to gate the `subclonal somatic` shortcut (see v2.0.1 bug
  #6) but was not written back out, making it impossible to audit which
  quality tier drove a given call downstream.

### Fixed

- **Stale regression fixtures** — `data/expected_samples_outcome/*.fmi.sgz.txt`
  and `*.fmi.sgz.full.txt` were regenerated to include the new `status`
  column. The fixtures had not been updated when the column was added, so
  the regression suite (`README.md § Regression Tests`) was reporting
  spurious `FAIL` results even though the underlying classification logic
  (call strings, zygosity, log-odds, allele burden, confidence intervals)
  was unchanged and correct in all 4 samples. Verified via full field-level
  diff, not just byte-diff, before regenerating.

---

## [2.0.1] - 2026-03-10

### Fixed

#### `basicSGZ_v2.py`

- **Output columns** — removed `zygosity` from `*.basic.sgz.txt`; the column
  was computed internally but must not appear in the final TSV (output now has
  five columns: `mutation`, `pos`, `depth`, `frequency`, `germline/somatic`).
- **Frequency formatting** — `frequency` values are now always written with
  exactly two decimal places (e.g. `0.60`, `0.10`) by passing
  `float_format='%.2f'` to `pandas.DataFrame.to_csv`; previously trailing zeros
  were dropped by pandas (e.g. `0.6`).

#### `fmiSGZ_v2.py`

- **Observed alt-count rounding** — changed `round(depth × vaf)` (Python 3
  banker's rounding, rounds ties to even) to `int(depth × vaf + 0.5)`
  (round-half-up) to match the intended arithmetic.  Only affects variants
  where `depth × vaf` is exactly a half-integer (e.g. FBXW7 in sample2:
  `0.55 × 590 = 324.5` → `324` before, `325` after), correcting `logOR_G`
  from `−29.1` to `−28.7`.
- **`logOR_G = −inf` for extreme low-VAF variants** — when `binomtest` returns
  `p = 0.0` for the germline hypothesis due to floating-point underflow (e.g.
  NUP98, VAF = 0.02, depth = 1199), the log-odds is now computed in log-space
  using `scipy.special.logsumexp` over individual `binom.logpmf` tail values
  instead of returning `−∞`.

### Added

- `scipy.special.logsumexp` import in `fmiSGZ_v2.py` (required by the
  log-space `logOR_G` fix above).
- Regression test suite under `data/samples/` and
  `data/expected_samples_outcome/` covering 4 samples × 3 output files = 12
  reference comparisons; documented in `README.md § Regression Tests`.

---

## [2.0.0] - 2025-01-01

### Added

#### Both scripts

- `compute_tmb()` function — counts somatic calls and divides by effective
  panel size in Mb.
- `--panel-size MB` CLI argument — TMB output (`*.tmb.txt`) is produced only
  when this flag is provided.
- `*.tmb.txt` output file with columns: `sample`, `n_somatic`,
  `panel_size_mb`, `TMB`.

#### `basicSGZ_v2.py` (previously undocumented — new public interface)

- Comprehensive module-level docstring covering algorithm, output format,
  usage, and dependencies.
- NumPy-style docstrings for all public functions: `read_pathology_purity`,
  `compute_pvalues`, `run_sgz`, `compute_tmb`, `parse_args`, `main`.
- Type annotations on all function signatures (`typing.Optional`,
  `typing.Union`) for Python 3.8 compatibility.
- Structured logging via `logging.getLogger(__name__)` at module level;
  info-level messages for purity value, variant count, classification summary,
  and TMB.
- `zygosity` column in `*.basic.sgz.txt` output — previously computed
  internally but silently dropped.
- `--version` CLI argument.
- Warning log message when the purity value in `-f` cannot be parsed.

#### `fmiSGZ_v2.py`

- `compute_tmb()` function; counts all calls whose string contains `"somatic"`
  (`"somatic"`, `"probable somatic"`, `"subclonal somatic"`).
- `--panel-size` argument in `_arg_parser()`.
- TMB computation block in `main()` with `*.tmb.txt` output.
- `panel_size_mb` parameter documented in `main()` docstring.

#### CNA converters (new scripts)

- `cnvkit_to_cna_model.py` — converts CNVkit `.cns` called-segment files to
  the fmiSGZ CNA model format.  Handles both tab-separated and space-padded
  CNVkit output.  Converts 0-based BED coordinates to 1-based inclusive
  coordinates.  Marks missing `baf` and `cn1` fields as `NA`.
- `gatk_to_cna_model.py` — converts GATK4 `ModelSegments` `*.modelFinal.seg`
  files to the fmiSGZ CNA model format.  Skips SAM-style `@` header lines.
  Derives integer total copy number and minor-allele copy number from the
  posterior log2-ratio and minor-allele fraction using the tumour purity +
  ploidy mixture model.  Back-computes `mafPred` from called integer copy
  numbers.  Supports non-diploid tumours via `--ploidy`.

### Fixed

#### Both scripts

- Output file stem derivation: replaced fragile `fname[:-4]` string slice with
  `os.path.splitext()` so any file extension is handled correctly.

#### `basicSGZ_v2.py`

- `compute_pvalues`: `depth` was passed as `float` to
  `scipy.stats.binomtest(n=…)`, which requires an integer; values are now cast
  to `int` before calling `binomtest`.
- Pandas `SettingWithCopyWarning`: added `.copy()` after the chrX filter to
  ensure the working DataFrame is not a chained slice.

#### `fmiSGZ_v2.py`

- `logger` scope bug: `logger` was only defined inside the `__main__` guard,
  causing `NameError` when `main()` was called programmatically; moved to
  module level as `logger = logging.getLogger(__name__)`.
- `__main__` block was creating a root logger via `logging.getLogger()`
  (no argument); corrected to use the module-level `logger`.

### Changed

- Replaced deprecated `scipy.stats.binom_test` with `scipy.stats.binomtest`
  (requires scipy ≥ 1.7) in both classifiers.

---

## [1.0.0] - 2020-01-01

### Added

- Initial implementation of `fmiSGZ` (Foundation Medicine Somatic/Germline
  Zygosity classifier).
- Copy-number-aware four-hypothesis binomial exact test framework (G1, G2, S1,
  S2).
- Priority-ordered decision tree producing call strings: `germline`,
  `probable germline`, `somatic`, `probable somatic`, `subclonal somatic`,
  `ambiguous_*`, `nocall_*`.
- Zygosity assignment (`het`, `homozygous`, `homoDel`, `not in tumor`).
- Bootstrap allele-burden estimation (1 000 resamples, 95 % CI).
- Python 2 implementation using deprecated `scipy.stats.binom_test`.

---

[2.0.2]: https://github.com/your-org/SGZ_SP_TMB/compare/v2.0.1...v2.0.2
[2.0.1]: https://github.com/your-org/SGZ_SP_TMB/compare/v2.0.0...v2.0.1
[2.0.0]: https://github.com/your-org/SGZ_SP_TMB/compare/v1.0.0...v2.0.0
[1.0.0]: https://github.com/your-org/SGZ_SP_TMB/releases/tag/v1.0.0
