# Revised WMHS paper

**When the reference hides uncertainty: a data-based simulation audit**

This independent revision audits the recorded justification for a 50-cell
reporting rule. The interval-comparison revision starts from `bf6b825`, independently
assesses the new critique, and adds a fixed-design comparison of five intervals
at 30 and 50 cells: 32,000 simulated studies and 160,000 interval evaluations.
It leads with the measured reference-versus-parameter coverage discrepancy,
compresses the synthetic external-control illustration, and distinguishes
nominal calibration from the weaker historical acceptance criterion.
The original experiment implementations and saved precision results are unchanged.

- `main.pdf`: revised anonymous manuscript.
- `main.tex`, `sections/`, `refs.bib`, `figures/`: editable LaTeX sources and figures.
- `source.zip`: portable LaTeX/Overleaf source bundle, with `main.tex` at its root.
- `INCOMING_BRANCH_REVIEW.md`: current independent assessment of Diego's repaired-interval results and Claude's alternate manuscript.
- `REPAIR_CRITIQUE_REVIEW.md`: preceding independent assessment, experimental findings and remaining limitations.
- `REPAIR_FOLLOWUP_DESIGN.md`, `repair_followup.py`: fixed design and implementation of the new interval comparison.
- `PRECISION_CRITIQUE_REVIEW.md`: preceding assessment of Monte Carlo reversals, interval scope and supporting evidence.
- `ANALYTICAL_ARTIFACTS.md`: inventory of local diagnostic sources excluded from the submission.
- `STRUCTURE_REVIEW.md`: preceding assessment of venue fit, structure, novelty and bootstrap claims.
- `THIRD_CRITIQUE_RESPONSE.md`: preceding independent assessment and experimental follow-up.
- `FOLLOWUP_DESIGN.md`: analysis settings fixed before the larger run.
- `precision_followup.py`: 160,000-study precision and paired interval comparison.
- `INCOMING_PAPER_REVIEW.md`: assessment of the pulled paper, adopted findings and corrected overclaims.
- `SECOND_CRITIQUE_RESPONSE.md`: assessment of the preceding critique and fixed-pool revision.
- `CRITIQUE_RESPONSE.md`: historical assessment of the first critique.
- `diagnostics/`, `diagnose_low_counts.py`: exploratory null/sparsity diagnostic and provenance.
- `evidence/review_record.json`: review-safe documentary statement and numerical context.
- `evidence/case_provenance.json`: local exact-source mapping; excluded from the review bundle.
- `REVISION_NOTES.md`: substantive corrections, supporting evidence and remaining limits.
- `build_verification.json`: exact PDF checksum, pagination and build checks.
- `results_audit.json`: independently recomputed counts and result-table checksums.
- `results_manifest.json`: 38 original pinned input tables, four external-control tables and five
  control-panel audit inputs.

The current build has **9 pages of main text**, **2 pages of references**
(starting on page 10), and **3 pages of appendix**, for **14 pages overall**.
Earlier versions are preserved in `archive/expanded_18page_source.zip` and
`archive/condensed_14page_source.zip`; the first critique revision is preserved
in `archive/first_critique_source.zip`. The pre-precision revision is in
`archive/pre_precision_source.zip`. These archives are historical local
artifacts and are excluded from the current anonymous source bundle. The required responsible-use statement is wholly
within the main text. The official, unmodified NeurIPS 2026 style replaces the
2024 compatibility shim. Its conference download URL and checksum are recorded
in `style_provenance.json`. The current WMHS call allows nine main-text pages:
https://wmhs-neurips.github.io/WMHS/

## Build

From the repository root, with TeX Live (`pdflatex`, `bibtex`) and Python with
`pypdf` installed:

```sh
papers/whms_bode/build.sh
```

An alternate Python interpreter can be supplied explicitly:

```sh
PYTHON_BIN=/path/to/python papers/whms_bode/build.sh
```

The build fails for LaTeX errors, unresolved citations/references, overfull
boxes, more than nine main-text pages, a stale PDF, missing source/figure
files, a changed official style, or detected author-identifying strings.
The bibliography deliberately starts on a fresh page, making the page-limit
check unambiguous. The check does not replace human anonymity review.

## Verify and regenerate from saved results

Install `requirements-review.txt` into a separate Python environment. This is
an artifact-reading and rendering environment, **not** a claim to reproduce
the original simulation environments; those are recorded in result sidecars.
The first editorial revision used only saved results. The follow-up revision
adds an explicitly exploratory 6,400-draw diagnostic of the original interval
under the null and alternative. The latest revision replays the same draws,
adds coverage against the fixed eligible-pool mean, and computes source-pool
variance summaries. It does not select or validate a replacement
threshold. The latest follow-up adds 5,000 fresh studies per setting and a
paired Welch interval comparator, with settings fixed in `FOLLOWUP_DESIGN.md`.
The incoming-paper integration adds a synthetic external-control
comparison: all 240 primary and nine balance summary rows were independently
reproduced. Original experiment code and saved tables remain unchanged.

```sh
PYTHONPATH=. python -m pytest papers/whms_bode/test_results.py papers/whms_bode/test_precision_followup.py papers/whms_bode/test_repair_followup.py papers/whms_bode/test_incoming_interval.py tests/test_external_control_demo.py tests/test_retained_control_audit.py -q
python papers/whms_bode/make_tables.py --check
python papers/whms_bode/audit_results.py
python papers/whms_bode/make_tables.py
python papers/whms_bode/make_diagnostic_table.py
python papers/whms_bode/make_fixed_pool_tables.py
python papers/whms_bode/make_external_control_table.py
python papers/whms_bode/make_figures.py
python papers/whms_bode/make_precision_results.py
python papers/whms_bode/make_repair_results.py
papers/whms_bode/build.sh
python papers/whms_bode/package_source.py
```

`make_figures.py` reads version-pinned tables and reproduces the historical controlled-calibration figure, now retained only as a local analytical artifact.
The controlled calibration figure deduplicates grid views only after asserting
that rates at shared counts are identical. The learned-estimator table and extended estimator checks remain local
analytical artifacts, excluded from the submission bundle; their purposes are
listed in `ANALYTICAL_ARTIFACTS.md`. The manuscript includes a compressed
learned-generator check, with actual denominators and synthetic training provenance.
The redundant normal-mean interval-width demonstration is no longer compiled.

Verification completed for this revision:

- 45 paper/precision/repair/incoming-result tests, 35 external-control harness tests and 21
  retained-control audit tests passed,
  including target-dependent rankings and the corrected comparison denominators.
- All 30 numerical cells in the retained learned-estimator analysis table match
  the saved results.
- In the earlier revision, 55 existing paper-number/submission tests passed against the original
  manuscript. These corroborate retained results but do not validate the new
  manuscript by themselves; they did not catch several corrected issues.
- Build and anonymity checks passed with no errors, undefined references or
  overfull boxes.
- All 14 rendered pages were visually inspected, including figures and tables.

The temporary review environment used while editing was
`/tmp/whms_bode_env/bin/python`; it need not exist on another machine.
The source ZIP can be compiled without the repository's data tables because
all generated LaTeX and figure assets are included. Regenerating figures or
running the claim tests requires this repository and its pinned results.

The new diagnostic can be rerun with local Lee raw data available:

```sh
PYTHONPATH=. python papers/whms_bode/diagnose_low_counts.py
python papers/whms_bode/make_diagnostic_table.py
python papers/whms_bode/make_fixed_pool_tables.py
```

It exports only aggregate condition summaries, not cell-level data. Four counts
and one controlled-run seed stream were chosen to investigate the reported
anomaly; the run is exploratory, conditional on the existing cohorts. All 16
alternative settings match the saved first-seed summaries. The null and
first-crossing checks are included in the dedicated tests. The source ZIP
contains paired diagnostic rates, Monte Carlo precision, source-pool summaries,
candidate reversals, synthetic external-control aggregates and the review-safe
documentary record in addition to the
standalone LaTeX inputs. The exact internal filenames and source mapping are
excluded; the build and packaging checks scan for these identifiers.

The fixed target is the mean contrast in each eligible two-patient empirical
pool, computed from all source cells before either synthetic arm is sampled.
It is fixed conditional on that pool. Pair selection changes the pool; this is
not patient-population inference or coverage against an all-patient pooled mean.
`make_fixed_pool_tables.py` regenerates the paired table, Wilson intervals,
reversal table and descriptive source-pool signal-to-noise calculations.

## Precision follow-up

With local Lee data available, run:

```sh
PYTHONPATH=. python papers/whms_bode/precision_followup.py
python papers/whms_bode/make_precision_results.py
```

The follow-up samples the exact marginal mature-count law of the original
generator using five fresh seed streams; it does not reproduce the older
seed-by-seed outputs. The 32 settings give 160,000 simulated studies and
320,000 interval evaluations. All are finite. Exports contain aggregate
coverage, null/effect rejection, seed rates and paired differences only.
The maximum binomial Monte Carlo SE is 0.71 percentage points. The
`precision_provenance.json` sidecar hashes the fixed design, implementation,
source-pool summaries and result files. The three precision CSVs ship in the
anonymous ZIP, alongside the scientific figure showing null and effect
curves together. The main text now contains four tables and one figure.

## Reporting-band and interval-repair follow-up

With the same local source data available:

```sh
PYTHONPATH=. python papers/whms_bode/repair_followup.py
python papers/whms_bode/make_repair_results.py
```

The design compares percentile B=200, percentile B=2,000, BCa B=2,000,
bootstrap-t B=2,000 and Welch on identical study samples. The bootstrap methods
share resamples. Coverage is scored against the fixed eligible-pool effect;
no drawn-reference benchmark is substituted for this target. Five fresh seed
streams produce 2,000 studies per cohort/pool/count/shift setting. All 160,000
interval evaluations return finite intervals. The maximum binomial Monte Carlo
SE is 1.12 percentage points. Four aggregate CSVs include Wilson intervals,
paired differences, per-seed rates and numerical diagnostics. Their provenance
sidecar hashes the design, implementation and results. No cell-level data are
exported.

At 30 cells, the original interval has 78.10–89.05% alternative coverage and
14.20–40.45% null rejection. Increasing the bootstrap budget at 50 cells reduces
null rejection by only 0.5–1.0 percentage points. Bootstrap-t improves the
SMC/reference null rejection from 8.45% to 6.05% [5.09, 7.18], but alternative
detection decreases from 76.3% to 65.9%. No method meets the joint criteria
at either tested count; this study does not select a replacement interval or
validate the entire reporting range. BCa adjusted endpoints also reach poorly
resolved bootstrap tails in a substantial fraction of studies; these finite-budget
results are not a general impossibility claim about skew-aware methods.

BCa endpoints are checked against SciPy's independent multisample jackknife
implementation on continuous, discrete, skewed and sparse samples. Further
checks cover the original percentile/Welch routines, studentized endpoint
construction, affine transformations, arm exchange, degeneracy, exact null
coverage/rejection accounting, numerical exports and provenance hashes.

## Selective integration of incoming branches

After `9c6afbc`, the incoming diagnosis at `ae8c109` and alternate manuscript
at `0176c335` were independently reviewed. The paper adds source-shape/arm-size
sensitivity, a qualified patient-t detection/availability comparison, and the
exact-assignment-balance implication of the existing estimator identity. It keeps
the fixed-pool result in the main text and does not import the incoming repair
pass labels, causal skew shares, or blanket patient-method abstention claims.

All 64 diagnosis rates were reproduced from local source data (128,000 interval
evaluations). All 1,440 repair seed rows were reaggregated and checked against
the saved summaries. The source bundle includes the two aggregate tables used
in the paper and an interpretation README. The full local exports and pinned
code/result hashes are in `diagnostics/incoming_*`; no raw cells are exported.
`audit_incoming_interval.py` regenerates and verifies these exports from the
pinned Git objects. `reproduce_incoming_diagnosis.py` reruns the diagnosis using
an extracted incoming source tree and local Lee data; its command-line arguments
are documented in the script and `INCOMING_BRANCH_REVIEW.md`.

Validation now includes 101 relevant checks. The nine-page main text and
three-page appendix are retained. `INCOMING_BRANCH_REVIEW.md` explains each
accepted and rejected change, including differences between interval targets,
implementations, Monte Carlo criteria and resampling units.
