# Revised WMHS paper

**Auditing a single-cell reporting rule: coverage, false positives and target choice**

This independent revision starts from repository commit `b503971` and addresses
the review of `WHMS_BIO_PAPER.pdf`. After pulling through `a336594`, it adopts
the verified synthetic external-control example from the incoming paper while
preserving the later fixed-pool and false-positive analyses. This integration
does not modify the other manuscript, experiment implementations or saved results.

- `main.pdf`: revised anonymous manuscript.
- `main.tex`, `sections/`, `refs.bib`, `figures/`: editable LaTeX sources and figures.
- `source.zip`: portable LaTeX/Overleaf source bundle, with `main.tex` at its root.
- `THIRD_CRITIQUE_RESPONSE.md`: current independent assessment and experimental follow-up.
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
- `results_manifest.json`: 38 original pinned input tables plus four external-control tables.

The condensed build has **8 pages of main text**, **2 pages of references**
(starting on page 9), and **3 pages of appendix**, for **13 pages overall**.
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
PYTHONPATH=. python -m pytest papers/whms_bode/test_results.py papers/whms_bode/test_precision_followup.py tests/test_external_control_demo.py -q
python papers/whms_bode/make_tables.py --check
python papers/whms_bode/audit_results.py
python papers/whms_bode/make_tables.py
python papers/whms_bode/make_diagnostic_table.py
python papers/whms_bode/make_fixed_pool_tables.py
python papers/whms_bode/make_external_control_table.py
python papers/whms_bode/make_figures.py
python papers/whms_bode/make_precision_results.py
papers/whms_bode/build.sh
python papers/whms_bode/package_source.py
```

`make_figures.py` reads version-pinned tables and reproduces the retained controlled-calibration figure.
The controlled calibration figure deduplicates grid views only after asserting
that rates at shared counts are identical. The learned-estimator table and extended estimator checks remain analytical
artifacts, included as noncompiled source files in the bundle. The manuscript
retains only the reference and interval-control examples needed for its argument.

Verification completed for this revision:

- 29 paper/precision tests and 35 external-control harness tests passed,
  including target-dependent rankings and the corrected comparison denominators.
- All 30 numerical cells in the retained learned-estimator analysis table match
  the saved results.
- In the earlier revision, 55 existing paper-number/submission tests passed against the original
  manuscript. These corroborate retained results but do not validate the new
  manuscript by themselves; they did not catch several corrected issues.
- Build and anonymity checks passed with no errors, undefined references or
  overfull boxes.
- All 13 rendered pages were visually inspected, including figures and tables.

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
anonymous ZIP, alongside a second scientific figure showing null and effect
curves together. The main text now contains three tables and two figures.
