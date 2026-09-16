# Revised WMHS paper

**When a recovery curve fails to validate a decision rule: a single-cell calibration case study**

This independent revision starts from repository commit `b503971` and addresses
the review of `WHMS_BIO_PAPER.pdf`. The original `paper/wmhs` manuscript,
experiment implementations, configurations and saved results are unchanged.

- `main.pdf`: revised anonymous manuscript.
- `main.tex`, `sections/`, `refs.bib`, `figures/`: editable LaTeX sources and figures.
- `source.zip`: portable LaTeX/Overleaf source bundle, with `main.tex` at its root.
- `CRITIQUE_RESPONSE.md`: point-by-point assessment of the follow-up critique.
- `diagnostics/`, `diagnose_low_counts.py`: exploratory null/sparsity diagnostic and provenance.
- `evidence/case_provenance.json`: exact motivating-document excerpts and checksums.
- `REVISION_NOTES.md`: substantive corrections, supporting evidence and remaining limits.
- `build_verification.json`: exact PDF checksum, pagination and build checks.
- `results_audit.json`: independently recomputed counts and result-table checksums.
- `results_manifest.json`: the same pinned input tables as the original submission.

The condensed build has **9 pages of main text**, **2 pages of references**
(starting on page 10), and **3 pages of appendix**, for **14 pages overall**.
Earlier versions are preserved in `archive/expanded_18page_source.zip` and
`archive/condensed_14page_source.zip`. The required responsible-use statement is wholly
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
under the null and alternative. It does not select or validate a replacement
threshold. Original experiment code and saved tables remain unchanged.

```sh
python -m pytest papers/whms_bode/test_results.py -q
python papers/whms_bode/make_tables.py --check
python papers/whms_bode/audit_results.py
python papers/whms_bode/make_tables.py
python papers/whms_bode/make_diagnostic_table.py
python papers/whms_bode/make_figures.py
papers/whms_bode/build.sh
python papers/whms_bode/package_source.py
```

`make_figures.py` reads version-pinned tables and reproduces the retained controlled-calibration figure.
The controlled calibration figure deduplicates grid views only after asserting
that rates at shared counts are identical. All ten learned-estimator rows are
regenerated from the primary result table by `make_tables.py`; no hand-entered
generator-family maxima are used there.

Verification completed for this revision:

- 17 dedicated claim tests passed, including the actual tolerance-switch
  count, overlapping exclusion categories and direct performance at 50 cells.
- All 30 numerical cells in the regenerated learned-estimator table match
  the saved results.
- 55 existing paper-number/submission tests passed against the original
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
```

It exports only aggregate condition summaries, not cell-level data. Four counts
and one controlled-run seed stream were chosen to investigate the reported
anomaly; the run is exploratory, conditional on the existing cohorts. All 16
alternative settings match the saved first-seed summaries. The null and
first-crossing checks are included in the dedicated tests. The source ZIP
contains the diagnostic CSV, candidate-reversal summary and documentary evidence
in addition to the standalone LaTeX inputs.
