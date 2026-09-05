# A1 — the reference-recipe recovery experiment, and what it does not license

**Tables:** `results/2026-09-05_736e0df/pseudobulk_recovery{,_gap}.parquet`
**Code:** `src/harness/pseudobulk_recovery.py`, `src/harness/jobs/run_pseudobulk_recovery.py`

> **Read this before quoting the tables.** Their sidecars were written before the
> correction below and frame the result as justifying a W1 rebuild. **It does
> not.** The numbers are unaffected; the interpretation in
> `what_this_answers` is wrong. Any re-run of the job emits the corrected text.

## What was measured

Two references built from the same Lee/SMC cells and the same 770 genes, by
**averaging over cells**:

| | recipe |
|---|---|
| `linear` | `mean(CP10K)` |
| `log1p`  | `expm1(mean(log1p(CP10K)))` |

Then pseudobulk at fractions we set, deconvolved with both, scored against the
realised draw. Two legs: `cross` (SMC reference → KUL3 pseudobulk) and `within`
(SMC reference → held-out SMC pseudobulk, the batch confound check).

| leg | method | linear | log1p | gap |
|---|---|---|---|---|
| cross | nnls | 0.945 | 0.779 | **+0.166** |
| cross | nusvr | 0.983 | 0.857 | **+0.126** |
| within | nnls | 0.983 | 0.930 | +0.053 |
| within | nusvr | 0.986 | 0.959 | +0.027 |

Batch costs 0.003 (`within`-linear 0.986 vs `cross`-linear 0.983). The recipe
costs 0.126–0.166, and costs *more* cross-cohort — 0.027 within against 0.126
across for nu-SVR.

The log recipe also returns 24–86 exactly-zero epithelial fractions per 200
samples where linear returns 0–9.

## The correction, 2026-09-05

**This does not justify rebuilding the committed S matrices**, and the tables'
sidecars say it does.

`run_full_reference` does no cell-averaging. It accumulates **one pseudo-cell
per cell type carrying that type's summed counts**
(`run_full_reference.py:314`), then CP10K-normalises and `log1p`s that single
row. With one row per type there is no within-type averaging, hence no Jensen
gap, and `expm1` inverts the profile **exactly** — measured at 4.7e-06, which is
float32 noise, and pinned by
`tests/test_bulk_deconvolution.py::test_expm1_is_an_exact_inverse_for_the_committed_construction`.

So the Stage 4 run, which passed `--linearise-reference`, was already on the
correct linear scale. Its gate failure at 0.462 and 0.479 is **not** a scale
artifact, and a rebuild would emit `expm1` of what already exists.

## What the result does say

A reference built by averaging over individual cells must be built on the linear
scale, and the penalty compounds under cross-cohort transfer. That is a real
finding about how to build a reference. It is not a finding about this one.

`build_signature_sparse(profile_scale="linear")` exists for exactly that case.
It is redundant for the current pseudo-cell construction, and six tests pin that
marker selection does not move with it.

## Standing caveats

- Instrument-level throughout. Nothing here is a result about colorectal cancer.
- Pseudobulk has no ambient contamination and no library preparation, so the
  absolute correlations (0.94–0.99) do not transfer to real bulk. **The gap
  transfers; the absolutes do not.**
- `r_non_epithelial` and `r_epithelial` are identical in every row. With six
  compartments summing to 1 that is an algebraic identity, not two independent
  confirmations.
