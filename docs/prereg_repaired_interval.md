# Pre-registration — diagnosing and attempting to repair the within-patient interval

**Status:** registered before any result table for this question was produced.
**Branch:** `wmhs/repaired-interval` (from `wmhs/build-fixes`).
**Owner:** W2. **Date registered:** 2026-09-15.

Written to answer one reviewer criticism of the NeurIPS workshop submission:

> "The paper diagnoses without prescribing. It declines to propose or test a
> repaired interval. Candidates would include a studentised or bias-corrected
> bootstrap, a patient-level or hierarchical bootstrap, or a pseudobulk-based
> interval. The inflated null rejection is striking, but the root cause is only
> partly explained: the all-zero tumour samples explain the 5-cell case, not the
> 12-30% rate at 50 cells or the 9.5% rate at 800 cells in pooled SMC."

A fix that fails is a result. A fix that was never tried is not. The candidate
list, the outcome measures and the falsifiers below are fixed **now**, so that
whatever comes back is reportable rather than selectable.

---

## 0. Which interval this is about, and which it is not

Two different bootstraps live in this repository. They are not the same object
and conflating them would make this exercise worthless, so the distinction is
stated once, here, and every table produced under this registration carries an
`interval_object` column naming which one it is.

| | **THIS registration** | **Not this** |
|---|---|---|
| Code | `src/harness/interval.py::within_patient_intrinsic_ci` | `src/reference/interval_calibration.py` |
| Resampling unit | **cells**, within one synthetic patient | **patients**, within a cohort |
| Estimand | one patient's intrinsic term | the cohort's mean effect |
| Typical n | 5 – 800 mature cells | 10 – 44 patients |
| Known result | the inflated null rejection the reviewer asks about | the percentile interval is `z·sqrt((n−1)/n)/t(n−1)` times the width it claims; Student-t is calibrated at every n measured |

`interval_calibration.py`'s closed form (`width_ratio`,
`expected_false_positive_rate`) is an **arithmetic** statement about a
percentile bootstrap of a mean of `n` exchangeable values. It therefore also
applies, as a floor, to the cell-level bootstrap at cell count `n`. But that
floor is **0.8 percentage points at n = 50** (`width_ratio(50) = 0.966`,
`expected_false_positive_rate(50) = 5.8%`) and **essentially nothing at
n = 800** (`width_ratio(800) = 0.9993`). It cannot be the explanation for
12–30% at 50 cells or 9.5% at 800, and this registration does not offer it as
one. It is entered as the measured floor `P3` below, precisely so that the gap
between it and the observed rate is visible as a number.

## 1. The object, stated as algebra

With `weighting="normal"` — what the sweep and the committed diagnostic run —
`src/estimator/kitagawa.py::decompose` gives

```
intrinsic_hat  =  f_n · (mean_t − mean_n)
```

and `within_patient_intrinsic_ci` holds `f_n` and `f_t` **fixed** across
resamples. `f_t` does not enter the normal-weighted intrinsic term at all. The
interval is therefore an affine image of a two-sample percentile bootstrap of a
difference of means, with a scale factor of `f_n = 0.40` that is constant across
the whole sweep.

Under the null `shift = 1.0`, `pseudobulk._apply_shift` returns early — an exact
no-op, not a resample — so the mature cells of both arms are **i.i.d. draws,
with replacement, from the same finite pool `F`** of the held-out patients'
mature cells (`_draw_cells` uses `replace=True`, and both arms share one
`eligible` mask). Hence:

* the parametric intrinsic term is exactly zero;
* `intrinsic_hat` is exactly unbiased for it;
* **null rejection is exactly the non-coverage of 0 by a two-sample percentile
  bootstrap of a scaled difference of means**, with `n_n = round(0.40 · 2000)
  = 800` held fixed and `n_t` the swept quantity.

Everything below is an experiment on that statement.

## 2. The baseline being explained

From `papers/whms_bode/diagnostics/low_count_diagnostic.csv`, produced by
`papers/whms_bode/diagnose_low_counts.py` at `R = 200` replicates, one seed
stream, inner `B = 200`. Null (`shift = 1.0`) rejection against a nominal 5%:

| cohort / pool | 5 | 50 | 100 | 800 |
|---|---|---|---|---|
| SMC / pooled | 63.0% | 28.5% | 18.0% | **9.5%** |
| SMC / reference | 37.5% | **12.0%** | 7.0% | 7.0% |
| KUL3 / pooled | 67.0% | **30.0%** | 19.0% | 4.5% |
| KUL3 / reference | 61.5% | 15.5% | 16.0% | 6.5% |

The reviewer's "12–30% at 50 cells" is the range across those four cells; the
"9.5% at 800" is SMC/pooled. At `R = 200` the binomial MCSE is up to 3.5 pp,
which is why this registration re-measures at `R = 2000` before attributing
anything: at 200 replicates, SMC/pooled 9.5% and KUL3/pooled 4.5% at 800 cells
are not separated from each other, let alone from nominal.

**Pre-registered expectation, recorded so it can be wrong:** the committed
grid's `R = 200` will prove to have overstated the spread across cohort/pool
cells, and the 800-cell rates will come in closer together and closer to
nominal than the table above suggests. This is a prediction, not a result, and
it is registered here so that confirming it cannot be presented as a discovery.

## 3. Task 1 — the diagnosis. Causes tested by construction

Each cause is tested by **substituting the resampling population `F`**, or by
**changing the design**, and changing nothing else. Not by argument.

### 3.1 Pool families (substituted for `F`)

| id | pool | removes | keeps |
|---|---|---|---|
| `P0` | `empirical` — `F` itself | nothing (baseline) | everything |
| `P1` | `zeros_removed_mean_matched` — `F` restricted to strictly positive values, rescaled by a constant so its mean equals `mean(F)` | the zero atom | the mean, the positive-part shape |
| `P2` | `skew_matched_no_zeros` — continuous strictly-positive distribution matched to `mean(F)`, `var(F)`, `skew(F)` | the zero atom **and** discreteness | mean, variance, skewness |
| `P3` | `gaussian_matched` — `Normal(mean(F), var(F))` | zeros, discreteness **and** skew | mean, variance |

`P3` is the floor: what survives there is the `z`-vs-`t` and plug-in-sd term of
§0, and it is checked against `expected_false_positive_rate` at the effective n
as a consistency check on the whole apparatus (falsifier F5).

### 3.2 Design factor (the reviewer's third candidate)

| id | design | |
|---|---|---|
| `D_frac` | `fixed_fraction` — `n_n = 800` always, `n_t` swept | what the committed sweep does |
| `D_bal` | `balanced` — `n_n = n_t` = the swept count | cells per arm held fixed instead of the fraction |

### 3.3 Attribution, and how it is reported

At **n_t = 50** and **n_t = 800**, in percentage points:

```
zero inflation             = R(P0) − R(P1)
skew, given no zeros       = R(P2) − R(P3)
non-Gaussian shape, total  = R(P1) − R(P3)
fixed-fraction design      = R(D_frac) − R(D_bal)   (at matched pool)
floor (z/t, plug-in sd)    = R(P3) − 5%
residual, unexplained      = R(P0) − 5% − (sum of the above)
```

Both the **sequential** (order-dependent) attribution and the **single-factor**
deltas are reported. If they disagree, the causes are not additive, and that is
reported in those words rather than resolved by silently picking an order.

Both pool readings (`pooled`, `reference`, per `calibration_gap.POOLS`) and both
cohorts (SMC = GSE132465, 10 patients; KUL3 = GSE144735, 6 patients) are run and
reported separately — CLAUDE.md invariant 4. The reviewer's two headline numbers
come from different readings and are not compared across them without saying so.

## 4. Task 2 — the candidate repairs. Frozen list

All applied to the **same generated samples**, from the same seeds, so that
differences between them are differences between methods and nothing else.

| id | method | notes |
|---|---|---|
| `R0` | `percentile_cells` | the status quo. Baseline, not a candidate. |
| `R1` | `studentised_cells` (bootstrap-t) | pivot `t* = (θ* − θ̂)/SE*`; `SE` the delta-method Welch SE of the scaled difference of means, recomputed on each resample. |
| `R2` | `bca_cells` | two-sample BCa; `z0` from the bootstrap draws, acceleration from a jackknife over the cells of both arms. |
| `R3` | `welch_t_cells` | the non-bootstrap reference: the cell-level analogue of the method `interval_calibration.py` found calibrated over patients. |
| `R4` | `patient_cluster_bootstrap` | resamples **patients** (as clusters of cells) within each arm, not cells. |
| `R5` | `pseudobulk_patient_t` | one value per (patient, arm) = that patient's mature-cell mean; Student-t over patients. The pseudobulk-based interval. |

`R4` and `R5` need more than one patient per arm. The committed sweep holds out
`n_held_out = 2`. Both are therefore run at **`n_held_out ∈ {2, 5}`**, so that
"this method fails" is separated from "two clusters is not a sample". A
patient-level interval built on 2 or 5 patients inherits exactly the small-sample
problem `interval_calibration.py` quantifies (`width_ratio(2) = 0.11`,
`width_ratio(5) = 0.62`). **If that is why `R4` and `R5` fail, it is the more
interesting answer and will be reported as the finding, not as a footnote.**

**The list is frozen here.** Any candidate added after a result table is seen
will carry `post_hoc = true` in the results and be named as such in prose.

## 5. Outcome measures

Per (cohort, pool, candidate, count, design, holdout size):

| measure | measured at | target |
|---|---|---|
| `null_rejection` | `shift = 1.0` | **≤ 5%** (nominal α) |
| `coverage` | `shift = 0.5` (the pre-registered detectable effect, `calibration.PREREGISTERED.detectable_shift`) | **≥ 90%** |
| `discrimination` | `shift = 0.5` | **≥ 80%** |
| `median_width` | both | reported, and as a ratio to `R0` |
| `abstention_rate` | both | reported |

Coverage at the detectable effect and discrimination are reported **in the same
table row** as every null-rejection number, so a candidate that buys calibration
by losing all power is visible as such and cannot be quoted without its cost.
This is the rule `interval_calibration.check_power_carries_its_own_calibration`
already enforces one layer up, applied here.

A candidate **passes at a count** iff, at that count, **all three** hold:

1. `null_rejection ≤ 0.05 + 2·MCSE` (not significantly above nominal),
2. `abstention_rate ≤ 0.10`,
3. `discrimination ≥ 0.80` at `shift = 0.5`.

Anything else is not a pass, whatever the null rate is.

## 6. Monte-Carlo precision

Every rate is reported with `MCSE = sqrt(p(1−p)/R)` **in the same table cell**.

`R = 2000` replicates per cell, drawn as **5 independent seed streams × 400
replicates**, so seed-to-seed spread is visible rather than assumed away. At
`p = 0.05`, MCSE = 0.49 pp; at `p = 0.095`, MCSE = 0.65 pp. The 5% vs 9.5%
contrast the reviewer names is therefore resolvable at better than 5 MCSE.

Inner bootstrap `B = 2000`, held identical across candidates. The committed
sweep's `B = 200` is too few for a BCa tail (its 2.5% adjusted quantile would
rest on ~5 draws) and is not used here; `B` is not a free parameter to be tuned
per candidate.

Seeds are fixed and recorded. Results land in `results/<date>_<sha7>/` with
`.meta.json` sidecars via `src/common/io.py::write_versioned_table`, which
refuses a dirty tree — so code is committed before any table is written
(CLAUDE.md invariant 10).

## 7. Falsifiers

Registered in advance. Each is a result, and will be reported in these words.

* **F1 — nothing works.** No candidate meets the §5 pass criterion at any tested
  count, in either cohort, under either pool reading. Reported as: *the interval
  cannot be repaired at these cell counts and cohort sizes*, and the paper
  prescribes nothing.
* **F2 — a hollow pass.** A candidate meets the null target only by
  **abstaining** (`abstention_rate > 0.10`) or by becoming **uselessly wide**
  (`median_width > 3 × R0`, or `discrimination < 0.80` at the detectable effect
  where `R0` reached it). Recorded as `verdict = "passes_by_abstention"` or
  `"passes_by_width"`, never as a pass.
* **F3 — the diagnosis does not close.** More than **40%** of the excess
  `R(P0) − 0.05` at either focal count is left unattributed by the causes in
  §3.3. Reported as: *the diagnosis leaves most of the rate unexplained*, with
  the residual quoted in percentage points.
* **F4 — non-monotone or count-local success.** A candidate that passes at a
  small count but not at a larger one, or is non-monotone in `n` beyond
  2·MCSE, is reported as such and is not described as a fix.
* **F5 — the apparatus disagrees with its own arithmetic.** If `R(P3)` departs
  from `expected_false_positive_rate` at the effective n by more than 3 MCSE,
  the Gaussian floor is not behaving as the closed form says it must and the
  whole diagnosis is suspect. This is checked and reported **before** any other
  number.

## 8. House rules this registration is bound by

* **Every guard gets a committed failing input.** Any check introduced here has
  an input in `tests/test_checks_can_fail.py` that forces it to fail.
* **The vectorised estimator path is tested against the real one.** The
  candidate intervals compute `f_n · (mean_t − mean_n)` directly, for speed. A
  test asserts this is elementwise identical to
  `kitagawa.decompose(...).intrinsic`. If W4 changes `decompose`, that test
  fails rather than these results silently diverging from the estimator they
  claim to calibrate.
* **`None` is not `0.0`** (invariant 1). A candidate that cannot form an
  interval abstains; abstention is counted in its own column and never folded
  into a rate, following `calibration.coverage_and_discrimination`.
* **Estimate per study, never pool across studies** (invariant 4). SMC and KUL3
  are run and reported separately.
* **Nothing under `paper/` or `papers/` is edited.** Those trees are read for
  the baseline numbers and integrated by the paper owner.

## 9. Data provenance

The pool `F` is the real Lee cohort's mature-cell `GUCA2A` counts, loaded by
`calibration_gap.load_cohort_arrays`. The raw matrices are not checked into this
worktree; they are read from the primary checkout via `BRP_DATA_DIR`, at the
byte sizes recorded in `data/manifest.csv`. Every table carries
`pool_source = "lee_raw"` in its sidecar. Had the matrices been unavailable, a
surrogate would have been used and labelled `pool_source = "surrogate"`, with no
claim that the percentages transfer — that contingency did not fire.
