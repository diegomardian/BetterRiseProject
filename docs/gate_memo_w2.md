# W2 gate memo — draft

**Status: DRAFT ON SYNTHETIC DATA.** Every number here comes from simulated
cells. Nothing in it is a result about colorectal cancer, and none of it should
be read into the gate until the same runs have been done on real cells. What it
does establish is that the machinery works and would produce an answer.

Date: 2026-08-16 · Owner: W2 · Reads against
[execution_plan.md §5](../execution_plan.md#5-week-5-gate)

---

## Summary

| Gate | Status | Basis |
|---|---|---|
| **G1** ambient correction | **not W2's** | W1's retention-vs-abundance statistic. The harness starts from corrected cells and cannot test the correction. |
| **G3** estimator recovers known ground truth | **preliminary PASS**, synthetic | Oracle arm recovers the known split within 1% wherever the mature compartment is non-empty; interval coverage 0.90–1.00 above 20 mature cells |
| **G4** <50% of patients below threshold | **cannot answer yet** | Needs real patient-level mature-cell counts. `gate_g4_verdict()` is implemented and waiting on data. |
| **G2** control tiers separate | **not yet run** | Needs W1's pilot |

---

## 1 · G3 — does the estimator recover known ground truth?

**Preliminary pass.** On synthetic cells the oracle arm — W4's `decompose()`
run on cell-level summary statistics — recovers the parametric split to within
1% at every mature fraction above zero.

Two guards make that statement worth something rather than circular:

- **The harness can reject.** `tests/test_harness.py` injects an estimator that
  always returns `0.0` and requires G3 to fail it. A harness that has never
  rejected anything is not evidence.
- **The oracle arm goes through the real estimator.** An earlier version read
  the harness's own realised truth back out and scored a perfect 1.0 for reasons
  having nothing to do with the estimator. Fixed; the arm now calls
  `decompose()`.

---

## 2 · Cutpoints — calibrated for the first time

The week-5 deliverable was blocked on an interval that did not exist. W4's
`attach_intrinsic_ci` broadcasts a **cohort-level** band onto every patient row —
identical for a patient with 800 mature cells and one with 21 — so coverage
against `n_cells_mature` was flat and no cutpoint could be derived.
`src/harness/interval.py` supplies a within-patient interval instead
([open_decisions #10](open_decisions.md) explains why that is not an invariant-5
violation).

```
n_cells_mature  coverage  discrimination  median_ci_width  verdict
             5     0.475           0.500            2.583  not_estimable
            20     0.925           1.000            1.994  ok
            30     0.950           1.000            1.708  ok
            40     0.900           1.000            1.470  ok
           100     0.925           1.000            1.018  ok
           800     1.000           1.000            0.596  ok
```

CI width runs 2.583 → 0.596, a **4.3× range**. The curve responds to cell count,
which is the whole requirement.

| cutpoint | provisional | calibrated (synthetic) |
|---|---|---|
| `ok` | 50 | 20 |
| `wide` | 20 | 20 |

### Do not swap `CUTPOINTS` on this evidence

Three reasons, all of which point the same way:

1. **Synthetic.** The reference is built from the same generative process that
   made the data.
2. **Discrimination never binds.** It sits at 1.000 everywhere above n=5 because
   the simulated effect is enormous next to the noise — intrinsic ≈ −12 against
   a CI width of ~2. So coverage alone drives the cutpoint. On real data, where
   the effect is smaller relative to noise, discrimination will bind and push the
   cutpoint **up**.
3. **Grid resolution.** The crossing lies between 5 and 20 mature cells and the
   grid does not resolve it. "Calibrated ok = 20" means *20 is the smallest
   tested count that passes*, not *20 is where it crosses*.

The honest reading is that the provisional 50 is not yet contradicted. Rerun on
real cells with a denser grid between 5 and 50 before touching `positivity.py`.

---

## 3 · Bake-off — which deconvolver

Eight samples, 800-gene signature, five cell types, synthetic.

| method | r | RMSE | RMSE (mature) |
|---|---|---|---|
| ν-SVR | 0.9999 | 0.0012 | 0.0011 |
| NNLS | 0.9971 | 0.0051 | 0.0033 |

ν-SVR wins, roughly 3× on the mature fraction. **Both are far better than
anything real data will give** — clean synthetic mixtures with a perfect
reference. CIBERSORTx and BayesPrism have not been run; per the staging decision
their absence is recorded rather than silently dropped.

---

## 4 · Signature width — §2.1 error #4, partly settled

ν-SVR, same samples, best *k* markers by specificity (not the panel — invariant 2
keeps target genes out of any reference, so the narrow arm is the best 11
**non-target** markers, which is a stronger test than 11 random genes).

| genes | RMSE | RMSE (mature) |
|---|---|---|
| 800 | 0.0011 | 0.0011 |
| 250 | 0.0012 | 0.0012 |
| 50 | 0.0012 | 0.0013 |
| 25 | 0.0014 | 0.0016 |
| **11** | **0.0052** | **0.0058** |

**What this supports:** 11 genes is materially worse — about 4.7× on the mature
fraction. The Executive Brief's proposal to deconvolve on the panel is wrong in
the direction claimed.

**What this does not support:** the "500–2000 genes" figure. On this data 25
genes already recovers nearly all of it, and everything from 50 up is flat.
High dimensionality earns its keep against reference mismatch, ambient
contamination and many more cell types — none of which this cohort has. Treat the
500-gene floor as unvalidated here rather than confirmed.

---

## 5 · Negative controls

| control | intrinsic, relative to target |
|---|---|
| housekeeping | **0.005** |
| permuted | 0.635 — *read against label-blind, not against zero* |

**Housekeeping** genes show essentially nothing, as they must.

**Permutation does not go to zero, and should not.** Silencing 40% of cells
moves the mean of any random subset, so after shuffling, Δ(per-cell mean) is the
whole-sample difference, diluted. The correct test is against a **label-blind**
reference — the same arithmetic with each arm's mean taken over all cells:

```
permuted     intrinsic  7.9316
label_blind  intrinsic  7.9100     agreement within 0.3%
```

Under shuffled labels the estimator returns precisely what an estimator ignoring
labels returns. It extracts nothing the labels were carrying. That is a sharper
claim than "small", and it is the one that would catch an estimator reading
batch, depth or patient identity — any of those would push the permuted arm
*away* from label-blind.

An earlier version shuffled labels *before* generating, which silenced a random
subset of cells — a real effect, not a null — and reported 23%. That was a
mis-specified control reading as a partial estimator failure.

The **compositional** term is not testable this way: the generator imposes the
mature fraction when it draws cells, so a shuffle preserves it by construction.

---

## 6 · The finding that survives every caveat

```
bulk over-reporting: 41/80 not-estimable rows got a confident number
median |intrinsic_hat| where the truth is UNDEFINED: 103.50
```

At zero mature cells the intrinsic term is undefined. The oracle arm knows,
because it counts cells. Deconvolution assigns a non-zero mature fraction to a
sample containing no mature cells, the division goes through, and bulk returns a
large confident number on **half** the rows where the honest answer is "not
estimable".

No amount of reference quality fixes this. Bulk cannot count cells, so it cannot
apply a positivity rule at all — and for the same reason it cannot form a
per-patient interval either. Three faces of one structural limit, and the reason
the third segment needs single-cell data to exist.

---

## 7 · What the harness cannot say

Unchanged from [the design spec](harness_design_spec.md) §3, restated because a
gate is exactly where it gets forgotten:

- It cannot test ambient correction. That is G1 and it is W1's.
- It cannot test whether the labels are biologically right. If they are wrong,
  the harness will confirm the estimator faithfully recovers the wrong thing.
- It cannot promise a real tumour behaves like the model.

A G3 pass means *the estimator recovers what we told it to recover*. It does not
mean the answer is true.

---

## 8 · What this memo needs before the gate

1. Every run above, on real cells. Blocked on
   [open_decisions #8](open_decisions.md) — the `lee_io` raw-count accessor,
   proposed on `w2/lee-raw-counts` and waiting on W4.
2. A denser calibration grid between 5 and 50 mature cells.
3. G4 numbers, which need real per-patient mature-cell counts.
4. G2, which needs W1's pilot.
5. W4's view on [open_decisions #9](open_decisions.md) — the `doubly_robust`
   weighting folds the interaction into both arms, which invariant 7 forbids.
   Unresolved, and it changes what the decomposition reports.
