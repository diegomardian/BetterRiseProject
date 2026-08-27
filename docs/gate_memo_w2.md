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
| **G1** ambient correction | **changed twice — premise and statistic; W2 has ratified the new statistic, see §0** | MA-vs-percentile per prereg amendment 2. The harness can say the criterion discriminates; it starts from corrected cells and still cannot test the correction itself. |
| **G3** estimator recovers known ground truth | **preliminary PASS**, synthetic | Oracle arm recovers the known split within 1% wherever the mature compartment is non-empty; interval coverage 0.90–1.00 above 20 mature cells |
| **G4** <50% of patients below threshold | **cannot answer yet** | Needs real patient-level mature-cell counts. `gate_g4_verdict()` is implemented and waiting on data. |
| **G2** control tiers separate | **not yet run** | Needs W1's pilot |

---

## 0 · G1 has changed twice, and this is the whole picture

Neither change is W2's criterion, but between them they reshape what the gate is
deciding, so they belong at the top rather than in a footnote. W1 asked
([#37](https://github.com/diegomardian/BetterRiseProject/issues/37)) that the two
be folded together so the gate meets one coherent G1 rather than two surprises
arriving separately. This is that.

### 0.1 · The premise changed on 2026-08-16 — W1's finding

W1 established ([open_decisions #11](open_decisions.md)) that **GSE178341 ships
no unfiltered droplets** — the deposit is post-`dropletUtils`, all 181 GSM
records carry no per-sample supplementary files, and the one unlinked
Broad-hosted lead turned out to have exactly 371,223 rows, the published post-QC
count. There are no empty droplets in any public source.

Consequences, in the order they bite:

- **CellBender cannot run at all.** It models every barcode including cell-free
  ones. Without empty droplets there is nothing for it to learn from.
- **SoupX runs only in degraded mode**, and DecontX replaces CellBender as the
  second method since it infers contamination from cluster structure instead.
- **G1 as written is no longer quite the question.** It asks whether ambient
  correction eliminates the intrinsic signal. We can no longer run the
  correction the criterion assumed, so the gate is choosing between *degraded
  correction* and *no correction*, not between two good methods.

That does not obviously fail G1 — the plate-based subset in the ICBI atlas has
essentially no soup, so an intrinsic signal surviving there remains strong
evidence it is not contamination, and that route is untouched by this. But the
pre-committed consequence for a G1 failure (pivot to snRNA-seq and spatial) is
now closer than it was on the 15th, and the gate should be told so rather than
discovering it in the room.

### 0.2 · The statistic changed on 2026-08-25 — and W2 has now checked it

W1's [prereg amendment 2](prereg_amendment_2_g1_tier_d.md) reports that decision
#17's statistic is arithmetically broken: apparent loss measured as a raw Δ of
per-cell means carries abundance inside it, and tier D — the half of G1 that
carries the falsification logic — holds one gene, over which a correlation is
undefined. It proposes the MA construction genome-wide with panel genes located
by percentile of M within their own abundance bin, and three thresholds committed
before any real M or A exists.

**W1 cannot ratify its own pre-registration.** The part of the team that holds
ground truth is the harness, so W2 checked it: `src/harness/g1_amendment.py`,
`tests/test_g1_amendment.py`, five simulated worlds with known answers.

**W2 ratifies the amendment, with one premise that should be written into it.**

**Confirmed — the old statistic is broken, and worse than reported.** On a world
where every gene loses exactly 30% and nothing depends on abundance anywhere,
#17's statistic returns **ρ = −0.997** while the proposed one returns **+0.006**.
Reproduced through W2's own construction, not by running W1's script.

W2 adds a second defect W1 did not report. `scipy.stats.spearmanr` on one
observation returns **nan**, and it does not raise. `abs(nan) > 0.5` is `False`,
so #17's rule — *fail if |ρ| > 0.5* — would not have errored on tier D. **It
would have reported PASS**, for the one tier whose whole job is to be able to
fail. That is the fifth instance of this repository's recurring defect, and the
first one located in the gate itself. `amendment2_verdict` refuses a non-finite
percentile rather than comparing it.

**Confirmed — the replacement can fire.** This is the check the amendment does
not contain, and the reason it matters is that showing the old statistic fails a
null establishes only that the old one is broken. On a world where loss is a
function of abundance and nothing else, G1a returns **|ρ| = 0.94**; on genuine
biology it returns **0.08**. It discriminates.

**Confirmed — the three thresholds discriminate.** 60 replicates per world,
seed 20260827,
[`results/2026-08-26_e77907b/g1_amendment_ratification.parquet`](../results/2026-08-26_e77907b/g1_amendment_ratification.parquet):

| world | truth | P(PASS) | owed |
|---|---|---|---|
| `broad_loss_tier_d_retained` | the project's claim is true | **1.000** | PASS |
| `isolated_tier_a_loss` | tier A gone, nothing else moves | **0.517** | PASS — see below |
| `uniform_loss` | everything loses 30%, no biology | 0.017 | FAIL |
| `pure_soup` | loss is abundance and nothing else | 0.017 | FAIL |
| `tiers_drift_together` | MS4A12 lost as hard as tier A | 0.000 | FAIL |

Under a pure null the gate passes **2.6%** of the time. Per threshold: tier A
5.1%, tier D **50.1%**, separation 21.2%.

### The caveat the gate should hear, because it is not in the amendment

**Threshold 2 asks MS4A12's within-bin percentile to be ≥ 0.50, and a gene that
is unchanged against an unchanged background sits at 0.50 by definition.** So on
`isolated_tier_a_loss` — a world where the project's claim is *true* — the gate
is a coin flip: P(PASS) = 0.517, and the 50.1% null rate above is the same fact
seen from the other side. Against a broadly-lost background MS4A12 sits at
0.890 ± 0.015 and P(PASS) = 1.000.

**G1's power is therefore not a property of G1.** It depends on there being a
broad loss background for the retained control to stand out against. That premise
is the project's own — differentiation-marker loss is broad — and is probably
satisfied. It was simply never stated, and if it fails, G1 rejects a true signal
half the time for reasons that have nothing to do with ambient RNA.

Amendment 2 made tier D **computable**. It did not make n = 1 **powerful**, and
those are different achievements. This is not a reason to reject it: the
statistic it replaces returned ρ ≈ −1 on a null and could not be evaluated on
tier D at all. It is a reason to write the premise into the pre-registration and
to read a G1 **FAIL** with the caveat attached.

**W2 has not proposed different numbers.** Adjusting a pre-committed threshold
during its ratification is precisely the move the amendment exists to avoid, and
`tests/test_g1_amendment.py` pins 0.20 / 0.50 / 0.30 so that W2 cannot drift them
later either.

### W2's answers to the four ratification questions on #37

1. **Loss becomes log₂ fold change, abundance becomes the MA average — yes.**
   Independently confirmed, and the rank construction has a property worth
   recording: within-bin percentiles are invariant to library-size
   normalisation, because a global rescale of one arm shifts every M by the same
   constant and leaves within-bin ranks untouched (max movement 0.016 under
   CP10K). A mean-based rule would not have been.
2. **Genome-wide with within-bin percentiles, for impossibility rather than
   convenience — yes.** With one gene you cannot compute a correlation but you
   can compute a percentile; that is a real argument and it survives checking.
3. **The three thresholds — yes, with the premise in §0.2 written into the
   pre-registration alongside them.**
4. **Tier B reported but not gate-bearing — yes.** n = 3 cannot support a
   threshold, and MLH1's M carries a different meaning from tier A's.

**What W2 still cannot do:** adjudicate G1's substance. The harness starts from
corrected cells and has no view on whether the correction worked. Everything
above is about whether the *criterion* can distinguish worlds — which is a
question the harness can answer, and the one that was open. What it still cannot
say is which world this cohort is in. Quantifying how much intrinsic signal
survives a given level of simulated ambient contamination would turn "degraded
correction" into a number; that is not built and not in week-5 scope, and it
remains the obvious next harness job if the gate wants it.

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
