# Harness design spec

**W2, week 1. Reviewed by all four before the generator was written.**

Status: **draft for review** · Owner: W2 · Reviewers: W1, W3, W4

This document defines what "known ground truth" means in this project. It is the
load-bearing assumption of the entire week-5 gate: if ground truth is defined
wrongly, G3 passes or fails for the wrong reason and nobody finds out until the
paper is being written. Everything numeric in here is pre-registered — fixed
before the attenuation sweep runs, so the cutpoints at week 5 are *derived* and
not chosen.

---

## 1 · What "known ground truth" means here

The harness builds matched (normal, tumour) pseudobulk pairs from real cells,
where we control the composition and the per-cell silencing. Two truths come out
of that, and **both are recorded on every sample**:

| | What it is | What recovery against it measures |
|---|---|---|
| **Parametric** | the composition and shift we asked for | estimator bias **plus** sampling noise |
| **Realised** | the empirical fractions and means in the cells actually drawn | estimator bias alone |

They differ because an integer number of cells cannot hit an arbitrary fraction,
because a cell type may be thin in the held-out pool, and because counts are
random. That difference is not a nuisance to be minimised — it is the reason
both are stored. **Report only one and a sampling artefact reads as estimator
bias, or the reverse.**

Implemented in [`src/harness/truth.py`](../src/harness/truth.py) as
`parametric_truth()` and `realised_truth()`, both attached to every
`GroundTruth`.

---

## 2 · The generative model

Cells come from **held-out patients** — the split is on `patient_id`, never on
cells. Invariant 5 is written about the bootstrap, but it holds here for the
same reason: cells within a patient are not independent draws.

1. Fix a normal composition `f_n` and a tumour composition `f_t` over cell types
   at the chosen rung. The mature fraction in `f_t` is the knob the attenuation
   sweep turns; it goes to zero deliberately.
2. Draw `N` cells per sample to match the composition. The realised mature count
   is what `classify_estimability()` sees.
3. Scale the mature-cell mean of each target gene by a known factor `s`, in the
   **tumour sample only**.
4. Sum counts across cells to make the pseudobulk pair.

### 2.1 How a multiplicative shift reaches integer counts

The shift is defined on the **mean**. Real cells give counts, so the factor is
realised by a mechanism that preserves the expectation exactly:

| `s` | mechanism | `E[new]` |
|---|---|---|
| `== 1` | untouched | `counts` |
| `< 1` | `Binomial(counts, s)` | `s · counts` |
| `> 1` | `counts + Poisson(counts·(s−1))` | `s · counts` |

`s = 1` is a genuine no-op and not a resample, so **the null is exactly null**,
not merely small. That is what lets the `s = 1.0` row of the sweep serve as a
control rather than as a second measurement.

### 2.2 The identity, on paper

With `m_n = μ`, `m_t = μ·s`, `Δf = f_t − f_n`, `Δm = μ(s−1)`:

```
normal weighting:   Δf·m_n  +  f_n·Δm  +  Δf·Δm   =  f_t·m_t − f_n·m_n   ✓
tumour weighting:   Δf·m_t  +  f_t·Δm  −  Δf·Δm   =  f_t·m_t − f_n·m_n   ✓
```

Both close exactly. `assert_identity_closes()` checks it on every generated
sample, because a truth that is not self-consistent cannot judge anything.

> **Correction to execution_plan.md §4.** The boxed definition there reads
> `compositional = Δ(mature fraction) × normal per-cell mean` and
> `intrinsic = tumour mature fraction × Δ(per-cell mean)`. Those two come from
> *different* weightings and do not sum to the total with a single interaction
> term. Each weighting closes on its own, as above, and
> [`kitagawa.decompose()`](../src/estimator/kitagawa.py) implements them
> coherently one at a time. The prose is a loose gloss on "the split is not
> unique" rather than a single triple. **No code changes — W4's implementation
> is right. Worth a line in the methods section so a reader does not try to add
> the boxed three together.**

---

## 3 · What the harness cannot test

Written down in week 1 so it is not argued about in week 5.

- **Ambient contamination.** That is G1 and it belongs to W1's
  retention-vs-abundance statistic. The harness starts from cells that have
  already been corrected; it cannot tell you the correction worked.
- **Whether the labels are biologically right.** It can tell you the estimator
  recovers a split defined by *these* labels. If the labels are wrong, the
  harness will confirm the estimator faithfully recovers the wrong thing.
- **Whether real tumours behave like the model.** Silencing in vivo is not a
  clean multiplicative factor applied uniformly to a cell population. The
  harness bounds estimator error, not biological realism.
- **The direction of the attenuation bias in real bulk.** It measures attenuation
  *on pseudobulk built from our own reference*, which is the optimistic case.
  Real bulk has platform differences the harness does not simulate.

A pass on G3 means "the estimator recovers what we told it to recover." It does
not mean the answer is true.

---

## 4 · Pre-registered numbers for cutpoint calibration

**Fixed now, before any sweep output is looked at.** This is the entire reason
the week-5 cutpoints are derived rather than chosen.

| Quantity | Value | Why this value |
|---|---|---|
| Detectable effect | `s = 0.5` | halving of per-cell output; middle of the published attenuation band (×0.6–0.8) and well inside what methylation silencing produces |
| Coverage target | 95% CI covers truth ≥ **90%** of the time | a nominal-95% interval that undercovers below 90% is not an interval |
| Discrimination target | CI excludes zero ≥ **80%** of the time at `s = 0.5` | conventional power, stated as a property of the interval rather than of a test |

Cutpoints follow mechanically:

| Cutpoint | Definition |
|---|---|
| `ok` | smallest `n_mature` with coverage ≥ 90% **and** discrimination ≥ 80% |
| `wide` | smallest `n_mature` with coverage ≥ 90% but discrimination < 80% |
| below `wide` | coverage fails, or CI width exceeds the parameter's plausible range → **not estimable** |

`calibrate_cutpoints()` returns a `Cutpoints` whose `source` field names the
sweep that produced it, so the number in the code traces back to a run.

**Both the provisional (50 / 20) and the calibrated values go in the gate memo.**
If they differ substantially, that is a finding about the provisional numbers.

---

## 5 · The attenuation sweep

- **Swept:** mature-cell fraction, log-spaced ≈ 0.001 → 0.5 (the x-axis); shift
  `s ∈ {1.0, 0.8, 0.5, 0.25}`.
- **Fixed:** cells per sample, sequencing depth, the reference S matrix, rung,
  axis.
- **Replicates:** 50–100 per grid point, each a fresh held-out patient draw with
  its own seed.
- **Output:** the `attenuation` table in
  [`src/harness/results.py`](../src/harness/results.py) — one row per (grid
  point, replicate), with both truths and the recovered terms.

The curve's *shape* is a result and is therefore asserted nowhere. What **is**
asserted: the `s = 1.0` row recovers zero within CI at every mature fraction. If
that fails, the harness is broken and nothing else in the sweep means anything.

---

## 6 · Negative controls

| Control | Expected | Implemented in |
|---|---|---|
| Within-patient label permutation | destroys **both** terms | `controls.permute_labels_within_patient` |
| Housekeeping genes | shows **neither** term | `controls.housekeeping_panel` |
| Broken estimator (returns 0.0 always) | **must fail G3** | `tests/test_harness.py` |

The permutation is within-patient so it destroys the label–expression
association while leaving composition and every batch effect intact. A
between-patient shuffle would scramble composition too and would not isolate
what we need isolated.

The third one is the important one. A harness that has never rejected anything
is not evidence that the estimator works — it is an absence of evidence that it
does not. So the test suite injects an estimator that always returns zero and
requires the harness to reject it.

---

## 7 · Harness self-validation, 2026-08-15

First end-to-end run. **Synthetic cells** — 12 patients, 5 cell types, 4200
cells, 601 genes, GUCA2A mature-restricted, patient-level scaling so patients
are not interchangeable draws. 28 grid points × 20 replicates × 2 arms.

```
null parametric intrinsic exactly zero : True
null sampling noise / real effect      : 0.020

median recovered / true intrinsic
                     bulk                 oracle
shift                0.25   0.50   0.80   0.25   0.50   0.80
frac_mature_tumour
0.01                1.085  1.091  0.856  0.995  0.999  1.024
0.02                1.042  1.037  1.210  1.004  1.009  0.989
0.05                1.034  1.034  1.125  1.003  0.995  1.003
0.10                1.001  1.036  1.075  1.001  1.003  1.007
0.20                1.008  1.021  1.076  1.001  1.001  1.002
0.40                1.000  1.002  0.987  1.000  1.001  0.986
```

**What this does support.** The pipeline runs end to end and the oracle arm
recovers the known split to within 1% wherever the mature compartment is not
empty. That is a preliminary **G3 pass on synthetic data**: W4's `decompose()`
is arithmetically sound and the harness can measure it. The null behaves — the
parametric intrinsic is exactly zero and sampling noise is 2% of a real effect.

**What this does NOT support.** The bulk arm also came back near 1.0, and that
must not be read as "bulk recovers the intrinsic term fine." The synthetic
cohort has no ambient contamination, no batch or platform effects, clean
cell-type separation, and — decisively — **the reference is built from the same
generative process that made the data**. Published attenuation comes almost
entirely from reference mismatch, which this setup does not have. The bulk
numbers here say the algebra is right, not that the method is usable. The real
attenuation curve needs real cells; expect it to look nothing like this.

**The finding that does survive the caveat**, because it is structural rather
than a matter of reference quality:

```
bulk over-reporting: 41/80 not-estimable rows got a confident number
median |intrinsic_hat| where the truth is UNDEFINED: 103.50
```

At zero mature cells the intrinsic term is undefined. The oracle arm knows,
because it counts cells. Deconvolution assigns a non-zero mature fraction to a
sample containing no mature cells, the division goes through, and bulk returns
a large confident number — on **half** the rows where the honest answer is
"not estimable". No amount of reference quality fixes this: bulk cannot count
cells, so it cannot apply a positivity rule at all.

This is the thesis in one table, and it is why the third segment needs
single-cell data to exist. `attenuation.bulk_overconfidence()` computes it.

---

## 8 · Open for review

1. **Is `s = 0.5` the right detectable effect?** It is defensible from the
   attenuation literature, but if W4 expects MLH1 silencing to be closer to
   complete, a smaller `s` would make the cutpoints stricter and G4 harder.
2. **Should the sweep vary sequencing depth**, or hold it fixed and treat depth
   as a separate one-dimensional check? Fixed is cheaper and the depth effect is
   probably second-order next to mature fraction.
3. **50 or 100 replicates per grid point?** 100 halves the Monte Carlo error on
   the coverage estimate at the cost of roughly a day of compute. At 50, a
   coverage estimate of 0.90 has a standard error of about 0.04, which is
   uncomfortably close to the threshold it is being compared against.
4. **W2's cohort** — see [open_decisions.md](open_decisions.md) #6. Needs W4's
   agreement before any download.
