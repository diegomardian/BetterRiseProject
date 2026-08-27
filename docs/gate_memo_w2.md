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
| **G4** <50% of patients below threshold | **cannot answer yet — but the cohort can support the question, see §9** | Needs real patient-level mature-cell counts. `gate_g4_verdict()` is implemented, now reports a Wilson interval and a `resolvable` flag, and is waiting on data. At 36 matched patients the effective decision line is 33%, not 50%; at SMC's 10 it is 10% and G4 is not answerable. |
| **Ambient sensitivity** (feeds G1) | **measured**, synthetic | §10 — at the 10% exclusion cap real terms retain 94% and a compositional-only world acquires an intrinsic term worth 4.6% of its compositional one. The artefact is one-directional. |
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
say is which world this cohort is in.

Quantifying how much intrinsic signal survives a given level of simulated
ambient contamination — the job §0 kept flagging as the obvious next one — **is
now built: see §10.** At the 10% exclusion cap, ambient manufactures an apparent
intrinsic term worth about 5% of a real one, and it **cannot** manufacture a
compositional term at all. "Degraded correction" now has a number attached to it.

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
3. G4 numbers, which need real per-patient mature-cell counts. **The gate's
   precision at the real cohort sizes is no longer open — see §9.** What is
   still needed is the counts themselves.
4. G2, which needs W1's pilot.
5. W4's view on [open_decisions #9](open_decisions.md) — the `doubly_robust`
   weighting folds the interaction into both arms, which invariant 7 forbids.
   Unresolved, and it changes what the decomposition reports.

---

## 9 · What the gate costs at the real cohort sizes — done, and it changes one verdict's status

Everything downstream of the plan was designed for roughly **60 patients**. The
cohorts are smaller: **GSE178341 has 36 matched of 62**, **SMC has 10 paired, not
23**. W2 promised W4 this re-costing twice and had not done it. It is done.

Two quantities are affected and they are not interchangeable.

### 9.1 · G4's verdict — a proportion against a pre-committed line

G4 asks whether fewer than 50% of patients are below the positivity threshold.
That is a proportion, so its precision is binomial in the number of *patients*,
and 50% is the point where a proportion is least precisely estimated.

`gate_g4_verdict()` now returns a Wilson interval on `fraction_below` and a
`resolvable` flag beside `passes`. **The 0.50 rule is unchanged** — moving a
pre-committed line after seeing what n does to it is the move this project
refuses everywhere else. What changes is that "PASS" and "PASS, and this cohort
could not have said otherwise" stop reading as the same sentence.

| cohort | n | clean PASS needs | effective line | rule says |
|---|---|---|---|---|
| SMC paired | 10 | ≤ **1** patient below | **10.0%** | 50% |
| GSE178341 matched | 36 | ≤ 12 | 33.3% | 50% |
| the plan's assumption | 60 | ≤ 22 | 36.7% | 50% |

**Decision #19's PASS survives.** It recorded 6 of 36 matched patients below
threshold — 16.7%, 95% CI [7.9%, 31.9%]. A defensible PASS at n = 36 allows up to
12, so that verdict clears with room to spare. This was the thing worth checking
and it came back clean.

**G4 is not answerable on SMC alone.** At n = 10 a defensible PASS needs 1 or
fewer of 10 patients below threshold; at 2 the interval already straddles the
line. Anything SMC says about G4 is a point estimate its own size cannot support,
and the gate should treat it as such rather than as a second opinion.

**One correction to decision #19's own table.** It reported the mixed population
as `32/62 = 51.6% FAIL`. Read with an interval that is [39.4%, 63.6%] — it
**contains the 50% line**, so the mixed verdict is *indeterminate*, not a clean
FAIL. Decision #19 is still right and for the reason it gave: mixing populations
reports a cohort-design fact as a positivity finding. But "mixing flips the gate"
overstates it. Mixing flips the *point estimate* and produces a verdict the
cohort cannot resolve, which is a weaker claim than the table implied.

**What more patients do and do not buy.** P(G4 returns PASS), exact binomial:

| true fraction below | n=10 | n=36 | n=60 |
|---|---|---|---|
| 0.35 | 0.751 | 0.954 | 0.988 |
| 0.45 | 0.504 | 0.670 | 0.742 |
| **0.50** | **0.377** | **0.434** | **0.449** |
| 0.60 | 0.166 | 0.083 | 0.044 |

More patients sharpen the verdict *away* from the line, not *on* it. If the truth
sits near 50%, no cohort this project can reach turns G4 into a decision
procedure — the pre-committed consequence is what makes it one.

### 9.2 · The cohort band on the decomposition terms

Measured through W4's own `bootstrap_over_patients`, not from a formula, because
the estimator nulls the intrinsic term below the positivity cutpoint and a
formula does not know that. 15 replicate cohorts per size, 300 bootstrap draws,
seed 4242:

| term | n=36 vs plan | n=10 vs plan |
|---|---|---|
| compositional | ×1.22 | ×2.05 |
| intrinsic | ×1.23 | ×2.02 |
| interaction | ×1.27 | ×2.15 |
| *1/√n prediction* | *×1.29* | *×2.45* |

All three widen together, close to and slightly under the 1/√n expectation. So
the shortfall from 60 to 36 costs about a quarter of the interval width — real,
survivable, and much smaller than the G4 effect in §9.1.

**A finding that did not survive replication, recorded because it nearly went in
the memo.** At 5 replicates the interaction term appeared to degrade faster than
the other two (×1.69 and ×3.50 rather than ×1.27 and ×2.15), which would have
been a genuinely interesting claim about the term invariant 7 keeps separate. At
15 replicates it is gone: the three terms are indistinguishable. The first run
was noise, and it was noise pointing at an interesting conclusion, which is the
kind this project has learned to re-run before quoting.

### 9.3 · What follows for the gate

- **G4 should be read off GSE178341 matched patients only**, which is already
  decision #19, and now for a second and independent reason: SMC cannot resolve
  it and pooling the two would import that.
- **A G4 verdict should be quoted with its interval.** The machinery does this
  now; the memo and any results table should carry `resolvable` alongside
  `passes`.
- **The gate should be told that the effective line is 33%, not 50%**, at the
  cohort it will actually decide on.

`src/harness/gate_cost.py`, `tests/test_gate_cost.py`.

---

## 10 · "Degraded correction" is now a number — the ambient sensitivity sweep

§0.1 says GSE178341 ships no unfiltered droplets, so CellBender cannot run and
the gate is choosing between *degraded correction* and *no correction*. §0 has
been calling that "degraded" without saying how much, and flagged the sweep as
the obvious next harness job. It is built.

**This does not say whether this cohort's correction worked.** The harness starts
from corrected cells and cannot know. It says what a given level of *residual*
contamination does to a decomposition whose truth is known — which is the
question that makes the choice decidable.

### The mechanism, and why the answer is not symmetric

Ambient RNA is the sample's own average expression, redistributed across every
barcode. So the soup in a **normal** sample is rich in mature-colonocyte
transcripts, and the soup in a **depleted tumour** is not.

Take a tumour whose loss is purely compositional — mature cells gone, survivors
untouched, true intrinsic term exactly zero. Contamination pulls each arm's
mature cells toward *their own arm's* soup, the normal arm's soup is richer, and
the mature-cell means separate **where the truth says they must not**. A
compositional-only world acquires an apparent intrinsic term.

### The numbers, at decision #16's 10% exclusion cap

A sample at exactly 10% is *kept*, so this is the worst case the cohort tolerates
by design rather than a hypothetical. 20 replicates, seed 20260827.

| world | term | truth | at 10% ambient |
|---|---|---|---|
| compositional only | compositional | −17.54 | −16.49 — **94% retained** |
| compositional only | **intrinsic** | **exactly 0** | **−0.79 — manufactured, 4.6% of the compositional term** |
| compositional only | interaction | exactly 0 | +0.60 — manufactured, 3.4% |
| intrinsic only | intrinsic | −15.18 | −14.27 — **94% retained** |
| intrinsic only | **compositional** | **exactly 0** | **exactly 0.0000 — nothing manufactured** |
| both | all three | — | 94–96% retained |

### Two findings, and the second is the useful one

**1 · Real terms are attenuated by about 6% at the 10% cap**, and about 18% at
30%. Contamination shrinks a true signal; it does not annihilate or flip one.
That is a mild cost and it is the same for all three terms.

**2 · The artefact is one-directional, and structurally so.** Contamination
manufactures **intrinsic** signal out of **compositional** truth, and *never the
reverse* — the compositional term stays at exactly 0.0000 at every level tested.
That is not luck. The compositional term is a function of the mature-fraction
*difference*, contamination moves means rather than fractions, so where the two
arms share a mature fraction the term is identically zero however much soup is
added. **Ambient can invent silencing. It cannot invent depletion.**

The manufactured intrinsic term grows smoothly: 1.6% of the compositional term at
2% ambient, 2.8% at 5%, **4.6% at 10%**, 6.3% at 15%, 12.1% at 30% — against a
1.2% floor from sampling noise alone.

### What this means for the gate

- **The 10% exclusion threshold is well placed.** At the cap, ambient
  manufactures an intrinsic term worth about 5% of a real one. That is small
  enough not to change a conclusion and large enough that it should be stated
  next to any intrinsic estimate, not buried.
- **The two known artefacts push in opposite directions.** Invariant 6's bulk
  attenuation shrinks the intrinsic term (toward "compositional", which is the
  prior hypothesis). Ambient inflates it (toward "intrinsic"). They do not
  compound, and a result that survives both is stronger than one tested against
  either.
- **G1's premise change is survivable at this level.** Choosing *no correction*
  over *degraded correction* costs roughly 5% of the intrinsic term in
  manufactured signal at the worst contamination the cohort admits. That is a
  number the gate can weigh against the risk of a correction nobody can validate.

### The bound, stated so it is not assumed away

Contamination here perturbs **expression, not labels**. Enough soup also pushes
cells across a maturity threshold, which would move the *fractions* as well as
the means, and that is W1's axis question rather than W2's. **Every number above
is therefore a lower bound on ambient's total damage.** Simulating the label
effect here would hide W1's uncertainty inside a W2 number, which is worse than
leaving the bound visible.

### One bug this nearly shipped, recorded because the shape recurs

The first version decided "the truth is zero here" by testing whether the clean
term was near zero. In a compositional-only world the *realised* intrinsic term
is zero **plus sampling noise** — about −0.05 against a compositional term of
−17.5 — so the test called it non-zero, formed a retention ratio against a
near-zero denominator, and reported a confident **2.38×** where the honest
statement is "the truth is zero and 0.79 appeared". The flag now comes from the
regime's **design** (`AmbientRegime.parametric_zero_terms`), and
`summarise_ambient` refuses a sweep that has lost it.

That is the same error as scoring coverage against realised rather than
parametric truth (§3 of the handoff), one module along, and it produced a number
that looked like a finding.

`src/harness/ambient_sensitivity.py`, `tests/test_ambient_sensitivity.py`.
