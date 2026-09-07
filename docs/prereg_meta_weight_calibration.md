# Pre-registration — the meta layer weights by a quantity with no finite mean

**Written:** 2026-09-07 · **Author:** W1 (Bode) · **Status:** proposed ·
**Concerns** `src/harness/meta.py` (W2-owned — see §8) and the committed
`results/2026-09-05_61ba221/` meta premise.

> **This is not blind, and pretending otherwise would be the defect it
> describes.** §2 lists exactly what was computed before this document existed
> and what it showed. The rule in §4 is fixed by a closed form that does not
> depend on any estimate, and §6 states in advance that applying it may
> overturn a committed verdict.

---

## 1 · What this is

`src/harness/meta.py` pools per-study estimates by inverse-variance weighting:
`w = 1/se²`, with `se` the standard error of the mean over patients
(`coexpression_meta.per_study_stats`, `sd(ddof=1)/sqrt(n)`). DerSimonian–Laird
then treats those SEs **as known constants**.

They are not known. They are estimated from as few as **three patients**, and

    w_hat / w_true  =  1 / (chi²(n-1) / (n-1))

so the weight is an *inverse* chi-square. Its expectation is

    E[w_hat] / w_true  =  (n - 1) / (n - 3)          [finite only for n >= 4]
    Var[w_hat] finite                                 [only for n >= 6]

**At n = 3 the inverse-variance weight has no finite mean.** `E[1/chi²(2)]`
diverges. This is not a small-sample approximation being loose; the quantity
the estimator averages does not have an average.

**And the floor sits one patient below that.** `MIN_PREMISE_PATIENTS = 3`
(`coexpression_silencing.py:105`) is documented as *"3 needed to tell a shift
from noise"* — a statement about whether the **estimate** is computable. Nothing
asked whether the **weight** derived from it is integrable. That is the defect,
in one line: **the floor was set by the estimability of the numerator and used
to admit a denominator.**

### The two things it is not

**It is NOT a claim that KRT8's verdict is wrong.** §2 records that the
observed `I² = 87.6%` sits at null `p = 0.015` — it survives calibration at the
committed floor. Whatever this document changes, it does not change that.

**It is NOT a licence to drop a study.** Khaliq_2022 is the most influential
point and it will be tempting to remove it. Post-hoc outlier removal is the
defect this repository catalogues. §4 fixes a floor by arithmetic that was
written down before the floor curve was computed, and §5 requires the **whole
curve** to be reported, not the floor that gives the nicest answer.

## 2 · What was already computed, before this document existed

Declared because the alternative is a fake blind. All of it is a re-read of
`results/2026-09-05_61ba221/` and `results/2026-09-05_3380d15/`; no new data.

**Design-side — depends on the eleven `n`s, not on any estimate:**

| | |
|---|---|
| `E[w_hat]/w_true = (n-1)/(n-3)` | verified against 4M draws to three decimals for n >= 4; at n = 3 the running mean wanders 9.7 → 22.2 → 13.0 → 16.5 across 10⁴ → 10⁷ draws and does not settle |
| null `I²` at the committed `n`s | median **0.269**, 90th pct 0.701, 95th pct 0.785 — under *perfect homogeneity* |
| `P(I² > 0.75 | homogeneous)` | **6.9%** |
| `P(Cochran's Q rejects at 0.05)` | **32.5%**, against a nominal 5% |
| the same eleven studies if each had n = 29 | median null `I²` = **0.000**, `P(I² > 0.75)` = 0.02% |

That last row is the whole point: **essentially all of the baseline `I²` is
manufactured by the small studies**, through the weight and not through the
biology.

**Data-side — post-hoc, and labelled as such wherever it is quoted:**

| | |
|---|---|
| Khaliq_2022 share of Cochran's Q (KRT8) | **57.1%** |
| `I²` without it | **60.7%** — below the ceiling |
| Khaliq fixed-effect weight, n = 3 | **20.4%** |
| Pelka fixed-effect weight, n = 29 | **10.0%** |
| Khaliq's own ACTB control shift | **−0.589**, the only one of eleven outside the ±0.5 tolerance |
| observed KRT8 `I² = 0.876` against the null | `p = 0.015` |
| observed ACTB `I² = 0.628` against the null | `p = 0.158` — the **84th percentile of homogeneity**, reported as a number suggesting substantial heterogeneity |

**One explanation was tried and the data refused it.** If Khaliq's KRT8 shift
were part of a study-wide detection-scale drift, per-study KRT8 and ACTB would
track each other. Across all eleven, Pearson `r = 0.486` (p = 0.130) —
**and `r = −0.176` (p = 0.627) with Khaliq removed.** The correlation *is*
Khaliq. Global scale drift is not the explanation and is not claimed.

## 3 · Why this was not caught

`meta.py`'s docstring says `MAX_I_SQUARED = 0.75` is *"Cochrane's conventional
'considerable heterogeneity' boundary, committed here rather than chosen after
seeing a value."* Pre-committed — correctly, and that part worked. But a
threshold committed in advance is still uncalibrated if nobody computes what it
does under the null **at the `n`s the design actually has**. Higgins' 25/50/75
are rules of thumb for a literature of trials with hundreds of subjects each.
At n = 3 the same number means something else.

**Same shape as §3a of `docs/HANDOFF.md`, one layer up.** There, an interval was
0.82× the width it claimed, by a factor that is a function of `n` alone. Here, a
weight is `(n-1)/(n-3)` times what it claims, also a function of `n` alone. Both
are arithmetic, both were invisible because nothing raised, and both make a
check report a verdict it has not earned. §3a's version made a check too eager
to find an effect; this one makes a gate too eager to refuse.

## 4 · The rule, fixed here, by arithmetic that does not see the answer

**The patient floor for entry into an inverse-variance meta-analysis is set by
the integrability of the weight, not the computability of the estimate.**

| floor | criterion | status |
|---|---|---|
| n >= 4 | `E[w_hat]` finite (`chi²` df > 2) | **the minimum defensible floor** |
| n >= 6 | `Var[w_hat]` finite (`chi²` df > 4) | **the pre-registered primary** |

**Primary is n >= 6.** A weight with a finite mean but infinite variance is
still a weight that can take any value; DL treats it as a constant, and the
whole apparatus of Q, tau² and I² assumes it is one. n >= 6 is the smallest
floor at which the weight has two moments, which is the least that assumption
needs. It is chosen from the chi-square df condition, and that condition was
written down in §1 before any floor curve was run.

**n >= 4 is reported alongside as the weaker sensitivity**, not as a fallback to
switch to. n >= 3 (the status quo) is reported so the change is visible.

**Nothing else changes.** Same statistic, same estimator, same tolerance, same
`MAX_I_SQUARED`, same studies otherwise. Only the entry floor moves, and it
moves for a reason stated as a closed form.

## 5 · What is reported — the curve, never a point

`config/labeling_axes.yaml` already requires this of the rung axis: *"A single
point estimate would present a modelling choice as a measurement."* The floor is
a modelling choice of exactly that kind, and it gets the same treatment.

For **each** floor in {3, 4, 6} and **each** control gene, report together:

1. `k`, and which studies were dropped by name
2. pooled estimate, CI, tau², `I²`, Q
3. **the null `I²` distribution at that floor's `n`s** — median and 95th
   percentile, from the same simulation, in the same row
4. the calibrated `p` for the observed `I²` against that null
5. the resulting `premise_verdict`

**A heterogeneity verdict may not be reported without its own null.** Rows 2 and
5 do not exist in a frame without rows 3 and 4, enforced the way
`check_power_carries_its_own_calibration` enforces its analogue in
`interval_calibration.py` — a function that refuses the frame, with the failing
input committed in `tests/test_checks_can_fail.py`.

## 6 · What this may do, said before it is run

**Raising the floor to n >= 4 removes Khaliq (n = 3) and Qi (n = 3). §2 records
that `I²` without Khaliq alone is 60.7%. So this rule will probably move KRT8
from UNRESOLVED to a readable pooled estimate, and that would change a
committed conclusion.**

Stating it here is the point. If it happens it must be read as **"the verdict
was not stable to a floor set by a criterion nobody had applied"** — which is a
statement about the instrument — and **not** as "KRT8's controls hold after
all," which the data does not support at any floor.

| branch | reading |
|---|---|
| KRT8 stays UNRESOLVED at n >= 6 | The heterogeneity is real and is not a weighting artefact. The committed conclusion stands and gains a calibration it did not have. |
| KRT8 becomes readable at n >= 6 **and** at n >= 4 | The verdict was carried by studies whose weights are not integrable. **The 13-study closure of the carcinoma question is weaker than `docs/NEXT_AVENUES.md` records it** and that document must be corrected. |
| The two floors disagree | Report both, claim neither. The instrument does not resolve it. |
| ACTB changes verdict in any direction | Larger than KRT8: ACTB is the control the premise is *built* on. Report prominently. |

**In every branch the pooled estimate remains gated by the premise**, and in
every branch the per-study table is the primary object. Nothing here licenses
reading detection deltas that were gated before.

## 7 · What this does not touch

The **adenoma** results (avenue A, §6h, §6j) do not go through `meta.py` at all
— they are single-cohort, Student-t over patients. This document has no bearing
on them. It concerns exactly one committed object: the 13-study carcinoma meta
premise.

## 8 · Ownership

`src/harness/meta.py` is **W2's** under CONTRIBUTING §2 and this document does
**not** modify it. The calibration lives in a new W1-owned module and the floor
is applied at the job layer, which is W1's. If the result argues for changing
`MIN_STUDIES` or `MAX_I_SQUARED` themselves, that is a W2 PR with two approvals
and is out of scope here.
