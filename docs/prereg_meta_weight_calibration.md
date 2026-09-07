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

---

## RESULT — computed 2026-09-07, after this document was committed in `a9dc874`

`results/2026-09-07_40515ce/`. Job `src/reference/jobs/meta_floor_sensitivity.py`,
a local re-read of the committed per-study table. No new data.

### Verdict: FLOOR-UNSTABLE — and §5's branch table could not express the outcome

| gene | floor | k | pooled [95% CI] | I² | null median I² | calibrated p | verdict | **cause** |
|---|---|---|---|---|---|---|---|---|
| ACTB | n≥3 | 11 | +0.152 [−0.013, +0.317] | 0.628 | 0.270 | 0.159 | HOLDS | within tolerance |
| ACTB | n≥4 | 9 | **+0.201 [+0.061, +0.341]** | 0.458 | 0.137 | 0.225 | HOLDS | within tolerance |
| ACTB | n≥6 | 6 | **+0.212 [+0.040, +0.383]** | 0.607 | 0.000 | 0.062 | HOLDS | within tolerance |
| KRT8 | n≥3 | 11 | −0.453 [−0.711, −0.194] | **0.876** | 0.270 | 0.016 | UNRESOLVED | **I² over the ceiling** |
| KRT8 | n≥4 | 9 | −0.412 [−0.590, −0.235] | 0.625 | 0.137 | 0.090 | UNRESOLVED | **homogeneous, straddles tolerance** |
| KRT8 | n≥6 | 6 | −0.482 [−0.685, −0.280] | 0.715 | 0.000 | **0.018** | UNRESOLVED | **homogeneous, straddles tolerance** |

**KRT8 is UNRESOLVED at all three floors and it is not the same UNRESOLVED.**
At the committed floor the studies are declared *not to be estimating a common
quantity* and the pooled value is called "substantively meaningless." At both
floors where the weight has a mean, the studies **are** homogeneous by the
project's own ceiling, the pooled value is readable, and what it says is that
KRT8 fell by 0.41–0.48 log2 — close enough to the ±0.5 tolerance that the
interval straddles it.

**So the committed finding does not survive, but not in the direction §6
predicted.** §6 guessed KRT8 would "probably move to a readable pooled
estimate." It moved off heterogeneity and landed on the tolerance instead. The
substantive change is real and it is this: *"the thirteen studies disagree"*
becomes *"the studies agree, and what they agree on is a control shift the
tolerance cannot clear."* The first is a statement that no reading is possible.
The second is a reading, and it is worse for the premise.

### §5's branch table had the defect it was written to catch

Its four branches are phrased in verdict **labels** — "stays UNRESOLVED",
"becomes readable". `premise_verdict` reaches UNRESOLVED down two routes, so
"stays UNRESOLVED" was satisfied by an outcome that means the opposite of what
that branch says it means. **The pre-registration's own falsifier could not
distinguish the result that occurred.**

The job's first draft shipped the same mistake in code: `read_verdict` compared
labels and printed **STABLE — "the heterogeneity is a property of the data, not
of the weighting"**, which is precisely backwards. `verdict_cause` now compares
the route; the input that forces the failure is
`test_read_verdict_does_not_call_a_changed_cause_stable`.

**Recorded rather than quietly fixed, because it is the twentieth instance of
this repository's thesis and the first found inside a pre-registration.** A
three-valued verdict whose states are multiply-caused is a check that cannot
fail on the distinction that matters. §5 of any future prereg must name the
cause, not the label.

### The ceiling disagrees with its own null, and it does so in both directions

| k | null median I² | null 95th pct | P(I² > 0.75) | Q rejects at 0.05 |
|---|---|---|---|---|
| 11 (n≥3) | **0.270** | 0.787 | 6.97% | **32.5%** |
| 9 (n≥4) | 0.137 | 0.694 | 2.78% | 20.1% |
| 6 (n≥6) | **0.000** | 0.630 | 0.97% | 9.7% |

All three rows are under **exact homogeneity** — between-study variance is zero
by construction, so every value is manufactured by estimating the weights.

**At k=11 the fixed 0.75 is too strict**: homogeneity alone gives a median I² of
0.270 and clears the ceiling 7% of the time, so a reader seeing "I² = 62.8%,
within the ceiling" for ACTB is looking at the 84th percentile of no
heterogeneity at all.

**At k=6 it is too lax, and this one bites.** KRT8's I² = 0.715 sits *below* the
ceiling — reported as homogeneous — while its calibrated **p = 0.018**, because
at six studies of n≥6 the null median is 0.000. `ceiling_and_null_agree` is
`False` on exactly that row. The one number does opposite things at the two
`k`s, and nothing in `meta.py` knew, because `MAX_I_SQUARED` is Higgins' rule of
thumb for trials with hundreds of subjects each.

**Cochran's Q is worse and is quoted in the committed detail string.** "Q =
80.92 on 10 df" appears in `results/2026-09-05_61ba221/`'s verdict text. At
these patient counts Q rejects **32.5%** of the time against a nominal 5%.

### The closed form, on the real studies

| study | n | weight share | patient share | share of Q | E[ŵ]/w |
|---|---|---|---|---|---|
| Khaliq_2022 | **3** | **20.4%** | 2.5% | **57.1%** | **no finite mean** |
| Chen_2024 | 20 | 33.3% | 16.9% | 16.7% | 1.12× |
| Pelka_2021 | **29** | **10.0%** | 24.6% | 5.5% | 1.08× |
| Joanito_2022 | 4 | 1.6% | 3.4% | 4.8% | 3.00× |
| Lee_2020 | 15 | 3.3% | 12.7% | 0.0% | 1.17× |

**A three-patient study carries twice the weight of a twenty-nine-patient one**
and contributes 57% of the heterogeneity that terminated the reading. Its
per-patient SD on KRT8 is 0.159, the smallest of the eleven — its three patients
came in at −1.235, −1.300 and −0.997. That is not an error in that study; it is
three draws agreeing, which at 2 degrees of freedom happens, and inverse-variance
weighting squares the consequence.

**No study was dropped by name, at any point.** The floors are the chi-square df
conditions and they were fixed in §4 before the curve was run.

### What this changes, and what it does not

**Changes.** `docs/NEXT_AVENUES.md` listed "more carcinoma single-cell data" as
explicitly not worth doing because *"the 13-study result closed that."* The
closure was an UNRESOLVED-by-heterogeneity verdict, and that verdict holds only
at a floor where the weights are not integrable. **The carcinoma question is
closed by a different and stronger route** — §2 of `docs/HANDOFF.md`, five
independent ones — but not by this one, and NEXT_AVENUES is corrected.

**Does not change.** ACTB HOLDS at every floor; the premise's own control is not
disturbed. The **adenoma** results (avenue A, §6h, §6j) never touch `meta.py`
and are untouched here. And at no floor does KRT8 hold: the premise is not
rescued, it is refused for a better-stated reason.

**One observation flagged post-hoc, because §5 did not pre-specify it.** ACTB's
pooled interval includes zero at the committed floor (+0.152 [−0.013, +0.317])
and excludes it at both others (+0.201, +0.212). So the two controls move in
**opposite** directions — ACTB up ~0.21, KRT8 down ~0.48 — which is not a global
detection-scale drift, since drift moves them together. §2 already recorded that
the across-study ACTB/KRT8 correlation is `r = −0.176` without Khaliq. Both
readings say the same thing and neither was pre-registered; they are hypotheses
for a future design, not results.
