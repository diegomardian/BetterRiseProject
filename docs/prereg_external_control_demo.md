# Pre-registration — does an observed-data reference ever cause HARM? A plasmode external control arm

**Written:** 2026-09-15 · **Author:** W2 (harness) · **Status:** proposed ·
**Branch:** `wmhs/external-control-demo` · **Answers** the reviewer objection
recorded in §1 against `paper/wmhs/sections/full/blind.tex` and
`refdesign.tex`.

> **This document is committed before any result table exists, and before the
> module that produces one exists.** The git history is the evidence: this
> commit precedes both `src/harness/external_control_demo.py` and the first
> `results/<date>_<sha>/external_control_*.parquet`. Every number in §9 is an
> *algebraic expectation* computed from the design constants in §4 — the
> weighted averages written out by hand — not an observed outcome. They are
> written down so that a result agreeing with them cannot be mistaken for a
> result chosen to agree with them, and so that a result *disagreeing* with
> them is a falsification rather than a bug to be quietly fixed.

---

## 1 · The objection, in the reviewer's words

> *"The equality audit only protects against something when a study uses an
> observed-data reference in place of theta as ground truth. The paper never
> shows that happening, in its own pipeline or in the seven audited
> repositories."*

The objection is correct about the evidence. `docs/prevalence_audit_result.md`
found **0 of 7**; six compute the reference from generating parameters, which is
the practice the paper recommends. So the paper proposes a safeguard against a
mistake it has never demonstrated anyone making, and never demonstrated causing
harm. `blind.tex` is already careful about this — *"We therefore claim the
failure is possible and cheap to exclude, not that it is widespread"* — but
"possible" is doing a lot of work when the demonstrated consequence is zero.

**What is missing is a harm case.** Not a proof that a check is uninformative
(§blind has that) but a workflow where a competent analyst, following ordinary
practice, is *led to a wrong conclusion they would act on* by validating against
an observed-data reference, and where the audit catches it.

### 1.1 The reviewer's load-bearing assumption

The objection assumes `theta` is always available, so that reaching for an
observed-data reference is a gratuitous choice. **In plasmode and
resampling-based designs it is not available**, and those are the standard
settings for synthetic and external control arms:

- the cohort is built by resampling real patient records, so there is no
  parametric truth — only an effect **induced** on real covariates;
- when the induced effect is modified by covariates, the **marginal** estimand
  is a functional of the realised sample (which records were drawn, in what
  mix), not a number in the config file;
- our own audit's `PamelaShaw/PlasmodeSimulation` case computed exactly such an
  induced marginal effect (`psi0.sim`) and was classified **NO** only on
  criteria 3a and 3b — it is computed from the fitted outcome model, on the
  *original* cohort, with treatment forced to 1 and to 0. **A variant that does
  not hold it out would be YES.**

That is the gap this experiment walks into. An observed-data reference is not an
error there; it is often the only thing to hand.

---

## 2 · The workflow, and why a reviewer should recognise it

A single-arm trial of a new agent. No randomised comparator is feasible, so the
comparison is an **external control arm** drawn from a disease registry —
ICH E10 §2.4 territory, and the design behind a growing share of oncology
submissions. Before touching the real trial, the trial statistician validates
the planned estimator by **plasmode simulation on the registry**: resample real
registry records, induce a treatment effect on the resampled trial arm, and
check that the estimator recovers it.

The five steps, each one standard:

1. **Registry pool.** `N_POOL` real control records. Each carries a covariate
   stratum `g` (line of therapy: treatment-naive / one prior line / heavily
   pretreated), its own **real** outcome `Y0_i` (months of progression-free
   survival), and a latent per-patient responsiveness `r_i` fixed to the record.
2. **Enrolment.** `n_trial` records are drawn **without replacement** under
   stratum-specific eligibility weights. New agents enter in refractory
   populations, so the trial arm is enriched for heavily pretreated patients
   relative to the registry. This is the case-mix shift, and it is the defining
   nuisance of external control arms.
3. **Induction.** The trial arm's outcome becomes `Y1_i = Y0_i + delta_g(i) *
   r_i`. The record keeps its real `Y0_i`. The benefit is **effect-modified**:
   large in treatment-naive patients, near zero in heavily pretreated ones. Also
   standard — a new agent's benefit is usually largest where disease biology is
   least eroded.
4. **External controls.** `n_ctrl = CONTROL_RATIO * n_trial` further records are
   drawn uniformly from the remaining pool, keeping their real `Y0_i`. Registries
   are large, so the control arm is several times the trial arm.
5. **Analysis file.** Trial arm (`treated=1`, outcome `Y1`) stacked on external
   controls (`treated=0`, outcome `Y0`). Arm membership depends on stratum, so
   the strata confound — which is exactly why anybody weights.

**The estimand is the ATT**: the average benefit **in the enrolled trial
population**. That is what the label would claim, what the confirmatory trial
would be powered on, and what ICH E10 and every external-control guidance
document says the comparison is for. It is not the effect in the registry.

### 2.1 Why there is no `theta`

The marginal benefit in the enrolled population is
`mean over the trial arm of delta_g(i) * r_i`. It depends on **which records
enrolment happened to draw** and on their **latent responsiveness**. No scalar
in the configuration equals it. The analyst genuinely must construct a reference
from the draw. §5 gives the two ways to do that; they are not equivalent, and
the whole experiment is the gap between them.

---

## 3 · The defect, and why the analyst would not otherwise notice it

**The defect: the estimator targets the pooled analysis population instead of
the enrolled trial population.** Concretely, saturated IPW written with ATE
weights `1/e` and `1/(1-e)` where the estimand requires ATT weights `1` and
`e/(1-e)`. One line. The same defect is reachable by an entirely different
route — g-computation that standardises the stratum contrasts over **all rows in
the analysis file** rather than over the trial arm — and both are in the sweep,
because the point is that this is a way of thinking, not a typo.

Because the external control arm is `CONTROL_RATIO` times the trial arm, the
pooled case mix is dominated by the **registry**. So the estimator silently
answers *"what would this agent do in the registry population?"* and the analyst
reads the answer as *"what does this agent do in the patients we enrolled?"*.
Under effect modification those are different numbers, and here the first is
nearly twice the second.

**Four reasons a competent analyst does not catch it:**

1. **The number looks right.** It is tightly estimated, it has the right sign,
   it is the right order of magnitude, and it does not move with cohort size
   except to tighten.
2. **Positivity is fine.** Stratum propensities are 0.083 / 0.217 / 0.571
   (§9). No stratum is empty in either arm, no weight is extreme, no trimming
   rule fires. The diagnostic that exists to catch external-control trouble
   stays quiet.
3. **The covariate-balance table passes.** This is the sharp one. Standard
   practice is to report standardised mean differences between arms *after*
   weighting. ATE weights balance both arms to the **pooled** mix, so the
   post-weighting SMDs are ~0 and the balance table is clean. The diagnostic
   checks balance *between arms*, not balance *against the enrolled
   population*, and the defect lives entirely in the latter. **We will compute
   this table and report it** (§7, table 3): if it is not clean, that is
   reportable evidence against our own realism claim.
4. **And the validation confirms it** — which is §5, and the point of the
   experiment.

### 3.1 Why the wrong conclusion matters

Pre-specified clinically meaningful difference: **`MCID = 1.5` months**. The
design constants put the pooled-population answer above it and the enrolled-
population truth below it (§9). The analyst concludes the agent clears the bar
in the population studied, and would act on that: file, or power a confirmatory
trial, on an effect roughly twice what the enrolled patients would get.

---

## 4 · The design constants

Fixed here, before the module exists. Any change is an amendment below.

| constant | value | why |
|---|---|---|
| `N_POOL` | 50,000 | registry-scale, and large enough that drawing arms without replacement does not deplete a stratum |
| strata | 3 | treatment-naive / one prior line / heavily pretreated |
| `POOL_MIX` | (0.55, 0.30, 0.15) | registries are dominated by earlier-line patients |
| `TRIAL_MIX` | (0.15, 0.25, 0.60) | a new agent enters in refractory disease |
| `POOL_MEANS` | (12.0, 8.0, 4.0) months | PFS falls with line of therapy |
| `POOL_SD` | 3.0 | within-stratum spread |
| `DELTA` | (3.0, 1.5, 0.25) months | benefit concentrated where biology is intact |
| responsiveness `r_i` | Gamma(shape 4, scale 0.25); mean 1, sd 0.5 | fixed to the record; makes the marginal truth a property of *which records were drawn*, so no parameter equals it |
| `CONTROL_RATIO` | 3.0 | external control arms are larger than the trial |
| `MCID` | 1.5 months | the decision threshold |
| `n_trial` sweep | 50, 100, 200, 400, 800, 1600 | external control arms are small |
| replicates | 200 per cell | matches `trial_recovery` |
| `POOL_SEED` | 20260915 | the pool is built once and is the same object in every replicate — the plasmode mechanic |
| equality tolerance | `1e-10` in months | explicitly chosen in outcome units, as §blind requires |

The outcome is continuous and uncensored. That is a **stated simplification**:
it removes non-collapsibility and censoring-induced departure (already covered
in §blind's survival arm) so that the case-mix mechanism is isolated. It is a
limitation, recorded in §10.

**The pool is synthetic.** No real registry is available to this repository
(`data/` is empty; see `data/README.md`). What is preserved is the plasmode
*mechanic* — a fixed finite pool, resampled; each record keeps its own outcome;
only the effect is induced — not the provenance of the numbers. This is a
limitation and is recorded in §10, not argued away.

---

## 5 · The references. This is the experiment.

Five references, computed on the same draw. Two are the ones an analyst plausibly
writes; three exist to attribute *which property of a reference does the work*.

| id | definition | what it is |
|---|---|---|
| `obs-pooled` | `sum_g (n_g/n) * (Ybar_1g - Ybar_0g)` over the whole analysis file | **THE MISTAKE.** Observed arm contrast, standardised to the pooled analysis population. Exactly §blind's `T(D)`. |
| `obs-trial` | the same contrast standardised to the **trial arm's** realised mix | observed-data, right population |
| `po-trial` | `mean over the trial arm of (Y1_i - Y0_i)` | **THE HONEST REFERENCE.** The realised individual causal effects among the enrolled, using both potential outcomes. Available in plasmode because we constructed `Y1` ourselves; not reproducible by any estimator, which sees one arm per record. This is `PlasmodeSimulation`'s criterion 3b, kept. |
| `po-pooled` | `mean over ALL analysis records of (Y1_i - Y0_i)` | **CONTROL.** Breaks the equality but keeps the wrong population. |
| `induced-trial` | `sum_g (n_g^trial/n^trial) * delta_g`, ignoring `r_i` | secondary anchor: the parameter-only construction, which is *not* the realised truth |

`obs-pooled` is not a strawman definition. It is "the average treatment effect
in my cohort", computed the way the phrase reads, with `df.groupby('stratum')`
and `len(df)` — the same functional `trial_recovery.py` ships as
`_standardised_effect` and the same one §blind analyses.

**The mechanism under test.** The analyst's target population is wrong in *two*
places at once — in the estimator's weights, and in the reference's weights —
because it is one belief expressed twice. The reference therefore moves with the
defect, and the validation confirms the error instead of catching it. The
audit's job is to notice that the estimator and the reference are the same
functional, and so that the validation could not have failed.

---

## 6 · The estimators

All seven are ordinary external-control practice.

| id | targets | status |
|---|---|---|
| `att-standardisation` | enrolled population | **correct** |
| `att-ipw` | enrolled population | **correct**; saturated ATT weights, algebraically equal to the above |
| `att-matching-all` | enrolled population | **correct**; every stratum-matched control, the most common external-control method |
| `ate-ipw` | pooled population | **THE DEFECT** |
| `ate-standardisation` | pooled population | **THE DEFECT**, reached by the other route |
| `ols-stratum-dummies` | variance-weighted population | a third population again; included because it is what most people would actually type |
| `unadjusted` | neither | the arm the curve can catch |

A draw returns `NaN` for an estimator whenever a stratum is missing an arm, or a
saturated propensity hits 0 or 1. Non-finite outputs are excluded pairwise and
**counted in the table**, per §blind.

### 6.1 The secondary arm that attacks our own equality claim

Saturated-discrete design is what makes the equality exact. So a secondary sweep
re-runs `ate-ipw` with the propensity fitted by **unpenalised logistic
regression on the stratum dummies** — degenerate in intent, and in arithmetic
only up to the solver's convergence tolerance. If the equality evaporates under
the implementation an analyst would actually reach for, that weakens the harm
case and must be reported (F3c).

---

## 7 · The primary outcome

Three tables under `results/<date>_<sha>/`.

**Table 1 — `external_control_primary`.** One row per
`(estimator, reference, n_trial)`:

- `max_residual` — `max |theta_hat - T(D)|` over valid pairs. The maximum, not
  the median (§blind).
- `equality_flagged` — `max_residual <= 1e-10`. The audit's verdict.
- `bias`, `rmse` — against that reference.
- `ratio_median`, `ratio_q25`, `ratio_q75` — the recovery curve the analyst
  looks at.
- `n_valid`, `n_excluded`.

**Table 2 — `external_control_shift`.** The case-mix sweep for F3d.
`TRIAL_MIX` is interpolated from `POOL_MIX` (`lambda = 0`, no shift) to the §4
value (`lambda = 1`), at `lambda in {0, 0.25, 0.5, 0.75, 1.0}`, at the largest
cohort size. Carries the ATE/ATT gap, the equality residual, and whether the
`MCID` decision flips. **`lambda = 0` is the experiment's own null**: with no
case-mix shift the defect must vanish and the two references must agree.

**Table 3 — `external_control_balance`.** Post-weighting standardised mean
differences between arms under ATT and ATE weights, plus the SMD of each
weighted pseudo-population against the **enrolled** case mix. Tests §3's
claim (3) that the diagnostic an analyst actually runs stays quiet.

### 7.1 The primary claim, as four propositions

- **P1 — the validation passes, and passes *best* for the defective arm.**
  Against `obs-pooled`: `ate-ipw` and `ate-standardisation` have
  `max_residual <= 1e-10` at every `n_trial`, `ratio_median` within
  `[0.99, 1.01]`, and the **smallest RMSE of all seven estimators** at every
  `n_trial`. The correct estimators show an apparent error that **does not
  shrink with `n`**. The analyst does not merely fail to catch the defect — the
  validation *selects the defective estimator over the correct one*.
- **P2 — the honest reference inverts the ranking.** Against `po-trial`:
  `att-*` bias → 0 and RMSE falls with `n`; `ate-*` carries a bias that does not
  shrink.
- **P3 — the wrong conclusion is actionable.** The defective estimator's median
  estimate exceeds `MCID = 1.5` at every `n_trial`; `po-trial` is below it at
  every `n_trial`.
- **P4 — breaking the equality is NOT sufficient.** Against `po-pooled`, the
  equality is broken (`max_residual > 1e-10`) and the defect is **not**
  revealed: `ate-*` is approximately unbiased there. The reference must also
  name the right target population. Without P4 the demonstration would license
  "any non-reproduced reference will do", which is false.

---

## 8 · THE FALSIFIERS

Fixed before the run. **Not to be weakened after seeing results.** Any deviation
is an amendment appended below with a UTC timestamp, leaving the original text.

**F1 — the validation does not pass, so there is nothing to be misled by.**
*If, against `obs-pooled`, the defective arm fails any of: (a)
`max_residual <= 1e-10` at every `n_trial`; (b) `ratio_median` in
`[0.99, 1.01]` at every `n_trial`; (c) strictly smallest RMSE among all seven
estimators at every `n_trial` — then the analyst had a visible warning and this
is not a harm case.*

**F2 — the honest reference buys nothing.** *If, against `po-trial`, the
defective arm's `|bias|` is below `0.2` months (about 20% of the ATT) at
`n_trial >= 400`, OR if `|bias|` at `n_trial = 1600` is below half its value at
`n_trial = 100` (i.e. it shrinks like an ordinary sampling error), then the
honest reference does not reveal the defect and the audit buys nothing.*

**F3 — the workflow is contrived.** Five separate charges, each with its own
pre-specified check. **F3 fires if any of a–e fires.**

- **F3a — nobody writes those weights.** *Fires if the defect is only reachable
  through `ate-ipw`. Checked by requiring `ate-standardisation` — plain
  "average the stratum contrasts over my analysis file" — to exhibit the same
  equality and the same bias. If it does not, the defect is a weight-formula
  curiosity rather than a way of thinking.*
- **F3b — nobody defines truth that way.** *Fires if `obs-pooled` is not
  bit-identical, as a functional, to the reference the paper itself analyses.
  Checked by asserting in the test suite that `obs-pooled` on a trial-shaped
  frame equals `trial_recovery._standardised_effect` on the same frame.*
- **F3c — the equality is an artefact of the saturated design.** *Fires if the
  secondary logistic-propensity arm (§6.1) does not flag equality at the
  pre-specified `1e-10`, AND does not flag it at a tolerance an analyst would
  plausibly choose either (`np.allclose` default, `atol = 1e-8`). If the
  equality only exists for hand-written saturated weights, the harm case
  narrows to that implementation and we say so.*
- **F3d — the case-mix shift is implausibly large.** *Fires if the `MCID`
  decision flips only at `lambda = 1.0` — i.e. only at the most extreme shift we
  wrote down. External control arms routinely show large baseline imbalance, so
  a flip that requires the maximum shift and nothing less is a designed
  coincidence, not a realistic one. Reported as the smallest `lambda` at which
  the flip occurs.*
- **F3e — a standard diagnostic already catches it.** *Fires if Table 3 shows
  the post-ATE-weighting between-arm SMD exceeding `0.1` in any stratum — the
  conventional balance threshold. If the routine balance table would have caught
  this, the analyst had an independent warning and the audit is redundant.*

**F4 — the demonstration cannot fail.** *Fires if, at `lambda = 0` (no case-mix
shift), the defect does NOT vanish: specifically if `|ATE - ATT| > 0.05` months
or if the `ate-*` bias against `po-trial` exceeds `0.1` months there. The
experiment must have a null in which the alleged harm is absent, or it is not
measuring case-mix shift and is measuring something else.*

**F5 — the honest reference is secretly the parametric one.** *Fires if
`po-trial` and `induced-trial` agree to within `1e-10` at every `n_trial`. If
they do, the latent responsiveness `r_i` is not doing its job, the realised
truth is a function of parameters plus realised counts alone, and the reviewer's
assumption that `theta` is available survives largely intact. This would not
falsify P1–P3 but it substantially weakens §1.1, and must be reported in the
first paragraph of the result.*

---

## 9 · What we expect, written down before the run

Computed by hand from §4, from the weighted averages alone. **No code has been
run.**

Pooled case mix at `CONTROL_RATIO = 3`:
`(TRIAL_MIX + 3 * POOL_MIX) / 4 = (0.4500, 0.2875, 0.2625)`.

| quantity | expected | how |
|---|---|---|
| ATT (enrolled population) | **0.975** months | `0.15*3.0 + 0.25*1.5 + 0.60*0.25` |
| ATE (pooled population) | **1.847** months | `0.4500*3.0 + 0.2875*1.5 + 0.2625*0.25` |
| overstatement | **×1.89** | |
| `MCID` | 1.5 | ATE above, ATT below → decision flips |
| `e(g) = P(treated \| g)` | 0.083 / 0.217 / 0.571 | no positivity failure, no trimming |
| `max_residual`, `ate-*` vs `obs-pooled` | `0.0` or ~`1e-15` | same functional; `ate-standardisation` bitwise, `ate-ipw` to rounding |
| `max_residual`, `att-*` vs `obs-trial` | ~`1e-15` | the same identity at the right population |
| apparent bias, `att-*` vs `obs-pooled` | ≈ **−0.87** months, flat in `n` | `ATT − ATE`; the correct estimator is *rejected* |
| bias, `ate-*` vs `po-trial` | ≈ **+0.87** months, flat in `n` | the defect, revealed |
| bias, `ate-*` vs `po-pooled` | ≈ **0**, shrinking | P4: equality broken, defect hidden |

If the run contradicts these, the prediction was wrong and that is the result.

---

## 10 · Limitations, recorded in advance

1. **The registry pool is synthetic.** The plasmode mechanic is preserved; the
   provenance of the records is not real. Any claim of the form "this happens in
   real registry data" is out of scope.
2. **Continuous, uncensored outcome.** Non-collapsibility and censoring are
   deliberately excluded (§4).
3. **Discrete strata.** The equality is exact because the design is saturated.
   §6.1 is the probe; whatever it returns is the honest extent of the claim.
4. **This is a constructed workflow, not an observed one.** It demonstrates that
   the harm is *reachable by ordinary practice*. It does not raise the audit's
   measured prevalence above 0 of 7, and the paper must not read it as doing so.
5. **The honest reference is itself constructed.** `po-trial` requires an
   explicit choice of target population. §11 is about what that means for the
   paper's claim.

---

## 11 · What this refines in the paper, if it works

The reviewer's assumption is that `theta` is available. If P1–P4 hold, the
finding is sharper than "sometimes it isn't": it is that in this class of design
the analyst must **assemble** a marginal truth, and assembling it requires
naming a target population — the very decision the defect corrupts. The
safeguard is then not "use `theta` instead"; it is "**construct the reference
from information the estimator cannot consume, and state the population it is
constructed for**". P4 is what forces the second clause. If P4 fails, the second
clause is unsupported and must not be written.

---

## 12 · Deliverables

- `src/harness/external_control_demo.py`, CLI matching `trial_recovery.main`.
- `results/<date>_<sha>/external_control_{primary,shift,balance}.parquet` with
  `.meta.json` sidecars.
- `tests/test_external_control_demo.py`, containing for every guard the input
  that **forces it to fail**, per `tests/test_checks_can_fail.py`.
- A `## RESULT` section appended below, stating which falsifiers fired.

---

## RESULT

**Run:** `results/2026-09-15_6b4a31a/` · seed `20260915` · 200 replicates per
cell · 6 cohort sizes · 8 estimators · 5 references · 5 case-mix shifts ·
50 balance replicates. Tables:
`external_control_{primary,shift,balance,falsifiers}.parquet`.

**No falsifier fired.** F1, F2, F3a, F3b, F3c, F3d, F3e, F4, F5 all clear, each
with the number that decided it recorded in
`external_control_falsifiers.parquet`. P1–P4 all hold. The demonstration exists.

Reading F5 first, as §8 requires: `max |median(po-trial) − median(induced-trial)|
= 0.0098` months. The realised truth is genuinely not the parameter-plus-counts
construction, so F5 does not fire — **but the gap is about 1% of the ATT**, and
§13 below says what that costs §1.1.

### 1 · P1 — the validation passes, and it passes *best* for the defective arm

Against `obs-pooled`, the reference the analyst wrote:

| estimator | max residual | RMSE | ratio median (IQR) | median estimate |
|---|---|---|---|---|
| `ate-standardisation` | **0.0 bitwise**, every n | 0.0 | 1.000000 (0.000000) | 1.85–1.93 |
| `ate-ipw` | 4.7e−15 … 6.0e−15 | 1.6e−15 | 1.000000 (0.000000) | 1.85–1.93 |
| `ate-ipw-logistic` | 6.1e−14 … 5.4e−13 | 1.7e−13 | 1.000000 (0.000000) | 1.85–1.93 |
| `att-standardisation` (**correct**) | 1.11 … 2.38 | **0.876 … 1.015, not shrinking** | 0.512–0.552 | 0.94–1.03 |
| `ols-stratum-dummies` | 0.69 … 1.87 | 0.529 … 0.669 | 0.694–0.721 | 1.28–1.38 |
| `unadjusted` | 4.72 … 7.06 | 4.31 … 4.37 | −1.25 to −1.37 | −2.42 |

The defective arm has the smallest RMSE of all seven primary estimators **at
every cohort size**. The correct estimator looks 48% attenuated and looks not to
improve with n. An analyst running this bake-off does not merely fail to catch
the defect: **the validation instructs them to discard the correct estimator and
deploy the broken one.**

`ate-ipw-logistic` matters here. Hand-counted saturated weights give bitwise
zero, which a suspicious analyst might notice. The fitted-propensity
implementation — what anyone actually writes — gives a residual of 1.7e−13
months, which does not read as "something is wrong". It reads as "my estimator
recovers the truth to numerical precision".

### 2 · P3 — the conclusion is wrong in a way that would be acted on

Median defective estimate **1.82–1.93 months** at every cohort size, against a
pre-specified MCID of **1.5**. Median `po-trial` — the benefit the enrolled
patients actually received — is **0.97–1.00 months** at every cohort size.
Overstatement **×1.89**. The agent clears the bar on the analysis and misses it
in the population studied. This is a filing decision, or the effect size a
confirmatory trial is powered on.

### 3 · P2 — the counterfactual: the honest reference reveals it

Same replicates, same estimators, reference swapped to `po-trial`:

| estimator | bias at n=50 → n=1600 | RMSE at n=50 → n=1600 | ratio median |
|---|---|---|---|
| `ate-ipw` | +0.844 → **+0.885** (flat) | 1.082 → **0.892** (floors) | **1.88**, not 1 |
| `att-standardisation` | −0.015 → +0.013 | 0.608 → **0.105** (falls as n^−1/2) | 0.98 → 1.02 |

**The ranking inverts completely.** The recovery curve that was pinned at exactly
1 now sits at 1.88 and stays there; the estimator that looked 48% attenuated is
unbiased and converges. The same check, run against a reference the estimator
cannot reproduce, gives the opposite and correct answer.

### 4 · P4 — and breaking the equality is NOT enough

The control arm is the reason this is not just "use any other reference".
Against `po-pooled` — potential outcomes, so no estimator reproduces it, but the
*pooled* population:

| estimator | max residual | bias at n=1600 |
|---|---|---|
| `ate-ipw` | 0.36 … 2.04 (**equality broken**) | **+0.003** — defect invisible |
| `att-standardisation` (**correct**) | 1.18 … 2.68 | **−0.869** — looks broken |

Breaking the equality bought nothing. It reproduces the original wrong answer
with a different reference. **The reference must name the right target
population**; not being reproduced is necessary and not sufficient. §11's second
clause is supported.

### 5 · The equality screen is a property of the PAIR, not a verdict on the estimator

`equality_flagged`, all six cohort sizes:

| estimator | `obs-pooled` | `obs-trial` | `po-trial` | `po-pooled` |
|---|---|---|---|---|
| `ate-*` (defective) | **True** | False | False | False |
| `att-*` (correct) | False | **True** | False | False |
| `ols`, `unadjusted` | False | False | False | False |

The three correct estimators are exactly as "blind" against `obs-trial` as the
defective ones are against `obs-pooled`. **The screen does not identify the
defective estimator.** What it identifies is that a particular validation could
not have failed. That is the honest claim, it is the one `blind.tex` already
makes, and this table is the sharpest version of it we have.

### 6 · F4 — the demonstration's own null, which it passes

At `shift = 0` the trial enrols the registry's own case mix. `ate-ipw` bias
against `po-trial` is **0.0017 months** — the harm is gone. And the equality is
still there: max residual against `obs-pooled` is 4.4e−15, still flagged.
**Equality is not itself harm.** The screen flags a validation that could not
fail, which is true at `shift = 0` as well; whether that costs anything depends
on the design. This is a real qualification on how the paper should sell the
audit, and it came out of the pre-registered null rather than out of review.

### 7 · F3d — the decision flips well before the extreme

| shift | trial mix | `ate-ipw` median | `po-trial` | bias | decision flips |
|---|---|---|---|---|---|
| 0.00 | (0.55, 0.30, 0.15) | 2.14 | 2.13 | 0.002 | no |
| 0.25 | (0.45, 0.29, 0.26) | 2.07 | 1.85 | 0.23 | no (both clear MCID) |
| 0.50 | (0.35, 0.28, 0.38) | 2.00 | 1.57 | 0.43 | no (both clear MCID) |
| 0.75 | (0.25, 0.26, 0.49) | 1.95 | 1.29 | 0.67 | **yes** |
| 1.00 | (0.15, 0.25, 0.60) | 1.89 | 1.00 | 0.88 | **yes** |

Smallest flip at `shift = 0.75`. The bias is monotone and substantial from
`shift = 0.25`. F3d does not fire.

### 8 · F3e — the diagnostic a protocol actually runs stays silent

Post-weighting standardised mean differences, 50 replicates at n=1600:

| weighting | worst \|SMD\| **between arms** | worst \|SMD\| **vs the enrolled population** |
|---|---|---|
| ATE (defective) | **5.6e−16** — passes perfectly | **0.964** |
| ATT (correct) | 4.5e−16 | 2.3e−16 |
| unweighted | 1.13 — fails, as expected | 0.964 |

The balance table an external-control protocol reports is **immaculate** under
the defective weights, because ATE weighting balances the two arms against each
other — against the pooled mix. The pseudo-population is a full standardised
difference away from the patients enrolled, and that column is not in anyone's
Table 1. Positivity is also clean throughout (stratum propensities 0.083 /
0.217 / 0.571; zero exclusions in all 48,000 primary draws).

### 9 · §9's predictions against the run

| quantity | predicted | observed |
|---|---|---|
| ATT | 0.975 | 0.97–1.00 |
| pooled ATE | 1.847 | 1.85–1.93 |
| overstatement | ×1.89 | ×1.89 |
| `ate-standardisation` residual vs `obs-pooled` | 0 or ~1e−15 | **0.0 bitwise** |
| `att-*` apparent bias vs `obs-pooled` | ≈ −0.87, flat | −0.876 … −1.015, flat |
| `ate-*` bias vs `po-trial` | ≈ +0.87, flat | +0.844 … +0.921, flat |
| `ate-*` bias vs `po-pooled` | ≈ 0 | +0.003 at n=1600 |

Nothing was back-filled and nothing missed.

### 10 · Our own verdict on F3 (contrivedness)

No F3 sub-charge fired. Our honest reading, stated more harshly than the
falsifiers require:

**What survives scrutiny.** The workflow is one a reviewer will recognise: an
external control arm from a registry, validated by plasmode resampling, with
case-mix shift and effect modification. Every component is ordinary. The defect
is reachable by two independent routes — an IPW weight slip and plain "average
over my analysis file" — so it is a way of thinking, not a typo (F3a). The
mistaken reference is bit-identical to the functional the repository already
shipped and the paper already analyses, asserted in the test suite (F3b). The
equality survives a fitted propensity (F3c). The decision flips at three
quarters of the shift, not only at the extreme (F3d). Every routine safeguard —
positivity, weight extremity, post-weighting balance — passes cleanly (F3e).

**Three places a determined reviewer can still push.**

1. **The winner's recovery ratio is exactly 1.000000 with zero IQR**, not "near
   1 and tightening". A curve with no dispersion at any n is itself odd, and an
   alert analyst might ask why. Our answer is `ate-ipw-logistic`: the
   implementation anyone actually writes gives 1.7e−13 months, which reads as a
   triumph rather than an alarm. But this is a real caveat and the paper should
   not claim the curve is indistinguishable from a healthy one — it is
   distinguishable, by a reader who already knows to look, which is the reader
   the audit is for.
2. **The analyst must get the target population wrong twice** — once in the
   estimator, once in the reference. That is the mechanism, not a coincidence:
   it is one belief about "the population" expressed in two places, and neither
   place writes it down as a parameter; it is the implicit denominator, `len(df)`
   in one and the weight formula in the other. `trial_recovery.py`'s own
   docstring already describes this as "the natural first thing to build when
   the simulator and the estimator are written by the same person in the same
   afternoon". We think this is realistic. We accept that it is the load-bearing
   assumption of the whole demonstration, and that a reviewer who rejects it
   rejects the result.
3. **The registry pool is synthetic** (§10.1). The plasmode mechanic is
   preserved; the provenance is not real.

We do **not** think the scenario is contrived. We do think its realism rests
entirely on point 2, and the paper should say so in that form rather than
claiming more.

### 11 · What this does and does not establish

It establishes that the harm is **reachable by ordinary practice**, with a
measured consequence: a ×1.89 overstatement that crosses a decision threshold,
produced by a validation that ranks the broken estimator first.

It does **not** raise the audit's measured prevalence above **0 of 7**. This is a
constructed workflow, not an observed one. The paper must not cite it as
evidence that anybody has done this.

### 12 · The correction this forces on the paper's own remedy

Two, both from pre-registered controls rather than from review:

- **P4.** "Use a reference the estimator does not reproduce" is insufficient.
  `po-pooled` is not reproduced and is just as misleading. The remedy is
  *"construct the reference from information the estimator cannot consume, **and
  state the target population it is constructed for**"*.
- **F4.** At zero case-mix shift the equality persists and the harm does not.
  Equality flags a validation that could not have failed; it does not flag a
  wrong answer. The audit is a statement about evidence, never about validity —
  which is what `blind.tex` already says, now with a case where the distinction
  has a price tag and a case where it does not.

### 13 · What it costs §1.1

F5 did not fire, but only by 0.0098 months — about 1% of the ATT. The latent
responsiveness makes the realised truth formally non-parametric and practically
almost parametric. So the strongest form of §1.1 — *"`theta` does not exist
here"* — is **not** what the run supports. The supported form is weaker and more
useful:

> The marginal estimand cannot be read off the configuration. It must be
> **assembled**, and assembling it requires naming a target population — which is
> exactly the decision the defect corrupts. An analyst who assembles it from the
> generating parameters *and states the population* catches this. An analyst who
> reads "the average effect in my cohort" off the analysis file does not.

That is a narrower claim than the pre-registration hoped for, and it is the one
the paper should make.
