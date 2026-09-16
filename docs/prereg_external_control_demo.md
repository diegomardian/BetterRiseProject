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
