# Pre-registration — the residual audit under a LEARNED generative patient model

**Written:** 2026-09-13 · **Author:** W2 (harness) · **Status:** proposed ·
**Branch:** `wmhs/learned-generator` · **Tests** the prediction made and
explicitly not run in `paper/wmhs/sections/full/refdesign.tex`, final paragraph.

> **This document is committed before any result table exists.** The git history
> is the evidence: this commit precedes the first
> `results/<date>_<sha>/trial_learned_generator_*.parquet`. Nothing in §7 has
> been computed on anything at the time of writing. What §8 records is an
> *algebraic expectation* derived from the estimators' functional forms before
> the run, not an observed outcome — it is written down precisely so that a
> reader can see it was not back-filled, and so that a result agreeing with it
> is not mistaken for a result chosen to agree with it.

---

## 1 · The claim under test, in the authors' own words

`refdesign.tex` closes with a prediction the paper flags as unrun:

> *"First, the equality is likely to become **approximate**: a learned generator
> and a re-fit outcome model share most but not all of their inputs, so the
> residual would sit near zero rather than at it — the regime §blind shows the
> boolean screen cannot see, and the one the residual scale ρ exists for. A
> practitioner who runs only `allclose` would find nothing. Second, step 1 of
> the audit inverts in difficulty: pre-specifying T(D) is trivial for a
> stratified simulator and is the substantive modelling decision for a
> generative one."*

Decomposed into three testable pieces:

| | claim | tested here |
|---|---|---|
| **(a)** | with a learned generator the equality degrades from exact to approximate; the residual sits near zero rather than at it | **yes, primary** |
| **(b)** | therefore `np.allclose` at default tolerance returns True and sees nothing, while ρ remains non-zero and does see it | **yes, primary** |
| **(c)** | pre-specifying `T(D)` inverts from trivial to being the substantive modelling decision | **partially** — §3 carries two admissible references and §7.4 reports what each costs; this is evidence for (c), not a test of it |

**(a) and (b) may be false. This document fixes in advance what would show
that, and §6 is not to be weakened after the numbers are seen.**

There is a fourth thing being tested, which the paper does not state as a claim
but does state as a *mechanism*, and which §6's F3 is aimed at:

> *"a learned generator and a re-fit outcome model share most but not all of
> their inputs"*

That sentence asserts the residual is a property of the **generator's
learnedness**. §8 records why we expect it to be a property of the
**estimator's functional form** instead, and F3 is the falsifier that
distinguishes them.

---

## 2 · The generator class, and its hyperparameters

### 2.1 What stays fixed

The ground truth remains `src/harness/trial_recovery.simulate_trial`,
**imported unchanged**: three covariate strata with means
`(10.0, 14.0, 25.0)`, prevalences `(0.5, 0.3, 0.2)`, stratum-dependent
treatment probabilities `(0.25, 0.50, 0.75)` so assignment is confounded,
a homogeneous additive effect `theta = 3.0`, and outcome SD `4.0`. A *training
cohort* of `n` patients is drawn from it. Every generator below is fitted to
that training cohort and then asked for a *synthetic cohort* of the same `n`.

This two-stage shape — real draw → fit → synthetic draw → estimate — is the
virtual-control-arm wiring the paper is about. The estimator never sees the
training cohort.

### 2.2 The generators

`torch` is **not pinned in any `env/*.yml`** (checked: `w1_reference`,
`w2_harness`, `w3_bulk`, `w4_estimator` pin `numpy`, `pandas`, `scipy`,
`scikit-learn=1.5.1` and no deep-learning stack). Per the standing rule that
results run on the pinned environment, **no VAE and no torch**. The learned
generators are sklearn, which is enough: a conditional Gaussian mixture over
`(covariates, treatment, outcome)` is a legitimate learned generative model of
tabular patient data and is what a practitioner without a GPU actually reaches
for.

| name | class | learned? | capacity knob |
|---|---|---|---|
| `parametric-oracle` | re-draw from `simulate_trial` at the same `n`, `theta`, fresh stream | **no — control** | — |
| `parametric-plugin` | MLE plug-in: `pi_g` from counts, `p_g` from within-stratum treated share, `mu_gd`/`sigma_gd` from cell mean and SD, Gaussian outcomes | **saturated** | — |
| `gmm-joint-K1` | `GaussianMixture(n_components=1)` on `[onehot(stratum)(3), treated(1), outcome(1)]`, `.sample()`, discretised back | yes | K=1 |
| `gmm-joint-K2` | as above | yes | K=2 |
| `gmm-joint-K4` | as above | yes | K=4 |
| `gmm-joint-K8` | as above | yes | K=8 |
| `mlp-conditional-gaussian` | `MLPRegressor` for `E[Y|g,a]`, Gaussian noise at the residual SD; `(stratum, treated)` resampled from the empirical joint | yes | `hidden_layer_sizes=(16,)` |

**Fixed hyperparameters**, chosen before the run and not tuned against any
residual:

- `GaussianMixture(covariance_type="full", reg_covar=1e-6, tol=1e-3,
  max_iter=200, n_init=1, init_params="kmeans", random_state=<derived>)` —
  sklearn defaults except `random_state`, which is derived from the seed
  stream. **`tol=1e-3` is sklearn's default and is load-bearing for §8**; it is
  left at the default deliberately, because the question is what a practitioner
  running default sklearn gets.
- `MLPRegressor(hidden_layer_sizes=(16,), solver="lbfgs", max_iter=300,
  random_state=<derived>)`. `lbfgs` over the default `adam` on runtime grounds
  alone (0.60s vs 2.46s per fit at n=5000, measured before the run).
- Discretisation for `gmm-joint-*`: `stratum = argmax` over the three one-hot
  columns, `treated = 1` if that column `> 0.5`. Sampled rows are otherwise
  taken as drawn. This is the standard "fit a GMM to a table and sample rows"
  recipe and its crudeness is part of what is being measured, not a defect to
  be hidden.

### 2.3 The saturation check — reported, never a filter

A learned generator that memorises the training cohort's cell means **is** the
parametric generator wearing different vocabulary, and its residual would be
zero for an uninteresting reason. So every generator reports, as a column on
its own row:

```
max_g,d | mu_hat_gd(generator) - ybar_gd(training cohort) |
```

`generator_is_saturated = (that maximum < 1e-9)`.

This is a **reported row, not an exclusion**. A generator that turns out to be
the parametric one in disguise is more convincing sitting in the table labelled
as such than quietly dropped from it. `parametric-plugin` is expected to be
saturated **by construction** and is included precisely to calibrate what
saturation looks like in this column.

---

## 3 · `T(D)` — the two admissible references, and why

The paper (§blind) defines the realised reference as an **observed-data**
functional of the draw: `T(D) = h(S(D))`, with `S` denoting summaries of the
draw, and is explicit that it "is an observed contrast, not the finite-sample
mean of paired potential-outcome differences". Under a learned generator there
are two things "the draw" can mean, and choosing between them is exactly the
inversion claim (c) is about. Both are pre-specified; the primary is fixed now.

### 3.1 `T_draw` — PRIMARY

```
T_draw(D~) = sum_g (n~_g / n~) * (ybar~_g1 - ybar~_g0)
```

computed on the **synthetic cohort** `D~` by the *same function*
`trial_recovery._standardised_effect` that the parametric arm uses. Imported,
not re-implemented.

**Why this is the right primary.** It is the paper's own §blind definition
applied to the draw the estimator actually sees. The estimator is handed `D~`
and nothing else; an audit screening `theta_hat` against a reference the
estimator had no access to is testing something other than functional reuse.
It is also the reference a practitioner can actually compute, which is the
condition §blind imposes on the whole audit.

### 3.2 `T_model` — SECONDARY, pre-specified

The fitted generator's own implied population standardised effect:

- `parametric-oracle`: `= theta = 3.0` exactly, in closed form.
- `parametric-plugin`: `sum_g pi_hat_g * (mu_hat_g1 - mu_hat_g0)`, closed form.
- `mlp-conditional-gaussian`: `sum_g pi_hat_g * (m_hat(g,1) - m_hat(g,0))` from
  the fitted regressor, closed form.
- `gmm-joint-K*`: **no closed form survives the discretisation step**, so it is
  estimated by an independent Monte-Carlo draw of `n_oracle = 200_000` from the
  same fitted generator on a dedicated seed stream, and **its Monte-Carlo
  standard error is reported as a column beside it**.

**Why it is secondary and why it is here anyway.** `T_model` is the analogue of
the latent reference in §blind's censoring paragraph: it asks a different
question (has the estimator recovered the generator's *population* effect)
and answers it with error of its own. For the GMM arms that error is roughly
`0.02` in outcome units and therefore **cannot resolve residuals below about
`2e-2`** — which is itself the quantitative content of claim (c): under a
learned generator the model-based reference is not free, it is an estimation
problem with its own floor, three to eleven orders of magnitude above the
residuals the observed-data reference resolves. That number is a deliverable,
not a nuisance.

**No third reference.** The residual against `theta_requested = 3.0` is the
recovery curve, which is the thing being audited, and is reported as
`difference_vs_requested` for context only — never as a residual.

---

## 4 · The estimators

All are re-fit **on the synthetic cohort**, which is the natural wiring the
paper says practitioners will use. Four are imported from `trial_recovery` /
`trial_blindness` unchanged so the numbers are comparable with the published
tables rather than a second implementation of the same idea.

| name | what it is | relation to `T_draw`'s functional |
|---|---|---|
| `gcomp-saturated` | standardisation over realised strata using cell means | **the same functional**, imported `gcomp_from_generator` |
| `ipw-saturated` | IPW, propensity `= phat_g` per stratum | algebraically the same functional |
| `gcomp-gmm-outcome-K1` | outcome model = per-cell `GaussianMixture(1)`, `mu_hat_gd = sum_k pi_k mu_k`; standardised by realised `n~_g/n~` | **same functional up to EM convergence** — see §8 |
| `gcomp-gmm-outcome-K4` | as above, `K=4` | same, with more EM to converge |
| `gcomp-mlp-outcome` | outcome model = `MLPRegressor` on `[stratum dummies, treated, treated x stratum]`; `mu_hat(g,d)` predicted | saturated in principle, **optimisation error in practice** |
| `ols-stratum-dummies` | OLS on treatment + additive stratum dummies | **NOT the functional** — variance-weighted, imported |
| `unadjusted` | difference in arm means | **NOT the functional**, and confounded |

### 4.1 The control that makes the screen mean anything

`ols-stratum-dummies` is the **negative control estimator** and it carries the
experiment's internal validity. It is consistent for `theta` under this
homogeneous effect — it is *accurate* — and it is not the reference's
functional, targeting the variance-weighted contrast instead
(Frisch–Waugh, `trial_blindness.varweighted_effect`).

**If `ols-stratum-dummies` does not show a residual clearly larger than the
near-zero arms, the screen is measuring nothing and no other row in this
experiment may be read.** The published parametric value is `0.0663` at
n=5,000, so "clearly larger" is pre-specified as **`> 1e-2` at every cohort
size**. `unadjusted` (published recovery ratio ≈ 2.4, residual ≈ 4.99) is the
second, coarser control.

This is the distinction between *residual near zero because the functional is
shared* and *residual near zero because the estimator is simply accurate*. A
design without a consistent-but-different-functional arm cannot tell them
apart.

---

## 5 · The sweep

| axis | values |
|---|---|
| cohort size `n` | `100, 200, 500, 1000, 2000, 5000` — matched to `trial_recovery` / `trial_blindness` so the parametric control is directly comparable |
| seed streams | `6` (`--seeds`), base seed `20260913` (`--seed`) |
| replicates per seed | `30` (`--replicates`) → **180 draws per cell** |
| generators | the 7 of §2.2 |
| estimators | the 7 of §4 |
| references | `T_draw` (primary), `T_model` (secondary) |

Grid = 7 × 7 × 6 × 6 × 30 = **52,920 estimates**, each on its own synthetic
cohort. RNG streams are `np.random.default_rng([seed, arm_code, ...keys])` in a
written-down key order, following `trial_blindness.SEED_KEY_ORDER`; the training
draw, the generator fit, the synthetic draw, the `T_model` oracle draw and each
estimator get **separate streams**, so that adding an estimator cannot change
another estimator's numbers.

`--quick` runs a smoke grid (2 sizes, 2 seeds, 3 replicates) which **is not a
result** and is labelled so in its sidecar.

### 5.1 Non-finite outputs

A synthetic cohort can lose a stratum-arm cell — a positivity failure — and
`_standardised_effect` returns `NaN` there, as does every estimator that needs
both arms. These are **counted and reported in their own columns**
(`n_nonfinite_reference`, `n_nonfinite_estimate`, `n_valid_pairs`) and excluded
from the maximum, per §blind: *"report excluded draws and do not compute the
maximum if none remain."*

**A cell with zero valid pairs reports `max_residual = NaN` and
`screen_verdict = "no valid pairs"`, never `True`.** `np.allclose([], [])` is
`True`, and a screen that returns green because everything was filtered out is
the exact failure this repository's test file exists to catch. There is a
forcing test for it.

---

## 6 · THE FALSIFIER

Fixed before the run. **Not to be weakened after seeing results.** Any
deviation is an amendment appended below with a UTC timestamp, stating what
changed and why, leaving the original text in place.

**F1 — nothing degraded.** *If, for the learned generators, the natural wiring
(`gcomp-saturated` — the re-fit outcome model that is the reference's
functional) returns a residual that is **exactly `0.0`**, bitwise, then
prediction (a) is FALSE as stated: the equality did not degrade to approximate
merely because the generator became learned.*

**F2 — the boolean screen works fine.** *If the residual for the learned
generators' natural wiring is large enough that `np.allclose(residual, 0.0)` is
**False** at numpy's default tolerance, then prediction (b) is FALSE: a
practitioner running only `allclose` would find exactly what there is to find,
and the claimed blind regime is not where this lands.*

**F3 — the stated mechanism is false.** *If the residual under `T_draw` is the
**same** for `parametric-oracle` and for the learned generators at a matched
estimator, then the mechanism `refdesign.tex` gives for (a) — "a learned
generator and a re-fit outcome model share most but not all of their inputs" —
is FALSE. The residual would then be a property of the (estimator, reference)
pair, and the generator's learnedness would not enter it at all.*

**F4 — capacity is irrelevant.** *If `max_residual` under `T_draw` is flat in
`K` across `gmm-joint-K1/K2/K4/K8` at a matched estimator, then "vary capacity
to see where exact becomes approximate" does not describe this audit, and the
capacity axis belongs to `T_model` rather than to the screen.*

**F5 — the learned generator is the parametric one.** *If a generator's
`generator_is_saturated` column is `True` (§2.3), its rows are not evidence
about learned generators and must not be quoted as such, however its residual
comes out.*

**F6 — the control estimator collapses.** *If `ols-stratum-dummies` does not
exceed `1e-2` at every cohort size, the screen has no demonstrated ability to
separate shared-functional from accurate, and **no other row in the experiment
may be read**. This falsifies the experiment, not the paper.*

**F7 — ρ says nothing extra.** *Claim (b) has two halves. If ρ is zero, or
non-finite, or fails to order the arms in the same direction the residual does,
then the half of (b) asserting that "the residual scale ρ remains nonzero and
does see it" is unsupported even where the `allclose` half holds.*

---

## 7 · The primary outcome

One table, `trial_learned_generator_primary`, one row per
`(generator, estimator, reference, n_patients)`, carrying:

1. `max_residual` — `max |theta_hat - T(D~)|` over the 180 valid draws. The
   maximum, not the median: §blind is explicit that "a zero median can conceal
   nonzero residuals, which is why the equality screen uses the maximum".
2. `is_exactly_zero` — `max_residual == 0.0`, bitwise. This is F1's column.
3. **`allclose_residual_vs_zero`** — `np.allclose(max_residual, 0.0)`.
   Default `rtol=1e-5`, `atol=1e-8`; since the second argument is `0.0` the
   `rtol` term vanishes and **this reduces to a pure `atol = 1e-8` screen**.
4. **`allclose_estimate_vs_reference`** — `np.allclose(estimate, reference)`
   elementwise over valid pairs. Here the threshold is
   `atol + rtol*|T| = 1e-8 + 1e-5*|T| ≈ 3.0e-5` at `theta = 3`.
5. **`flip_atol`** — the tolerance at which screen (3) flips True→False. For a
   single maximum this is `max_residual` itself: the screen is True iff
   `atol >= max_residual`.
6. **`flip_atol_paired`** — the tolerance at which screen (4) flips:
   `max_i(|r_i| - rtol*|T_i|)`, clipped below at 0.
7. `rho` — `median|theta_hat - T| / median|T - theta_requested|`, the
   definition reused verbatim from `trial_blindness.information_ratio`, plus
   `generator_noise_share = 1/(1+rho)`.
8. `n_valid_pairs`, `n_nonfinite_reference`, `n_nonfinite_estimate`.
9. `generator_is_saturated`, `max_cellmean_deviation` (§2.3).
10. `median_recovery_ratio`, `difference_vs_requested` — context, not residuals.

**Items 3–6 are reported side by side deliberately.** `np.allclose(r, 0.0)` and
`np.allclose(estimate, reference)` are both things a practitioner types, they
differ by about three orders of magnitude at `theta = 3`, and which one gets
typed decides what the screen can see. That gap is a result and it belongs in
the primary table rather than in a footnote.

### 7.1 The parametric control in the same table

`parametric-oracle` and `parametric-plugin` are rows of the **same** table, at
the same cohort sizes and seeds, so the contrast between exact and approximate
equality is *visible in one place rather than asserted across two documents*.

### 7.2 Cohort-size trend

Reported as a derived column and in `trial_learned_generator_trend`: the sign
and magnitude of `max_residual` regressed on `log n` per
`(generator, estimator, reference)`, classified `shrinks` / `grows` / `flat`
(|slope| within 10% of the residual's own scale ⇒ `flat`). This decides whether
the problem gets better or worse with more data, which is the question a
practitioner scaling up a synthetic cohort actually has.

### 7.3 Capacity

`trial_learned_generator_capacity`: `max_residual` and `T_model` error against
`K ∈ {1,2,4,8}` at fixed `n`, so a reader can see where — and whether — exact
becomes approximate as the generator gains capacity.

### 7.4 What each reference cost

Carried in the sidecar and in §3.2's MC-error column: the resolution floor of
`T_model` per generator, against the resolution floor of `T_draw` (machine
precision). This is the evidence offered for claim (c).

---

## 8 · What we expect, written down before the run

Stated so it cannot be back-filled, and so that agreement with it is visibly
not the product of choosing an analysis that agrees.

**We expect F1 and F3 to FIRE, and therefore prediction (a) to be false as
stated.** The reason is algebraic and was available before running anything:
`T_draw(D~)` and `theta_hat(D~)` are **both functionals of the same synthetic
cohort `D~`**. The estimator is handed `D~`; whether its output coincides with
`T_draw(D~)` is a property of *its own functional form*, and the procedure that
produced `D~` cannot enter. `gcomp-saturated` computes
`sum_g (n~_g/n~)(ybar~_g1 - ybar~_g0)`, which is `T_draw` evaluated by the same
code path, so the residual should be bitwise `0.0` for a learned generator
exactly as for a parametric one. The paper's mechanism — "share most but not
all of their inputs" — describes the relationship between the *fitted
generator* and the *outcome model*, but the screen does not compare those two
objects. It compares two functionals of one cohort.

**We expect approximate equality to appear anyway, from somewhere else.** At
any EM fixed point a Gaussian mixture satisfies
`sum_k pi_k mu_k = ybar` **exactly**: the M-step gives
`mu_k = sum_i r_ik y_i / sum_i r_ik` and `pi_k = sum_i r_ik / n`, so
`sum_k pi_k mu_k = (1/n) sum_i y_i sum_k r_ik = ybar`, because responsibilities
sum to one. A GMM-based outcome model's implied cell mean therefore **is** the
sample cell mean, up to how far EM was allowed to run. So
`gcomp-gmm-outcome-K1/K4` should land near zero but not at it, with the floor
set by `GaussianMixture(tol=1e-3)` and `max_iter` — **numerical constants, not
statistical properties of the generator.**

That is the same shape as the finding already in `trial_blindness.py`, where
`ipw_logistic_l2`'s residual floor is `lbfgs`'s convergence tolerance rather
than any modelling decision. If it reproduces here, **the honest replacement
for the paper's prediction (a) is: approximate equality under a learned
generator arrives from a solver's convergence constant, not from the generator
sharing some-but-not-all inputs with the estimator** — and a practitioner can
move it by four orders of magnitude by passing `tol=1e-8`, without changing a
single statistical decision.

**We expect (b) to hold only for that one wiring.** `gcomp-gmm-outcome-*` at
default `tol` should sit above `atol=1e-8` and below `3e-5` for at least some
cells — `allclose(r, 0.0)` False, `allclose(estimate, reference)` True — which
would make the answer to (b) *depend on which of the two lines the practitioner
typed*. `gcomp-mlp-outcome` should be clearly non-zero and caught by both.

**We expect F4 to fire** under `T_draw` and **not** under `T_model`: capacity
changes the synthetic cohort's distribution, hence `T_model`, and cannot change
whether two functionals of that cohort agree.

These expectations are **predictions, not results.** They are as falsifiable as
the paper's and are reported as wrong if they are wrong.

---

## 9 · Deliverables

- `src/harness/trial_learned_generator.py` — module + CLI matching
  `trial_recovery.main` (`--seed`, `--replicates`, `--seeds`, `--results-dir`,
  `--allow-dirty`, plus `--quick`).
- `results/<date>_<sha>/trial_learned_generator_primary.parquet` (+ `_trend`,
  `_capacity`), each with a `.meta.json` provenance sidecar.
- Tests in `tests/test_learned_generator_can_fail.py`, in the style of
  `tests/test_checks_can_fail.py`: **every guard gets a committed input that
  forces it to FAIL**, plus a positive control.
- The sidecar records the **pinned** environment (`scikit-learn=1.5.1`,
  `numpy=1.26.4`, `pandas=2.2.2`) alongside the **live interpreter's** actual
  versions. They differ on this machine and the gap is recorded rather than
  papered over, because a residual floor set by a solver's default tolerance is
  exactly the kind of number that moves between library versions.

---

## RESULT

**Not run.** This section is filled in after the tables exist, in the same
commit that adds them.

---

## Amendments

*(none)*
