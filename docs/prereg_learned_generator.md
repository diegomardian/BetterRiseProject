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

**Run 2026-09-13** at `19e016b`, tables in `results/2026-09-13_19e016b/` —
`trial_learned_generator_primary.parquet` (840 rows), `_trend.parquet` (140),
`_capacity.parquet` (480), each with a `.meta.json` sidecar. 7 generators x 10
estimators x 6 cohort sizes x 6 seeds x 30 replicates = **180 draws per cell**,
52,920 estimates. Seed `20260913`.

**Headline: prediction (a) is FALSE AS STATED, and prediction (b) fails in a way
worse than "the boolean screen is blind" — the residual scale rho, which the
paper offers as the remedy, is blind in exactly the same case.**

### 1 · F1 and F3 FIRED — the natural wiring is bitwise zero, learned or not

`gcomp-saturated` is the wiring `refdesign.tex` describes: compute the effect on
the cohort you just generated. Its residual against `T_draw`:

```
max |theta_hat - T_draw| = 0.0     exactly, in 42 of 42 cells
                                   30 of 30 learned-generator cells
                                   is_exactly_zero = True everywhere
```

Not `1e-15`. **Bitwise zero**, at every cohort size, for `gmm-joint-K1/K2/K4/K8`
and `mlp-conditional-gaussian` exactly as for `parametric-oracle`. The
parametric control and the learned arms are the same number in the same table:
`0` and `0`.

**F1 fired**: the equality did not degrade to approximate. **F3 fired**: the
residual is identical for parametric and learned generators, so the generator's
learnedness is not the operative variable.

**F2 did not fire.** `np.allclose(residual, 0.0)` is `True` on the natural
wiring — but because the equality is *exact*, not because it is approximate. The
paper is right that a practitioner running `allclose` finds nothing, and wrong
about why.

### 2 · F7 FIRED — both halves of (b) fail together, and this is the finding

Prediction (b) has two halves: `allclose` returns True and sees nothing, **while
the residual scale rho remains nonzero and does see it.** The second half is the
remedy the paper offers for the first.

```
rho for gcomp-saturated, over every generator and every cohort size:
    unique values = [0.0]
```

**rho is identically zero, for the same reason `allclose` is True: the residual
is identically zero.** `rho = median|theta_hat - T| / median|T - theta|` has a
zero numerator, so the scale statistic is blind in precisely the case it is
offered to rescue. There is no regime here in which `allclose` misses something
that rho catches.

This is a failure of the **proposed remedy**, not only of the predicted
mechanism, and it is the most useful thing in this experiment: a reader who
accepts the paper's audit and adds rho as the safeguard has added nothing.

### 3 · The blind band is real, and only PENALISATION reaches it

41 of 420 `T_draw` cells land strictly between the two screens. **Every one of
them is a Ridge arm.** Nothing else in the experiment is anywhere near the band:
every other estimator is either machine-precision or caught by both screens.

```
                               max |theta_hat - T_draw|, over all cohort sizes
gcomp-saturated                0                      <- invisible to both
ipw-saturated                  1.62e-14               <- invisible to both
gcomp-gmm-outcome-K1           9.33e-15               <- invisible to both
gcomp-gmm-outcome-K4           2.22e-14               <- invisible to both
gcomp-ridge-outcome-a1e-5      8.65e-05               <- BLIND BAND (31 cells)
gcomp-ridge-outcome-a1e-3      8.64e-03               <- BLIND BAND (10 cells)
gcomp-ridge-outcome-a1e-1      0.780                  <- caught by both
gcomp-mlp-outcome              3.30                   <- caught by both
ols-stratum-dummies            1.71                   <- caught by both
unadjusted                     11.3                   <- caught by both
```

**The two flip tolerances, side by side.** `np.allclose(a,b)` tests
`|a-b| <= atol + rtol*|b|`, so the second argument is the one `rtol` scales:

| screen | threshold | at theta=3 |
|---|---|---|
| `np.allclose(residual, 0.0)` | `atol + rtol*0` = **`atol` alone** | **`1e-8`** |
| `np.allclose(estimate, reference)` | `atol + rtol*abs(T)` | **`~3.0e-5`** |

The band between them is a factor of **~3,000**, and it is not hypothetical: at
`n=5000`, `gcomp-ridge-outcome-a1e-5` on `gmm-joint-K4` has residual
`9.567e-08`, `flip_atol = 9.567e-08` and `flip_atol_paired = 0.0` — the first
screen fires, the second never does at any cohort size. **Which of the two lines
a practitioner types decides whether the reuse is reported or invisible.**

**F8 did not fire**: the band is reachable, so the verdict on (b) is a measured
verdict and not "not tested".

### 4 · Both proposed mechanisms for (a) are dead, for different reasons

**The paper's mechanism — "a learned generator and a re-fit outcome model share
most but not all of their inputs" — fails on the algebra.** The screen never
compares the fitted generator with the outcome model. It compares `T(D~)` and
`theta_hat(D~)`, both functionals of the *same* synthetic cohort. How `D~` was
produced cannot enter. §1's bitwise zeros are that argument measured.

**This document's own mechanism — an EM convergence floor — fails on the
measurement.** §8 predicted `gcomp-gmm-outcome-*` would sit near zero with a
floor set by `GaussianMixture(tol=...)`. Measured:

```
GaussianMixture(K=4, tol=1e-3 )  ->  |sum_k pi_k mu_k - ybar| = 3.6e-15
GaussianMixture(K=4, tol=1e-10)  ->  |sum_k pi_k mu_k - ybar| = 5.3e-15
```

Four orders of magnitude of tolerance, no movement. The identity holds after
**any M-step** — it needs only that responsibilities sum to one — so convergence
is irrelevant and the mixture outcome model is degenerate to machine precision
(`9.33e-15` at K=1, `2.22e-14` at K=4). Recorded as wrong in Amendment 1, before
the grid ran, and asserted as a parametrised test.

Both mechanisms are stated here because a reader holding the paper's claim will
propose one of the two.

### 5 · The replacement claim, as the data support it

**Approximate equality is a property of REGULARISED ESTIMATION, not of the
generator being learned.** The only arms in the blind band are penalised, their
residual is a closed form in the penalty — `ybar_gd * n_gd / (n_gd + alpha)` on
an orthogonal cell design — and it tracks `alpha` across three decades:
`8.65e-05`, `8.64e-03`, `0.780` at `alpha = 1e-5, 1e-3, 1e-1`. The generator is
irrelevant to it: the same three magnitudes appear under `parametric-oracle`.

This is the **second independent instance** of the phenomenon in this
repository. `trial_blindness.ipw_logistic_l2` already records that an
L2-penalised logistic propensity's residual is "a function of a regularisation
constant" and that turning the penalty off leaves a floor at the solver's `tol`.
Arriving again, in a different estimator class, against a different reference,
under a learned generator, makes it general rather than an artefact of one
sklearn default.

**So the paper's framing — that this matters *more* for learned generative
models — is wrong as stated.** What the data support instead: *the residual
screen's blind band is entered by penalised and iteratively-fitted estimators,
whose distance from the reference is set by a regularisation constant or a
convergence tolerance rather than by any statistical decision; this is
independent of whether the generator was learned.*

### 6 · Capacity, cohort size, and the saturation column

**F4 FIRED — capacity does not move the residual under `T_draw`.**

```
max |theta_hat - T_draw| by mixture components K:   K=1      K=2      K=4      K=8
gcomp-saturated                                       0        0        0        0
gcomp-gmm-outcome-K1                           6.66e-15 7.55e-15 9.33e-15 7.99e-15
ipw-saturated                                  1.07e-14 1.62e-14 9.77e-15 1.33e-14
```

Flat, as the algebra requires: capacity changes the synthetic cohort's
*distribution*, and the screen compares two functionals *of* that cohort. Under
`T_model` capacity does move things (`4.07 -> 2.86 -> 2.72 -> 2.96`), because
that reference is a statement about the generator. **This is the cleanest
statement in the experiment that learnedness is not the operative variable.**

**Cohort size: the residual shrinks, and that makes detection WORSE.**

```
estimator                  trend                 log10 slope   n=100      n=5000
gcomp-saturated            at-machine-precision          --        0           0
gcomp-ridge-outcome-a1e-5  shrinks                    -1.73  8.65e-05    1.03e-07
gcomp-ridge-outcome-a1e-3  shrinks                    -1.73  8.64e-03    1.03e-05
ols-stratum-dummies        shrinks                   -0.818     1.71       0.118
```

The penalised arms fall like `n^-1.7`. Tracking one cell of
`gcomp-ridge-outcome-a1e-5` on `gmm-joint-K4` across the sweep:

```
n=100   5.54e-05   caught-by-both
n=200   3.97e-05   caught-by-both
n=500   2.55e-06   BLIND BAND
n=1000  8.24e-07   BLIND BAND
n=2000  3.09e-07   BLIND BAND
n=5000  9.57e-08   BLIND BAND
```

**A penalised estimator that both screens catch at n=100 crosses into the blind
band at n=500 and stays there.** More data does not ease the problem; it moves
the estimator further into the screen's blind spot, because the estimator
converges onto the reference's functional as the penalty is swamped. This was
not anticipated in §8.

**Saturation (prereg 2.3), reported and not filtered.** `parametric-plugin` is
saturated by construction at exactly `0.0`, which calibrates the column. **No
learned generator is saturated** — deviations `9.60` (K=1), `4.65` (K=2), `2.56`
(K=4), `2.28` (K=8), `0.032` (MLP) — so §1's bitwise zeros are not the parametric
generator in disguise. **F5 did not fire for any learned arm.** Capacity does
bring the generators closer to the training cell means, monotonically, and the
`T_draw` residual stays at `0` throughout.

**F6 did not fire.** The negative control `ols-stratum-dummies` has minimum
`max_residual = 0.0334` over every cell, above its pre-specified `1e-2` floor,
so the screen separates shared-functional from merely-accurate and the rows
above are readable.

### 7 · Evidence for claim (c)

`T_model` costs what §3.2 said it would. Its Monte-Carlo standard error is
`0.0196`-`0.0266` for the mixture generators against `0` for the closed-form
ones, so it **cannot resolve residuals below ~2e-2** — six orders of magnitude
coarser than `T_draw`, which resolves `0` from `1e-14`. Worse, it does not
discriminate estimators at all: `max_residual` against `T_model` is `4.07` for
*eight* of the ten estimators, because the generator's own bias dominates every
estimator's departure. **Choosing the model-based reference makes the screen
measure the generator instead of the estimator.** Claim (c) is supported, in the
specific sense that the reference choice is consequential and the cheap-looking
one is the wrong one.

### 8 · Exclusions and environment

Non-finite outputs are counted, not dropped silently: 41 draws lost the
reference to a positivity failure (an empty stratum-arm cell), and
`gcomp-gmm-outcome-K1/K4` additionally abstained on 109 and 417 draws where a
cell held fewer records than components. `gcomp-saturated` retains `7,517` valid
pairs. No cell reached zero valid pairs; the `no-valid-pairs` abstention is
exercised by its forcing test instead.

The interpreter is **newer than `env/*.yml` pins** and the gap is in every
sidecar: `scikit-learn` 1.5.1 pinned / **1.7.2 live**, `numpy` 1.26.4 / **2.3.5**,
`pandas` 2.2.2 / **2.3.3**, `scipy` 1.13.1 / **1.16.3**. This matters here more
than usual: §4 and §5 turn on solver defaults, and those move between versions.
The bitwise-zero results in §1 are algebraic and do not.

### 9 · What this does not say

It does not say the paper's audit is wrong. The residual screen works: it
separated the degenerate estimators from the informative ones at every cohort
size, and the negative control behaved. It says that **(a) is false as stated,
(b)'s remedy fails with it, and the phenomenon the paper attributes to learned
generators belongs to regularised estimation instead.** One generator family,
one outcome type, one effect structure, `theta = 3`; the algebraic results in §1
and §2 generalise, the magnitudes in §3 and §6 are this design's.

---

## Amendments

### Amendment 1 — a positive control for the blind band

**Timestamp:** 2026-09-13T00:00:00Z (UTC) · **Before any result table exists.**
**Status of the original text: unchanged. No falsifier in §6 is weakened,
reworded or removed.**

**What prompted it.** A micro-check of one library identity — not the
experiment — run while building the module:

```
GaussianMixture(K, tol=1e-3 ).fit(y)  ->  |sum_k pi_k mu_k - ybar| = 3.6e-15
GaussianMixture(K, tol=1e-10).fit(y)  ->  |sum_k pi_k mu_k - ybar| = 5.3e-15
```

**§8's stated expectation is wrong and is recorded here as wrong.** §8 predicted
that `gcomp-gmm-outcome-*` would sit near zero with a floor set by
`GaussianMixture(tol=...)`. It does not. The identity `sum_k pi_k mu_k = ybar`
holds exactly after **any M-step**, not merely at convergence: it requires only
that the responsibilities sum to one, which they do at every iteration. So
`tol` is not a floor, and the GMM-based outcome model is degenerate to machine
precision rather than approximately.

**Why that is a problem for the experiment as pre-registered.** With §8's
mechanism gone, every estimator in §4 is expected to be either exact to machine
precision (`gcomp-saturated`, `ipw-saturated`, `gcomp-gmm-outcome-*`) or clearly
non-zero (`gcomp-mlp-outcome`, `ols-stratum-dummies`, `unadjusted`, all `>=1e-2`
in the published parametric tables). **Nothing in the pre-registered set can
land between `atol = 1e-8` and `1e-8 + 1e-5*|T| ~ 3e-5`** — the band in which
`np.allclose(r, 0.0)` is False while `np.allclose(estimate, reference)` is True,
which is the regime prediction (b) is *about*.

An experiment in which F2 fires **because no arm could have landed in the band**
is an experiment that cannot come out against us. That is the paper's own
complaint, pointed at this design, and it has to be fixed before the run rather
than explained afterwards.

**The addition — three arms, one knob, analytically transparent.**

| name | outcome model | implied cell mean |
|---|---|---|
| `gcomp-ridge-outcome-a1e-5` | `Ridge(alpha=1e-5, fit_intercept=False)` on saturated cell indicators | `ybar_gd * n_gd / (n_gd + alpha)` |
| `gcomp-ridge-outcome-a1e-3` | `alpha = 1e-3` | as above |
| `gcomp-ridge-outcome-a1e-1` | `alpha = 1e-1` | as above |

The design matrix is one indicator per `(stratum, treated)` cell, so it is
orthogonal and the penalised solution is the shrunken cell mean in closed form.
The residual against `T_draw` is therefore a smooth, *chosen* function of
`alpha`, and `alpha` in `1e-5 .. 1e-1` brackets the blind band at these cohort
sizes (bias per cell `~ ybar * alpha / n_gd`, and `ybar ~ 25`).

**What these arms are for, stated in advance so they are not over-read.** They
are a **positive control for the screen's blind band**, exactly as
`tests/test_checks_can_fail.py` demands a forcing input for every guard. They
establish that the band is reachable and that the two `allclose` spellings
genuinely disagree inside it. They are **not** evidence for the paper's
prediction (a): what puts an estimator in the band is a regularisation constant
the analyst typed, not the generator being learned. If the only arms landing in
the blind band are the Ridge ones, that is reported as the finding, and (a)
remains false as stated.

**Added falsifier F8, for the new arms only.** *If no `gcomp-ridge-outcome-*`
arm lands strictly inside the blind band (`allclose(r,0.0)` False **and**
`allclose(estimate, reference)` True) at any cohort size, then this experiment
has no demonstrated ability to exhibit the regime prediction (b) describes, and
its verdict on (b) is "not tested" rather than "false".*

