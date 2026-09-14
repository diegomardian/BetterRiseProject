# Prevalence audit: does anyone else wire the estimator to its own reference?

**Search date:** 2026-09-13 · **Branch:** `wmhs/prevalence-audit`
**Pre-registration:** [`prereg_prevalence_audit.md`](prereg_prevalence_audit.md), committed at `e7434db` before the first search
**Queue:** [`prevalence_audit_queue.md`](prevalence_audit_queue.md), committed at `15b0680` before the first repository was opened
**Machine-readable tables:** `results/2026-09-13_5c888f7/prevalence_audit_verdicts.parquet` and `…_screening.parquet`, with sidecars

---

## Headline

**We examined 7 public repositories meeting pre-registered criteria. All 7 were
classifiable. 0 were YES, 7 were NO, 0 were UNCLEAR.**

This is a null result and we report it as one.

The stronger statement is the one worth putting in the paper. Six of the seven
were NO at the first step of the rule: the validation reference is the
generating parameter, or a counterfactual quantity the estimator never sees.
That is not a near miss — it is the practice the paper's own analysis
recommends, found in every corner of the queue that had anything to find. The
seventh reached the functional comparison and came out NO there.

**We found no evidence that anyone else in this sample does what our pipeline
did.** The paper must say so.

---

## The pre-registered rule, in brief

Locate the code that (a) generates a draw `D`, (b) applies the estimator to it,
and (c) compares the result to a reference `T`. Then triage `T` by provenance:

| `T` is… | verdict |
|---|---|
| **3(a)** a function of the generator's parameters only | **NO** |
| **3(b)** computed from records the estimator never consumes — both potential outcomes, latent event times, held-out data | **NO** |
| **3(c)** computed from the same observed records the estimator consumes | go to step 4 |
| **3(d)** indeterminable from the code | **UNCLEAR** |

Step 4 asks whether `T` and `θ̂` are the *same functional* — equal identically
on every draw in the support, by algebra. Check input records, then weighting,
then conditioning set, then nuisance estimation (cross-fitting breaks it), then
the observation process; stop at the first difference. **YES** needs
affirmative evidence: a shared code path, a textual identity, or a written
derivation.

Full text in [§3 of the pre-registration](prereg_prevalence_audit.md).

**One amendment, A-1, dated 2026-09-13.** Six of ten pre-registered query
strings returned zero results — `gh search repos` ANDs over name, description
and topics and does not index READMEs, so a six-word natural-language query
matches nothing. A-1 relaxed those six to their indexable cores. It was
prompted by the search tool, not by a repository, and was written before any
repository's code was opened, so no verdict required re-application. The
classification rule was **not** amended at any point.

---

## The table

Verdict basis codes are from [prereg §3](prereg_prevalence_audit.md). `3a` =
reference is a generating parameter; `3b` = reference uses data the estimator
cannot see; `4.4` = estimator and reference differ in nuisance estimation.

| # | repository | commit | domain | estimator (file:line) | reference (file:line) | prov. | verdict |
|---|---|---|---|---|---|---|---|
| 5 | [brycewang-stanford/Auto-Empirical-Research-Skills](https://github.com/brycewang-stanford/Auto-Empirical-Research-Skills) | `d31f17f` | benchmark grading agent-written causal pipelines | cross-fitted partialling-out DML — `benchmark/lib/dml.py:103-119` | full-sample OLS coefficient on `d` — `benchmark/lib/dml.py:89-93`, used at `check_benchmark.py:442` | 3c | **NO** (4.4) |
| 7 | [rje42/causl](https://github.com/rje42/causl) | `0659ed8` | R package: specify, simulate from, and fit causal models | `fit_causl()` MLE of the frugal causal parameter — `tests/testthat/test-fit.R:9` | `pars$y$beta = c(0, 0.5)`, the generating vector — `test-fit.R:4`, compared at `:18` | 3a | **NO** (3a) |
| 8 | [Zhengxian-Fan/target-trial-emulation](https://github.com/Zhengxian-Fan/target-trial-emulation) | `35983ac` | target trial emulation on CPRD + plasmode benchmarking arm | PSM+Cox, IPTW+Cox, TMLE, pooled by Rubin's rules — `estimate_treatment_effects.py:73-152` | `CONFIG['theta_treat'] = -0.5` — `simulate.py:14`, entering the hazard at `:68` | 3a | **NO** (3a) |
| 10 | [natsousa/ecasim](https://github.com/natsousa/ecasim) | `53a201e` | external control arm, propensity-score methods | Cox for the treatment term, unweighted / IPTW / matched — `R/analyze.R` | `log_hr_true <- log(0.7)` — `R/simulate.R:61` | 3a | **NO** (3a) |
| 18 | [PamelaShaw/PlasmodeSimulation](https://github.com/PamelaShaw/PlasmodeSimulation) | `7ae26ee` | plasmode simulation for causal inference (arXiv:2504.11740) | 11 ATE estimators per replicate — `00_plasmode_utils.R:257-263` | `psi0.sim = ey1 - ey0`, standardised over the **original** cohort with treatment forced — `00_plasmode_utils.R:278-286` | 3a + 3b | **NO** (3a, 3b) |
| 24 | [brycewang-stanford/StatsPAI](https://github.com/brycewang-stanford/StatsPAI) | `dd4640e` | benchmark of causal estimators across designs | per-design point estimate — `smart/benchmark.py:238-251` | `spec['true_effect'] = 0.5`, the generator's own `effect` — `benchmark.py:48-104`, `utils/dgp.py:119` | 3a | **NO** (3a) |
| 25 | [igerber/diff-diff](https://github.com/igerber/diff-diff) | `ff03d08` | difference-in-differences library | Callaway–Sant'Anna ATT(g,t) and aggregations — `tests/test_csdid_ported.py` | DGP constants: ATT = 1, 3 with 1.5× dose response — `tests/helpers/csdid_dgp.py:283, 310-323` | 3a | **NO** (3a) |

Full justifications, quoting the code, are in the parquet's `justification`
column. Nineteen further queue entries were walked and screened out; each is
recorded with the criterion it failed in
`prevalence_audit_screening.parquet`.

---

## The one case that reached step 4, and the false positive it produced

Only **one** of the seven computes its reference from the realised simulated
draw at all. That is the branch where a YES could live, and it is worth setting
out in full because the pre-registered procedure changed the answer.

`Auto-Empirical-Research-Skills` grades a candidate's cross-fitted DML estimate
against a reference the checker recomputes from the same CSV:

```python
def true_theta(rows: list[dict]) -> float:
    """Fully-controlled OLS coefficient on d (recovers THETA by construction)."""
    X = [[1.0, _num(r, "d")] + [_num(r, k) for k in CONTROLS] for r in rows]
    y = [_num(r, "y") for r in rows]
    return lalonde.ols(X, y)[1]                     # benchmark/lib/dml.py:89-93
```

Running the shipped code gives

```
true_theta (full OLS)  = 1.4999999999999607
dml_theta  (cross-fit) = 1.5000000000000004
max|dml - true_theta|  = 3.975e-14
```

A residual at machine precision. Under a boolean `allclose` screen this reads
as functional reuse. **It is not.** The paper's own instruction — *"Only an
algebraic argument establishes an identity"* — is what resolves it.

The generator writes `y = THETA * d + _g(x)` with `_g` exactly linear in the
controls and **no noise** (`benchmark/lib/dml.py:48-58`). OLS is linear in the
response and `_g` lies exactly in the span of the controls, so for *any* subset
fit the coefficient vector is exactly `THETA·β_d + β_g`. The residuals
therefore satisfy `r_y = THETA · r_d` **pointwise**, in-fold and out-of-fold
alike, and every ratio `Σ r_d r_y / Σ r_d²` returns `THETA` regardless of how
the folds are cut. The equality is a property of this noiseless,
exactly-linear data file — not of the two functionals.

Adding Gaussian outcome noise to the same rows and changing nothing else
separates them at once:

```
with outcome noise:  true_theta = 1.599471   dml_theta = 1.474678   residual = 0.124793
```

Cross-fitting versus a full-sample fit is exactly the step-4.4 difference, and
it is the same phenomenon the paper reports in its own trial sweep. **Verdict
NO**, on the pre-registered basis, with a constructed counterexample draw as
the evidence §3.4 asks for.

Two things follow that are worth more than the verdict.

**First, this is the false positive the pre-registration was written to catch.**
An audit that ran `np.allclose` and stopped would have reported a YES here and
been wrong. The prevalence question cannot be answered by a screen; it needs
the algebra, per repository.

**Second, the repository draws the distinction the paper argues for, by name.**
Its golds titled `*-recovers-true` compare the candidate to `true_theta`; those
titled `honest-*` compare the candidate to the checker's own run of the *same*
estimator, and are described in-repo as ensuring "no fabricated numbers"
(`benchmark/tasks/dml-recovery.toml:54-62`). The second family *is* estimator
and reference sharing a functional — deliberately, with a stated
anti-fabrication purpose, under a different name. That is precisely §blind's
"equality identifies redundancy", implemented as a separate check rather than
folded into the recovery claim. `igerber/diff-diff` separates its tiers the
same way: recovery against the DGP constant, and agreement against R's output,
labelled apart in the file's own docstring (`tests/test_csdid_ported.py:4-8`).

A reader who took an `honest-*` gold as *the* validation would reach YES on
that gold. The pre-registration classifies the primary experiment, which is the
recovery gold, and we record the alternative reading rather than bury it.

---

## What the null looks like from the other side

The queue also produced a clean set of statements about how public code
actually defines its truth. These are worth the paper's space because they make
the recommendation concrete:

- **From the generating parameter.** `causl` checks `out$pars$y$beta` against
  the `pars` list written down before any draw existed, with a chi-square
  critical value — a check that can fail. `StatsPAI` computes
  `bias = point - spec['true_effect']`. `diff-diff` asserts
  `abs(att_e2 - 3.0) < 1.0`. `ecasim` and `Zhengxian-Fan` pin `log(0.7)` and
  `theta_treat = -0.5` in the generator's configuration.
- **From counterfactuals the estimator cannot see.** `PamelaShaw` standardises
  the fitted outcome model over the *original* cohort with treatment forced to
  1 and to 0, computed once before the replication loop, and compares via an
  explicit `bias = colMeans(est - psi0)`. Two screened-out repositories
  (`Netflix-Skunkworks/oci-agent`, `lfiaschi/bayesian-autoresearcher`) use
  `mean(mu1 - mu0)` from shipped potential-outcome columns.
- **The generator infrastructure does the same.** `cran/Plasmode` hands
  downstream users `RR = mean(p_1)/mean(p_0)` computed from the true
  coefficient vector on the sampled design matrix with exposure set to 1 and 0
  (`R/PlasmodeBin.R:174-180`). It never touches the realised outcomes `ynew`,
  so a downstream estimator *cannot* be the same functional of it. The most
  widely used plasmode generator makes the failure mode hard to commit by
  accident.

One honest asymmetry, noted because it cuts the other way: `StatsPAI`'s
heterogeneous mode draws unit effects as `effect + N(0, 0.3)` yet keeps `0.5`
as the reference, so its reported "bias" mixes estimator error with the
realised-versus-requested gap. That is the *other* failure the paper describes
— conflating `θ̂ − T` with `T − θ` — and it appears in code, unprompted. It is
not the failure this audit was looking for, and we did not pre-register a rule
for it, so we report it as an observation and not as a count.

---

## Counts

| quantity | n |
|---|---|
| queue built (positions) | 41 |
| queue walked (positions 1–26) | 26 |
| screened, not examined | 19 |
| **examined and classified** | **7** |
| classifiable (not UNCLEAR) | 7 |
| **YES** | **0** |
| **NO** | **7** |
| **UNCLEAR** | **0** |
| reference computed from the realised draw (3c, reaching step 4) | **1** |
| code executed for confirmatory numerics | 1 |

Stopping point: queue position 26. The target was 5 (floor 3); the audit
continued past it because the remaining screening was cheap, and stopped when
the queue's yield of eligible repositories thinned. No eligible entry in
positions 1–26 was skipped.

**No prevalence rate is quoted**, per prereg §8. Ten queries are not a random
sample of practice and 7/41 is not a denominator anyone should divide by.

---

## Coverage limits

These are part of the result, not a disclaimer attached to it.

1. **The domain the paper is actually about returned no public code at all.**
   Q9 — `synthetic control arm virtual control arm simulation validation
   reproducible code` — yielded zero examinable artifacts. Its top hit,
   [arXiv:2507.16048](https://arxiv.org/abs/2507.16048), *Evaluating
   virtual-control-augmented trials for reproducing treatment effect from
   original RCTs*, has no code-availability statement; the full text was
   checked. The remainder were vendor and editorial pages. Q1b returned two
   repositories, one of which is an HTML brochure. **Industrial
   virtual-control-arm pipelines are not public.** The population the paper's
   claim concerns is largely unobservable by this method, and no amount of
   further searching recovers it. This is the single most important limit here.

2. **Only one repository reached the branch where a YES can live.** Six were
   decided at step 3 by a reference that is a constant or a counterfactual. A
   rule that decides most cases in one step is a good rule, but it means the
   audit's power to *detect* functional reuse rests on a sample of one at the
   step that matters. A null over n=1 in the operative branch is weak evidence
   of absence, and the paper should not claim more.

3. **Single rater.** One reader applied the rule; there is no second rater and
   no agreement statistic. The rule is written for independent reproduction,
   which is the mitigation and not a substitute.

4. **Reading, not running.** Six of seven verdicts rest on code reading. Only
   the AERS case was executed (it is pure-stdlib and deterministic).
   `PamelaShaw` cannot be run at all — the underlying cohort is not shareable.

5. **Not a random sample.** GitHub's relevance ranking is opaque and favours
   well-described repositories. The round-robin queue prevents cherry-picking;
   it does not deliver representativeness.

6. **Language and platform.** The examined set is 4 Python, 3 R, all on GitHub,
   all English-described. No Stata, SAS, or Julia reached examination; the one
   Stata repository in the queue (`bldestavola/TTE-Short-Course`) was teaching
   material with no DGP.

7. **Absence of code is not absence of the practice.** Nineteen entries were
   screened out, and several were screened out *because* their validation code
   was absent — `mariannawicks-sys` ships placeholder scripts containing only
   their own filenames; `kathoffman` ships simulated demo data whose generator
   is not in the repository. What those pipelines compare against is unknown
   and unknowable from the artifact. The audit measures published, readable
   code, which is a biased view of practice in a direction we cannot sign.

---

## What this does to the paper's claim

The paper's central analysis is algebraic and this audit does not touch it. If
`θ̂(D) = h(S(D)) = T(D)` then the residual is zero and the recovery curve
measures the generator; that holds whether the wiring is common or unique.

What the audit bears on is the paper's *implicit* scope claim — that the
failure mode matters beyond one lab. On that:

- It **does not support** the claim empirically. We looked, under a rule fixed
  in advance, and found nobody doing it.
- It **does not undermine** the paper either. Every repository that had a
  reference at all used one the algebra cannot collapse, which is direct
  confirmation that the recommended practice is the practised one and that the
  recommendation is coherent.
- It **leaves open** the case the paper actually cares about, and sharpens why.
  §matrix already predicts that a *learned* generative patient model makes this
  failure the default wiring rather than an oversight, because the realised
  effect must itself be estimated from the sampled cohort and the natural
  estimator for it is the estimator under test. Every repository we examined
  has a closed-form generator — the favourable case, by the paper's own
  argument. **The audit cannot reach the case the paper predicts is dangerous,
  because that code is not public.** The two limits point at the same gap: the
  learned-generator, virtual-control-arm pipelines are exactly the ones that
  are not on GitHub.

The honest sentence for the paper is that the failure mode is demonstrated, its
prevalence in public code is zero out of seven under a pre-registered rule, and
the setting where it is predicted to arise by default remains unaudited for
want of public artifacts. A reader deciding whether to run the screen on their
own pipeline should be told that this is cheap insurance against something we
did not observe elsewhere, not a correction to a widespread error.

---

## On fairness to the work examined

Every repository in the table does something this audit judges correct, and
several do it with more care than the audit required — `causl`'s chi-square
check against the generating vector, `PamelaShaw`'s standardisation computed
outside the replication loop, `AERS` and `diff-diff` naming their reproduction
checks apart from their recovery checks. Nothing here describes anyone's
intent, competence or care; each verdict describes what a named file does at a
named line, and the underlying code was read before it was characterised.
