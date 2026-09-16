# Assessment of the Monte Carlo / supporting-evidence critique

Pulled before editing; the branch was current at `02c7405`. Reviewed the source
graph, saved grid rates, precision outputs, estimator implementation, learned
generator design and the cited methodological literature. The revisions below
use existing evidence; no new simulation was run.

## 1. Local reversals: demote the empirical claim

**The concern is valid; the assertion that all underlying curves are monotone
is not established either.** Recomputed the cited counts and Wilson intervals
from the pinned SMC/reference dense-grid results. Every listed setting has 200
valid replicates. Examples at each selected reversal are:

| Seed | Counts | Estimated detection | Marginal 95% Wilson intervals |
| --- | --- | --- | --- |
| 3 | 45 → 50 | 82.0% → 76.0% | [76.09, 86.71]% → [69.63, 81.39]% |
| 6 | 50 → 60 | 80.0% → 79.5% | [73.91, 84.95]% → [73.37, 84.51]% |
| 7 | 50 → 60 | 81.0% → 77.5% | [75.00, 85.83]% → [71.23, 82.74]% |

The marginal Monte Carlo standard errors here are about 2.7–3.0 percentage
points; 3.5 points is the grid-wide maximum, not each point's error. These
marginal intervals are descriptive precision checks, not a simultaneous
monotonicity test. Selection of the first passing count also matters. The
observed crossings do not establish a real decrease, but neither this evidence
nor the review establishes that noise is the entire explanation.

Removed the standalone reversal subsection and the introductory claim that
the grid demonstrates failure above a boundary. The grid section retains a
brief description of the noisy selections and states explicitly that the
5,000-study follow-up does not test counts 45, 60 or 70. Section 4 retains the
logical requirement to assess the reporting range. Existing reversal exports
remain inspectable, clearly classified as descriptive diagnostics.

The higher-precision pooled-draw detection trough is a **different result**:
5,000 studies at 50/100/800 cells, with null rejection and estimated-variance
diagnostics reported beside detection. Retained that evidence and made the
conclusion explicit about which result it describes. No unrun 40–110 study is
implied.

## 2. Interval scope and bootstrap literature

**A fair request for specificity, partly already addressed.** The current
abstract already said “this interval substitution,” and the conclusion named
Welch. It did not conclude that every interval repair fails. Nonetheless,
“substituting a Welch interval does not repair the rule” is clearer, so that
wording now appears in the abstract and the table caption names Welch too.
BCa, bootstrap-t/studentized, and Poisson- or negative-binomial-based intervals
are explicitly untested in the limitations.

Added [DiCiccio and Efron (1996)](https://doi.org/10.1214/ss/1032280214) on
bootstrap confidence-interval accuracy and bias/skewness correction. The paper
now states why the simulation is useful: it quantifies the actual reporting
rule's operating characteristics and a paired comparator, rather than
rediscovering that simple intervals can fail for skewed distributions. Welch
is informative as an analytically computed comparator without bootstrap
resampling noise, not as a claim to the best method for these counts. The
failure of either method does not establish how BCa, bootstrap-t or a fitted
count model would perform.

## 3. SBC: restore and distinguish the related work

**Accepted.** Restored citations to [Talts et al.](https://arxiv.org/abs/1804.06788)
and [Modrák et al.](https://arxiv.org/abs/2211.02383). SBC tests Bayesian
computation using prior-predictive calibration; test-quantity choice affects
sensitivity. The present audit evaluates frequentist coverage and rejection
at specified effects conditional on empirical pools, and documents an actual
recovery-to-reporting-rule inference. It is not a new SBC method or a
replacement for SBC. Removed citations in an earlier structural edit were
an omission worth correcting.

## 4. ECA: surface the construction, retain the diagnostic controls

**Accepted the transparency change, not the claim that the example has no
use.** Its opening now states that estimator and observed-data reference
intentionally share a population mistake, that exact recovery may alert an
analyst, and that it is a pedagogical stress test. Removed the buried duplicate
qualification from the appendix.

The review overlooks an informative existing control: the pooled
potential-outcome reference breaks exact equality yet still favors the
estimator of the wrong population effect. Thus detecting a zero residual alone
does not address the population-alignment issue. This remains an illustration
of an established estimand distinction, not evidence of how often analysts
make the error. Sample-size RMSE wording now explicitly compares the smallest
and largest sizes; it does not imply monotone improvement at all sizes.

## 5. Figure order and the transferable result

**The figure-order criticism is stale.** Figure 1 in the reviewed build already
shows `precision_rejection.pdf`; the historical benchmark is Figure 2. Retained
that order. The latter caption now starts by calling it a historical selection
diagnostic and points readers to the direct interval assessment. Adding a
fixed-pool curve over the full historical grid would require data not supplied
by the four-count precision study; no such curve was fabricated.

**Accepted promotion of the shared-reference calculation.** It has its own
subsection and an abstract statement: under a halving, the drawn benchmark
removes three quarters of the reference-arm contribution to the estimation
error variance. This is not three quarters of total variance, nor a universal
coverage-inflation theorem. The 5.8–6.8 percentage-point difference is reported
specifically at 800 cells; the five-cell reversal in direction remains stated.

## 6. Learned-generator evidence and venue fit

**A concise restoration adds useful scope; a full study would dilute the
paper again.** Restored a compact supporting appendix check and a main-text
pointer. It identifies four Gaussian mixtures, an MLP conditional Gaussian
and two parametric controls, trained/tested using synthetic three-stratum
cohorts. Actual counts are 75,580 evaluations on 7,558 generated cohorts,
with 74,726 finite estimate–reference pairs. The review quoted scheduled counts
(75,600 / 7,560); the appendix distinguishes them and reports all exclusions.
The zero residual for saturated G-computation persists across all 42
settings because the estimator and reference are identical on each evaluable
dataset. Nonsaturated estimators provide a contrasting check.

This is evidence about what a learned generator can change in a reference
check. It is not a learned clinical trajectory model, a clinical realism test,
or an empirical replication of the single-cell reporting error. Venue fit
therefore remains a limitation, not something resolved by adding “MLP.”

## 7. Significance, readability and source maintenance

The evidence is still one empirical self-audit. Zero affirmative cases in a
small, nonrandom external search neither establish prevalence nor establish
that the wider failure mode is rare. Kept that result and its limits; did not
replace it with a popularity claim. The practical contribution remains an
inspectable instance where recovery was taken as evidence for a reporting rule
whose operating characteristics subsequently failed direct evaluation.

Cut the normal-mean quarter-width interval demonstration and its main-text
reference: valid arithmetic, but redundant after the actual interval study.
Compressed historical-grid and estimator-identity prose to accommodate the
learned check without enlarging the appendix. Standardized authored manuscript
prose to US spelling. Expanded the inactive-source inventory in
`ANALYTICAL_ARTIFACTS.md`; local generated diagnostics remain available and
checked, while inactive extended analyses are removed from the submission ZIP.
Namespaced the extended analysis labels to prevent collisions if reused.

## Verification

- 86 existing tests passed (30 paper/precision, 35 external-control, 21
  retained-control audit); the learned-result denominator and identity checks
  are included in the paper tests.
- Built with the unchanged official NeurIPS 2026 style; no LaTeX errors,
  undefined references, overfull boxes or detected identifying strings.
- Final layout: 9 main-text pages, 2 reference pages and 3 appendix pages.
- Rebuilt the anonymous source ZIP, compiled it in isolation, compared all
  page text, and visually checked all final PDF pages.
