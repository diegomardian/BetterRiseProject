# WMHS: proposed additional runs

Status: all four experiment paths were implemented and exercised on 2026-09-13.
The boundary fixes, pinned submission inputs, recursive full-paper verification,
CI compilation gate, clean control, calibration sensitivity and independent
lifelines comparison are implemented. Experiment 1's controlled union-grid
pilot completed, but its exact eight-seed, 200-replicate final run is deferred:
the pilot projects roughly 200 CPU-hours, and the current working tree contains
uncommitted manuscript work. Do not turn its one-replicate pilot into a claim.

The completed scratch runs are versioned under `/private/tmp` and their
sidecars truthfully record a dirty tree. They are verification evidence, not
submission inputs. Adopt results only after committing the producing code and
rerunning from that clean SHA. The shorter paper remains unchanged.

- Controlled-grid pilot: both cohorts, both draw pools, one seed, one outer
  replicate and 10 inner draws; 47.95 seconds.
- Clean control: 2,000 replicates at each of the null and non-null effects.
- Calibration sensitivity: fixed 45 SMC and 15 KUL3 holdout pairs at 200,
  1,000 and 5,000 inner draws, plus separate leave-one-patient-out runs;
  74.32 seconds.
- Independent Cox comparison: eight bespoke/lifelines pairs, all valid and
  within the fixed `1e-5` tolerance (maximum difference about `1.4e-7`).

## Before any new experiment

Repair the review's reproduced boundary cases: undefined survival groups must
not receive a passing verdict; a zero requested effect must retain difference
metrics without dividing by zero; RMST must distinguish an unidentified tail
from a survival curve that has already reached zero. Add focused regression
cases for those behaviours and for mixtures of valid and invalid outputs.
Use neutral audit outcomes (numerical equality, departure, insufficient valid
pairs), not claims of independence or validity.

Make paper verification read the full paper's actual input graph, including
sections/full/. Pin the result files used for the submission, require their
presence in its verification job, and check manuscript compilation in CI.
Keep these engineering changes separate from new scientific results.

## 1. Isolate grid and binning effects

Question: how much of the 100/90/70 candidate difference comes from evaluated
counts, binning, random draws and simulation budget?

- Keep SMC and KUL3 separate. Preserve the existing target, axis, weighting,
  draw-pool alternatives and acceptance criteria.
- Generate a union of the existing count settings. Key random streams to the
  actual configuration, not its position in a grid, so shared settings reuse
  identical draws. Use the same replicate and inner-bootstrap budgets.
- Compare subsets under fixed bin edges. Separately compare the existing
  adaptive bins, and report per-count rates without binning where supported.
- Record attempted draws, valid intervals, abstentions, coverage,
  discrimination, and crossing brackets. Do not turn a missing crossing or
  coincident cutpoints into an endorsed threshold.

Deliverable: a controlled comparison table and one crossing plot. Report Monte
Carlo uncertainty conditional on the cohort separately from patient sampling
uncertainty. If candidates still depend on configuration, retain that result;
there is no requirement to recover a particular cutpoint.

## 2. Show how residual equality relates to performance

Question: can readers distinguish reference reuse from successful calibration?

Use a prospectively specified small synthetic comparison containing a correctly
calibrated sample mean equal to its empirical reference, the same point estimate
with deliberately narrow intervals, and an estimator with a nonzero residual
and known bias. Include a valid estimator that departs from the reference.
Include null and non-null requested effects. Define the reference before
examining outputs and compute it independently where possible.

Deliverable: a compact table of residual, bias, RMSE, interval coverage and
valid-output rate. Show that neither residual class certifies quality. Do not
score legitimate equality as a scientific failure or use a favourable residual
as a substitute for coverage.

## 3. Check the calibration's two sources of uncertainty

Question: do interval Monte Carlo noise or particular patients change the
candidate around the observed crossing?

On fixed samples near the crossing, compare 200 inner-bootstrap draws with
larger budgets (for example 1,000 and 5,000). Hold the outer samples constant.
Balance representation of the 45 SMC and 15 KUL3 holdout pairs, then assess
leave-one-patient-out sensitivity separately. Keep patient omission distinct
from repeated simulation of an unchanged cohort.

Deliverable: candidate/bracket stability, conditional coverage and
discrimination, and an influence summary naming how many patients change the
conclusion. Stable simulation output is not a population confidence interval.
Start with the relevant count window rather than repeating all 1.7 million
non-null oracle replicates.

## 4. Verify one independent implementation and package the example

Question: does the audit remain usable outside the bespoke estimator code?

First run the existing Cox comparison with lifelines in an environment where
that dependency is required rather than silently skipped. Add a small,
independently implemented estimator/reference example with a fixed tolerance,
explicit invalid-output handling, and requested-parameter performance metrics.
Package the existing comparison as a reproducible command or notebook.

Deliverable: an independently checked comparison and instructions another
researcher can execute in a fresh environment. Do not claim prevalence across
clinical pipelines from this single example.

## Execution order and scope

Complete the boundary and verification fixes, then prioritize the controlled
grid comparison and the clean-control table. Follow with calibration sensitivity
and independent reproduction. Time a small pilot before choosing final budgets;
retain fixed seeds and versioned summaries, not unnecessary simulated matrices.
No GPU or new biological dataset is required for these proposed experiments.

Keep the manuscript focused on validation. Do not reopen the CRC mechanism
search, add unrelated biological endpoints, or increase replicate counts merely
to repeat an algebraic identity. Any new results should replace weaker material
or directly resolve a stated limitation rather than lengthen the appendix by
default.
