# Fixed design before the precision/interval follow-up

This analysis responds to review after commit ccb6264. It is an exploratory
follow-up on the same empirical patient pools, not a preregistered clinical
study or validation on new patients. The settings below were fixed before
running the new simulations.

- Four cohort/pool combinations: SMC and KUL3; pooled or reference tissue.
- Mature-cell counts 5, 50, 100, 800; reference arm 800; f_N=0.4.
- Effect multiplier 0.5 and null multiplier 1.0.
- 5,000 outer studies per combination (160,000 total), from five new seed
  streams 202609160 through 202609164, 1,000 studies each.
- Uniform selection of an eligible two-patient pair, independent replacement
  draws of mature-cell target counts, then binomial thinning of tumour counts.
  This is the exact marginal law of the original generator for the intrinsic
  contrast; non-mature cells do not enter that statistic. It is a fresh Monte
  Carlo run, not a bitwise replay of earlier draws.
- Compare the existing 200-resample percentile interval with a two-sided
  Welch-Satterthwaite t interval on exactly the same samples. Multiply the
  difference-of-means interval by 0.4. Welch is a transparent comparator,
  not an assumed repair for zero-inflated, skewed counts. No interval selection
  or new cell-count threshold is fit to these outcomes.
- Use fixed eligible-pool truth for coverage, and retain drawn-reference
  inclusion to isolate the scoring difference. Report null rejection,
  alternative rejection, signed rejection, interval width, missing intervals,
  and paired differences between methods with Monte Carlo standard errors.
- Diagnose non-monotonicity using both effect and null curves, all-zero samples,
  and the ratio of sample-based standard error to the exact pool-conditional
  standard deviation. No causal mechanism is claimed from a single summary.
- Report Wilson intervals and seed-specific rates. More draws reduce Monte
  Carlo error; they add no independent patients.
- Validate the vectorized percentile calculation against the original interval
  routine and Welch calculations against scipy.stats.ttest_ind. Verify the
  new pool summaries against the existing 120 pair summaries before running.
- Preserve the previous diagnostic outputs and historical grid selections.
  Report any unfavorable result from the comparator rather than trying a
  series of methods until one passes.
