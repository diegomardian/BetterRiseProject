# Fixed design: reporting-band and skew-aware interval follow-up

Written after review of bf6b825 and before producing this study's results.
This is an exploratory follow-up on previously examined empirical pools, not
an independent patient validation or a preregistered clinical study.

- Use the same 120 eligible two-patient pools, checked against pinned summaries.
- Counts: 30 (an interior point of the flag-wide band) and 50 (reporting boundary).
  One interior point does not validate the whole 20–49 range.
- Four cohort/pool conditions; shifts 0.5 and 1; 800 reference mature cells;
  reference fraction 0.4; independent cell draws with replacement and binomial
  thinning, using the exact marginal law of the original contrast.
- Five fresh seed streams 202609170–202609174, 400 studies per stream: 2,000
  studies per setting, 16 settings, 32,000 studies and 160,000 interval evaluations.
- Compare percentile B=200, percentile B=2,000, BCa B=2,000,
  bootstrap-t B=2,000, and Welch on every study. The bootstrap methods share
  the same 2,000 draws; the B=200 comparator uses the first 200 of those draws
  in each arm. This has the original B=200 sampling law, not its old random stream.
- BCa: midpoint ties for the bias correction; two-sample jackknife influence
  contributions -(X-Xbar)/n_N and (Y-Ybar)/n_T for acceleration. Follow SciPy's
  two-sample construction, checked against its independent jackknife implementation.
  Undefined bias/acceleration, a nonpositive quantile-transform denominator,
  or reversed/nonfinite endpoints produces an unavailable interval, not a fallback.
- Bootstrap-t: t*=(estimate* - estimate)/SE*, with unbiased within-arm variances;
  endpoints estimate - reversed-quantiles(t*)*SE. If observed SE or any bootstrap
  SE is zero/nonfinite, return unavailable. Count this explicitly.
- Score coverage against the fixed eligible-pool effect, not the drawn-reference
  benchmark. Report coverage and rejection at both shifts, abstention, width,
  Wilson intervals, per-seed rates, paired changes and Monte Carlo SEs.
  At shift 1, valid-interval coverage plus null rejection must equal one exactly.
- Distinguish nominal 95% coverage from the historical minimum of 90% at the
  alternative. Report detection separately. No method or cutoff is selected
  as a replacement from these same studies, even if a point estimate passes.
- Validate against the original percentile and Welch routines; check BCa
  endpoints against scipy.stats.bootstrap on continuous, discrete/skewed and
  sparse samples. Independently check the studentized formula, affine scaling,
  and reversal under exchanging arms. Verify output counts and provenance.
- Report all five methods and both counts, including unfavorable results.
- Preserve original precision code and outputs. Keep only aggregate exports.
