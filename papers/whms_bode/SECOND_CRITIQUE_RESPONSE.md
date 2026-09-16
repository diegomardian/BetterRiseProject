> Subsequent integration of the pulled external-control example is reviewed in
> `INCOMING_PAPER_REVIEW.md`; the fixed-pool results below remain in the paper.

# Independent assessment of the second critique

The strongest recommendation was to evaluate the fixed source-pool target on
the existing diagnostic draws. That improves the evidence, and the revised
paper now does it. I also accepted the recommendations to reduce textbook
material, promote threshold reversals, enlarge the figure and remove internal
filenames from the review package. I did not accept the characterization of
the existing benchmark as inherently invalid, the interpretation of the cited
source-cell counts, or a categorical claim that the paper is out of scope.

## 1. The target comparison was missing: accepted, with a necessary definition

The drawn-reference score was a real property of the implementation, not a
meaningless quantity. But it did not answer coverage of the generator's fixed
mean effect. Describing that distinction without measuring the alternative
left avoidable uncertainty in the argument.

The generator selects a two-patient source pool before drawing cells. Its
correct conditional target is therefore

    theta_P(s) = f_N (s - 1) mean(GUCA2A in all eligible mature source cells in P).

The mean is computed before either synthetic arm is drawn. It is fixed under
cell resampling conditional on P. Patient-pair selection changes P, and the
reported rate averages conditional-pool coverage over those pairs. Using a
single all-patient pooled mean would change the estimand without changing the
actual two-patient sampling law. The paper makes this distinction explicit;
it does not rename conditional-pool coverage as new-patient coverage.

The original 6,400 configurations were replayed with the same seeds, generator
and interval. Every previously recorded aggregate is unchanged. The 16
alternative conditions also match the original first-seed controlled summaries.
Only the evaluation target was added.

At the 50-cell reporting boundary:

| Cohort / pool | Drawn-reference inclusion | Fixed-pool coverage | 95% Wilson interval for fixed-pool coverage |
|---|---:|---:|---:|
| SMC / pooled | 91.5% | 86.5% | 81.1–90.6% |
| SMC / reference | 92.5% | 91.0% | 86.2–94.2% |
| KUL3 / pooled | 87.0% | 87.0% | 81.6–91.0% |
| KUL3 / reference | 86.5% | 84.5% | 78.8–88.9% |

The first three uncertainty intervals include the design's 90% coverage
requirement; only the KUL3/reference upper endpoint is below 90%. The paper
therefore does not claim statistically established failure of that criterion
in all four cases. None meets both empirical coverage and discrimination
requirements in this seed stream, because discrimination is 52.5%, 78.0%,
54.0% and 55.0%, respectively.

At 800 cells, drawn-reference inclusion is 99–100%, whereas fixed-pool coverage
is 90.5–94%. The paired decreases are 6–9.5 percentage points. Both tables and
Monte Carlo standard errors are reported.

There is also an exact explanation. The fixed-target error contains the full
reference-arm sampling error; subtracting the drawn benchmark instead changes
that term's coefficient from 1 to s. Conditional error variances differ by
f_N^2 (1-s^2) sigma_P^2 / n_N. At s=0.5, three quarters of the reference-arm
variance contribution is removed from the error being scored, even though the
bootstrap still resamples both arms. This clarifies the very high large-count
benchmark inclusion without pretending the variance formula alone determines
finite-sample percentile coverage.

## 2. Too much main-text space for known identities: accepted

Removed the standalone equality section and moved a much shorter set of trial,
censoring and learned-estimator checks into Appendix C. The two machine-precision
residual tables are omitted from the submission; the full numerical source
artifacts remain in the working folder. A short mean/interval example remains
in the appendix. Main-text space now goes to the paired target experiment,
its variance explanation, the null results and the cohort comparison.

## 3. Small patient sample and unexplained cohort gap: partly accepted

The ten- and six-patient evidence base remains small. Cell resampling cannot
repair that, and no broader patient-population claim was added.

However, 662 versus 3,824 compares reference-only versus pooled-tissue cells
within SMC, not SMC versus KUL3. The critique conflated those comparisons.
The revised paper compares matching reference-only definitions across cohorts:
662 source cells in SMC and 844 in KUL3.

Across all 45 and 15 eligible pairs, respectively, median source-pool squared
coefficients of variation are 4.35 and 16.99, and median zero fractions are
39.3% and 62.1%. Under the implemented independent replacement sampling and
binomial thinning, an exact variance calculation gives median signal-to-standard-
deviation ratios of 3.02 and 1.53 at 50 tumour cells and 800 reference cells.
This is a concrete partial explanation for weaker KUL3 detection, despite its
larger source-cell pool. It does not attribute the biological origin of the
dispersion difference or fully explain the 70-versus-400 first-crossing ratio.
Both limits are stated directly.

## 4. A self-audit of a known type of mistake: accepted as a scope limit

The introduction now explicitly calls this a self-audit of an unpublished
project. The zero-of-seven external audit is retained as a limit on external
evidence. The claim to contribution is the measured target discrepancy,
interval behavior, and threshold instability in this case, not discovery of
a new general statistical principle or proof of prevalence elsewhere.

A documented negative result can still be useful. The critique is persuasive
about scope and presentation, but does not establish that the empirical work
has no value. The revisions concentrate on evidence that readers can inspect.

## 5. Venue fit: real concern, but not a categorical exclusion

The [official WMHS call](https://wmhs-neurips.github.io/WMHS/) explicitly includes
negative results, uncertainty quantification, calibration, abstention and
validation/evaluation work. Its central focus remains patient world models.
Thus this paper has a defensible evaluation connection and a weaker direct
application connection than a patient-trajectory or clinical-trial study.

The new title, *Validating simulation-based reporting rules: coverage, false
positives and unstable thresholds*, describes the actual evaluation content.
The abstract immediately names the colorectal single-cell setting. I did not
relabel it as a clinical-trial experiment or invent a deployed-world-model
validation to make it appear closer to the venue. Editorial changes cannot
guarantee that a program committee will find the connection sufficient.

## Other presentation and review-package changes

- Rewrote the abstract to lead with fixed-pool coverage and false positives.
  Removed repeated identity/validity disclaimers from the main text.
- Made the reversal analysis its own subsection and a three-row table, stating
  the first passing, later failing and all-subsequent-passing counts.
- Enlarged Figure 1 vertically, increased text and line sizes, moved the legend
  outside the panels, labelled the operating targets, and annotated the
  nonmonotone pooled KUL3 exclusion curve.
- Replaced internal filenames with Artifact A/B in the manuscript. The review
  artifact retains the original validation sentence and relevant numerical
  context without its file paths. Exact local mapping and checksums remain
  outside the anonymous source ZIP. Build and packaging checks now detect the
  two internal filenames as well as author/repository identifiers.

The final revision is **8 main-text pages, 2 reference pages and 3 appendix
pages (13 total)**. Original experiment implementations and pinned result tables
are unchanged. The added fixed-target scores, source-pool summaries, precision
calculations and analytic checks are confined to `papers/whms_bode`.
