# Independent assessment of the third review

## Judgment

The strongest criticisms are actionable: the headline rates needed tighter
Monte Carlo precision, the paper should compare an alternative interval on
the same samples, and the non-monotone detection curve needed analysis. I
implemented those changes rather than accepting the review's predicted
outcomes. The result is a more informative negative finding: the coverage
criterion is not violated everywhere, and a standard alternative helps but
does not repair the 50-cell rule.

The score and the promised increase to 7–7.5 are opinions, not verifiable
acceptance probabilities. A broader conceptual advance or demonstrated
clinical-world-model performance cannot be created by changing the prose.

## New study and what it found

`FOLLOWUP_DESIGN.md` fixed the counts, effects, seed streams, comparator and
summaries before running this follow-up. Five fresh streams give **5,000
studies per setting**, across two cohorts, two tissue pools, four mature-cell
counts and null/halving effects: **160,000 studies and 320,000 intervals**.
All intervals are finite. This samples the same marginal mature-cell
experiment, not the same old random draws. Every source-pair mean and variance
matches the previous analysis; the vectorized percentile routine matches the
original implementation on actual-pool checks. Welch calculations agree with
SciPy's independently implemented unequal-variance test interval.

| Cohort / pool | Percentile coverage at halving | Percentile null rejection | Welch coverage at halving | Welch null rejection |
|---|---:|---:|---:|---:|
| SMC / pooled | 87.14% | 25.38% | 90.16% | 23.46% |
| SMC / reference | 91.02% | 9.80% | 92.52% | 8.86% |
| KUL3 / pooled | 84.66% | 32.24% | 89.70% | 30.78% |
| KUL3 / reference | 86.30% | 19.02% | 88.14% | 18.76% |

These are 50-cell results. The new manuscript table includes Wilson intervals
for coverage, null rejection and detection, rather than only point estimates.
Three original-interval coverage upper endpoints fall below 90%; the
SMC/reference lower endpoint is above 90%. Both methods' null-rejection lower
endpoints exceed 5%, and both methods' detection upper endpoints fall below
80% in all four settings. Thus more simulation **did not turn every coverage
estimate into a failure**, contrary to the review's implication.

The largest binomial Monte Carlo SE falls from 3.54 to 0.71 percentage points.
These are conditional simulation uncertainties; 5,000 resamples do not add
patients or establish transportability. Original results remain preserved.

## Which criticisms I accepted, qualified or rejected

**Monte Carlo precision: accepted; promised conclusion rejected.** The
200-study diagnostic was useful for finding problems but too imprecise for
several coverage comparisons. The larger independent run is worthwhile. The
review also conflates outer simulation replication with inner bootstrap
resampling: increasing the latter does not increase the number of independent
studies used to estimate coverage. Both budgets are now explicit.

**Alternative interval: accepted as a paired diagnostic, not as an assumed
fix.** I chose the transparent Welch-Satterthwaite interval for the difference
of means, multiplied by the fixed reference fraction. It uses the same target
and samples and avoids another bootstrap tuning parameter. Its normal/t
approximation is not guaranteed for sparse skewed data; that is why we measure
its performance. It improves coverage by 1.5–5.0 points but leaves too many
false positives. This tells readers that simply substituting this interval
is insufficient. It does not imply no other interval can work, or establish
that an increased threshold alone would solve the problem. The formula was
checked against the [NIST/SEMATECH primary reference](https://www.itl.nist.gov/div898/handbook/prc/section3/prc31.htm).

**Non-monotone detection: accepted.** In pooled KUL3, the larger run gives
49.0%, 44.0% and 57.7% detection at 50, 100 and 800 cells. The corresponding
null rejections are 32.2%, 20.3% and 6.5%. The trough persists, although the
old single-seed 33.5% estimate at 100 cells overstated its depth. At 50 versus
100 cells, the median sample-SE/exact-SD ratio rises from 0.66 to 0.86, and
the fraction below one half falls from 29.1% to 7.6%. All-zero tumour samples
occur in only 0.32% and 0%, respectively. The evidence supports reduced
spurious rejection as sampled variability becomes better represented, followed
by increasing signal-to-noise; it does not prove a complete causal decomposition
of the curve. The new figure shows null and alternative together. The review's
phrase “5-cell dip” is also imprecise: the five-cell detection is spuriously
high, and the trough occurs later.

**Presentation: accepted.** The title is now “Auditing a single-cell reporting
rule: coverage, false positives and target choice.” The abstract leads with
the practical lesson, the evaluated rule, and the paired interval finding.
Direct fixed-pool/null evaluation comes before grid history. The main text
has three tables and two figures, replacing the earlier six tables and one
figure (the review counted seven tables). Repeated sensitivity and reversal
tables were removed while retaining their results in prose and saved sources.
Extended learned-generator/censoring checks remain accessible as analytical
source artifacts rather than taking appendix space needed for the new methods.
The appendix remains three pages.

**Limited novelty: substantially valid, but not repairable through stronger
claims.** The underlying statistical principles are established. The paper
now makes its contribution as a measured audit with a failed attempted repair
more explicit. The derivations explain this experiment rather than purporting
to invent estimand specification or interval calibration.

**Venue fit: a real limitation, not a categorical exclusion.** The current
[workshop scope](https://wmhs-neurips.github.io/WMHS/) expressly includes negative
results, uncertainty calibration and abstention, alongside its clinical-world-model
focus. The empirical single-cell evidence remains indirect for that focus.
The constructed external-control example establishes a possible consequential
population mismatch; its exact self-reference match is arithmetic, but its
bias against the intended target and the ranking reversal are useful measured
illustrations. I retained it without upgrading it to clinical validation.

**Self-audit/prevalence: valid scope limit, retained.** There is no new positive
external audit case, and no amount of resampling changes that. The paper still
states this openly. Marker disjointness rules out direct target-gene inclusion
in the labelling score; it should not be taken as proof that all depth or
annotation problems disappear.

**Minor factual qualifications.** The source pools are not all at most 844
cells: that number concerns KUL3 reference tissue; pooled mature-cell totals
are 3,824 and 2,475. The earlier headline already rested on several operating
criteria, not a claim that every coverage condition failed. The revised text
makes that distinction easier to see.

## Remaining limitations

Only one alternative interval was evaluated, fixed before the run. No new
threshold was selected and no independent patient dataset was added. The
larger simulation budget settles the stated conditional precision questions,
not clinical generalization, method optimality, or an acceptance prediction.

## Verification

- 29 paper and follow-up tests plus 35 external-control harness tests pass.
- New checks cover the original percentile implementation, independent Welch
  calculations, all denominators and seed aggregation, paired comparisons,
  provenance, generated tables and the revised headline claims.
- The previous 42 pinned result tables are unchanged; new aggregate CSVs and
  a provenance sidecar are separate outputs.
- The final PDF has **8 main-text pages, 2 reference pages and 3 appendix
  pages**. All 13 rendered pages were inspected. Build and anonymity checks
  pass with no errors, undefined references or overfull boxes.
- The portable source bundle independently compiles to 13 pages with identical
  extracted page text. Exact checks are recorded in `build_verification.json`.

This revision is based on `ccb6264`. Additional remote control-panel audit
commits arrived during this review and are not incorporated in this local
revision; they are separate from the interval follow-up evaluated here.
