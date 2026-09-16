# Independent review of the repaired-interval and Claude branches

Reviewed against `wmhs/build-fixes` at `9c6afbc`, after fetching origin.
The user's branch spellings resolve to `origin/wmhs/repaired-interval`
(`ae8c109`, results produced by `effd3b2`) and
`origin/wmhs/claude-revision` (`0176c335`, alternate paper in
`papers/wmhs_claude`). Neither branch was merged. The paper being edited is
`papers/whms_bode`; the other manuscript and incoming implementations are intact.

## What was worth adding

**Source substitutions give evidence beyond interval substitution.** The incoming
null diagnosis evaluates four distribution families and two arm-size designs at
50 and 800 cells, in four cohort/pool conditions. Each rate uses five streams
of 400 studies and 2,000 bootstrap resamples. At 50 tumor cells:

| Condition | Empirical, 800 reference cells | Remove zeros and rescale mean | Gaussian, same mean/variance | Empirical, 50 reference cells |
|---|---:|---:|---:|---:|
| SMC / pooled | 22.75% | 8.80% | 5.45% | 8.95% |
| SMC / reference | 9.70% | 7.10% | 5.45% | 7.45% |
| KUL3 / pooled | 31.80% | 12.65% | 5.45% | 11.45% |
| KUL3 / reference | 16.95% | 10.95% | 5.45% | 8.75% |

These results support sensitivity to source shape and arm imbalance. They do
**not** identify separate causal shares. Removing zeros also changes variance,
skewness and higher moments; replacing a distribution with a gamma or Gaussian
changes its shape; reducing reference-arm size changes precision as well as
balance. The balanced study evaluates the null, not detection at the alternative.
It supplies no recommendation to discard reference observations. The Gaussian
rows use the same standardized random draws across conditions, so their identical
rates are not four independent replications. At 800 cells the designs coincide;
the empirical rates range from 4.55% to 6.40%.

**The equal-arm symmetry explanation is useful, after correcting its proof.**
For independent identically distributed arm means A and B, exchangeability gives
A-B the same distribution as B-A. This proves symmetry, even with a skewed source.
A zero third moment alone does not prove symmetry. For source third central moment
mu3, the unscaled difference of independent means has third central moment
mu3*(1/n_T^2 - 1/n_N^2); equal sizes cancel it, but symmetry follows from the
stronger exchangeability argument. Symmetry does not itself guarantee percentile
coverage. The main text now connects the mechanism with this distinction.

**Changing to patients as the resampling unit has a measurable cost here.**
At 50 cells and two held-out patients, paired patient-t has 2.77–6.29% null
rejection but only 4.21–8.02% detection among returned intervals. Alternative
abstention ranges from 0.35% to 34.55%. These are useful evidence against an
unqualified instruction to “just bootstrap patients.” However, this method equally
weights the patients present in both arms, unlike the original cell-weighted pool
contrast. Its saved coverage is against the sampled-reference benchmark. The
paper reports the rejection/detection/availability tradeoff and target change,
not a fixed-pool coverage comparison or a verdict about patient inference generally.

**Exact assignment balance changes which estimators coincide.** Claude's assignment
section contains a valid algebraic consequence of the existing OLS weighting
identity. Exact 1:1 allocation within each stratum makes OLS's normalized
n_g*p_g*(1-p_g) weights equal n_g/n; the unadjusted difference also has those
weights. The appendix now says so. Balance in expectation does not ensure this
realized identity. A numerical fixture checks the three estimators directly.

## What was not imported, and why

1. **Repair-pass labels are not fixed-pool validation.** In the incoming driver,
   `truth = sample.truth.parametric[...]["intrinsic"]` is the sampled-reference
   quantity at the alternative. At the null it is zero, so null rejection remains
   interpretable. At shift 0.5 its inclusion rate cannot replace fixed-pool
   coverage. Furthermore, `verdict()` does not use alternative coverage at all:
   it checks null rejection against 5% + 2*MCSE, availability, width and detection.
   Failing to distinguish a rate from 5% is not an equivalence demonstration, and
   this tolerance changes with simulation budget. The main paper's existing
   fixed-pool comparison is retained. Exported column names explicitly call the
   incoming quantity `draw_reference_inclusion` and label verdicts as source
   classifications.

2. **“Eleven of 144 meet the null target” misstates the counts.** Eleven rows have
   the source's joint `pass` verdict; 46 satisfy its null-tolerance criterion.
   The single additional non-baseline pass at 100 cells is not a validated
   reporting range. The claimed “1.2 Monte Carlo SE gap” is not a paired comparison:
   saved marginal counts do not identify the within-study covariance of methods.
   The claimed full factorial also overcounts: at holdout size five, only the
   percentile baseline and two patient-level methods were evaluated.

3. **“All patient-level methods abstain” is false.** With two held-out patients
   at 50 cells, null abstention is 0.4%/0.75% in SMC pooled/reference and 33.25%
   in both KUL3 readings. At five held-out patients and 50 cells it is zero.
   Failures differ by method and setting: excess rejection, low detection and/or
   limited availability. Three reachable bootstrap compositions per arm do not
   themselves trigger the implementation's abstention condition, which requires
   fewer than two represented patients (or degenerate variance for patient-t).
   No blanket claim about the number of clusters being the sole cause is warranted.

4. **The attribution does not establish “58–73% caused by skew.”** That percentage
   is a ratio involving a gamma substitution, not an isolated intervention on
   skewness. Matching mean and skew does not control other standardized moments.
   Some gamma skews are clipped; smaller skew alone does not guarantee smaller
   rejection, so the claimed lower-bound interpretation is unsupported. The
   original code's attribution combines overlapping comparisons and is flagged
   `over_explained` in all four 50-cell conditions. Its residual is not identically
   zero, contrary to the alternate appendix. The main paper uses directly measured
   contrasts without presenting their sum as a partition of causes.

5. **BCa and bootstrap-t do not merely fix a z-versus-t term.** Their purposes and
   constructions go beyond replacing a critical value. The source BCa uses pooled
   delete-one differences with n-1 denominators and a strict below-point bias
   correction; our implementation follows SciPy's multisample influence formula
   and midpoint ties. The source studentized method drops invalid resample pivots;
   ours marks an interval unavailable if any pivot is invalid. These are different
   implementations as well as different outer draws, so the two comparisons are
   not pooled. The existing paper's code and results remain unchanged. The source's
   Gaussian “closed form” is an approximation involving population-variance Welch
   degrees of freedom and normal bootstrap tails; agreement within Monte Carlo
   error does not make it an exact finite-resample formula or universal floor.

6. **The alternate 50-cell table mixes experiments.** Its percentile entries are
   diagnosis rates (22.75%, 9.70%, 31.80%, 16.95%), while the candidate verdicts
   come from the separate generated-sample repair study, whose baseline rates
   are 23.05%, 8.90%, 31.75%, 17.35%. Both studies are meaningful, but this cannot
   be presented as a same-draw table. The imported exports keep them distinct.

7. **The older 800-cell figure is not proved to be “just Monte Carlo noise.”**
   The new experiment changes B from 200 to 2,000 and expands the streams. Its
   estimates are more precise, but the observed difference is not isolated to
   one explanation. Our paper already reports the separate 5,000-study precision
   run and retains its uncertainty rather than replacing it opportunistically.

8. **Do not move the core fixed-pool evidence to the appendix.** The current paper
   leads with the measured reference/parameter discrepancy, which is more specific
   than the alternate abstract's general validation advice. The alternate also
   restores “changing the interval alone does not repair the rule,” broader than
   the tested methods warrant. The current wording, explicit nominal-versus-
   historical standards, and compact ECA illustration are kept. Much of Claude's
   other wording is already inherited from earlier shared revisions.

9. **The restored clean-control table adds length, not needed evidence.** Giving
   the same point estimate a much narrower interval illustrates a known distinction
   already made directly by the audited intervals. The extended assignment table
   and large historical-grid figure also remain local analytical material. Only
   the concise assignment-design consequence is adopted. This keeps the paper
   focused and avoids growing the appendix for redundant demonstrations.

## Verification and reproducibility

- Read the producing code, dated design/amendment, five result tables and sidecars.
- Reconstructed pooled repair counts from all 1,440 seed rows; checked conditional
  rates, complementarity at the null, abstentions and the 144-row summary.
- Reproduced all 64 diagnosis rates from the original local source data with the
  pinned incoming code, all five streams and B=2,000. This verifies 128,000
  simulated interval evaluations, including coincident designs and shared
  Gaussian streams; it does not add independent evidence by repeating a run.
- Incoming algorithm and relevant guard checks: 86 passed, 2 skipped; seven
  repository-wide checks were excluded from the partial extracted source tree.
  An initial broad invocation failed on absent history/results/handoff files in
  that extraction, not on the interval calculations. Those failures are not
  described as a clean full-branch test suite.
- Five new integration checks verify count reaggregation, source hashes, stated
  ranges, availability, nonadditive attribution, symmetry and assignment balance.
  Full paper checks and PDF/source-bundle verification accompany this revision.

`audit_incoming_interval.py` reads the pinned Git objects and regenerates the
aggregate exports with hashes. The anonymous source bundle includes diagnosis
rates, the repair summary and a README explaining its target/denominator limits.
Local analysis additionally retains seed rows, pooled counts, attribution output
and provenance. No cell-level data are added to the repository.

To reproduce the diagnosis from local Lee data, extract `src/` from `ae8c109`
into a temporary directory, then run `reproduce_incoming_diagnosis.py` with
`--source-dir` set to that directory, `--raw-dir` to `data/raw/lee`, and `--work-dir`
to a scratch directory. This retains raw cached arrays only in that scratch
location and verifies the saved aggregate rates. Scientific conclusions depend
on the conditional source-pool design; this is not new patient validation.
