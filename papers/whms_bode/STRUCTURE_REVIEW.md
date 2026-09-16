# Independent assessment of the venue/structure critique

Reviewed the current manuscript, its source graph, the paired interval implementation,
the synthetic external-control implementation and the pinned results. Pulled through
`7c3c818` before editing. This is an editorial and explanatory revision; it adds no
simulation results and makes no claim that restructuring supplies missing evidence.

## Decisions on the six criticisms

1. **Venue fit: a real weakness, partly addressable by structure.** The previous
   version introduced the clinical example late, after a long single-cell audit.
   The [WMHS call](https://wmhs-neurips.github.io/WMHS/) includes external controls,
   validation, uncertainty, calibration and negative results, so the subject has a
   defensible connection. That does not make it a demonstrated patient world model.
   The new title and introduction state one thesis about validating the target,
   interval and reporting policy. Section 2 now leads with the synthetic trial
   example; the biological setup moves into the empirical audit in Section 3.
   The clinical example now explicitly defines the finite enrolled-cohort target,
   the two standardization weight sets, the reference controls and sample-size
   behavior. These additions explain existing code and results, rather than
   dressing the example up as real clinical evidence. Venue fit remains a judgment
   for reviewers; it is not resolved by vocabulary alone.

2. **Novelty: substantially correct, but it does not make the audit valueless.**
   Estimand specification and joint evaluation of operating characteristics are
   established. The introduction now describes the contribution as a reproducible
   worked audit, rather than a new validation principle. Its specific evidence is
   the documented recovery-to-rule inference, paired target comparisons, measured
   operating characteristics, and controlled threshold-selection analysis.
   The review's claim that ICH E9(R1) was already cited was inaccurate for this
   version. Added the [official 2019 guideline](https://database.ich.org/sites/default/files/E9-R1_Step4_Guideline_2019_1203.pdf)
   for its requirement to identify the population targeted by the clinical question
   (Section A.3.3), not as a source for every statistical lesson in this paper.

3. **One self-audit: correct and not repairable by editing.** The synthetic
   registry is an illustration, not an independent empirical replication.
   The external search still contributes no affirmative external case. The
   introduction now distinguishes the two evidential roles immediately, and the
   limitations explicitly retain the lack of prevalence and patient-world-model
   evidence. No new claim of generalization was introduced.

4. **Bootstrap budget: correct observation, narrower rebuttal than suggested.**
   The headline precision follow-up holds the percentile interval at B=200.
   The B=200/1000/5000 sweep concerns a separate candidate-selection experiment
   and cannot be presented as a sensitivity analysis of the headline rates.
   Section 3.3 now states that Welch intervals have no bootstrap resampling noise
   but still reject the null in 8.9–30.8% of the same simulated studies. Thus finite
   inner-bootstrap noise cannot be the sole explanation for the failure shared
   by both methods. This does **not** estimate what a higher B would do to the
   percentile interval, rule out an improvement, or establish a universal minimum
   bootstrap budget. The remaining limitation is stated beside the result.

5. **Coherence: valid concern, but the four-study description is partly stale.**
   The extended estimator residual matrix was already excluded from the compiled
   manuscript; it remains a noncompiled analytical artifact in the source bundle.
   The revised paper now has two central examples and one common question:
   does the validation support the claim made for the downstream rule? The main
   biological-control discussion is shortened to a scope paragraph; supporting
   details remain in the existing appendix. Section 4 organizes the practical
   consequences into target construction, joint rule evaluation, and separation
   of selection, assessment and transfer. It also explains why a legitimate
   finite-cohort reference can depend on the simulated records: the requirement
   is target alignment, not blanket independence from the dataset.

6. **Repeated hedging: partly correct.** Removed redundant general disclaimers
   from the introduction and conclusion, and replaced the long interval caveat
   with a direct statement that this substitution fails the joint requirements.
   The abstract now says “this interval substitution does not repair the rule.”
   It no longer extrapolates from Welch to all possible interval replacements.
   Retained caveats that change a number's meaning: fixed-pool versus new-patient
   inference, diagnostic pre-gate intervals, the historical grid's scoring target,
   and the fact that a comparison of medians is not a decision-error probability.
   These are necessary definitions, not removable rhetorical caution.

## Incoming changes assessed separately

The pull added useful grid and acronym definitions; these are retained, with
acronyms also expanded in the compiled supporting section. “Pre-registered grid”
is narrowed to the grid recorded in the original design, avoiding an unnecessary
claim about the status of that record.

It also restored a duplicate appendix block and sensitivity table. Removed them
again because the same numerical evidence is already reported and the restored
wording overstates what the checks establish. In particular, the restored
all-epithelial MS4A12 endpoint of -0.996 conflicts with the saved minimum
-0.999611 (rounded -1.000 in the retained text). Depth matching and lower-abundance
comparison genes do not exclude every measurement artifact; bulk depletion does
not identify a within-cell mechanism; repository chronology cannot establish
what evidence was known outside the repository. The accurate, qualified account
and the sensitivity results remain available in the manuscript.

## Validation and scope

- 86 existing tests passed: 30 paper/precision, 35 external-control and 21
  retained-control audit tests. No experiment implementation or saved result was changed.
- Manuscript builds with the unchanged official NeurIPS 2026 style; no LaTeX
  errors, unresolved references, overfull boxes or detected identifying strings.
- 9 main-text pages, 2 reference pages and 3 appendix pages: 14 pages total.
- Rebuilt anonymous source ZIP; verified it compiles independently and matches
  all manuscript page text. All final PDF pages were visually inspected.

The remaining high-value empirical gaps are unchanged: an independent external
case, testing the actual percentile interval over bootstrap budgets at the
reporting boundary, and validation on patient-level clinical simulations.
The present revision improves the argument and its limits without claiming to
have filled those gaps.
