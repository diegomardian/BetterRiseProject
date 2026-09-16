# Revision notes

## Scope

The revised paper is an independent copy under `papers/whms_bode`. Its starting
point is commit `b503971` and the attached 16-page PDF. The original manuscript,
code, preregistrations, configuration, and saved experimental results were not
changed. This revision corrects and reorganizes the interpretation of existing
evidence; it does not present newly run simulations or improved clinical
performance.

## Scientific changes

### 1. Center the demonstrated calibration problem

The new title, abstract, introduction and conclusion lead with the failure of
the original decision-rule justification. Section 2 now gives the controlled
union-grid experiment before the more limited historical sweeps. Those
historical results are summarized in Appendix B, including their original
replicate budgets and the reasons they cannot isolate a resolution effect.

The main figure is regenerated from the controlled **per-count** results on
two cohorts and two pools. It displays eight-seed medians and min–max bands,
with the 50-cell threshold and the two operating targets marked explicitly.
The bands are conditional simulation variability, not population confidence
intervals. Grid views are deduplicated after verifying that every shared
count has identical rates and denominators.

The text now reports behavior at the actual proposed threshold. At 50 cells,
SMC/reference discrimination is 75.5–81.0%, with both criteria met in only two
of eight seeds. KUL3/reference discrimination is 45.0–55.0%. This directly
supports the critique of a common 50-cell rule without implying that every
SMC configuration fails at that count or that 70 is a biological constant.

### 2. Make the draw-dependent target explicit

Code inspection found an important qualification absent from the old paper.
`src/harness/pseudobulk.py` computes the normal mean from the **drawn** mature
normal cells and passes it to `parametric_truth`; `src/harness/truth.py` then
uses `f_N * (s - 1) * mean_normal` for the normal-weighted intrinsic reference.
Thus the stored target depends on the sampled baseline. It is not a fixed
population parameter based on the mean of the eligible source pool.

New Equation 2 states this target. The methods, calibration section,
conclusion and limitations explain that the saved “coverage” rates measure
interval inclusion of this draw-dependent benchmark under the existing
specification. The bootstrap resamples both arms while fixing fractions.
These rates must not be presented as ordinary coverage of a fixed population
parameter or uncertainty across a population of patients.

The original calculation is preserved rather than silently replacing its
target and treating old results as if they came from a different experiment.
A fixed-pool or patient-population calibration would require a separately
specified experiment. The current evidence still challenges the recorded
justification and reports discrimination at the specified shift, but does not
supply a validated replacement threshold.

### 3. Separate equivalence, target accuracy and decision-rule performance

The equality audit is now a check on the scope of the comparison. Equality
with a reference can hold for a valid estimator; nonzero residuals do not
certify validity. The paper explicitly says that equality does not cause the
calibration failure. The clean control is retained because it demonstrates
this distinction concretely.

The trial examples retain assignment balance, OLS weighting, cross-fitting,
censoring and RMST. They establish possible equalities in stylised simulations,
not prevalence in deployed virtual arms. The full residual matrix is now in the main text; the residual-ratio
definition and implementation cross-check details remain in the appendix.

The practical recommendation is to declare the target, reference functional,
observation process, tolerance and valid-pair denominator, then evaluate target
error, interval behavior and abstention separately. The manuscript recommends
independent evaluation of a selected rule, but does not claim that such
validation was performed in the existing experiments.

### 4. Correct the learned-generator section

Several corrections go beyond phrasing:

| Item | Earlier paper | Revised, checked against saved tables |
|---|---|---|
| Expanded design | 52,920 estimates despite ten estimators | 75,600 scheduled evaluations on 7,560 cohorts |
| Generator failures | Not enumerated in the text | 75,580 recorded evaluations on 7,558 cohorts; the two missing cohorts are in the parametric-plugin arm at n=100 |
| Finite estimate–draw-reference pairs | Not enumerated | 74,726, with 854 excluded |
| Exclusion categories | Not reconciled | 410 nonfinite references, 731 nonfinite estimates, overlapping in 287 evaluations |
| Mixture-family ridge maximum at alpha=1e-5 | 8.6e-5 | 7.9e-5, rounded from the actual mixture-family maximum |
| Displayed learned estimators | Selected rows, two mixture outcomes merged | All ten estimators, including the MLP outcome model and all three ridge penalties |
| Tolerance transitions from n=100 to n=5000 | Fourteen, based on an approximate constant 3e-5 threshold | Thirteen using the saved actual Boolean comparison results |
| Mixture model-reference MC standard errors | 0.020–0.027 | 0.019–0.027 after rounding the saved extrema |

The regenerated table uses maxima within each stated generator family, not
values copied from a maximum over all families. It includes the MLP-outcome
maximum of 3.3 on mixture-generated data, which was absent from the old table.

The text no longer describes rho=0 at exact equality as a failure of a proposed
remedy: that is its expected value when its denominator is positive. Nor does
it infer that shrinking ridge residuals imply worsening validity as sample
size increases. The n^-1.75 description was a fit to maxima over a finite grid,
not a derived asymptotic rate, and is removed from the claims.

The mixture mean identity is explained using its M-step equations. The fitted
generator's population effect is treated as a different target; a common
reference across estimators is not a defect. The decomposition into estimator
sampling departure and generator target discrepancy replaces the claim that
this comparison measures only generator bias. Numerical Monte Carlo
uncertainty is tied to the finite simulation budget, not treated as an
irreducible precision limit. Parametric and MLP references are distinguished
from Monte Carlo mixture references.

### 5. Bound the external audit

The seven recorded NO verdicts are retained. “Already the norm” and claims to
bound prevalence are removed. The paper states that the sample is nonrandom,
has one reader, includes only one executed repository, and did not retrieve
examinable code in the query specific to virtual control arms.

The noiseless case is explained consistently with the trial-design results:
distinct functions may agree on a restricted data-generating support. Adding
noise changes that support; the counterexample shows they are not identical
on general noisy data, rather than disproving agreement under the shipped
design. The original classification is not retroactively changed.

### 6. Preserve negative biological results and clarify raw diagnostics

The failed biological controls, annotation dependence, patient-level limits,
non-identifiability and uncertainty of mechanism remain explicit. The new
paper does not claim that a biological mechanism or virtual patient model has
been validated.

The recovery sweep includes raw zero-cell arithmetic before the reporting
gate, which substitutes an empty-arm mean with zero in the historical harness.
Appendix B now explains why the resulting finite ratio is a computational
diagnostic, not an intrinsic effect in an absent cell population. Reporting
and interval construction still abstain in those cases.

The original dataset publication is now cited: Lee et al., Nature Genetics
52, 594–603 (2020), DOI 10.1038/s41588-020-0636-z. The text describes the
analyzed paired subsets rather than implying the entire SMC accession has
ten patients. Citation details and accession association were checked against
GEO and the publisher.

## Presentation and reproducibility

- The official NeurIPS 2026 package, downloaded from the conference-linked
  author kit, replaces the 2024 compatibility shim. It is not modified.
- The condensed main text, including responsible use and the two promoted
  tables, fits on nine pages; references start on page 10. No font or margin
  compression was used.
- The visible malformed bold command is removed by rewriting its paragraph.
- Two figures are regenerated from the same pinned inputs: controlled
  calibration and trial recovery. Redundant recovery-sweep and historical-grid
  plots were removed in condensation; their numerical findings and limitations
  are summarized in the appendix.
- The incomplete-looking frame around the audit formula is replaced with
  ordinary displayed mathematics and code.
- A dedicated build gate checks freshness, page count, cross-references,
  overfull boxes, basic anonymity and the style checksum.
- Dedicated result tests check the revised claims rather than relying only
  on tests that reference the old source tree.

## Evidence that remains limited

The revision improves accuracy and the argument; it does not create evidence
that is missing. No fixed-population coverage experiment, held-out validation
of a newly selected threshold, representative external prevalence audit, or
validation on a deployed learned patient model is claimed. Patient-level
uncertainty remains unquantified in the threshold study. These are boundaries
of the results, not problems that can be fixed by stronger wording.

## Condensation to the nine-page main-text limit

The first independent revision had eight main-text pages, two reference pages
and eight appendix pages (18 total). The condensed version has nine main-text
pages, two reference pages and three appendix pages (14 total). The original
18-page portable sources are preserved in `archive/expanded_18page_source.zip`.

- Promoted the reference-by-design comparison and the complete ten-estimator
  learned-generator table into Section 4, where readers need them to assess
  the claims. Added interpretation of the learned-table exceptions.
- Removed a second trial-results table that duplicated the retained figure,
  recovery-sweep and historical-grid plots, and the tangential synthetic
  benchmark. Historical calibration and recovery accounting remain in prose.
- Consolidated supporting checks into four appendix sections: simulation
  settings; calibration details and limits; trial and generator checks; and
  external audit and research context. Compressed repeated explanations and
  the LLM-use disclosure.
- Retained fixed simulation settings, valid-pair denominators, overlapping
  exclusion counts, empty-cell behavior, tolerance definitions, reference
  dependence, biological limitations and the external audit's scope limits.
  No new experiments or stronger evidential claims were introduced.

This reduction uses content selection and organization, with the official
style, font sizes and margins unchanged. The final PDF and portable source
bundle were rebuilt and checked after condensation.
