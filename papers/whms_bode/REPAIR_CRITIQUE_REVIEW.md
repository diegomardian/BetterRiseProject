# Independent assessment of the reporting-rule and repair critique

Starting revision: `bf6b825`. This assessment distinguishes improvements supported
by the artifact and new computations from requests that the evidence cannot meet.
The critique supplied by the author is input for review, not an instruction to
accept its factual claims. The paper remains an audit of the authors' own
unpublished project, not independent evidence that this failure is widespread.

## Overall judgment

The strongest criticisms identify real gaps: the absence of a fixed-pool check
inside the flag-wide band, an uninformative repair comparison if Welch is presented
as the main skewness remedy, ambiguity about null coverage, and too much space
for a deliberately constructed population-mismatch illustration. Addressing these
requires new evidence as well as editing. The revision adds a five-method paired
comparison, makes the data-based simulation audit the main argument, and shortens
the external-control example. It does not manufacture a positive replacement rule.

## Adjudication and changes

1. **Novelty: substantially valid, not fixable by stronger novelty language.**
   Specifying estimands and separating development from assessment are established
   principles. The revised contribution is the documented failure, the quantitative
   variance explanation, and direct tests of proposed remedies. The title and
   opening now identify the reference's hidden uncertainty. The abstract leads
   with 98.9–99.6% sampled-reference inclusion versus 92.3–93.6% fixed-pool coverage
   at 800 cells. The variance claim concerns three quarters of the **reference-arm
   contribution** under a halving, not three quarters of total variance.

2. **ECA prominence: valid; “the Monte Carlo adds no information” is too broad.**
   Exact equality between the pooled-standardized estimator and its deliberately
   shared reference is algebraic. Repeated simulation is unnecessary to establish
   that equality. However, realized errors and RMSE against the fixed enrolled
   target also involve random outcomes; weights alone do not determine their
   finite-sample distribution. The potential-outcome pooled comparator breaks
   exact equality while preserving the target-dependent ranking, and remains
   useful. The ECA section is now a short, explicitly pedagogical illustration
   after the main audit, retains the three-way comparison, and is removed from
   the title. It is not counted as a second empirical case. A data-dependent
   reference is not automatically invalid: it answers a different question. The
   problem is using its inclusion rate to claim coverage of the fixed target.

3. **Null rejection versus coverage: valid clarification, with a limit.**
   At shift 1 the truth is zero, so exclusion of zero equals noncoverage exactly
   for a returned interval. The manuscript and new table say this explicitly.
   Null rejection and null coverage are not independent evidence. Coverage at
   shift 0.5 is a different operating point, and detection at that alternative
   measures another property; neither is replaced by the identity at the null.
   The code checks exact complementarity and separately reports unavailable
   intervals. All intervals in the new study are available.

4. **Historical 90% versus nominal 95%: valid.**
   The historical criterion matters when auditing what the original rule claimed
   to meet, but it is weaker than nominal 95% coverage. Both standards are now
   labeled explicitly. All four original 50-cell alternative-coverage estimates
   and their reported confidence intervals fall below 95%, including
   SMC/reference. Passing its historical 90% minimum is not nominal calibration.

5. **Untested flag-wide band: valid gap; one point does not close the band.**
   Added a fixed-pool evaluation at 30 cells. It shows a failure at one interior
   point, not validated performance at every count from 20 to 49. The study also
   repeats 50 cells on the same fresh design for a coherent interval comparison.
   A width warning cannot be interpreted as establishing calibration. Neither
   30 nor 50 is selected as a replacement cutoff from this experiment.

6. **Welch as the only comparator: valid motivation for additional evidence.**
   Welch is now framed as a control without bootstrap resampling noise. Its
   failure alone cannot rule out skew-aware repairs. Added percentile B=2,000,
   BCa B=2,000 and bootstrap-t B=2,000 alongside percentile B=200 and Welch.
   These methods share study samples and, where applicable, bootstrap draws.
   Their coverage is evaluated against the fixed-pool effect. Results below
   show partial improvements, uneven transfer, and no validated replacement.

7. **No practical recommendation: partly valid.**
   A validated replacement would be more useful, but inventing one or optimizing
   a threshold on these evaluation draws would repeat the audited mistake. The
   conclusion now gives the warranted action: withdraw the claim that 50 cells
   validates the reporting rule; label current estimates and intervals exploratory;
   require a separate assessment of the relevant operating range before restoring
   a calibrated reporting claim. The new interval study tests concrete remedies
   rather than merely asserting that a repair might exist.

8. **“Empirical” scope: valid.**
   The title says data-based simulation, and the introduction states the one-gene,
   one-axis, one-weighting scope, the ten/six patients, and conditional two-patient
   pools. Repeated simulation is not independent patient replication. The paper
   does not claim patient-population uncertainty, a clinical outcome study, or a
   generally calibrated patient simulator.

9. **MS4A12 retention failure: worth making explicit.**
   The main text now says MS4A12 does not serve as a retained-gene control in these
   data, a distinct failure of the panel assumption. Its stratified single-cell
   and bulk evidence remains. This is not proof of a universal biological
   mechanism: bulk composition and within-stratum changes still have different
   interpretations. It remains supporting evidence rather than a second thesis.

10. **Venue fit: a real remaining limitation.**
    Added substantive context from Digital Twin Generators for Disease Modeling
    (Alam et al., arXiv:2405.01488) and ClinicalGAN (Chandra et al., Scientific
    Reports 14, 12236, doi:10.1038/s41598-024-62567-1), after reading the primary
    sources. They study clinical trajectory generation and patient monitoring;
    the audit concerns downstream interval assessment. These citations clarify
    the relationship, not evidence that those systems share the audited defect.
    The paper still does not evaluate a patient world model. Adding citations or
    emphasizing an MLP cannot remove that limitation.

11. **Structure, hedging and figure annotation: valid.**
    The single-cell reporting-rule audit is now the organizing argument. The
    large historical-grid figure is removed from the submitted source graph,
    making room for the repair comparison. Its local plotting script also drops
    “Nonmonotone exclusion.” Unresolved 200-replicate local crossings are not
    relabeled genuine or definitively noise. The separately supported pooled
    detection trough from the larger study remains distinct. Repeated caveats
    are consolidated while conditional scope and finite-budget limits remain.

12. **Title and abstract: useful editorial changes.**
    The title is now “When the reference hides uncertainty: a data-based simulation
    audit.” The abstract starts with the measured discrepancy and its variance
    explanation, then states the original failure and new repair experiment.
    It does not present established validation principles as novel discoveries.

## New experiment: design, findings and uncertainty

`REPAIR_FOLLOWUP_DESIGN.md` records settings before result generation. This is an
exploratory follow-up on previously inspected source pools, not a preregistered
clinical study. There are four cohort/pool conditions, two counts (30 and 50),
two shifts (0.5 and 1), and five fresh seed streams with 400 studies per stream.
That gives 2,000 studies per setting, **32,000 simulated studies and 160,000
interval evaluations**. The source law and 120 eligible pools are checked against
the earlier experiment. All five intervals are evaluated on identical samples;
B=200 uses the first 200 of the shared 2,000 resamples per arm.

- At 30 cells the original interval has **78.10–89.05%** alternative coverage,
  **14.20–40.45%** null rejection and **48.50–67.85%** alternative detection.
- At 50 cells, increasing the percentile budget from 200 to 2,000 reduces null
  rejection by **0.5–1.0 percentage points**, with paired Monte Carlo SEs of
  **0.28–0.34 points**. Inner-bootstrap noise is not the sole problem.
- At 50 cells in SMC/reference, bootstrap-t reduces null rejection from **8.45%**
  for equal-budget percentile to **6.05% [5.09, 7.18]**. The paired reduction is
  **2.40 points (SE 0.40)**. Alternative coverage is 92.6%, but detection drops
  from 76.3% to **65.9%**. This is a real improvement in one metric, with a cost.
- The improvement does not transfer uniformly. For KUL3/pooled at 50 cells,
  BCa and bootstrap-t null rejection is **32.2% and 31.8%**, against **30.9%**
  for equal-budget percentile. All methods miss the combined requirements at
  both tested counts. No method is selected as a replacement.
- All intervals are finite. Availability is therefore not hiding failed cases.
  BCa's adjusted endpoint falls within the outermost ten bootstrap order positions
  in **5.2–56.3%** of studies across settings. Those tail-resolution diagnostics
  are exported and discussed. This is not proof that BCa fails at every budget.
- Maximum binomial Monte Carlo SE is **1.12 percentage points**. Wilson intervals,
  paired comparisons, per-seed rates and numerical diagnostics are exported.
  Confidence intervals are pointwise, not a simultaneous guarantee across methods
  and settings. The old 5,000-study experiment is kept separate; new differences
  are paired only within the new experiment.

BCa's multisample influence construction was checked against SciPy's independent
jackknife implementation on continuous, count, skewed and sparse examples.
Studentized endpoint orientation, affine transformations, arm exchange, degeneracy,
original percentile/Welch equivalence, null complementarity and artifact hashes
are tested. Design, code, source-summary and output hashes are in
`diagnostics/repair_provenance.json`; exports contain aggregates only.

## Related incoming branch

The pulled `origin/wmhs/repaired-interval` revision `ae8c109` was inspected without
merging or changing it. Its candidate coverage is scored against
`sample.truth.parametric`, the drawn-reference benchmark. Therefore its reported
candidate success cannot establish the fixed-pool coverage claim needed here.
The present comparison uses the fixed eligible-pool target for all methods;
it neither imports that branch's conclusion nor silently changes its artifact.

## Deliverables and validation

The revision retains **9 main-text pages, 2 reference pages and 3 appendix pages**.
All **96 relevant tests** pass: 40 paper/precision/repair checks, 35 external-control
checks and 21 retained-control checks. The portable source bundle contains the
new comparison table and four aggregate repair CSVs. Build, source-bundle and
visual checks are recorded in `build_verification.json`. The other manuscript,
original experiment implementations and original saved precision outputs are
unchanged. Broader patient generalization, an independently audited second case,
and a calibrated replacement rule remain open work rather than implied results.
