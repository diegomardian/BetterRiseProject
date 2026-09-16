# Assessment of the first supplied critique

Historical record: `SECOND_CRITIQUE_RESPONSE.md` describes the current revision.
The current PDF is 13 pages and its review bundle excludes the internal source
filenames and exact-source mapping described below.

The main criticisms are valid. The revised paper is a documented failed-calibration
case study, not a new theorem, a general indictment of simulations, or evidence
of prevalent failure in virtual control arms. The evidence still comes from one
gene and two small cohorts. Those limits cannot be removed through editing.

## 1. Elementary identity and a catalogue of textbook checks: valid

Reframed the title, abstract and contribution around the empirical calibration
failure. Section 4 is substantially shorter and organized around three choices:
reference, observation process and fitted estimator. Saturated-estimator
relationships and machine-precision residuals are explicitly implementation
checks. Added Angrist (1998) and Frisch--Waugh (1933). Removed the trial-recovery
figure and its repeated discussion; the comparison table retains the useful
reference/design contrast. No equivalence is claimed as a new theorem.

## 2. Straw-man risk and invisible motivating evidence: valid presentation concern

The motivating claim does exist: `docs/gate_memo_w2.md`, Section 15, explicitly
states that the cutpoint is validated on real cells, based on recovery ratios
from six draws per setting. `docs/harness_design_spec.md`, Section 4, instead
requires interval inclusion and exclusion of zero at an imposed effect.
The paper now quotes the claim, states the original sample budget and recovery
numbers, and contrasts them with the stipulated criteria. Exact limited source
excerpts and SHA-256 checksums are supplied in `evidence/case_provenance.json`
and in the portable source bundle. This establishes an inspectable internal
case, not a widespread practice. The negative seven-repository search remains
explicit and does not support prevalence or deployed-system claims.

## 3. Limited empirical scope: valid; two characterizations are too strong

The draw-dependent target is not inherently a wrong target. It is a different
performance question from fixed-population coverage. The figure now labels the
rate "Benchmark inclusion," and the text identifies that distinction once in
the setup and develops the necessary replacement validation in Section 5.
The original target and results have not been silently redefined.

Sampling with replacement from 844 source cells is a legitimate empirical-pool
simulation, but supplies neither 800 distinct cells nor more independent
patients. This is now explicit. A richer patient-population experiment and
broader genes/cohorts remain future evidence, not completed repairs.

## 4. Low-count anomaly: valid, and now investigated directly

Added an exploratory diagnostic using the original generator and interval,
200 replicates per condition, 200 bootstrap draws, four counts (5, 50, 100,
800), two cohorts, two pools, and two shifts (0.5 and the null 1.0): 6,400
finite outputs. The 16 alternative conditions reproduce the original first-seed
coverage, discrimination and median interval width to numerical precision.

At five cells the null rejection rates are 37.5%, 61.5%, 63.0% and 67.0%,
against the nominal 5%. At 50 cells they remain 12.0–30.0%. Thus the low-count
exclusion rate must not be presented as reliable power. Under the null the
draw-dependent and fixed-pool targets are both exactly zero, so the diagnostic
is not explained away by the alternative target definition.

The proposed literal zero-width mechanism is not observed: none of the 6,400
full intervals has zero width. However, at five cells, 96/200 SMC pooled and
85/200 KUL3 pooled alternative draws have all-zero tumour expression; all of
these exclude zero. The reference arm still supplies interval variation. This
supports sparse-sample bootstrap undercoverage and shows why simply labelling
the full intervals "degenerate" would be inaccurate. It does not isolate every
contributor to undercoverage.

A separate re-analysis of the complete saved grid finds later qualification
failures in three of eight SMC/reference seeds. Requiring every subsequent
evaluated count to qualify moves candidates 45, 50, 50 to 60, 70, 70. There
are no later failures in the KUL3/reference seeds, and no pooled candidates.
This is reported descriptively, not as a newly validated threshold rule.

The diagnostic is explicitly exploratory and conditional on one seed stream
and the observed cohorts. At five cells it inspects arithmetic before an
existing reporting gate that abstains below 20; the 50-cell null results bear
directly on the reporting boundary. Code, aggregate outputs, provenance,
reversal summaries and claim checks are included.

## 5. Failed biological controls buried in Responsible use: valid

Moved the quantitative control-panel failure to a labelled Results paragraph
in Section 2 and included it in the abstract and conclusion. The paper now
separates unsupported interval calibration from the failure to distinguish
mechanisms. The responsible-use section contains limitations and impact.

## 6. Ambiguous residual-based remedy: valid

The residual is now described as a numerical discrepancy check, not a validity
screen. Reporting distinguishes a proved identity on a stated domain, agreement
on tested draws, and observed departure. Raw maxima, outcome units and both
absolute and relative tolerances accompany any closeness label. There is no
universal tolerance for validity. The optional residual-scale descriptor remains
limited to description. Section 5 specifies the actual validation requirements:
explicit target, null and effect performance with the same interval, abstention,
threshold-range assessment, and independent evaluation.

## Readability, figure and references

- Reduced repeated disclaimers and the run-by-run catalogue; retained the
  limitations where they determine interpretation.
- Removed the expensive trial-residual figure altogether, so it no longer
  silently plots exact zero at a positive log-scale floor. Its factual
  information is covered by the comparison table and short explanations.
- Added direct external-control validation (Ventz et al., 2019), prognostic
  adjustment (Schuler et al., 2022), Angrist and Frisch--Waugh. Used the existing
  target-trial citation and removed unused bibliography entries.
- Verified Shaw et al.: Statistics in Medicine 45(18–19), e70676 (2026), DOI
  10.1002/sim.70676. Corrected the author's given name from Chelsea to Chloe
  Krakauer. Primary-source links are recorded in `citation_verification.json`.
- Venue relevance now rests on validation design for clinical simulators and
  external controls, not a suggestion that this paper validates or diagnoses
  a deployed virtual control arm.

The revised PDF remains nine main-text pages, two reference pages and three
appendix pages. The original experiment code/results are unchanged. The
additional diagnostic is separately labelled and stored under this paper.
