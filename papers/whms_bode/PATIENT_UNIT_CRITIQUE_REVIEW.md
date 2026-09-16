# Independent assessment: inference unit, clarity and presentation

Starting paper: `ce272ca` on `wmhs/build-fixes`. Dimension scores were not used.
The review is valuable on presentation and on asking us to actually score the
patient-t procedure against fixed targets. Its proposed replacement thesis is
not supported by the experiment it quotes or by the new re-scoring.

## Central recommendation: test it, but do not adopt its conclusion

The review calls 2.77–6.29% null rejection with 4.21–8.02% detection a calibrated
patient-level repair and attributes the original failure to resampling cells
instead of patients. There are four different questions here:

1. Is patient sampling required to support inference over a population of patients?
   Yes. That is the intended scope of project invariant 5. The paper now states
   this directly. Cell resampling cannot provide independent patient replication.
2. Does that make cell resampling invalid for a conditional finite-cell-pool
   target under a generator that samples cells independently? No. Its inferential
   target and sampling law are explicitly conditional. The audited interval can
   fail that narrower task too. Changing inferential units does not explain away
   its measured conditional coverage failure.
3. Does the incoming patient-t experiment use all ten SMC or six KUL3 patients per
   simulated study? No. The quoted rows use **two held-out patients**, with some
   source pools containing mature cells from only one. The extension uses five.
   These results cannot establish the power of a study analyzing the full cohorts.
4. Does near-5% null rejection prove that patient-t is calibrated and solves the
   reporting rule? No. Coverage, target alignment, detection and availability still
   require evaluation. Even the quoted range includes conservative rejection and
   a KUL3/pooled point estimate above nominal. The source labels for the quoted
   two-patient rows are `passes_by_width` (SMC) and `passes_by_abstention` (KUL3),
   not uniformly `calibrated_but_underpowered`. Some five-patient rows have the
   latter label; others fail the null criterion.

## New evidence: exact same-draw re-scoring

The design and implementation were committed as `102801f` before results were
produced. The experiment replays the pinned incoming generator and unchanged
patient-t routine at 50 and 800 tumor mature cells, two and five held-out patients,
both effects and four cohort/pool conditions. Five original streams of 400 studies
give **32 settings and 64,000 studies**. All original availability, rejection and
sampled-reference inclusion counts reproduce exactly.

Before drawing either arm, we compute:

- The original fixed **cell-weighted** source-pool contrast.
- The fixed **equal-patient-weighted** contrast among selected source patients
  with defined mature-cell means.

Both are conditional finite-pool targets. Neither estimates uncertainty over a
new patient population. The interval equally weights the patients represented in
both simulated arms; missing representation can change that weighting or cause
abstention. We record availability rather than substituting zero. Changing the
scoring target cannot change rejection or detection for the same intervals.

At **50 cells, two held-out patients**, alternative coverage against the fixed
cell-weighted target is **96.03–97.84%** (equal-patient target: **96.10–97.63%**).
Thus the previously reported good coverage/low detection tradeoff survives the
scoring correction. Detection remains **4.21–8.02%** among returned intervals,
with **0.35–34.55%** alternative abstention. This positive clarification belongs
in the main text: it is no longer deferred on grounds of a missing fixed-target
check.

But **five held-out patients at 50 cells** give KUL3 null rejection of **18.1%**
(pooled) and **8.7%** (reference-only). Both pointwise Wilson lower limits exceed
5%. The proposition that patient-level resampling restores calibration generally
is therefore false for the tested method/design. At 800 cells, detection in
five-patient settings is **7.30–33.95%**, still far below the historical 80%
requirement. These findings bound this procedure's behavior, not the achievable
power of every patient-level method or of the full cohorts.

The manuscript promotes this result as a target/coverage/detection comparison,
not as proof that all failure was caused by the resampling unit. Aggregate exports
include both fixed coverages, original-reference inclusion, Wilson intervals,
availability, unconditional rejection/detection, source-patient counts, widths,
seed rates and provenance. Versioned parquet files are also saved under
`results/2026-09-16_102801f/`. No individual cells are exported.

## Other substantive criticisms

**Novelty:** substantially fair as a limitation, but not a reason to invent a
new theoretical claim. The contribution remains a documented failure, its exact
variance explanation, operating-characteristic comparisons and attempted repairs.
The broad phrase “do not use a reference sharing randomness” is itself too strong:
a realized-cohort estimand necessarily shares data with its estimator. The issue
is whether the reference represents the claimed target. Congeniality and inverse
crime are relevant context, not exact synonyms for every shared-reference failure.

**Venue fit:** still a real limitation. The paper does not evaluate a clinical
world model. We do not relabel the generator or imply that adding references
solves this. The clinical bridge is shorter and explicitly about target/population
alignment and downstream decision-rule validation. A venue change is a separate
author decision; this revision retains the requested workshop format.

**Generalization:** the narrow empirical scope is correct. One gene, labeling
axis, weighting and two source cohorts do not establish prevalence or patient
transfer. The assertion that only the algebra can transfer is too categorical:
the reproducible evaluation procedure can be reused, but the numerical failure
rates cannot be transported without new evidence. Repeated sampling improves
conditional numerical precision, not independent biological support.

## Presentation changes that improve the paper

- Added an analytic mechanism schematic: the same draw and same reported interval
  are assessed using error `f_N(e_T-e_N)` or `f_N(e_T-s e_N)`. At a halving, only
  the reference-arm variance contribution shrinks to one quarter. A schematic
  avoids presenting an arbitrary illustrative error density or one interval width
  as though it represented the empirical distribution of varying intervals.
- Expanded the null/detection plot to all four cohort/pool conditions. The earlier
  pooled-only plot was chosen to discuss the trough, but complete presentation is
  more balanced and makes the reference-only trajectory visible.
- Replaced the 20-row bare-number repair table in the manuscript with a three-panel
  comparison displaying pointwise 95% Wilson intervals for both counts and all
  methods. The numerical table remains a reproducible local artifact; exact
  aggregates accompany the source bundle. The figure distinguishes nominal 95%
  coverage from the historical 90% minimum.
- Shortened the external-control example to one main-text paragraph and retained
  its target definitions, controls and quantitative comparison in the appendix.
- Consolidated terminology: “warning range,” “proposed cutoff,” and an explicit
  definition of the evaluated grid. The historical grid discussion is shorter;
  shared-draw and Monte Carlo limitations remain. Internal documentary labels are
  introduced as validation/design records rather than unexplained “artifacts.”
- Added a sentence at the new repair comparison explaining why its 50-cell numbers
  differ slightly from the earlier precision table: fresh draws, separate runs.
- Expressed the historical grid's maximum MC SE as **3.5 percentage points**.
- Corrected “median signal-to-standard-deviation ratios” to **plug-in ratios using
  median pool moments**. The old arithmetic is unchanged; the statistic's label
  now describes what was computed.
- Restored Oaxaca and Blinder alongside Kitagawa, and scCODA/Milo/propeller in the
  composition context. Removed the unused DA-seq bibliography entry rather than
  inserting a decorative citation. The biological premise is not claimed to be
  validated by any of these methods.
- Renamed the local historical-figure generator to `make_historical_figures.py`
  and removed it from the submission regeneration command list. The figure is
  retained as a documented historical artifact, not a submitted figure.

The new material replaces repetitive prose and auxiliary grid material. The
patient-resampling distinction is stated where the result appears, with broader
limits consolidated in the scope and responsible-use sections.

## Citation checks and source reading

The decomposition context was checked against the original
[Oaxaca paper](https://inequality.stanford.edu/sites/default/files/media/_media/pdf/Classic_Media/Oaxaca_1973_Discrimination%20and%20Prejudice.pdf)
and [Blinder's publication record](https://www.princeton.edu/~blinder/articles.htm).
The restored biological references describe compositional/differential-abundance
methods: [scCODA](https://www.nature.com/articles/s41467-021-27150-6),
[Milo](https://www.nature.com/articles/s41587-021-01033-z), and
[propeller](https://pubmed.ncbi.nlm.nih.gov/36005887/).
They contextualize the biological question; they are not evidence for the audited
interval's operating characteristics.

## Validation

The new target tests check unequal cell-versus-patient weighting, coincident
null targets, exact replay of incoming counts, seed aggregation, availability,
reported rates and versioned-result provenance. Existing numerical checks are
retained. Prose assertions were updated where descriptions were shortened or
renamed; the underlying rate/denominator assertions were not weakened.
Build, visual inspection, anonymous source packaging and independent compilation
are recorded in `build_verification.json`. No dimension scores informed edits.
