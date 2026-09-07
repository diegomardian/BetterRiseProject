# Pre-registration — D2: bulk GUCA2A and survival after purity adjustment

**Status:** locked 2026-09-07, after the outcome-blind D2 feasibility `PASS`
and before any D2 survival model, outcome table, or survival-derived result was
fitted or inspected under this contract.

## 1. Question and boundary

D2 asks a prognosis question, not a mechanism question: among one
primary-tumour RNA sample per TCGA participant, is bulk **GUCA2A** associated
with time to outcome after adjustment for tumour purity and the locked clinical
covariates?

The D2 feasibility gate in
[`d2_feasibility_gate.md`](d2_feasibility_gate.md) is outcome-blind. Its
`PASS` licenses this pre-specification only; it does not license a model fit.
This document is separate from the Stage 4 analysis, which remains
`not_prespecified` for D2.

The conditional association is the estimand. Bulk GUCA2A can reflect tumour
purity and normal-mucosa content as well as biology within malignant cells.
Consequently, no D2 result can establish cell-intrinsic GUCA2A silencing,
selective loss of a GUCA2A-high population, causality, or treatment benefit.

## 2. Fixed inputs and analysis unit

- **Expression:** the committed TCGA log2-CPM matrix used by the D2 feasibility
  gate. `GUCA2A` must map one-to-one through
  `config/gene_index/gene_index_1.0.0.map.tsv`; an absent or ambiguous mapping
  stops the analysis.
- **Unit:** one primary-tumour RNA sample per TCGA participant, selected by the
  committed sample manifest. Aliquots and multiple samples never increase
  participant count.
- **Outcomes:** TCGA-CDR progression-free interval (PFI), disease-specific
  survival (DSS), and overall survival (OS), using the project’s committed
  curation and endpoint definitions. No alternate endpoint is introduced.
- **Purity:** ABSOLUTE copy-number-derived purity is the primary adjustment;
  `estimate_affy_extrapolated` is a separately reported sensitivity analysis.
  The two sources are never coalesced. This uses the committed purity table and
  the existing purity reader; purity is required because the predictor is bulk
  expression.

The analysis population for each endpoint is its complete-case intersection of
the above inputs and the fixed covariates below. The result must report the
participant and event count after every exclusion, including the reason for
each missingness exclusion.

## 3. Endpoints and directional prediction

PFI and DSS are primary endpoints under invariant 9; **PFI is the lead
primary endpoint**. OS is secondary and may not be used for a headline claim,
because COAD overall survival includes non-cancer death.

Before any model is fitted, the directional prediction is:

> Higher bulk GUCA2A is associated with lower hazard, conditional on purity and
> the fixed covariates (hazard ratio per +1 log2-CPM is less than 1).

DSS remains primary because it is less contaminated by non-cancer death, but
its smaller event count makes it a pre-specified consistency endpoint rather
than an independent headline result. A null or imprecise DSS interval is not
relabelled as evidence of equivalence.

## 4. Fixed model specification

For each endpoint, fit a Cox proportional-hazards model with Efron ties and
complete-case analysis. GUCA2A enters continuously, per one log2-CPM, with no
median split, optimal cut point, subgroup search, transformation search, or
outcome-derived threshold.

All models include the locked covariates from
[`config/covariate_set.yaml`](../config/covariate_set.yaml): AJCC stage, age
(linear), sex, MSI status, and ABSOLUTE purity. PFI and OS also include site;
DSS omits site exactly as the locked `expression_models` endpoint override
requires. `conflicting` MSI remains missing rather than becoming a third level.

`project` and `plate` are pre-specified strata. Project is already a model
stratum; plate is required for an expression predictor and is treated as a
stratum rather than 29 fixed effects or a post-hoc random-effect choice. The
reported model table includes the number of participants, events, and nonempty
event-contributing strata.

Two estimates are reported for transparency:

1. a descriptive clinical-adjusted model that omits purity; and
2. the primary purity-adjusted model above.

Only the second is the D2 test. The first can show the scale of compositional
confounding but cannot upgrade or rescue an adjusted result.

## 5. Estimability and diagnostics

For the lead PFI model, the analysis requires at least 10 events per
non-stratum regression degree of freedom after the GUCA2A term is included. If
that condition fails, D2 is **NOT ESTIMABLE** for its lead endpoint and no
directional survival conclusion is made. DSS is reported with its exact event
count and events-per-degree-of-freedom because its reduced covariate set was
pre-specified precisely for its lower event count; it cannot independently
establish a result.

Schoenfeld-residual proportional-hazards diagnostics are reported for every
term. A violation affecting GUCA2A or purity means the corresponding endpoint
is labelled **PH VIOLATED** and supplies no Cox-based directional conclusion.
There will be no data-driven time split, alternate transformation, or model
replacement to recover a result.

## 6. Decision rules and falsifiers

The lead PFI result supports the pre-specified association only when the
primary ABSOLUTE-purity-adjusted GUCA2A hazard ratio is below 1 and its two-sided
95% confidence interval excludes 1, with the estimability and PH conditions
met.

- If the PFI interval includes 1, D2 reports **no supported association**.
- If the PFI hazard ratio is at least 1, D2 reports a direction contrary to the
  prediction, not a rescued alternative story.
- If the DSS point estimate is in the opposite direction to PFI, D2 reports
  **endpoint discordance** and makes no general survival claim. An imprecise
  DSS interval alone is reported as such, not treated as a contradiction.
- OS is descriptive secondary evidence only; it cannot reverse or establish the
  primary conclusion.

Repeat the full endpoint set with `estimate_affy_extrapolated` purity, keeping
every other decision fixed. If the adjusted PFI GUCA2A direction reverses or
the 95% confidence-interval decision changes (excludes versus includes 1), the
result is labelled **PURITY-SOURCE SENSITIVE** and is not described as a robust
purity-adjusted association.

No p-value, subgroup, stage-specific model, interaction, or survival-derived
cut point outside these rules changes the D2 verdict. Such work, if ever
desired, is exploratory and separately labelled.

## 6a. Lock conditions and stated limitations

Added on review before the lock in §7. The first is an artifact requirement;
the second fixes the multiplicity strategy; the third is a limitation reported
rather than hidden or re-tuned after fitting.

**Purity missingness is not missing-at-random, and complete-case drops the
wrong participants.** `config/covariate_set.yaml` records ABSOLUTE purity
coverage at **0.894**; combined with MSI (0.959) and stage (0.968) against the
gate's 624 participants, complete-case lands near **518**. ABSOLUTE fails
preferentially on **low-purity** tumours — and purity is the confounder this
design exists to adjust for, so the 10.6% removed are the participants where
compositional confounding is largest. §6's `estimate_affy_extrapolated` re-run
is framed as a purity-**source** sensitivity and is not the same thing.
**Required addition to §7's artifact:** report the GUCA2A distribution and
participant count in the dropped set beside the retained set, per endpoint and
per purity source. That is outcome-blind and it is the only way the bias
becomes visible.

**Multiplicity is handled by structure, stated here in advance.**
Three endpoints and two purity sources is up to six intervals, and a reader will
count them. The hierarchy in §3 — PFI lead, DSS a pre-specified consistency
endpoint, OS descriptive only — together with §6 treating the purity re-run as a
labelled sensitivity rather than a second test, is a legitimate alternative to
an alpha correction. **It is only legitimate if it is stated as the strategy in
advance**, which this sentence does. No correction is applied and none is owed,
because exactly one interval — lead PFI, ABSOLUTE-adjusted — can support the
claim.

**The estimability rule is a rule of thumb with no calibration behind it.** "At
least 10 events per non-stratum regression degree of freedom" is a convention,
not a property of this design measured at this n. `docs/HANDOFF.md` §3a records
an interval in this repository that was 0.82× the width it claimed by a closed
form nobody had computed, and `docs/prereg_meta_weight_calibration.md` records a
75% heterogeneity ceiling that turned out to be too strict at one k and too lax
at another. **Ten-events-per-df is the same class of number.** It is kept
because it was fixed before execution and because replacing it now would be
retuning; but it is recorded here as uncalibrated, and no power figure is quoted
anywhere in this document — a power number in this repository must carry the
false-positive rate of its own method (`interval_calibration.check_power_carries_
its_own_calibration`), and none has been computed for a Cox model at this n.

**One thing that is a choice rather than a gap.** The repository's own bulk
evidence points at CDX2 rather than GUCA2A: the pre-registered CIMP screen
(`results/2026-09-05_9203809/`) found GUCA2A falling **less** than CDX2,
+0.544 [+0.219, +0.878], and `docs/NEXT_AVENUES.md` §C records that C2's
instinct to target CDX2 first is corroborated more strongly than its proposal
claimed. D2 is defined as the GUCA2A question and stays that way. A reader
should know the stronger bulk signal is in a different gene.

## 6b. The order the labels resolve in

§6 lists its labels as independent rules and a run must emit one verdict, so
the order is fixed here — **before any model has been fitted**, which is the
only time it can be fixed honestly.

1. **NOT ESTIMABLE** — the §5 events-per-degree-of-freedom condition on the
   lead endpoint. Nothing downstream is read.
2. **PH VIOLATED** — a Schoenfeld violation for GUCA2A or purity *in the lead
   endpoint's primary model*. §5 labels "the corresponding endpoint", so a
   violation in DSS or OS does not silence PFI; OS in particular is descriptive
   under §6 and cannot reverse the primary conclusion.
3. **PURITY-SOURCE SENSITIVE** — direction reversed or the interval decision
   changed between ABSOLUTE and ESTIMATE.
4. **ENDPOINT DISCORDANCE** — DSS and PFI point estimates in opposite
   directions. The verdict row also records whether the DSS interval includes
   1, because §6 distinguishes a contradiction from an imprecise DSS interval
   and a reader cannot tell them apart from the label alone.
5. Then the lead result: **SUPPORTED ASSOCIATION** (HR < 1, interval excludes
   1), else **NO SUPPORTED ASSOCIATION** (interval includes 1), else
   **DIRECTION CONTRARY** (interval excludes 1 and lies above it).

**The last two are ordered the way §6 writes them, and that ordering is the
point.** §6 reads "if the PFI interval includes 1, no supported association"
*before* "if the hazard ratio is at least 1, direction contrary" — so a null
interval whose point estimate sits a little above 1 is a null result, not a
contrary finding. Testing the sign first would turn an absence of evidence into
a reported direction, which is the one error this document cannot afford to
make in that direction.

## 7. Required artifact and lock

Before execution, this proposed document must be locked as its own D2 change,
with implementation tests that assert the endpoint roles, continuous predictor,
primary purity source, sensitivity purity source, covariate set, and strata.
The execution artifact must contain the coefficient and interval for every
endpoint/model, diagnostic status, complete-case attrition, event counts,
events per degree of freedom, and the decision labels above. It must not report
a Kaplan–Meier plot based on a data-chosen GUCA2A threshold.
