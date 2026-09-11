# Pre-registration — the Chen adenoma decomposition by lesion subtype

**Status: WRITTEN AND UNLOCKED, 2026-09-10.** Nothing here licenses a run. The
gate conditions in §9 are unmet, the labels have not been joined to any
decomposition output, and no subgroup estimate of any kind has been computed or
inspected. Locking requires §9 discharged in writing and the amendment recorded.

This is `docs/NEXT_AVENUES.md` roadmap item 3. It supersedes the AD-versus-SSL
framing there, and it does **not** supersede
[`wes_subtype_plan.md`](wes_subtype_plan.md), which remains the route for the
genotype label and is now scoped to one term rather than both (§4).

---

## 1 · Question and boundary

Does the **already-computed** avenue-A decomposition
(`docs/HANDOFF.md` §6h, `lineage` rung, 44 patients) differ between polyp
lesions that an independent modality assigns to different classes?

It is a **within-study subgroup heterogeneity analysis**. It is not:

- a replication of avenue A — same cohort, same cells, same estimator;
- a claim about all conventional adenomas or all serrated lesions — the arms are
  the lesions Chen deposited and cBioPortal labels, nothing more;
- a mechanism result. A subtype difference is an association between a
  morphological or genotypic class and a decomposition term.

Invariant 4 forbids pooling this with any other cohort. Invariant 7 requires the
interaction term to be reported separately and never folded into either arm.

## 2 · Fixed inputs and the analysis unit

| input | fixed value |
|---|---|
| transcript side | the committed avenue-A `lineage`-rung Chen adenoma decomposition, 44 patients. **No re-estimation, no re-tuning, no rung sweep.** |
| label side | cBioPortal study `crc_hta11_htan_2021`, sample clinical attributes and mutation MAF, both open access |
| join key | exact equality of `scRNA_biospecimen_id` to a cBioPortal `sampleId`. **Participant-level joins are refused**, as in `wes_subtype_plan.md` |
| crosswalk | `results/2026-09-09_0bf9734/wes_subtype_provenance_crosswalk.parquet` |
| unit of inference | **the patient** (invariant 5). Lesions are not independent units and are never bootstrapped over |
| seed and sha | recorded in the result sidecar (invariant 10) |

The label snapshot is pinned: the cBioPortal responses are written to
`data/` with sha256 in `data/manifest.csv` before the first join, and the
analysis reads the snapshot, not the live API.

## 3 · Labels, their provenance, and the arms

Two label families, declared separately under invariant 11.

**Family P — histopathology.** `POLYP_TYPE` ∈ {`AD`, `SER`}.
`label_provenance`: pathologist diagnosis from tissue morphology.
`claim_provenance`: single-cell transcript decomposition. No overlap of
derivation: the label is not computed from the transcripts.

**Family G — somatic genotype.** Truncating `APC` (nonsense, frameshift, splice,
start/stop-loss) versus `BRAF` V600E, positive-only arms.
`label_provenance`: FFPE whole-exome somatic calls.
`claim_provenance`: as above.

**Arm sizes, fixed here before any outcome is read** (`VERIFIED`, against the
crosswalk and the public label tables; the derivation is
`docs/DATA_HUNT_2026-09-10.md` §1):

**Two lesion denominators, and they are not the same set.** *Labelled* lesions
carry the label, including those belonging to the conflicted patient excluded
by rule 1. *Analysable* lesions are the ones that actually enter the arm. The
patient column is always post-exclusion. Both lesion columns are stated because
an arm table that gives one of them beside a post-exclusion patient count reads
as though they share a denominator, and the next person to recompute an arm is
then wrong by one lesion without being able to tell which number was meant.

| family | arm | labelled lesions | analysable lesions | **patients** |
|---|---|---|---|---|
| P | `AD` | 16 | 15 | **13** |
| P | `SER` | 11 | 10 | **9** |
| P | excluded — conflicting lesions within patient | 2 | — | 1 |
| P | no label | 25 | — | 21 |
| G | truncating `APC`, `BRAF`-negative | 6 | 6 | 5 |
| G | `BRAF` V600E, `APC`-truncating-negative | 5 | 5 | 4 |
| G | neither callable | 2 | — | 2 |

Family G's two columns agree: no lesion carries both markers, so rule 1 excludes
nobody there, and its one multi-lesion patient is concordant under rule 2.

**Verified against the pinned snapshot**, gate 4, `results/2026-09-10_461d20b/`:
all six figures reproduce, as do the conflicted patient and the three
concordant multi-lesion patients.

**The `SER` arm is mostly not sessile serrated lesions**, and that is fixed here
rather than discovered afterwards. **State the denominator:** of the **11
labelled** `SER` lesions, 4 are sessile serrated lesions and 7 are hyperplastic
polyps (3 microvesicular, 2 goblet cell-rich, 2 unqualified). **The analytic set
is 10 lesions, not 11** — one SSL belongs to the conflicted patient excluded by
rule 1 — so the arm that is actually analysed is **7 hyperplastic polyps and 3
sessile serrated lesions, and those 3 lesions come from only 2 patients.** Hyperplastic polyps are not the serrated precursor the
AD-versus-SSL question is about. The arm is therefore named `SER`, never `SSL`,
throughout; an SSL-only arm is 3 analytic lesions over 2 patients and is **not
analysable** — it is reported as a count and nothing else.

Family P covers **23 of the 44** avenue-A patients. DIS/VAL splits 9/4 for `AD`
and 6/3 for `SER`. The 21 uncovered patients are reported as attrition with
their `POLYP_TYPE` values (`Unknown`, or no cBioPortal record), and the
decomposition values of the covered and uncovered sets are compared
**descriptively** — the D2 §6a lesson is that a dropped set can differ
systematically in the direction of the outcome, and it is cheaper to look than
to be asked later.

### Rules that must be fixed now because they are choices

1. **Conflicting patients are excluded, not assigned.** A patient with both an
   `AD` and a `SER` labelled lesion (`HTA11_6801`) enters neither arm and is
   counted in the attrition table. Assigning them by any rule — larger lesion,
   first sample, worse dysplasia — is a modelling choice presented as a
   measurement.
2. **Concordant multi-lesion patients contribute one value.** Three patients
   (`HTA11_6818`, `HTA11_8622`, `HTA11_866`) carry two agreeing lesions. The
   patient-level decomposition value is what avenue A already emits per patient;
   lesions are not averaged into it and not entered twice.
3. **`Unknown` is not a third arm.** It is missing, and it is reported as
   missing.
4. **Case is normalised before any field is read.** `Not Stated` and
   `Not stated` are one category.

### `ADVANCED` and `ATYPIA` are excluded as primary stratifiers

Both are **34% and 23% uninformative** over the 35 covered lesions, with
missingness that tracks the arm — `ADVANCED` is uninformative for 5 of 11 `SER`
against 1 of 16 `AD`. Conditioning on them drops serrated lesions
preferentially. They may appear only as a **declared sensitivity** with the
complete-case attrition reported by arm, never as a primary contrast, and never
as a compositional stratifier at all (§4).

## 4 · The two terms have different standing, and this is the design

**This section is the reason the pre-registration exists**, and it must be
settled before the run rather than argued in a discussion section.

Invariant 11 tests **provenance** — whether a label descends from the
transcripts whose programme is then claimed. Family P passes. What invariant 11
does not test is whether the label and the endpoint are **measurements of the
same tissue property**.

A pathologist separates `AD` from `SER` by reading crypt architecture and the
maturation of epithelium toward the luminal surface. The **compositional** term
is the change in mature-cell fraction. These are largely one property in two
modalities. A compositional difference between P arms is therefore at risk of
being definitional — lesions selected on architecture, then measured for
architecture. The **intrinsic** term is per-cell output within mature cells,
which no pathologist grades, and is not exposed to this.

| term | family P (morphology) | family G (genotype) |
|---|---|---|
| **intrinsic** | **primary test** | supporting, underpowered |
| **compositional** | **descriptive only** — confounded with the selection criterion | clean in kind, **descriptive only** — 5 versus 4 patients is below every floor here (`HANDOFF` §3a) |
| **interaction** | reported separately (invariant 7), descriptive | descriptive |

**So the pre-committed position is: there is no compositional test available on
any label this project can currently obtain.** Family G is the cleaner of the
two and is worth the Synapse certification it is blocked on, but it buys a
better *descriptive* contrast, not a test. Any later report that presents a
compositional arm difference as a finding is in breach of this section.

## 5 · The primary contrast, and the interval

**Primary.** The difference in the patient-level **intrinsic** term between
family P arms, `AD` minus `SER`, at the `lineage` rung, on the corrected
detection scale (`src/reference/detection_scale.py`).

**The interval is a Student-t interval on the difference, not the percentile
bootstrap.** `docs/HANDOFF.md` §3a: at n≈10 the project's percentile bootstrap
is 0.82–0.91× the width it claims and excludes zero ~7% of the time under a true
null. At 13 and 9 that defect is squarely in range, and a two-arm contrast
compounds it. The precedent is
`src/reference/jobs/mlh1_positive_control.py`, which already does this.

### Two things §5 left under-determined, fixed 2026-09-10 before the run

Recorded here rather than settled at the keyboard, because both would otherwise
be choices made with the table in view.

1. **The between-arm interval is Welch's t, not one-sample Student-t.** §5's
   requirement is a t-interval rather than the percentile bootstrap, and that
   stands. But `student_t_interval` is one-sample, and this is a two-arm
   difference at 13 against 9. The project's own between-arm method is Welch —
   `src/reference/interval_calibration.py:481`, *"Welch rather than Student
   because the arms differ in size and in cell count per patient, so their
   variances differ by construction"* — which is this situation exactly.
   Welch is used, and named in the result.

2. **All three weightings are reported and none is primary.** `normal`,
   `tumour` and `doubly_robust`, never folded together, as
   `prereg_adenoma_decomposition` §164 already requires. §7's branches are
   evaluated **on each**, and whether they agree is part of the report.
   Choosing one weighting now, with no pre-registered basis, would be selecting
   an answer; choosing one after the run would be worse.

**Genes.** The panel's target and control roles as frozen, scored as whole tiers
and **every pair reported**, not the targets against a housekeeping comparator
alone — the §6d correction, where reporting only `GUCA2A − X` left the claim
about identity markers with no row behind it.

**`None` is not `0.0`.** A patient whose arm has too few mature cells to ask is
`not_estimable` with a reason, in both arms, and the per-arm not-estimable count
is reported beside every estimate (invariant 1).

## 6 · Directional prediction, made before the run and taken from the literature

The prediction is not ours and predates the data, which is what makes this
confirmatory rather than a subgroup hunt. Bashir et al., *Hum Pathol* 2019
(PubMed 30716341): GUCA2A is lost in conventional adenomas, serrated adenomas
and MSI tumours alike, while GUCY2C is near-eliminated in serrated lesions
specifically, attributed to loss of CDX2.

> **Predicted: the GUCA2A intrinsic term does not separate `AD` from `SER`.
> CDX2 falls in the `SER` arm and not in the `AD` arm.**

This is a prediction about CDX2 behaving differently between arms while GUCA2A
does not — the opposite shape from a generic "serrated lesions are more
different" expectation, and therefore falsifiable in a way that matters.

It is also consistent with what is already committed: §6d found CDX2
indistinguishable from housekeeping in a Chen cohort whose polyps are majority
conventional adenoma, which is the `AD`-arm half of this prediction, already
observed. **That half is therefore not evidence for the prediction** and is
excluded from any claim of confirmation. Only the `SER` arm is new information.

## 7 · Falsifiers and pre-committed consequences

| outcome | consequence, committed now |
|---|---|
| CDX2 intrinsic difference excludes zero in the predicted direction, GUCA2A's contains zero | The published prediction replicates in a third modality pairing. Reported as a **tier-level, within-study heterogeneity result** in one cohort. Not a mechanism claim. |
| Neither contrast excludes zero | Item 3 returns **NOT SEPARABLE AT THIS RESOLUTION** and closes. 13 versus 9 patients is the reason to expect this; it is written here so it is not later read as a biological null. |
| GUCA2A's intrinsic difference excludes zero and CDX2's does not | Contradicts the published prediction. Reported as such, **without** a post-hoc mechanism, and the `SER` arm's composition is inspected for the obvious artefact (§3's hyperplastic majority). |
| The arms differ in the compositional term only | **Not reportable as a finding** (§4). Stated as expected under the selection confound. |
| Fewer than 8 patients survive `not_estimable` in either arm | **NOT ESTIMABLE**, no interval reported, on any statistic. |

The floor in the last row is fixed now because it is the number that decides
whether anything is reported at all, and choosing it after seeing the
estimability counts is exactly the defect `HANDOFF` §3 catalogues.

## 8 · What this cannot say, at any outcome

- Nothing about survivorship. A GUCA2A-high population preferentially destroyed
  is not transcript-detectable here or anywhere in this project.
- Nothing about carcinoma. Five routes terminate there (`HANDOFF` §2).
- Nothing about adenomas outside this deposit. One cohort, one centre,
  23 patients.
- Nothing causal about `APC` or `BRAF`. Family G is an association at n=5 and 4.
- It does not repair avenue A's stated qualifier that it is one cohort. A
  subgroup of one cohort is still one cohort.

## 9 · Gate conditions — six discharged 2026-09-10, the seventh waived

Conditions 1–6 are discharged and recorded in
`docs/chen_subtype_data_gate.md`; 7 carries a recorded waiver. **§9 is
discharged: this document is ready to lock.**

1. **DISCHARGED** · **Data-use terms** for HTAN/cBioPortal reviewed and recorded, in the form
   `docs/` used for the UniToPatho CC-BY review. A public API is not a licence.
2. **DISCHARGED** · **Label provenance confirmed at source**, not inferred from an attribute
   name: that `POLYP_TYPE`/`POLYP_SUBTYPE` are pathologist diagnoses and not
   anything derived from the deposited transcripts. Chen et al. *Cell* 2021
   methods, plus the HTAN biospecimen data model.
3. **DISCHARGED** · **Snapshot pinned** — the cBioPortal responses written under `data/` with
   sha256 in `data/manifest.csv`, and the analysis reading the snapshot.
4. **DISCHARGED** · **An `id_provenance` check** that the 27 labelled lesions are the same
   physical specimens as the scRNA biospecimens, on the standard
   `wes_subtype_plan.md` §1 rule: exact equality, suffix families not trusted.
5. **DISCHARGED** · **`src/common/label_provenance.py` declarations written and passing** for
   both families, before anything is read.
6. **DISCHARGED** · **The attrition table produced first**, covered versus uncovered patients
   compared on the decomposition values, and inspected before the arm contrast
   is computed. If the uncovered 21 differ systematically, that is recorded and
   the contrast is reported as conditional on coverage.
7. ~~**W2 review of the interval choice**~~ — **WAIVED by the project owner,
   2026-09-10.** The condition read: *"since §5 departs from the project's
   default estimator. `src/harness/` is W2-owned and `CONTRIBUTING` §2–3 route
   changes there through a PR — and `HANDOFF` §6k records that three such
   changes have already landed without one. This one does not."*

   It now does. The waiver is recorded rather than the condition deleted,
   because the reason it was written has not changed: this is the fourth
   harness-adjacent change to proceed without the review CONTRIBUTING requires,
   and a reader who finds a §5 interval they disagree with should be able to see
   that no second pair of eyes was on it. No code under `src/harness/` was
   altered under this waiver.

## 10 · Standing

**Ready to lock.** §9 is discharged. On the locking commit, the arm
sizes in §3 and the floor in §7 are frozen as written, and only then may a
decomposition value be joined to a label.

**It does not outrank the write-up.** The WMHS deadline is 15 September 2026.


## RESULT — RUN 2026-09-10. Weighting-unstable; no target-specific conclusion.

`results/2026-09-10_dddd34f/`, clean tree, design locked at `3a6d446` before the
join. Both arms clear §7's floor: 13 and 9 estimable patients, **zero
not-estimable in either arm**, as gate 6 said in advance.

### The predicted half holds; the informative half does not

**GUCA2A's intrinsic difference contains zero under all three weightings** —
doubly robust −1.263 [−5.397, +2.872], normal −2.378 [−7.927, +3.172], tumour
−0.148 [−2.967, +2.671]. That is the prediction's first clause, and §6 already
warned it is the half observed in §6d, so it is **not evidence**.

**CDX2 excludes zero under `tumour` only**, +0.523 [+0.020, +1.026] — a lower
bound 0.02 from zero — and contains zero under doubly robust +0.428 [−0.118,
+0.975] and normal +0.334 [−0.279, +0.946]. The `SER`-arm clause, the only
informative one, is **not supported**.

### The §7 branch is not single-valued

| weighting | branch |
|---|---|
| doubly robust | **NOT SEPARABLE AT THIS RESOLUTION** |
| normal | **NOT SEPARABLE AT THIS RESOLUTION** |
| tumour | PREDICTION REPLICATES |

Two of three weightings close item 3; one does not. **The pre-committed
consequence cannot be applied**, and this is exactly what fixing "all three,
none primary" at lock was for: choosing `tumour` would have reported a
replication, and choosing either other would have closed the item, on the same
data.

### The pairs say the separation is not where the design looks

§5 required every pair for the §6d reason, and it earns it here. Of the 15 pair
contrasts, **every one that excludes zero involves `ACTB` or `KRT8`** —
`ACTB − CDX2` −6.854, `ACTB − GUCA2A` −5.163, `CDX2 − KRT8` +8.660,
`GUCA2A − KRT8` +6.969, and so on. `ACTB − KRT8` itself contains zero: the two
controls move together. Every pair among `CDX2`, `EPCAM`, `GUCA2A` and `MS4A12`
contains zero.

Per gene, the largest and most consistent arm differences are the **controls** —
`ACTB` −6.426 and `KRT8` −8.232 under doubly robust, both excluding zero on all
three weightings — while the targets do not move. A contrast whose housekeeping
and identity genes separate the arms more strongly than its targets is not
measuring something specific to the target programme.

### Reading, and what is not claimed

**The recorded status is `DISAGREES ACROSS WEIGHTINGS`, and the headline is
"weighting-unstable / inconclusive at this resolution; no target-specific
conclusion."**

**It is not "NOT SEPARABLE AT THIS RESOLUTION", and the correction matters.**
That was written here first, on the grounds that two of three weightings return
it. **§5 fixes all three weightings with none primary and pre-commits no
majority-vote rule**, so counting branches is an unregistered decision rule
adopted after seeing which way the count fell — the exact move the "none
primary" clause exists to forbid. The artifact was right and the prose was
wrong. Item 3 is closed as *inconclusive*, not as separated-or-not.

13 against 9 patients was named in §7 as the reason to expect a weak result, so
this is **not a biological null** — it is a resolution limit, and the arms are
what this deposit has.

**POST-HOC, and outside the locked design**
(`results/2026-09-10_4456e4d/chen_subtype_control_diagnostic.parquet`): the arms
do differ in measurement properties. `frac_mature_tumour` is the largest
separation in the panel — median **0.291** in `AD` against **0.720** in `SER`,
standardised difference **−1.677** — while `frac_mature_normal` is flat at
**0.043**. `SER` also carries ~30% more depth (9,290 against 7,171, −0.673) and
about twice the cells (282 against 139, −0.430).

That is §4's stated confound, measured: the arms differ most in mature-cell
fraction *in the lesion*, which is the property the pathologist reads to assign
the label. It is a description and no interval is reported. **It does not
establish that the depth or resolution difference caused the control movement**,
and §8 bars the interesting readings regardless.

**Family G is closed on the floor, not on access**
(`results/2026-09-10_4456e4d/chen_subtype_family_g_floor.parquet`): 5 and 4
estimable patients against §7's floor of 8, so **NOT ESTIMABLE**, no interval
on any statistic. The calls came from the open cBioPortal MAF. The Synapse
certification, the 403 on `syn23520239` and the Level 3 download ACL were never
the binding constraint on this analysis — the arm sizes were.

The compositional term was computed and is **not reported as a finding** under
§4, whatever it shows.
