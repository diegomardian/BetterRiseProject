# What is left to try, reviewed against the data

**Written 2026-09-05.** Every claim below was checked against the cached obs or a
committed table before being written down. Where a proposal did not survive that
check, it says so and why — a list of avenues is only useful if the dead ones
are marked.

Read `docs/HANDOFF.md` §2 first for what is already established.

---

## Tier 1 — in hand, and one of them is bigger than path C

### 1a. The decomposition on Chen_2021. **PRE-REGISTERED 2026-09-06 — `docs/prereg_adenoma_decomposition.md`.**

> **This heading used to read "at `best4`", and that was wrong.**
> `config/labeling_axes.yaml` is frozen and requires the split reported as a
> CURVE across rungs: *"A single point estimate would present a modelling
> choice as a measurement."* A `best4`-only reading is exactly that point
> estimate. The pre-registration runs every rung the cohort supports,
> `epithelial` included — it is degenerate by design and a curve whose lower
> bound is missing cannot show it is a curve.

The decomposition was abandoned because its ratio collapses:

    i/c = (f_N / Δf) × (m_T/m_N − 1)

As a gene's surviving per-cell mean → 0 the bracket → −1, the ratio becomes
`−(f_N/Δf)` — a property of the cell fractions, identical for every gene on the
same labels. In carcinoma that constant is −5.85 and five panel genes were
indistinguishable by it.

**In adenoma it does not collapse.** Measured on
`results/2026-09-05_d869bdd/icbi_adenoma.parquet`:

| gene | m_T/m_N (lineage) | bracket |
|---|---|---|
| KRT8 | 0.946 | −0.054 |
| ACTB | 0.834 | −0.166 |
| CDX2 | 0.791 | −0.209 |
| EPCAM | 0.737 | −0.263 |
| MS4A12 | 0.383 | −0.617 |
| GUCA2A | 0.374 | −0.626 |

An order of magnitude of spread, and nothing near −1. **The estimand the project
was built around may be identifiable here**, which is a larger claim than the
coexpression reading can make. It is the same run: the cells are already loaded.

*Temper it with this:* GUCA2A (0.374) and MS4A12 (0.383) are again
indistinguishable, while CDX2 (0.791) sits well above both. Same two blocks the
corrected specificity reading found — terminal differentiation down, intestinal
identity retained — and the same limit: the decomposition would separate GUCA2A
from housekeeping and from CDX2, and not from MS4A12. `m_T/m_N` is a ratio, so
unlike a detection delta it IS comparable across genes; that is why this table
survived the scale correction unchanged.

**Which makes the comparator choice the whole design.** Run against housekeeping
alone it will look like a clean gene-specific result. The identity markers are
what stop that, and they must be scored in the same run rather than added
afterwards. Identifiable is not the same as gene-specific.

**And score every pair, not the target against each.** The first specificity
table reported only `GUCA2A − X`, so the claim that intestinal identity is
retained — a statement about where CDX2 sits relative to the CONTROLS — had no
row behind it and was read off CDX2's delta looking small beside GUCA2A's. On
the raw detection scale that comparison reverses the answer. Full account in
`docs/HANDOFF.md` §6d; the corrected table is
`results/2026-09-05_9c43f4f/adenoma_specificity.parquet`.

#### Re-checked 2026-09-06, against the same table. Three things.

**The bracket table reproduces and 1a's core claim stands.** The quoted numbers
are a ratio of COHORT means; the median and geometric mean of the per-patient
ratios agree with them to within ~0.08 (GUCA2A 0.374 / 0.392 / 0.316; MS4A12
0.383 / 0.467 / 0.398; KRT8 0.946 / 0.951 / 0.943). *The arithmetic mean of
per-patient ratios does NOT* — it reads MS4A12 at 0.841 — because a ratio with a
small denominator explodes and MS4A12's normal-arm mean tops out at 8.0 CP10K
against GUCA2A's 38.7. That is an estimator artefact, not a finding, and it is
recorded here because it is the obvious first thing to compute and it is wrong.

**Estimability is fine at `best4` and worse at `lineage`, which is the opposite
of what the cell counts suggest.**

| rung | patients | median mature cells/arm | GUCA2A usable | MS4A12 usable |
|---|---|---|---|---|
| `lineage` | 44 | 255 | **41/44** | **41/44** |
| `best4` | 20 | 30 | **20/20** | **19/20** |

A patient is unusable when `m_N` or `m_T` is exactly zero — `not_estimable`
under invariant 1, never zero silencing. At `best4` the cells are fewer but they
are BEST4+ absorptive cells, which is where these markers actually live, so the
gene is more reliably present in the population being asked about. **This is
good news for 1a**: the rung it targets is the one that loses almost nobody.

**But the "temper it" note above is a `lineage` statement and does not carry to
`best4`.** At `best4` the two targets are **GUCA2A 0.517 against MS4A12 0.385** —
further apart than at `lineage` (0.374 / 0.383), not indistinguishable. Whether
that separation is real needs an interval, and an interval needs the per-patient
run. So: do not carry "GUCA2A and MS4A12 are again indistinguishable" into the
`best4` design, and do not carry the opposite either. It is open, and 1a is what
would settle it.

**What 1a still needs, which is not on this page:** `decompose()` takes
`frac_mature_normal` / `frac_mature_tumour`, and **no committed table carries
them.** `icbi_adenoma.parquet` has mature cell COUNTS with no denominator. So 1a
is not a re-read of committed tables the way the specificity correction was — it
is emit-the-fractions, pre-register, cluster run.

**And the interval it reports must not be the project's usual one.** At `best4`,
n=20, the percentile bootstrap over patients is 0.913× the correct width and
excludes zero **7.1%** of the time under a true null — see `docs/HANDOFF.md`
§3a. Use the Student-t interval, as `src/reference/jobs/mlh1_positive_control.py`
does.

### 1b. The atlas's own annotations. **PARTLY WRONG — two of three are unusable.**

- `SOLO_doublet_status` is **`singlet` for all 4,264,929 cells.** It is a
  constant. A cross-check against it cannot fail, which is this repository's
  own signature defect; the atlas already removed doublets, so there is nothing
  to cross-check. **Do not use it.**
- `microsatellite_status` is **94% missing on Chen_2021** — 7,999 MSS and 2,845
  MSI-H against 166,435 unannotated. MSI stratification of the ADENOMA reading
  is not available. It *is* available on Pelka (130,351 MSS / 110,279 MSI-H),
  so the pre-registered subgroup contrast belongs to the carcinoma cohort.
- `n_genes` / `total_counts` / `pct_counts_mito` are real and per-cell, and the
  QC path already computes its own from raw counts. Comparing the two is a
  genuine cross-check and the only one of the three that survives.

### 1c. Zheng_2022's within-patient gradient. **RAN 2026-09-08 — descriptive only, n = 3.**

Three patients carry normal → polyp → carcinoma. Descriptive only, and worth
reporting as such: it is timing evidence carcinoma cannot give, and it does not
need a premise to resolve across patients because the comparison is within one.
Report beside Chen_2021, never pooled with it (`MIN_STUDIES` is 3).

**"Descriptive only" was the right call and here is the number behind it**
(added 2026-09-06). At n=3 the project's percentile bootstrap over patients is
**0.372× the correct width** and excludes zero **25.1%** of the time under a
true null — one interval in four, on data with nothing in it. See
`docs/HANDOFF.md` §3a; the rate is `P(|t(n−1)| > z·sqrt((n−1)/n))`, a function
of n alone. **No interval from three patients may be reported here at all**, on
any statistic. Descriptive means descriptive: report the three trajectories and
let them be three trajectories.

**Result.** All three pre-specified patients retained measurable epithelial
cells at normal, polyp, and carcinoma after QC. GUCA2A detection fell from
normal to polyp in each trajectory (P1 0.579 → 0.040; P2 0.703 → 0.027;
P3 0.642 → 0.040), and MS4A12 did too. CDX2 did not follow one direction.
P1/P2 normal samples were much shallower than their polyp samples and P3 has
only 150 post-QC polyp epithelial cells, so the table is timing context only.
`results/2026-09-08_f6b22f6/zheng_three_stage_trajectories.parquet` and its
availability companion are the complete artifact. **Do not pool, test, or fit
these trajectories.**

---

## Tier 1+ — the missing positive control, and it outranks everything else

### The MLH1 reading is not primarily a biology test. It is the control this whole instrument has never had.

The atlas carries `MLH1_promoter_methylation_status` on **240,630 of Pelka's
340,686 cells** (70.6%) — 76,015 `meth`, 164,615 `no_meth`, and **zero patients
carry more than one value**, confirming it is a patient-level assay annotation
rather than anything derived from expression. That independence is what makes it
usable.

**Two corrections to the proposal.** It is 240,630 of 340,686, not all of Pelka.
And **Chen_2021 has no MLH1 annotation at all** — every cell is null — so this
cannot be folded into path C. It is a carcinoma-cohort reading.

**And the framing should change.** Proposed as the tier-B intrinsic contrast, it
is nearly tautological: promoter methylation silences MLH1, that is textbook,
and confirming it teaches little biology.

What it actually supplies is the thing this project has never had. **Every
negative and unresolved result the coexpression instrument has produced —
UNRESOLVED on three cohorts, UNRESOLVED at 13 studies, not-specific on adenoma —
rests on an instrument whose ability to detect a *known* silencing event has
never been demonstrated.**

MLH1 in methylated patients is a silencing event we know occurred, established
by an assay that is not transcription. So:

> Within the mature cells of MLH1-methylated patients, is MLH1 detection lower
> than in non-methylated patients?

- **Detected** → the instrument can see silencing when silencing is there. Every
  null it has returned becomes evidence rather than absence of evidence, and the
  13-study negative gets much stronger.
- **Not detected** → the instrument cannot see silencing at the panel's
  abundance, and **every null it has produced is uninformative**. That is a
  finding about the method, and it would reframe most of this project.

Either branch is worth more than another cohort. The cost is adding one gene to
`GENE_ROLES` and a stratum split, on data already on the cluster.

*Known risk, stated first:* the pre-registered MLH1 contrast
(`docs/prereg_g2_mlh1.md`) died on the detection floor at CP10K means of 0.039.
Detection at ≥1 UMI is more sensitive than a mean, which is the reason to
retry — but a feasibility check on mature-cell counts per stratum comes before
the reading, not after.

**Do the scale correction first, and the reason is about MLH1 rather than
tidiness.** At 0.039 CP10K and Pelka's median depth of 11,286 UMI, MLH1 sits at
a **4.3% detection rate** — the far end of the abundance range from GUCA2A's
44%. Its result will be read against this panel, and a cross-gene comparison on
the raw detection scale is not one the statistic supports (`docs/HANDOFF.md`
§6d). A positive control read on a non-comparable scale is uninterpretable in
exactly the way it exists to resolve. That correction is done:
`src/reference/detection_scale.py`.

### DONE, 2026-09-06 — built, pre-registered, and the DiD did not survive sizing

**`docs/prereg_g2_mlh1_within_stratum.md` is the design;
`src/reference/jobs/mlh1_positive_control.sh` is the run.** What is left is one
`qsub` against the cluster's atlas.

**The difference-in-differences recommended above is not available, and the
feasibility check is why.** It was the right correction — a between-patient
comparison is not what the rest of this project does — but only **29 of 62**
Pelka patients survive the pipeline's own filters, and the arms at `lineage` are:

| stratum | patients scored | median mature cells/arm |
|---|---|---|
| `mlh1_methylated` | **10** | 262 |
| `mlh1_intact_mmrd` | **4** | 127 |
| `mmr_proficient` | 14 | 182 |

**Four is the number the original prereg also reached after depth matching**, by
the GSE178341 route rather than this one. Two independent pipelines agreeing on
four makes it a property of the cohort, not of anybody's filters. And the
dilution cannot be fixed by stratifying, because the stratum you would stratify
into has four patients in it.

Measured rather than asserted (`results/2026-09-06_a0483ae/mlh1_two_sample_power.parquet`,
Welch, τ=0.2): at **75%** silencing the pre-registered DiD detects it **60.0%**
of the time. It is close to a coin flip at the effect size it exists to catch.

**So the reading is within-stratum, and it is powered:** 10 methylated patients,
~262 mature cells each, MLH1 at ~3.2% detection → about 8 positive cells per
patient per arm. **99.3%** power at 75% silencing, **73.7%** at 50%. The n=19
unmethylated arm is reported as secondary and CONFOUNDED (it mixes methylation
with MSI status); the n=4 arm is reported as UNDERPOWERED and carries no verdict
in either direction.

**Two things the sizing turned up that were not on anyone's list.**

*The atlas annotation and the week-0 clinical strata agree exactly* — 22 `meth`
against 22 `mlh1_methylated` on all 62 patients, no crossings. Two independent
derivations of the arm the reading is about, and a disagreement would have meant
the reading was not about the arm the prereg named.

*The interval was wrong before the gene was.* See `docs/HANDOFF.md` §3a: the
percentile bootstrap this project uses everywhere is **0.82× the width it claims
at n=10 and 0.53× at n=4**, by a closed form containing no data. The MLH1
reading reports a Student-t interval for that reason, and the measurement was
committed before the design was written so that the choice could not be a free
parameter.

*Provenance:* the earlier sizing figures in this section (23,256 meth / 54,623
no_meth cells, ~2,200 per arm for a 50% effect) came from the cached obs and
were **not from a committed table**. They are superseded by the numbers above,
which are.

---

## Reviewed 2026-09-06 — three proposals, checked against the repo

Three avenues (A, B, C) proposed after the MLH1 positive control returned
UNINTERPRETABLE. Reviewed to this file's own standard: **claims checkable
against the cached obs or a committed table were checked, and the ones that
could not be are marked as such.** Nothing external was verified from here.

### A. The decomposition on the adenoma already in hand — CONFIRMED as the flagship

**The premise is right and the arithmetic behind it holds** (see §1a's re-check).
This is the project's original deliverable, on the one substrate where the
collapse does not fire, and it needs no MLH1, no premise resolution and no
instrument sensitivity. It should be next.

**The input contract is as described, with two additions.** `decompose_cohort`
requires ten columns, not four: the six identifiers plus the four statistics.
`icbi_adenoma.parquet` carries `mean_normal`/`mean_tumour` (as
`cp10k_normal`/`cp10k_tumour`), `n_cells_mature` (as `n_tumour` — correct, the
tumour arm is what `classify_estimability` reads), `patient_id`, `study_id`,
`granularity_rung`; `labeling_axis` is a constant (`stem_pole`). **Only
`frac_mature_normal`/`frac_mature_tumour` are genuinely missing.** Confirmed.

**Three gaps the plan does not name, all in the compositional term — which is
the half the decomposition exists to produce.**

**1. The denominator is a decision, not a lookup, and it is load-bearing.**
Every patient in the MLH1 run reported *exactly* 25.0% of epithelial cells as
`unresolved_depth` — 351/1405, 404/1614, 1103/4411, thirty of thirty. That is
not data, it is `DEPTH_QUANTILE = 0.25`: the depth target is the 25th percentile
so a quarter of cells fall below it by construction. So a quarter of the
epithelium is excluded from any denominator, and **whether that quarter splits
evenly between the two arms is unknown and unchecked.** It is applied in
`assign_labels`, per patient, over both arms pooled and *before* depth matching.
If the arms differ in depth pre-matching, the exclusions concentrate in the
shallower one and the mature *fraction* differs between arms for a purely
technical reason — which is the compositional term measuring sequencing depth.
**Emit the per-arm unresolved share and pre-commit a tolerance on it.** This is
the one new risk in A and it is not "none".

*What is reassuring:* `depth_ratio` runs 0.97–1.05 across all 44 adenoma
patients. That is measured **after** matching, so it does not settle the
question — but it does mean matching is not being asked to close a large gap.

**2. The mean and the fraction would be computed on different populations.**
`rows_for_patient` selects mature cells, depth-matches *those*, and computes
`cp10k_*` on the survivors. The epithelial denominator a fraction needs is not
matched and cannot be, since matching is defined on the mature set. So
`mean_*` is post-matching and `frac_mature_*` is pre-matching unless the design
says otherwise. **Pre-commit which**, and report it on the row — this is
precisely the kind of population mismatch `build_decomposition_summary` already
records learning the hard way.

**3. Student-t in the schema slot collides with a settled decision.**
`docs/open_decisions.md` #10 puts `bootstrap_over_patients`'s **percentile**
band in the schema's `ci_low`/`ci_high` — W2 proposed, W4 confirmed 2026-08-22 —
and `src/estimator/` is W4's under CONTRIBUTING §2. So A may not simply swap the
interval. **Carry the Student-t interval BESIDE the schema band, in its own
columns**: additive, needs nobody's approval, changes no frozen decision, and
lets a reader see both. Note what the schema band is worth here — at `best4`,
n=20, the percentile bootstrap is a **7.1%** test, not a 5% one (§3a).

**4. The comparator rule is right and should be pre-committed as stated.** Score
housekeeping and the identity markers in the same run. §1a's re-check adds a
number to it: at `best4` the two targets are 0.517 against 0.385, *not*
indistinguishable as they are at `lineage`. Whether that separation is real is
exactly what A would settle, so the prereg must not assume it in either
direction.

### B. More polyp substrate — the ICBI claim is exact, and B1 is not a new idea

> #### RUN 2026-09-07. The instrument objection is dead; a REFERENCE-ARM objection replaced it.
>
> `results/2026-09-07_5ce00b8/becker_feasibility{,_by_arm}.parquet`.
>
> **The physical risk B1 was gated on does not exist.** §3 feared that GUCA2A
> and MS4A12, being cytoplasmic, would be lost by nuclear sampling. In Becker's
> `healthy_donor` arm, on the detection scale against Chen's mature-normal
> baseline, they are the **two best-preserved genes on the panel** — GUCA2A
> −0.258 (82% of Chen's rate) and MS4A12 −0.478, against ACTB −1.246, EPCAM
> −1.477, KRT8 −2.171 and CDX2 **−3.182**, the worst. The worst-hit gene is a
> nuclear transcription factor, which is the opposite of what the cytoplasmic
> story predicts. **snRNA-seq is not the problem.**
>
> **The reference arm is.** `normal` in this deposit is a FAP donor's uninvolved
> mucosa, not healthy colon — `becker_io.DISEASE_STAGE_MAP` kept them apart as a
> fifth arm precisely so this could be seen. GUCA2A falls **+1.689 in log** from
> healthy donor to FAP-unaffected against a control floor of +0.514 to +1.186:
> **2.38× beyond what the controls explain**, where MS4A12 (0.92×) and CDX2
> (0.50×) sit inside that floor.
>
> **So the paired design's reference is contaminated by the phenomenon it is the
> reference for.** Each polyp is contrasted against that donor's own uninvolved
> mucosa; if the mucosa is already GUCA2A-depleted, the contrast understates the
> loss. This is an **estimand** problem: labelling does not fix it and more
> donors do not fix it. B1 was demoted for cohort size (4 paired donors); this
> is a second and independent reason, and it is **structural to FAP** — the
> germline *APC* hit is in the whole colon, so "unaffected" is pre-lesional
> everywhere. Sporadic adenoma does not carry it in the same form.
>
> **n = 2.** The healthy-donor arm is B001 and B004. Below every patient floor
> in this repository, post-hoc, and not pre-registered. **It is an observation
> that motivates a design and it may never carry an interval.**
>
> **The gate closed the same day. B1 is NOT LICENSED.**
> `results/2026-09-07_2305f23/`. On §3's own quantity — mature cells of the
> reference arm — GUCA2A clears the 0.10 floor (0.079 → 0.148) **by depth**: its
> enrichment is +0.669 against a control band of +0.517 to +0.739, between the
> two housekeeping genes. The label is not at fault: **MS4A12 outruns its
> controls in all three arms**, so the markers do find mature colonocytes.
> GUCA2A does not concentrate in them — not in FAP mucosa and not in healthy
> donor colon. **Nothing further on this deposit is pre-committed**; the prereg's
> RESULT has the full account, including why the `tumour` arm's apparent
> enrichment cannot be read.
>
> **One methodological finding worth carrying forward.** The label was built
> from the mature-colonocyte program, which GUCA2A belongs to. Invariant 2 bars
> a *target* from its own label; it does not bar genes **co-regulated with** the
> target, and that is nearly as circular. Proposed as an amendment to invariant
> 2 — `src/schema.py`-adjacent, so a PR with two approvals, not taken unilaterally.

#### The substrate criterion, stated sharply enough to search on

The requirement is no longer "more polyp data." It is:

**sporadic adenoma, paired with the same patient's normal mucosa, at nuclear
sensitivity or better.**

FAP is disqualified by the reference-arm argument above, not by size. Chen_2021
is exactly this substrate, which is why B1 was the replication and why losing it
costs. **No accession is named here on purpose** — this repository has spent
real time on two deposits whose format and contents were assumed rather than
checked (`docs/HANDOFF.md` §4), so a candidate is not an avenue until its
`--inspect` has run.


**"ICBI is exhausted" is VERIFIED, precisely.** Two studies carry any `polyp`
cells at all, out of 49:

| study | polyp cells | patients |
|---|---|---|
| `Chen_2021_Cell` | 93,913 | 94 |
| `Zheng_2022_Signal_Transduct_Target_Ther` | 13,045 | **3** |

So B2's n=3 is confirmed, and there is no third polyp cohort in the atlas.

**B1 is the substrate this project froze its third axis against in week 0, and
nobody has ever fetched it.** `config/labeling_axes.yaml` — frozen, PR plus two
approvals to change — names axis 3 as:

    chromatin:
      basis: "A different measurement — chromatin accessibility"
      transcript_based: false
      source: "Becker/Chang multiome"

with the caveat *"Not transcript-based, and therefore the strongest defence
against label leakage. Week 13+."* **That is Becker 2022.** So B1 is not a new
avenue competing with the others; it is the pre-registered week-13 substrate,
and the circularity objection every design here carries was answered in advance
by a dataset nobody has downloaded. That materially raises its ranking, and it
means C3 is not a separate item — it is what axis 3 always was.

**B1's stated risk is the right one and here is its number.** snRNA-seq on
cytoplasmic transcripts. The panel's baseline detection in the adenoma normal
arm, which snRNA would have to preserve:

| rung | GUCA2A | MS4A12 | CDX2 | ACTB |
|---|---|---|---|---|
| `lineage` | 0.437 | **0.363** | 0.822 | 0.984 |
| `best4` | 0.582 | **0.279** | 0.868 | 0.992 |

**MS4A12 at 0.279 is the floor of the panel.** A nuclear protocol that halves
cytoplasmic detection puts it near zero, and a gene that cannot be detected
cannot be a tier member. **So B1's feasibility gate is DETECTION OF THE PANEL,
not cell count** — the shape of `icbi_premise_feasibility.py`, run on the panel
before anything else. The scATAC arm does not share this risk and is the more
valuable half regardless.

**B3 could not be checked from here.** No accession in it is verifiable against
the cached obs. Treat the cell counts as unconfirmed until a read-only
feasibility check exists.

### C. The survivorship discriminators — consistent, and C2's target choice is corroborated

**C1 (segmented spatial) — LOOKUP COMPLETE 2026-09-07, and it is not runnable
on a stock catalogue panel.** `results/2026-09-07_8fa00f2/panel_coverage.parquet`,
job `src/reference/jobs/panel_coverage.py`, gene lists in
`data/manifest.csv` with sha256.

| gene | role | CosMx 1K | CosMx 6K | Xenium Colon v1 | Xenium Prime 5K |
|---|---|---|---|---|---|
| **GUCA2A** | **target** | absent | absent | **present** | absent |
| MS4A12 | identity | absent | absent | present | absent |
| CDX2 | identity | absent | present | absent | present |
| EPCAM | epithelial | present | present | absent | present |
| KRT8 | control | present | present | absent | absent |
| **ACTB** | **control** | absent | absent | absent | absent |

The Human Colon panel supplies the target but no control; Prime 5K supplies an
identity marker but no target or control. **No stock panel has C1's required
target, identity marker, and control.** C1 needs custom probes or a different
spatial design — a different cost and lead time from “order the colon panel.”

*The tidy explanation is wrong, so do not repeat it.* Absence does not track
abundance: `KRT8` is the **highest**-expressing gene in the panel (25.8 CP10K,
normal arm) and it is on both, while `ACTB` at 19.5 is on neither. What the
present genes have in common is that they are cell-typing markers, which is what
these panels are built for — the CosMx 1K's own description is "robust cell
typing." This project needs a housekeeping control (discriminates nothing, so a
typing panel omits it) and functional maturity markers. That reading is an
interpretation and is labelled one; the table above is the measurement.

The Visium caveat C1 carries is already this repo's (`docs/HANDOFF.md` §6e: a
spot is still a mixture). Ruling Pelka's own GeoMx out as a primary is right —
region-level, 3 samples.

**C2's instinct to target CDX2 before GUCA2A is right, and a committed table
supports it more strongly than the proposal claims.** The pre-registered CIMP
screen (`results/2026-09-05_9203809/`) returned NOT SPECIFIC with GUCA2A falling
**less** than CDX2, +0.544 [+0.219, +0.878]. So at bulk level CDX2 is the gene
with the real fall, and it is the better first methylation target on the
evidence rather than only on field convention.

*One thing to carry with it:* CDX2 behaves differently in the two substrates.
In carcinoma it falls (MLH1 run, methylated arm: −0.684 [−1.182, −0.186],
excluding zero); in adenoma it is indistinguishable from housekeeping on the
load-bearing scale (§6d). C2 is a TCGA carcinoma study, so it is the arm where
CDX2-down is established — consistent, but the two must not be quoted as one.

**C3 folds into B1**, as above — it is axis 3.

### Ranking after this review

1. **A**, unchanged — the flagship, in hand, and the only one that produces the
   README's actual deliverable. Fix the three gaps above in the prereg.
2. **B1**, raised, and **now pre-registered** —
   `docs/prereg_becker_replication.md`. The frozen axes file's own week-13
   substrate, unfetched, and the only thing here that answers the circularity
   objection outright. It is also the **one closing path** for avenue A's
   largest open item: the cross-gene statistic was post-hoc on Chen_2021 and is
   fixed in `ac7eca1`, so a second substrate makes it confirmatory.
   **Blocked on three things nobody has checked** — the accessions, the file
   format (reported as Seurat objects, not a GEO matrix), and disk, since both
   cluster quotas were near full after the ICBI fetch. Gate on panel detection
   under a nuclear protocol, per gene: MS4A12 at 0.363 is the floor and GUCA2A
   failing ends the replication.
3. **C2**, on CDX2 rather than GUCA2A, with the 450k already in hand.
4. C1, B3, B2 as supplementary — **and C1's cost changed on 2026-09-07**: the
   target is not on a stock CosMx panel at any plex, so it is custom probes, not
   a catalogue order. See §C above.

---

## Tier 2 — the survivorship discriminators

None of the above separates *silencing* from *GUCA2A-high cells having been
preferentially destroyed*. Nothing transcript-based can. These are the
measurements that could, in cost order:

| | what it buys | cost |
|---|---|---|
| Cell-type-resolved methylation (EPISCORE) | methylated survivors = silencing; unmethylated survivors = death-selected. The direct discriminator | weeks of new machinery |
| Becker 2022 scATAC (HTAN) | a chromatin maturity label — identity not defined by transcription, which answers the circularity objection outright | a research project |
| Spatial Xenium/CosMx | identity from anatomy; tests the CDX2 tumour-bud prediction | data hunt; still survivorship-limited |

The CIMP screen weakened the methylation prior but did not close it: GUCA2A
silencing need not be CIMP-tied.

---

## Tier 3 — real results that are not mechanism

### Marker → survival
Bulk GUCA2A/CDX2 against DSS/PFI, purity-adjusted, on committed TCGA. Answers
the CDX2-adjuvant question — *does the marker carry the prognostic signal* —
which is a legitimate result and mechanism-agnostic. Laptop-runnable.
**Needs its own pre-specification**; the Stage 4 lock excludes it deliberately
(`not_prespecified`), and a directional prediction reported after the fact is a
story.

### ~~Explaining the 13-study heterogeneity~~ **ANSWERED 2026-09-07 — and not by iCMS.**

**`docs/prereg_meta_weight_calibration.md` RESULT; `docs/HANDOFF.md` §6k;
tables `results/2026-09-07_40515ce/`. Do not queue this.**

The heterogeneity is not a subtype story. **57% of Cochran's Q is one study,
Khaliq_2022, at n = 3** — holding 20.4% of the fixed-effect weight against
Pelka's 10.0% at n = 29. At a patient floor where the inverse-variance weight
has a finite mean (n ≥ 4; `E[w]/w = (n-1)/(n-3)`, no mean below that) the
studies are **homogeneous**, I² falls to 0.63 and 0.71, and KRT8 stays
UNRESOLVED for the opposite reason: they agree the control fell 0.41–0.48 log2
and the ±0.5 tolerance straddles it.

**The iCMS version was not runnable anyway**, which is worth recording since
this page recommended it. The atlas obs carries no iCMS column, and the
covariates it does carry cannot support a meta-regression at k = 11:
`matrix_type`, `suspension_type` and `enrichment_cell_types` are **constants**
across all eleven — a covariate that cannot vary is a check that cannot fail —
`tissue_cell_state` is fresh for 9 of 11, and `platform` is 8 of 11 10x with
each remaining level at k = 1.

---

## What I would do, in order

0. ~~**Fix the cross-gene scale.**~~ **DONE**, 2026-09-05
   (`results/2026-09-05_9c43f4f/`). It had to come first because MLH1 lands at
   the far end of the abundance range and would have been read against this
   panel.
1. ~~**The MLH1 positive control.**~~ **RAN 2026-09-06 — UNINTERPRETABLE.**
   The premise does not hold in the methylated arm (ACTB +0.443 on log2
   expression, tolerance 0.5). Not a negative result; the pre-committed gate
   fired. **And it closed the question rather than answering it:** the
   instrument can only be validated where the premise holds, the premise has
   held only on adenoma, adenoma carries no MLH1 annotation, and Pelka is the
   only study in the 49-study atlas that does. `docs/HANDOFF.md` §6g.
2. ~~**The decomposition on Chen_2021** (1a).~~ **RAN 2026-09-06 — IDENTIFIABLE,
   and it is the project's deliverable.** The `i/c` ratio runs 0.33 to 2.03
   against carcinoma's collapse onto −5.85; the pre-registered falsifier did not
   fire. Tier-level, not gene-specific. **Quote `lineage` only — `best4` is
   retracted**, and note that this section's old "at `best4`" framing was wrong:
   the frozen axes file requires a rung CURVE, not a point.
   `docs/HANDOFF.md` §6h.
3. ~~**The Wnt mechanism test** (D1).~~ **RAN 2026-09-06 — TECHNICAL**, a clean
   negative. GUCA2A sits inside the control floor and it survives an
   over-conditioning objection raised after the run. §6i.
4. ~~**The DIS/VAL stability split.**~~ **RAN 2026-09-06 — AMBIGUOUS**, with
   zero sign reversals across three specimen collections. §6j.

### Where it stands at the end of 2026-09-07 — read this first

**Nothing is half-finished except three things, and all three are named below.**

**This block outranks every dated section under it.** Those sections each say
they supersede what came before them, which is true of their own date and no
later one — so when one of them says an item is open and this table does not,
this table is right. On 2026-09-07 a reader of the ranked list below was sent to
retry Crowell `242` and to chase the Xenium panel, both of which had closed
hours earlier. Those entries are now corrected in place, but the precedence rule
is what stops the next one.

| open | what it needs |
|---|---|
| **The write-up** | Eight days. ~15 results across 2026-09-06/07 in neither paper, and `neurips_2026.sty` still unvendored so `./build.sh` has never run. **This outranks everything.** |
| ~~**D2 survival**~~ | **CLOSED 2026-09-07 — null on every endpoint.** `results/2026-09-07_b42a2cd/`; PFI ABSOLUTE 0.9697 [0.8564, 1.0979], every interval containing 1, PH clean, 11.27 events/df. Locked, implemented, run. `docs/prereg_d2_survival.md` §8. **Do not queue it.** |
| ~~**Invariant 11**~~ | **CLOSED 2026-09-08.** The guard, forcing tests, and declared call sites are merged; `docs/invariant_11_proposal.md` records the reviewed result. |

**Closed on 2026-09-07 and not worth reopening:** Crowell (7 of 7, §6l),
Becker/B1 (not licensed), C1's panel lookup (no stock panel carries a target and
a control together), the meta-layer weights (§6k), and the carcinoma `best4`
denominator.

**No conditional substrate remains executable locally.** Lesion-level Wnt
closed on its fixed detection gate; its Becker per-donor inventory is durable.
Item 2b is a distinct H&E morphology-to-molecular-label data gate, not a new
reading of avenue A or a replacement for Crowell/A2, and MHIST has now closed
as one source for it.

### Roadmap amendment, 2026-09-07 — seven items, after the Crowell run

Supersedes the corresponding entries above where they conflict. Each item names
what it is *not* as tightly as what it is, because five of the original six
were written one class stronger than the data supports. Item 2b was added
2026-09-07 and is a data gate, not an avenue in hand.

**1 · Crowell — CLOSED 2026-09-07, seven of seven.** The 504 cleared on a
retry; `242.h5ad` is `results/2026-09-07_da1f46b/` and the aggregate is
`7336e02`. `docs/prereg_crowell_multisection.md` RESULT. The result is a
**patient-level target-versus-epithelial-domain pattern**: GUCA2A −1.307
[−1.872, −0.741] and MS4A12 −1.241 [−1.874, −0.609] exclude zero and are
negative in 7 of 7 blocks; CDX2 +0.327 [−0.151, +0.805] and EPCAM +0.234 do not.
Section 7's pass branch, at 3.01× avenue A's width.

**The seventh block is the argument for §3's fixed order.** CDX2's lower bound
ran −0.251, −0.156, −0.020, −0.151 across four, five, six and seven blocks — at
six it was 0.020 from clearing zero *opposite* the target, and 242's CDX2 came
in at −0.395 and pushed it back out. Stopping at six, which felt complete, would
have reported a discriminator on the edge of resolving.

It is **not** evidence separating compositional from intrinsic loss
(Amendment 2 — the only control-role gene is an epithelial keratin) and **not** a
per-cell silencing result (feasibility §7). Below-floor observations are bounds,
not effect-size points (Amendment 4). **Nothing here is open.**

**2 · A2 / ML — contingent infrastructure only.** The public MxIF product is 42
pixel NPZ regions with **no cell segmentation and no cell-by-marker table**. A
configured Synapse account received HTTP 403 for root `syn23520239` on
2026-09-09; the metadata-only job stopped rather than treating a partial tree as
an absence. **No biological analysis proceeds from the public product alone.**
Crypt-position assignment needs independently annotated
anatomical ground truth, patient-held-out evaluation, and a
morphology-or-geometry-derived label — and it still does not touch the per-cell
false-negative problem (`docs/HANDOFF.md` §6g).

**2b · H&E foundation-model ML — a separate, conditional molecular-prediction
avenue.** H&E offers an orthogonal morphology measurement, so its validity is
not limited by the transcript false-negative floor. It does **not** validate
the compositional-versus-intrinsic decomposition, repair Crowell's control
limit, or turn same-cohort morphology into an independent replication.

The proposed estimand is deliberately narrow: a frozen H&E encoder predicts a
pre-specified, independently assayed molecular endpoint for the **same polyp
specimen**, with an abstention rule fixed before the held-out evaluation. It is
not “another polyp classifier.” Polyp subtype or dysplasia labels obtained
from the very H&E slide being classified are diagnostic-label imitation, not an
independent biological claim; they may be used only to test engineering
performance, never as the project result.

**Data gate before a model, image download, or encoder choice.** The candidate
source must provide (1) licensed full-resolution H&E or WSI files; (2) a
specimen-exact crosswalk from each image to an independently generated
molecular call, not merely a participant-level join; (3) a fixed molecular
endpoint with its prevalence, missingness and complete-case count shown before
fitting; and (4) enough patients and sites for patient-held-out evaluation.
The reported COLON MAP targeted-sequencing and imaging collections establish a
specimen-exact metadata join. HTAN r7 metadata now maps all 18 in-scope exact
H&E--VCF biospecimens to one CRDC object each. One authenticated CRDC,
header-only pilot also shows that a selected H&E object is readable without
downloading a slide. That does **not** establish full-resolution release status,
usable clinical strata, or a callable endpoint; the collection remains a source
to inventory, not data in hand. Public polyp image datasets can benchmark an
encoder but cannot supply the molecular endpoint.

**Inventory record:** `docs/he_molecular_data_gate.md` records the HTA11
candidate screen. `src/reference/jobs/he_molecular_gate.py` enforces the exact
biospecimen join and reports attrition. The currently exported metadata returns
`NOT LICENSED`: resolution is not established, case diagnosis/site are missing,
and no endpoint is pre-specified. Candidate molecular access is `Synapse`; the
metadata-only authenticated check has now resolved all 18 candidate entities,
but it does not establish VCF-content download. Its versioned output can clear
only the pending metadata-access reason in the main gate. The one-image CRDC
pilot is deliberately not a gate-wide resolution result. The r7 imaging
dimensions conflict with the native TIFF headers, so r7 is used for its
crosswalk only. The durable bounded audit now covers all 18 exact objects:
12 parsed headers carry the same generic 72-DPI-equivalent scale and six lack
the required physical-resolution tags. It yields 0/18 usable physical-scale
measurements, so the full-resolution clause remains unmet.
`src/reference/jobs/he_molecular_synapse_access.py` performs that check.
Neither job selects an endpoint or licenses image download.

**If the gate passes.** All tiles from a patient and specimen remain in one
split. Encoder and molecular endpoint are fixed before fitting; calibration and
the abstention threshold are selected only on training/validation patients; an
untouched patient-held-out test reports coverage, selective error, calibration,
and the ordinary non-abstaining baseline. A site-held-out test is required when
more than one acquisition site exists. Failure of the crosswalk, prevalence, or
split requirement ends the avenue as **NO SUBSTRATE**, rather than replacing
the molecular endpoint with a morphology-derived diagnosis.

**3 · AD versus SSL — CLOSED 2026-09-10, inconclusive on the open pathology
labels; the genotype arm closed on the §7 floor. See the 2026-09-10 amendment
below. Original entry follows.** ~~WES provenance/callability gate in progress, not yet an
analysis.** The ICBI atlas carries `polyp` with no subtype and Crowell is
**entirely TVA**; A2's public pixel inventory remains too small and cannot make
a biological claim. However, the Chen collection may carry a modality-orthogonal
FFPE WES label for a subset of avenue A's 44 lineage patients. This is not yet
a usable AD/SSL substrate: participants can have multiple scRNA polyp sample
IDs, while VCF IDs span more than one naming family. Participant overlap and
identifier suffixes therefore do not prove the same lesion.

The sole live route is [`docs/wes_subtype_plan.md`](wes_subtype_plan.md): an
outcome-blind `id_provenance` chain must establish exact equality between a
specific polyp scRNA biospecimen and a VCF's documented assayed/originating
biospecimen. Only then may authenticated VCF-content access and header
capability be checked. The candidate arms are positive-only, molecular-signature
labels (`BRAF` V600E-positive versus truncating-`APC`-positive), not all AD
versus all SSL; dual, neither, and uncallable lesions remain visibly
unclassified. A precommitted arm-size/precision rule is required before any
variant record or transcriptomic subgroup outcome is read. Failure at either
gate closes **NO SUBSTRATE AT THIS RESOLUTION**.

**4 · Lesion-level Wnt — feasibility-gated, and the gate is partly answered.**
D1 was within-patient and explicitly conceded between-lesion variation. The
substrate now exists but is thin: **Crowell gives 2 of the 6 blocks that were
audited for sub-regions with multiple adenoma sub-regions** (`110` and `221`,
three each — from `crowell_subdomains_exploratory.parquet`; `242` arrived after
that audit and has not been read for sub-regions). Becker has **42 polyp lesions across 12
donors** — re-read from the series matrix 2026-09-07. **The 43 quoted here
before was sequencing rows, not lesions:** `A002-C-010` has two
technical-replicate GSM rows and counting GSMs turns one physical polyp into
two. **31 of those lesions sit in the 4 donors that carry a reference arm** —
A001 (7), A002 (10), A014 (6), A015 (8). `becker_io.tumour_lesion_counts`
collapses replicates by `sample_id` and keeps the row count beside the
biological count so the collapse stays auditable.
**The independent unit is still 4 donors, not 31 lesions** (invariant 5, Becker
Amendment 2), so this needs a lesion-within-donor pre-specification before any
analysis. **The inventory is now an artifact:**
`results/2026-09-07_682ad52/becker_lesion_inventory.parquet`, written from a
clean tree by `--lesion-inventory` off the 8.8 KB series matrix alone. Twelve
donors, 42 lesions, 43 sequencing rows, and A002 carries the replicate in the
open at 10 lesions against 11 rows. `paired` is True for A001, A002, A014 and
A015 and sums to 31.
Sub-regions and lesions are **not independent patients**; any analysis is
patient-clustered or mixed-model, with the Wnt score and adjustment set fixed
first, and read as association. Crowell's sub-domains stay `exploratory` under
its own pre-registration and cannot be borrowed.

**5 · D2 survival — CLOSED 2026-09-07, a pre-registered null.** The whole
sequence ran in one day: gate PASS, lock, implementation, run.
`results/2026-09-07_b42a2cd/` — **PFI ABSOLUTE 0.9697 [0.8564, 1.0979]** and
every other interval also containing 1, proportional hazards clean on every
term, 11.27 events per non-stratum df. The verdict label reads
PURITY-SOURCE SENSITIVE because §6's reversal rule fired on 0.9697 against
1.0090, which is two nulls straddling 1; §6c records why the rule was not
re-decided after the fact. §6a's missingness concern was real — ABSOLUTE drops
107 participants whose GUCA2A is *lower* than the retained set's — and the
ESTIMATE population, which loses only 48, returns the same null. Full reading in
the pre-registration's §8. The original framing follows.

**The gate PASSED 2026-09-07; the pre-specification was then locked and run.** `results/2026-09-07_1dc7624/`: n=624, IQR 2.9607, 312
below and 312 above the median — so the "GUCA2A is at floor in most tumours"
worry is refuted for bulk TCGA, and the gate licenses a pre-specification and
nothing else. `docs/prereg_d2_survival.md` is written and **unlocked**; its §6a
names three decisions it leaves open (purity missingness is not MAR and drops
~106 participants preferentially at low purity; multiplicity handled by
structure; ten-events-per-df uncalibrated). **What is open is the lock and the
implementation, not the gate.** The original text follows.

Assess GUCA2A's distribution and usable variation in the committed TCGA
matrices **without touching outcomes**; record the gate and its threshold; model
only if it passes. If the marker is at floor in most tumours there is no
variation to regress on and the prespec would be written for a study that cannot
produce a non-null. Same shape as B1's detection gate, and B1 is why this
sentence exists. `docs/d2_feasibility_gate.md` fixes the input, unit, and
outcome-blind pass conditions before the cluster run.

**6 · Carcinoma — a candidate resolvability diagnostic, not a fraction.**
`best4` returns exactly-zero mature fraction on **1,314 of 1,350
method-by-sample rows (97.3%)**; that is **641 of 675 tumours zero under both
methods**. Every row is marked `estimability='estimated'` with an empty reason
field. **Never present it as a biological mature-cell fraction.**

**The denominator audit ran 2026-09-07 (`a51f4af`) and the number above was
being quoted wrongly, here included.** 1,350 is method-by-sample rows, not
samples. The correct statements: **1,314/1,350 rows (97.3%)** and **641/675
tumours zero under both methods (95.0%)**, with 673/675 zero under at least one
— so exactly two tumours are non-zero under both. `docs/stage4_carcinoma_audit.md`
records it as a resolvability diagnostic. The only supported carcinoma
conclusion is "not estimable at this resolution". **Closed.**

**Policy, not an invariant.** *ML may provide segmentation, geometry,
anatomical measurements, or a prediction against an independently measured
endpoint; biological claims remain patient-level estimates with explicit
abstention rules.* That belongs in `CONTRIBUTING`, because the standing
invariants are enforced by assertions and this one is not testable. The
**testable** half — *no transcript-derived label may define a population and
then be used to claim the state of its own defining programme* — **is now
`CLAUDE.md` invariant 11**, merged 2026-09-07 by invariant 3's route: a
`shared/...` PR, W2's approval, and a guard with a forcing input for every
branch. `src/common/label_provenance.py`; both reference jobs declare their
provenance, call the guard before reading anything, and write the declaration
into their result sidecars. See `docs/invariant_11_proposal.md` RESULT.

### Still open, ranked

5. **The write-up.** Six results from 2026-09-06 are in neither paper and the
   WMHS deadline is **15 September 2026**. This outranks everything below it.
   `docs/HANDOFF.md` §6f.
6. ~~**B1's feasibility gate.**~~ **RAN 2026-09-07 — NOT LICENSED.**
   `docs/prereg_becker_replication.md` RESULT, table
   `results/2026-09-07_2305f23/`. The gate is **CLEARED BY DEPTH — NOT
   LICENSED**: the arm-threaded read passes only because Becker's depth exceeds
   Chen's, and §3's physical premise — that these genes are visible in nuclei at
   in-situ depth — is refuted rather than supported. **Do not queue it.**
7. ~~**C1's panel lookup.**~~ **CLOSED 2026-09-07 on both vendors (`a51f4af`),
   and the answer is sharper than the CosMx-only reading.** The panels are
   **complementary and disjoint** on what this design needs:
   Xenium Colon v1 (322) carries GUCA2A and MS4A12 and **no control at all**;
   CosMx 6K (6415) carries KRT8/EPCAM/CDX2 and **neither target**; Xenium Prime
   5K carries EPCAM/CDX2 and no KRT8; CosMx 1K carries KRT8/EPCAM only.
   **No stock panel carries a target and a control together**, so C1 needs
   custom probes or a different design. Parses verified clean (322/322,
   5001/5001) so the absences are real, and the panel files are checksum-pinned
   in `data/manifest.csv`. **Do not queue it.**
8. ~~**Zheng's gradient**~~ **RAN 2026-09-08 — descriptive only.** All three
   normal→polyp GUCA2A trajectories fall; CDX2 is mixed, and stage-specific
   depth/QC differences preclude a pooled or intrinsic reading. The committed
   table is `results/2026-09-08_f6b22f6/`. **Do not queue it.**
9. ~~Tier 3's heterogeneity explanation.~~ **RAN 2026-09-07 — the heterogeneity
   is one n=3 study and an uncalibrated ceiling, not a subtype.** It did not
   strengthen the result in hand; it changed the reason behind it, which §2 of
   `docs/HANDOFF.md` now carries. See §6k. **Do not queue it.**

10. ~~**The remaining half of §6k.**~~ **CLOSED 2026-09-08.** The fixed 75%
   I² ceiling is retired. `src/harness/meta.py` now requires the
   patient-count-matched null tail probability at fixed α=0.05, and the clean
   rerun is `results/2026-09-08_cf5b419/`. `docs/HANDOFF.md` §6k records the
   approval/process history and the floor-specific verdict causes.

~~Tier 2 only after the MLH1 control says whether the instrument can see
silencing at all.~~ **That gate is now permanently open, and not the way anyone
wanted.** The control cannot be run on available data (above), so waiting for it
is waiting for nothing. The reasoning it encoded still stands and now points
elsewhere: **do not spend weeks feeding a transcript instrument whose
sensitivity cannot be established.** Prefer the avenues that do not depend on it
— A, which needs no premise and no sensitivity, and B1/C2, which replace the
label rather than trusting it.

### Plan, 2026-09-07 — the corrected order of work

The roadmap amendment above records **what each item is**. This section records
**what to do next, in what order**, and three of its four entries correct a next
step the amendment states differently. Where they conflict, this section is the
later decision.

**The write-up still outranks all of it** ("Still open" item 5, WMHS deadline
15 September 2026). What follows is the technical queue behind it, not ahead of
it.

**Plan 1 · The meta layer — COMPLETE 2026-09-08.** The unused
`se_from_interval` helper is gone. `db487f0` replaces the fixed 75% I² ceiling
with the pre-specified patient-count-matched null at α = 0.05; a premise verdict
now refuses an omitted or invalid calibration rather than treating it as
homogeneity. `b0ef250` closes the remaining summary defect: every reported
floor's `verdict_cause` is compared, so matching n≥3 and n≥6 labels cannot hide
a different n≥4 route.

The clean-tree re-run is `results/2026-09-08_cf5b419/`. **KRT8 remains
UNRESOLVED at all three floors, but is floor-unstable by cause:** calibrated
heterogeneity at n≥3 (p=0.0156) and n≥6 (p=0.0176), versus a homogeneous
tolerance-straddle at n≥4 (p=0.0898). ACTB holds at every floor. This is a
resolvability result, not permission to select a floor or claim that controls
hold. The result sidecar is clean (`git_dirty: false`) and pins both code
implementation history and the 200,000-draw null simulation.

**Outstanding process debt, recorded 2026-09-08.** `src/harness/` is W2-owned
and CONTRIBUTING §2-3 route every change to it through a PR; §6k of
`docs/HANDOFF.md` named the requirement specifically — "a W2 PR with two
approvals". `db487f0`, `b0ef250` and `32ef54d` landed directly on
`submission/competitor-bench` with neither. **The code is verified and the
numbers reproduce; what is missing is the review, and it is still owed.** This
paragraph exists because the completion note that replaced the plan dropped the
requirement rather than discharging it.

**Item 4 · Lesion-level Wnt — CLOSED 2026-09-08.** The outcome-blind detection
gate ran on the four paired Becker donors and returned **NO SUBSTRATE AT THIS
snRNA RESOLUTION**: NOTUM and TCF7 failed one or both pre-specified floors in
three donors. `results/2026-09-08_a96e3b3/` and
`docs/prereg_becker_lesion_wnt_detection_gate.md` record the result. No
lesion-level Wnt pre-specification, reduced signature, or model is licensed.

**Item 2 / A2 · Authenticate before rewriting the request.** The repo records
that the three Synapse entities **refuse anonymous read**. An anonymous-read
refusal is not evidence of an access gate — it is the expected response to an
unauthenticated client, and it does not distinguish "controlled access" from
"log in." Before any further work on the request text in
`docs/a2_synapse_access_request.md`:

```
synapse login
python -m src.reference.jobs.a2_synapse_inventory --no-write
```

Only if an **authenticated** read also refuses does the access-request path
become real. Item 2's substantive ruling is untouched either way: no biological
analysis proceeds from the pixel product alone.

**Order, and why.** Plan 1 and item 4 are complete. The A2 login check remains
external access work, not a biological analysis; it is not needed to interpret
any result above.

**Unchanged by this plan:** 2b's four gate conditions and its NO SUBSTRATE exit;
AD-versus-SSL as a data hunt with no substrate in hand; and the closures of
Crowell, D2 and carcinoma.

### Amendment, 2026-09-10 — item 3 is unblocked, and §B's closure has a contested reason

**Supersedes item 3 above and the `Still open` table where they conflict.** The
external data hunt is `docs/DATA_HUNT_2026-09-10.md`; it carries the corrections
made to its own first draft.

**Item 3 no longer needs Synapse to be sized.** The Chen deposit is served open
access by cBioPortal as `crc_hta11_htan_2021`: an hg19 mutation MAF for 30
biospecimens and sample-level pathologist labels (`POLYP_TYPE`, `POLYP_SUBTYPE`,
`ATYPIA`, `ADVANCED`), keyed by HTAN biospecimen ID and therefore joinable
specimen-exactly to `results/2026-09-09_0bf9734/`. Joined against the 44-patient
lineage universe: **`AD` 13 patients versus `SER` 9** after excluding the one
patient carrying both, against **5 versus 4** for the genotype arms.

**RAN 2026-09-10 — locked at `3a6d446`, closed inconclusive. See the RESULT
section of [`prereg_chen_lesion_subtype.md`](prereg_chen_lesion_subtype.md) and
`docs/HANDOFF.md` §6n; the paragraph below is the pre-run framing and is kept
for the reasoning in §4, which held.**

**The outcome: weighting-unstable, no target-specific conclusion.** GUCA2A's
intrinsic difference contains zero on all three weightings. CDX2 excludes zero
under `tumour` only, on a lower bound 0.02 from zero, and the recorded status is
`DISAGREES ACROSS WEIGHTINGS` — §5 fixes all three with none primary and
registers no majority-vote rule, so two-against-one is not a verdict. Of the 15
pair contrasts every one that excludes zero involves `ACTB` or `KRT8`, while
`ACTB − KRT8` itself contains zero: the controls separate the arms more strongly
than the targets, which does not support a target-programme-specific reading.
13 versus 9 patients was named in §7 as the reason to expect a weak result.

**Family G is closed on the floor, not on access.** 5 and 4 estimable against
§7's floor of 8 — NOT ESTIMABLE, no interval, from the **open** cBioPortal MAF.
So `wes_subtype_plan.md` does **not** stay alive: the Synapse certification, the
`syn23520239` 403 and the Level 3 download ACL were never the binding constraint
on this analysis, and no access decision would have changed the arm sizes.
**Do not queue the Synapse chase for this question.**

**Original pre-run framing follows.** Its §4 is
the load-bearing part and it changes what item 3 can claim: **the intrinsic term
is a test and the compositional term is not.** A pathologist's AD/SER call reads
crypt architecture and surface maturation, which is the compositional endpoint
in another modality — invariant 11 passes on provenance and cannot see that.
The genotype label is clean in kind but sits at 5 versus 4 patients. So
`wes_subtype_plan.md` stays alive for the compositional term, descriptively,
and the pathology label carries the intrinsic one. Seven gate conditions are
unmet; §9 lists them. *(Superseded: six were discharged and the seventh waived
on 2026-09-10; the genotype arm then closed on the floor.)*

**§B's Becker closure rests on three reasons and one of them is contested.**
Blomain et al., *Cancer Biol Ther* 2020 ([PMC7515455](https://pmc.ncbi.nlm.nih.gov/articles/PMC7515455))
report guanylin protein **and** mRNA maintained under APC heterozygosity and
lost only after LOH — FAP normal mucosa expressing, adjacent polyps not, by
immunofluorescence and qRT-PCR, with *Apc*^min/+ mice agreeing by RNA-seq. That
contradicts §B's "structural to FAP" reference-arm argument. Their human arm is
n=2 and ours is n=2 healthy donors, so this is **not resolvable either way** and
B1 is **not** reopened: the four-paired-donor reason and the detection gate at
`results/2026-09-07_2305f23/` are untouched. Recorded so the closure rests on
the two reasons that survive.

### ML amendment, 2026-09-08/09 — engineering branch, separate from the biology

The HTA11 H&E molecular-prediction gate remains closed. Public polyp morphology
labels cannot replace its specimen-exact, independently assayed molecular
endpoint. They can support an explicitly scoped engineering benchmark, but never
a biological project result.

**MHIST source screen — CLOSED FOR ITEM 2b, 2026-09-09.** The durable
metadata-only inventory is `results/2026-09-08_ccfa39c/mhist_feasibility.parquet`.
Its 3,152 224x224 PNG tiles and four-field annotation CSV have neither
full-resolution/WSI information, specimen-exact molecular data, a molecular
endpoint, nor a patient/slide key. The supplied HP-versus-SSA majority label is
read from the classified tile, so it is diagnostic-label imitation. MHIST
therefore fails all four of item 2b's fixed conditions and is **NO SUBSTRATE FOR
THE 2b MOLECULAR-PREDICTION ESTIMAND**. This does not close HTAN's distinct,
still-gated H&E--VCF candidate route. The RUA also prohibits derivative works,
so even a separate tile-level encoder sanity check remains license-pending until
the publisher authorizes embeddings, model fitting, and publication of aggregate
derived outputs.

UniToPatho and the abstention survey remain optional engineering work, not
project-result routes or priorities. The full source record and its limits are
in `docs/ml_abstention_roadmap.md`.

**Release-7 early-lesion source — CLOSED AS INDEPENDENT REPLICATION,
2026-09-09.** The clean-tree source-identity audit is
`results/2026-09-09_e019955/release7_early_lesion_identity.parquet`: all 55
fixed Release-7 participants overlap the cached 106-participant Chen universe.
Release-7 therefore contributes no independent replication; the earlier
"provisionally set aside" wording below is superseded. This does not assert
anything about its expression values or license a same-cohort reanalysis.

**Grady spatial release watch — not an analysis.** `GSE342139` (Visium) and
`GSE341852` (Xenium) are named as reviewer-token deposits whose custom-50-gene
panel has not been publicly verified. Recheck only on public release. A Xenium
analysis becomes a fresh data-gate question only if the released panel contains
the pre-specified targets and controls, patient/section mapping is available,
and its new analysis is pre-registered. Until then it is not a substrate, an
ML task, or evidence that the current single-cohort cap can be lifted.

**Two corrections to the H&E header audit, 2026-09-08.** The 12 parsed rows in
`results/2026-09-08_3cc8f73/` all reported *exactly* 352.777778 microns per
pixel, which is 25400/72 — a TIFF writer's 72-dpi default, not a measurement.
Taken literally it implies slides 2.95 to 14.45 metres wide. `he_tiff_header.py`
now classifies a writer default or an implausible scale and returns **`None`**
for the scale rather than the placeholder, on invariant 1's reasoning. **The
committed CRDC table predates that classifier**; its `x_microns_per_pixel`
column must be read as unresolved, not as 352.78, until the job is re-run.
Separately, the deposit's own imaging metadata records `PhysicalSizeX` 0.25 µm
at 40x for all 18 crosswalk biospecimens, so the images are plausibly
full-resolution and the header simply never carried the fact. **That evidence
has not been wired into the gate and `image_resolution_metadata_available`
remains `False`** — promoting it is a gate decision, not a defect fix.

The three HTAN H&E result tables also predate their jobs' invariant 11
declarations, which were added 2026-09-08. The jobs now declare and guard; the
committed sidecars do not carry the declaration and are regenerated on the next
run against the source TSVs, which are not currently on disk.

The preceding Release-7 provisional status is superseded by the clean-tree
2026-09-09 artifact named above. It remains neither an H5AD-read license nor a
same-cohort reanalysis authorization.
