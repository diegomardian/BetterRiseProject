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

### 1c. Zheng_2022's within-patient gradient. **Correct as stated, n = 3.**

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

| open | what it needs |
|---|---|
| **The write-up** | Eight days. ~15 results across 2026-09-06/07 in neither paper, and `neurips_2026.sty` still unvendored so `./build.sh` has never run. **This outranks everything.** |
| **D2 survival** | `docs/prereg_d2_survival.md` is **proposed, not locked**. Its gate PASSED (n=624, IQR 2.96). §6a lists three decisions it leaves open. Needs the lock, then implementation. |
| **Invariant 11** | `docs/invariant_11_proposal.md` drafted with a guard design, unimplemented. Needs a PR and two approvals per invariant 3. |

**Closed on 2026-09-07 and not worth reopening:** Crowell (7 of 7, §6l),
Becker/B1 (not licensed), C1's panel lookup (no stock panel carries a target and
a control together), the meta-layer weights (§6k), and the carcinoma `best4`
denominator.

**The one avenue with an unexplored substrate** is lesion-level Wnt — item 4
below — and its audit came back thin: 2 of 6 Crowell blocks have multiple
adenoma sub-regions, and Becker's per-donor lesion counts were never committed.

### Roadmap amendment, 2026-09-07 — six items, after the Crowell run

Supersedes the corresponding entries above where they conflict. Each item names
what it is *not* as tightly as what it is, because five of the six were
originally written one class stronger than the data supports.

**1 · Crowell — retry `242`, and report six of seven if it stays blocked.**
`docs/prereg_crowell_multisection.md` RESULT. The result is a **patient-level
target-versus-epithelial-domain pattern**: GUCA2A −1.378 [−2.046, −0.709] and
MS4A12 −1.272 [−2.053, −0.492] over six blocks, CDX2 and EPCAM not clearing
zero. It is **not** evidence separating compositional from intrinsic loss
(Amendment 2 — the only control-role gene is an epithelial keratin) and **not** a
per-cell silencing result (feasibility §7). Below-floor observations are bounds,
not effect-size points (Amendment 4). `242.h5ad` returned HTTP 504 from every
Zenodo endpoint on 2026-09-07; §3's order was followed exactly, so a missing
seventh is a **named server failure, not a stopping rule**.

**2 · A2 / ML — contingent infrastructure only.** The public MxIF product is 42
pixel NPZ regions with **no cell segmentation and no cell-by-marker table**; the
three Synapse entities refuse anonymous read. **No biological analysis proceeds
from it alone.** Crypt-position assignment needs independently annotated
anatomical ground truth, patient-held-out evaluation, and a
morphology-or-geometry-derived label — and it still does not touch the per-cell
false-negative problem (`docs/HANDOFF.md` §6g).

**3 · AD versus SSL — a DATA HUNT, not a pre-registrable analysis.**
**No transcriptomic substrate exists in hand.** The ICBI atlas carries `polyp`
with no subtype, from two studies only; Crowell is **entirely TVA**. The sole
AD/SSL source is A2's pixel inventory — `config/a2_mxif_regions.csv`, 7 AD and 8
SSL patients, of which **only 8 have a labelled normal/tumour pair: 3 AD against
5 SSL** — and item 2 rules that product out for biological claims. At 3 patients
the AD arm sits exactly on `MIN_STUDIES` and 8.07× avenue A's width.
**Requirement: a cohort with paired reference, both lesion types, subtype
metadata, and a patient-level n that is not 3.**

**4 · Lesion-level Wnt — feasibility-gated, and the gate is partly answered.**
D1 was within-patient and explicitly conceded between-lesion variation. The
substrate now exists but is thin: **Crowell gives 2 of 6 blocks with multiple
adenoma sub-regions** (`110` and `221`, three each — audited from
`crowell_subdomains_exploratory.parquet`). Becker has 43 polyps across 12
donors, seven of them with ≥2, but **only four donors carry a reference arm and
those per-donor counts were never committed** — they exist only in a
`--inspect` log and need the 8.8 KB series matrix re-read to become an artifact.
Sub-regions and lesions are **not independent patients**; any analysis is
patient-clustered or mixed-model, with the Wnt score and adjustment set fixed
first, and read as association. Crowell's sub-domains stay `exploratory` under
its own pre-registration and cannot be borrowed.

**5 · D2 survival — a feasibility gate before the pre-specification.**
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
field. **Never present it as a biological mature-cell fraction.** The source
artifact's denominator and semantics must be audited before quoting either
number; the only supported carcinoma conclusion is "not estimable at this
resolution".

**Policy, not an invariant.** *ML produces segmentation, geometry or anatomical
measurements; biological claims are patient-level estimates with explicit
abstention rules.* That belongs in `CONTRIBUTING`, because the standing
invariants are enforced by assertions and this one is not testable. The
**testable** half — *no transcript-derived label may define a population and
then be used to claim the state of its own defining programme* — is a candidate
**invariant 11**, and it reaches `CLAUDE.md` only through invariant 3's route: a
PR, two approvals, and a guard that can fail. It is already applied in practice
(Becker RESULT, Crowell multisection Amendment 1 §2) and unenforced in code.

### Still open, ranked

5. **The write-up.** Six results from 2026-09-06 are in neither paper and the
   WMHS deadline is **15 September 2026**. This outranks everything below it.
   `docs/HANDOFF.md` §6f.
6. **B1's feasibility gate** — as **C1 evidence only**. Its paired cohort is
   four donors (Becker Amendment 2) so it cannot replicate avenue A, but the
   GUCA2A/MS4A12-in-nuclei detection table bears on whether in-situ platforms
   can see these genes at all. Data is on disk; `--inspect` is the next command.
7. **C1's panel lookup** — **HALF DONE 2026-09-07.** CosMx answered and it is a
   clean negative: **GUCA2A is on neither the 1K nor the 6K panel**, nor is
   ACTB. C1 is not runnable on a stock CosMx panel; it needs custom probes for
   the target itself. **Xenium is still open** — 10x rate-limited every request
   — and it is the cheapest thing left on this list. Tier 1 §C above has the
   table and the two URLs.
8. **Zheng's gradient** (1c), descriptive, and **no interval may be reported
   from it**: at n=3 the percentile bootstrap is 25.1% and 0.372× the correct
   width.
9. ~~Tier 3's heterogeneity explanation.~~ **RAN 2026-09-07 — the heterogeneity
   is one n=3 study and an uncalibrated ceiling, not a subtype.** It did not
   strengthen the result in hand; it changed the reason behind it, which §2 of
   `docs/HANDOFF.md` now carries. See §6k. **Do not queue it.**

10. **The remaining half of §6k, and it is W2's not W1's.** `MAX_I_SQUARED = 0.75`
   is Higgins' rule of thumb with no null behind it: at k = 11 homogeneity alone
   gives a median I² of 0.270, at k = 6 it gives 0.000, and Cochran's Q rejects
   **32.5%** of the time at the committed patient counts against a nominal 5%.
   `src/reference/meta_calibration.py` computes the null; wiring it into
   `src/harness/meta.py` in place of the fixed ceiling is a **PR with two
   approvals** under CONTRIBUTING §2 and was deliberately not done here.

~~Tier 2 only after the MLH1 control says whether the instrument can see
silencing at all.~~ **That gate is now permanently open, and not the way anyone
wanted.** The control cannot be run on available data (above), so waiting for it
is waiting for nothing. The reasoning it encoded still stands and now points
elsewhere: **do not spend weeks feeding a transcript instrument whose
sensitivity cannot be established.** Prefer the avenues that do not depend on it
— A, which needs no premise and no sensitivity, and B1/C2, which replace the
label rather than trusting it.
