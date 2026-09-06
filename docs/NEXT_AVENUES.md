# What is left to try, reviewed against the data

**Written 2026-09-05.** Every claim below was checked against the cached obs or a
committed table before being written down. Where a proposal did not survive that
check, it says so and why — a list of avenues is only useful if the dead ones
are marked.

Read `docs/HANDOFF.md` §2 first for what is already established.

---

## Tier 1 — in hand, and one of them is bigger than path C

### 1a. The decomposition at `best4` on Chen_2021. **VERIFIED, do this.**

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

**C1 (segmented spatial)** is the standing Tier 2 ranking unchanged, and the
Visium caveat it carries is already this repo's (`docs/HANDOFF.md` §6e: a spot
is still a mixture). Ruling Pelka's own GeoMx out as a primary is right —
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
2. **B1**, raised — the frozen axes file's own week-13 substrate, unfetched,
   and the only thing here that answers the circularity objection outright.
   Gate on panel detection under a nuclear protocol.
3. **C2**, on CDX2 rather than GUCA2A, with the 450k already in hand.
4. C1, B3, B2 as supplementary.

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

### Explaining the 13-study heterogeneity
KRT8's I² of 87.6% is the reason the meta premise did not resolve. If that
disagreement tracks iCMS subtype — and Joanito 2022, the iCMS paper, is one of
the 13 — then the negative gets an explanation rather than staying a shrug.
Cheap: the per-study estimates are already committed.

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
2. **The decomposition at `best4` on Chen_2021** (1a). **NEXT after MLH1.** The
   algebra has been re-checked and does not collapse, and estimability at
   `best4` is 20/20 and 19/20 — better than at `lineage`. It is NOT a re-read
   of committed tables: `frac_mature_*` was never emitted, so it needs the
   adenoma job extended and a cluster run. Report the Student-t interval, not
   the percentile bootstrap (7.1% at n=20).
3. **Zheng's gradient** (1c), descriptive, alongside.
4. Tier 3's heterogeneity explanation — cheap, and it strengthens a result
   already in hand.

~~Tier 2 only after the MLH1 control says whether the instrument can see
silencing at all.~~ **That gate is now permanently open, and not the way anyone
wanted.** The control cannot be run on available data (above), so waiting for it
is waiting for nothing. The reasoning it encoded still stands and now points
elsewhere: **do not spend weeks feeding a transcript instrument whose
sensitivity cannot be established.** Prefer the avenues that do not depend on it
— A, which needs no premise and no sensitivity, and B1/C2, which replace the
label rather than trusting it.
