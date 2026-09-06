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
1. ~~**The MLH1 positive control.**~~ **BUILT AND PRE-REGISTERED**, 2026-09-06.
   Needs `qsub src/reference/jobs/mlh1_positive_control.sh` and nothing else.
   Still the only item here that can change how every previous result is read.
   The DiD correction recommended above did not survive sizing — see above.
2. **The decomposition at `best4` on Chen_2021** (1a). Same run, and the algebra
   has already been checked to not collapse.
3. **Zheng's gradient** (1c), descriptive, alongside.
4. Tier 3's heterogeneity explanation — cheap, and it strengthens a result
   already in hand.

Tier 2 only after the MLH1 control says whether the instrument can see silencing
at all. Spending weeks on a new assay to feed an instrument of unknown
sensitivity is the wrong order.
