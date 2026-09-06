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
`results/2026-09-05_3a1af9f/icbi_adenoma.parquet`:

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
indistinguishable, while CDX2 (0.791) sits well above both. The same tiered
pattern the detection reading found — terminal differentiation down, intestinal
identity retained — and the same limit: the decomposition would separate GUCA2A
from housekeeping and from CDX2, and not from MS4A12.

**Which makes the comparator choice the whole design.** Run against housekeeping
alone it will look like a clean gene-specific result. The identity markers are
what stop that, and they must be scored in the same run rather than added
afterwards. Identifiable is not the same as gene-specific.

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
(`docs/prereg_g2_mlh1.md`) died on the detection floor at CP10K means of 0.04.
Detection at ≥1 UMI is more sensitive than a mean, which is the reason to
retry — but a feasibility check on mature-cell counts per stratum comes before
the reading, not after.

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

1. **The MLH1 positive control.** It is the cheapest item here and the only one
   that can change how every previous result is read.
2. **The decomposition at `best4` on Chen_2021** (1a). Same run, and the algebra
   has already been checked to not collapse.
3. **Zheng's gradient** (1c), descriptive, alongside.
4. Tier 3's heterogeneity explanation — cheap, and it strengthens a result
   already in hand.

Tier 2 only after the MLH1 control says whether the instrument can see silencing
at all. Spending weeks on a new assay to feed an instrument of unknown
sensitivity is the wrong order.
