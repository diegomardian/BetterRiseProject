# Pre-registration — can this instrument see silencing in cells that are present?

**Written 2026-09-10, before any Replogle file was downloaded or any count
read.** The abundance and depth anchors in §3 come from the project's own
committed inputs (`results/2026-09-06_765eb29/`), not from Replogle.

## 1 · The question, and why it is not a biology question

`HANDOFF` §6g records the gap: **absence of transcript is not silencing unless
the instrument can detect a gene that is off in cells that are still present.**
The MLH1 positive control was to settle it and cannot be run on any available
data, so every intrinsic reading in this project rests on an unmeasured
sensitivity.

Replogle et al. 2022 genome-scale Perturb-seq supplies what no colorectal
dataset does: **cells retained by construction**, DNA-level silencing, measured
on-target knockdown, and thousands of baseline abundances. It is not a
colorectal result and is not offered as one. It measures **the instrument**.

**Deliverable:** detection power as a function of baseline abundance,
sequencing depth, cell count and knockdown magnitude — a curve, not a pass/fail.

## 2 · Perturbed versus control is defined by guide assignment, never by expression

**This is the load-bearing design decision.** Defining a perturbed cell as one
with low expression of the targeted gene, then asking whether the targeted gene
reads low, is the circularity invariant 11 exists to forbid — and it would
manufacture near-perfect power.

`label_provenance`: **CRISPRi guide identity**, a DNA-level assignment read
from the guide capture, modality `genotype`.
`claim_provenance`: transcript detection of the targeted gene, modality
`transcript`.

Non-targeting-control cells are the comparison arm. A cell with no confident
guide assignment is **dropped and counted**, never assigned to either arm.

## 3 · Abundance and depth bins, fixed now, anchored to this project's regimes

Fixed from `results/2026-09-06_765eb29/`, `lineage` rung, normal arm, across
the 44 adenoma patients — **the distribution, not its maximum**:

| gene | 10% | 25% | **median** | 75% | 90% | max |
|---|---|---|---|---|---|---|
| GUCA2A | 1.59 | 3.03 | **5.48** | 9.53 | 15.06 | 38.68 |
| MS4A12 | 0.19 | 0.61 | **1.06** | 3.42 | 5.94 | 7.96 |
| CDX2 | 2.64 | 3.21 | **3.87** | 4.73 | 5.48 | 6.21 |
| ACTB | 10.40 | 13.58 | **16.62** | 22.11 | 26.93 | 30.89 |

**A correction the anchors required.** `DATA_HUNT_2026-09-10.md` §278 states
"GUCA2A sits at 38.7 CP10K in the adenoma normal arm". **38.682 is the maximum
over 44 patients; the median is 5.476 and the mean 7.797.** Anchoring here at
38.7 would have targeted a regime roughly seven times more favourable than the
one the project operates in, and would have made the instrument look better than
it is by construction. The hunt's figure is corrected there.

**Pre-committed abundance bins**, in CP10K, chosen to bracket those regimes:
`[0, 0.1)`, `[0.1, 0.5)`, `[0.5, 1.5)`, `[1.5, 4)`, `[4, 8)`, `[8, 20)`,
`[20, ∞)`. The MS4A12 and GUCA2A medians fall in distinct bins; MLH1's reported
0.039 falls in the first.

**Pre-committed depth points**, by controlled downsampling of Replogle UMIs to
the project's own **10th, 25th, 50th, 75th and 90th percentiles** of
`depth_normal` — 3,518 / 4,678 / 5,664 / 8,114 / 10,008, used as **3,500 /
4,700 / 5,700 / 8,100 / 10,000** median UMIs per cell. These are percentiles
across the five-point range, **not quartiles**. Downsampling is binomial
thinning at a fixed seed.

**Pre-committed cell counts per arm**: **40 / 85 / 255 / 535**, the project's
`n_normal` **10th, 25th, 50th and 75th percentiles** (41 / 84 / 255 / 533).
Again percentiles, not quartiles.

**No target gene is selected after seeing its result.** Every Replogle gene
meeting the §4 support requirement enters its abundance bin. Bins are reported
whether or not they are favourable, and a bin with too few genes is reported as
**not estimable**, never merged into a neighbour.

## 4 · The positive-control gate, stated so it can fail

**Runs first among post-lock artifacts** (§7), on undownsampled data at full
cell count. Everything below is numeric because a gate whose pass condition is
described rather than defined cannot fail, which is this project's own thesis.

**Dataset.** Replogle 2022 **K562 essential-scale** (day 6), raw single-cell
counts from the Figshare+ deposit.

**Target selection, structural and outcome-blind.** Every targeted gene
qualifies that has (a) its own transcript detected in **≥ 20%** of
non-targeting-control cells, (b) **≥ 50** cells carrying a confident guide
assignment for it, and (c) **≥ 2** distinct guides. No gene is named in advance
and none is added or dropped after its result is seen.

**Statistic.** Per qualifying target, the **on-target log2 fold change** of that
target's own transcript: CP10K-normalised pseudobulk of its guide-assigned cells
against the non-targeting-control pool.

**Pass condition, both clauses required.**

1. **Median on-target log2FC ≤ −1.0** across qualifying targets — at least 50%
   knockdown, well inside CRISPRi's usual 60–90%.
2. **≥ 80%** of qualifying targets individually show log2FC **< 0**.

**Fail consequence.** The reader is wrong and **no curve is reported**. A
sensitivity curve computed by a broken reader is a claim about this code, not
about single-cell RNA.

## 4a · "Power", defined

**Independent unit: the target gene.** Guides against one target are not
independent replicates and are pooled before any test.

**Per cell of the grid** — abundance bin × depth point × cell count × knockdown
bin — the procedure is, for each qualifying target in that cell:

1. Draw the specified number of cells from the guide-assigned arm and the same
   number from the non-targeting pool, without replacement.
2. Binomial-thin both to the depth point.
3. Compute the detection statistic and its **Student-t interval**
   (`src/reference/interval_calibration.py`, the project's estimator, and the
   §3a reason for not using the percentile bootstrap).
4. **Detected** = the interval excludes zero **in the silencing direction**. An
   interval excluding zero upward is not a detection and is counted separately.

**Repeats: 100** draws per target per grid cell, fixed seed.

**Power** = the fraction of the 100 repeats in which silencing is detected,
reported **per target** and summarised across targets in the cell by median and
interquartile range. A grid cell with **fewer than 5 qualifying targets** is
reported **not estimable** and is never merged into a neighbouring cell.

**The false-positive rate is measured in the same run**, by applying the whole
procedure to **two disjoint non-targeting-control groups**, where the true
effect is zero by construction. A power figure without it is uninterpretable —
`interval_calibration` already refuses that pairing elsewhere. The nominal rate
is 5%; the measured rate is reported beside every power estimate and **no cell
whose measured false-positive rate exceeds 10% may be read as a power result.**

**The decisive-failure condition, fixed now** (§5's first row): power **≤ 0.20**
for the abundance bin containing GUCA2A's median (5.48 CP10K → bin `[4, 8)`),
at the project's median depth (5,700) and median cell count (255), at knockdown
**≥ 70%**, with the measured false-positive rate at or below 10%.

## 5 · Interpretation is asymmetric, and this is fixed before the run

| outcome | what may be concluded |
|---|---|
| **Silencing is not detectable at GUCA2A's regime even under Replogle's favourable conditions** | **Decisive for the pipeline this project specifies, under optimistic whole-cell conditions.** The claim is "*this pipeline* cannot resolve silencing at this abundance *even given whole-cell Perturb-seq data*" — never the unqualified "the instrument cannot resolve silencing". Whole-cell Perturb-seq and tissue nuclei are **not strictly nested conditions**, so a whole-cell failure bounds the nuclei case by argument, not by containment. §6g's UNINTERPRETABLE becomes a bounded statement rather than an absence. |
| Silencing **is** detectable at that regime | **Generic capacity only.** It does **not** show the colorectal null is biological, and does not license reading any intrinsic term as silencing. Cell lines at high depth are the easy case. |

**Replogle is an optimistic bound and only its failure direction is
informative.** This is the same asymmetry `DATA_HUNT_2026-09-10.md` applies to
methylation, and it is written here so it is not softened afterwards.

## 6 · What this cannot say, at any outcome

- Nothing about colorectal biology. K562 and RPE1 are cell lines.
- Nothing about **nuclei**. Replogle is whole-cell; Becker and Crowell are
  nuclei, and `prereg_becker_replication` already closed on that distinction.
- Nothing about ambient contamination in tissue, which cell lines do not have.
- **CRISPRi knockdown is partial, typically 60–90%; promoter hypermethylation is
  complete.** The knockdown-magnitude axis is reported for this reason, and a
  curve at 70% knockdown is not a curve for complete silencing.
- It does not repair any existing result. No committed estimate changes.

## 7 · Sequencing — four metadata-only gates before lock, then the control gate

**An earlier draft of this section required the §4 positive-control gate before
lock. That was a contradiction and is corrected here:** the gate reads counts,
and any result read before lock can shape the design it is supposed to test.
The gates below are **metadata-only** — schema, checksums, licence, declarations
— and none of them reads an expression value.

**Pre-lock, metadata only:**

1. **Terms and licence** for the Replogle Figshare+ deposit reviewed and
   recorded, in the form `docs/chen_subtype_data_gate.md` used.
2. **Snapshot pinned** — file identifiers and sha256 in `data/manifest.csv`.
   Checksums are computed over bytes; no count is parsed.
3. **Guide-assignment field confirmed at source** from the deposit's own schema
   and documentation, not inferred from a column name, and the
   no-confident-guide rule implemented against that schema.
4. **`src/common/label_provenance.py` declaration written and passing** for the
   guide/transcript pairing.

**Then lock**, immediately on those four passing.

**Post-lock, in order:**

5. **The §4 positive-control gate**, as the **first** post-lock artifact. It may
   fail, and if it does the pipeline is wrong and nothing further runs.
6. **The §4a sensitivity curve**, **only** if 5 passes.

## 8 · Standing

**Unlocked.** When gates 1–4 of §7 pass, this document is locked by commit, and
the bins in §3, the gate in §4, the power definition in §4a and the asymmetry in
§5 are frozen as written. The control gate and the curve then run against a
frozen design.

**It does not outrank the write-up.** The WMHS deadline is 15 September 2026.
