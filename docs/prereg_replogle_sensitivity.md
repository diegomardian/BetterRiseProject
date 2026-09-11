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
the project's own quartiles: **3,500 / 4,700 / 5,700 / 8,100 / 10,000** median
UMIs per cell. Downsampling is binomial thinning at a fixed seed.

**Pre-committed cell counts per arm**: **50 / 85 / 255 / 535** cells, the
project's `n_normal` quartiles.

**No target gene is selected after seeing its result.** Every Replogle gene
meeting the §4 support requirement enters its abundance bin. Bins are reported
whether or not they are favourable, and a bin with too few genes is reported as
**not estimable**, never merged into a neighbour.

## 4 · Positive-control gate, before any curve

Under the **undownsampled** data, at full cell count, a known CRISPRi target
must show detectable on-target suppression against non-targeting controls, with
adequate guide and cell support. **If the gate fails, the pipeline is wrong and
no curve is reported** — a sensitivity curve computed with a broken reader would
be an instrument claim about this code, not about single-cell RNA.

The gate is outcome-blind with respect to the abundance curve: it asks only
whether known suppression is recovered where it must be.

## 5 · Interpretation is asymmetric, and this is fixed before the run

| outcome | what may be concluded |
|---|---|
| **Silencing is not detectable at GUCA2A's regime even under Replogle's favourable conditions** | **Decisive.** Supports "this instrument cannot resolve silencing at this abundance", and §6g's UNINTERPRETABLE becomes a bounded statement rather than an absence. |
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

## 7 · Gate conditions — all unmet, all required before lock

1. **Terms and licence** for the Replogle Figshare+ deposit reviewed and
   recorded, in the form `docs/chen_subtype_data_gate.md` used.
2. **Snapshot pinned** — file identifiers and sha256 in `data/manifest.csv`
   before any count is read.
3. **Guide-assignment field confirmed at source** from the deposit's own schema,
   not inferred from a column name, and the no-confident-guide rule implemented.
4. **`src/common/label_provenance.py` declaration written and passing** for the
   guide/transcript pairing, before anything is read.
5. **The positive-control gate of §4 run and passed**, recorded as its own
   artifact.

## 8 · Standing

**Unlocked.** When §7 is discharged this document is locked by commit, and the
bins in §3 and the asymmetry in §5 are frozen as written.

**It does not outrank the write-up.** The WMHS deadline is 15 September 2026.
