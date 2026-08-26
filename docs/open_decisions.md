# Open decisions

Things the scaffold could not decide on its own, recorded so they get decided
deliberately rather than discovered in week 9. Close each with a line in the
weekly meeting and a PR that updates this file.

Open questions of a scientific kind live in
[execution_plan.md §10](../execution_plan.md#10-open-questions). This file is
for decisions that block or shape code.

---

## CORRECTIONS — 2026-08-22

A verification pass over W1's own recommendations, against the other
workstreams' branches rather than against memory. Five were wrong. They are
corrected in place below; recorded here because several were circulated.

**Raised with the other workstreams as issues, not as pushes to their branches**
(CONTRIBUTING §2 puts `src/bulk/` and `src/estimator/` off-limits to W1):
[#7 for W3](https://github.com/diegomardian/BetterRiseProject/issues/7) ·
[#8 for W4](https://github.com/diegomardian/BetterRiseProject/issues/8) ·
[#9 for W2](https://github.com/diegomardian/BetterRiseProject/issues/9).

**1 · The gene index. W1 never emitted one into the repo.** The decision board
said "W1 emits it at `config/gene_index/gene_index_1.0.0.*`". That file does not
exist on any branch — it was written on the cluster and never committed. What
*does* exist is **`gene_index_0.9.0`, built by W3**, when the fallback written
into decision #2 fired as designed. It is correct: unversioned Ensembl key,
version in its own column, symbols mapped, 60,616 genes from the GDC gene model —
exactly decision #3. **W1 must not emit a competing 1.0.0.** See #2 below.

**2 · ~~The hg19 warning was overstated.~~ RETRACTED 2026-08-22 — the original
number was right.** This correction said "a silent 8% gene loss" was asserted
rather than measured and was "too pessimistic". W3 then measured it (`cc06981`):

```
W1 43,078 | W3 60,616 | both 39,236
W1-only 3,842 (8.9% of W1) | W3-only 21,380 (35.3% of W3)
```

**8.9%.** The estimate was accurate and the retraction of it was not. What *was*
right in this correction is the mechanism: the loss is GENCODE release drift
(v28 against v36), not the assembly, which is exactly what decision #3 predicted
and why keying on unversioned ENSG was correct. **All 23 panel genes are present
on both indexes**, so tiers A–D are safe whichever index is adopted.

W3 and W1 reached the same conclusion independently: **1.0.0 should be the
intersection (39,236 genes)**, with each arm keeping its full matrix for its own
work. #2 and #3 are arithmetic now, not opinion — they need ratifying, not
deciding.

**3 · W4's cut points are pooled, not per-sample.** #13 told W4 to switch "or the
term is zero". That is wrong. `classify_maturity` takes a quantile of the score
over the **whole input**, so Δ(mature fraction) is free to move. Per-sample
quantiles were W1's bug, not W4's. The real argument for reference-arm cuts is
different and is stated correctly in #13 now.

**4 · Emitting the unmatched 26 as `not_estimable` is not free.** Two things were
missed. The schema has **no reason field**, so "no normal arm" is indistinguishable
in the output from "too few mature cells", and those have opposite implications.
And `gate_g4_verdict` counts patients whose `n_cells_mature` is below threshold —
26 guaranteed-zero rows would push a **pre-committed gate criterion** toward
failure for a reason that has nothing to do with positivity. See #9.

**5 · Telling W4 to adopt a 50% mito cap was overreach.** 29.8% is the epithelial
median measured on **GSE178341**. W4's cohort is Lee. The transferable claim is
"20% is a lymphocyte number and colonic epithelium runs higher", not the number.
See #12.

---

## W1 status board — 2026-08-22

Ranked by what breaks if they stay open. Seven are live.

| # | Subject | State | Why it matters |
|---|---|---|---|
| **14** | Neither labelling axis is a clean maturity measure | OPEN | **Blocks composition.** Kappa read 2026-08-20: `stem_pole` **0.313** (fair), `best4` 0.03 (unusable), axis 2 negative by construction. Also found: `lineage` and `crypt_position` are the *same partition* on axis 1. Rerun the depth sweep — kappa was measured at its worst target. |
| **13** | W1 and W4 label cells differently | OPEN | Cohorts not comparable at the week-5 gate. W4's pooled cuts do *not* zero the term — corrected reasoning in #13. |
| **12** | Mitochondrial cap | ANSWERED for W1 (50) | W4 is on 20 (verified in `src/estimator/ingest.py`). Ask is *measure it on Lee*, not *adopt 50*. |
| **10** | Refined tier-B test | **DRAFTED** | Pre-registration written and committed: [docs/prereg_g2_mlh1.md](prereg_g2_mlh1.md). Needs team ratification, not more drafting. |
| **9** | The 26 unmatched patients | OPEN | 42% of the cohort. Emitting as `not_estimable` **can flip gate G4** — needs W2. Cohort artifact first. |
| **11** | Sorted samples | OPEN | Implemented, unratified. W4 must nod to the arm asymmetry. |
| **8** | No unfiltered droplets | ANSWERED | CellBender out; SoupX + DecontX in. W4 has the same exposure on Lee. |
| **2/3** | Shared gene index | **ANSWERED** | W3 measured the overlap: 39,236 genes in both, all 23 panel genes present. 1.0.0 = the intersection. Ratify, do not re-decide. |
| **16** | Ambient: measure not correct; exclude >10% | **PRE-COMMITTED** | Median 2.2%, but 9 of 84 samples above 10%. Threshold set before counting its cost. |
| **15** | CNV calling fails differentially by MMR status | **CONFIRMED** | 62-patient run: **15/15 MMRp separable, 15/20 MMRd**, monotone across four strata, 3x enrichment gap. Bias runs **along** the pre-registered contrast. |
| — | CNV reference design | NEEDS SIGN-OFF | Matched normal with 30% held out. **Corrected once** — see `src/reference/malignancy.py`. |

Not decisions, but outstanding and cheap: tell W3 the gene index is emitted at
`config/gene_index/gene_index_1.0.0.*` **and that this deposit is hg19 while TCGA
is GRCh38**; run `src/reference/jobs/pull_icbi_metadata.py`, which has never been
run and sizes the plate-based subset — G1's fallback after #8, and the only route
to a deeper-sequenced axis 1 for #14.

---

## 1 · MUC2 and TFF3 are both targets and labels — OPEN

**Raised:** scaffold, 2026-08-15 · **Owner:** W1 + W2 · **Needed by:** week 3
(W1 builds labels in weeks 3–4)

MUC2 and TFF3 sit in **tier E** (targets) and in **labelling axis 2**
(opposite lineage: MUC2, TFF3, SPDEF, ITLN1). Invariant 2 says target genes
never appear in labels or the reference matrix. Read across the whole panel,
axis 2 loses half its markers and the "agreement across structurally different
axes" argument weakens considerably.

The scaffold takes the **narrow** reading, without deciding the question:
`build_signature()` excludes the target set *for the run in question*, not the
whole panel. Consequence as implemented — **a run testing MUC2 or TFF3 may not
use axis 2.** `tests/test_freeze.py` pins the collision set at exactly
`{MUC2, TFF3}` so a future panel or axis edit that adds a new one fails loudly.

Options:

| | Approach | Cost |
|---|---|---|
| a | Keep the narrow reading. MUC2/TFF3 results carry axes 1, 3 only. | Two tier-E genes get thinner axis coverage. Cheapest. |
| b | Drop MUC2 and TFF3 from tier E. | Panel edit — 2 approvals. Tier E is exploratory, so the loss is small. |
| c | Rebuild axis 2 on SPDEF, ITLN1 and non-panel goblet markers. | Axis edit — 2 approvals, and axes were frozen first on purpose. |

**Recommendation: (a) now, revisit at the gate.** It costs nothing today and
the gate will show whether axis 2 is carrying weight.

---

## 2 · Who owns the shared gene index — FALLBACK FIRED, W3 BUILT IT

**Raised:** scaffold, 2026-08-15 · **Acted on:** W3, 2026-08-17 ·
**Owner:** W1 + W3 · **Confirm at:** the next weekly

W1 emits S matrices on a fixed gene index; W3 emits bulk on the same index;
integration is a join. The plan does not say who *produces* the index, and both
need it in week 1. See [config/gene_index/README.md](../config/gene_index/README.md).

**Superseded by events — see CORRECTIONS #1.** The fallback fired. W1's index
was written on the cluster and never committed, so W3 built
`gene_index_0.9.0` from the GDC gene model, correctly and on purpose at 0.9.0
("provisional, superseded the moment W1 emits theirs").

**Corrected recommendation: W1 does not emit a competing index.** The original
reasoning — "single-cell features are the more constrained set" — was about which
set constrains the *join*, and that is honoured by intersecting, not by W1 owning
the file. Deconvolution needs each gene in both matrices, so the operative index
is the intersection of the GDC gene model and this deposit's features.

`src/reference/jobs/check_gene_index.py` measures that intersection and the panel
coverage inside it. Then, in the meeting:

- intersection close to the single-cell feature count → **adopt 0.9.0 as-is and
  promote it to 1.0.0**, W3 owning the file. #2 closes.
- materially smaller → **1.0.0 IS the intersection**, emitted once, by whoever the
  team names.

Either way: one index. A panel gene missing from the intersection is a reason to
revisit the index, not to drop the gene.
**Original recommendation: W1 emits it in week 1** from the GSE178341 feature
table, with the written fallback *"if W1's ingest slips, W3 emits a provisional
index from the GDC gene model and W1 conforms."*

**The fallback fired.** W1's week-1 work landed — ingest guards, per-batch QC,
ambient estimator, labels — and `config/gene_index/` still holds only its
README. Nothing under `src/` handles Ensembl IDs at all; W1's pipeline is
symbol-keyed throughout.

W3 has built it: [`src/bulk/gene_index.py`](../src/bulk/gene_index.py), emitted
as **`gene_index_0.9.0`** — provisional, and versioned to say so. It is *not*
1.0.0.

### Measured, 2026-08-17 — neither index is the right shared index

W1's Ensembl IDs were pulled from the Broad-hosted
`colon10x_default_dDvec_ensgID.csv.gz` (188 KB, the file identified in §11) and
intersected with W3's index:

| | Genes |
|---|---|
| W1 (GSE178341 features) | 43,078 |
| W3 (index 0.9.0, GENCODE v36) | 60,616 |
| **On both** | **39,236** |
| W1 only — not in v36 | 3,842 · **8.9% of W1** |
| W3 only — bulk-only | 21,380 · 35.3% of W3 |

**All 23 panel genes are present on both.** Tiers A–D are safe either way, which
was the thing most worth checking.

Note the 8.9%: decision #3 predicted *"the usual way a join silently loses ~8% of
genes"* and that is almost exactly what version drift between W1's CellRanger
reference and GENCODE v36 costs. The prediction was right.

**So the question is not "0.9.0 or 1.0.0".** Adopting either wholesale throws
away genes the other arm measured — W3 conforming loses 21,380, W1 conforming
loses 3,842. **The shared index should be the intersection: 39,236 genes,
committed as `gene_index_1.0.0`.** Each arm keeps its own full matrix for its own
work; the shared index is what integration joins on, and 39,236 is far more than
the 500–2,000 markers deconvolution needs (§2.1 error #4).

To settle at the weekly:

1. **Agree the shared index is the intersection**, not either arm's native set.
2. W1 confirms the reference release behind `dDvec_ensgID` so the 3,842 can be
   attributed to version drift rather than to reference filtering.
3. Promote to `gene_index_1.0.0` in its own commit; retire 0.9.0.

---

## 3 · Symbol vs. Ensembl ID on the shared index — ADOPTED AS RECOMMENDED

**Raised:** scaffold, 2026-08-15 · **Implemented:** W3, 2026-08-17 ·
**Owner:** W1 + W3

Follows from #2 and is the usual way a join silently loses 8% of genes. TCGA
STAR counts arrive as Ensembl IDs with versions; the panel and both labelling
axes are written as symbols; single-cell references are usually symbols.

**Recommendation, now implemented: Ensembl IDs (unversioned) as the index,
symbols as a mapped column.** The version suffix is stripped into its own column
because the same gene carries a different suffix in a different GENCODE release,
so a versioned key drops those genes silently. Panel lookups resolve through
`gene_index_0.9.0.map.tsv`; `panel_resolution_report()` reports unmapped and
ambiguous panel genes in week 1.

**Ambiguous symbols are not auto-resolved.** A panel gene mapping to two Ensembl
IDs is a decision, not a tie-break. Tiers A–D are eleven genes — pin them by
hand and write down why.

---

## 4 · COAD and READ: pool or stratify — OPEN

**Raised:** execution_plan.md §4 (W3 gotchas) · **Owner:** W3 · **Needed by:** week 3

Different treatment patterns. The plan says decide now and write down the
reason. It is here so the reason gets written somewhere findable.

---

## 5 · CODEOWNERS handles and branch protection — OPEN

**Raised:** scaffold, 2026-08-15 · **Owner:** W2 · **Needed by:** week 1

[`.github/CODEOWNERS`](../.github/CODEOWNERS) has placeholder handles, so
GitHub currently ignores most of it and the "two approvals for frozen code" rule
is convention rather than mechanism. Fill in the handles and turn on branch
protection for `main`.

---

## 6 · Which cohort W2 ingests — OVERTAKEN BY EVENTS, see §8

**Raised:** W2 planning, 2026-08-15 · **Resolved:** 2026-08-15 by W4 shipping
`src/estimator/lee_io.py`, which loads **both** cohorts (SMC and KUL3) behind one
`load_lee_cohort(which=...)`. The split proposed below did not happen and no
longer needs to — nobody duplicated anything, because W4 got there first.

What replaces it is a narrower question about the *shape* of the artifact rather
than who downloads it. See **§8** below.

<details>
<summary>Original proposal, kept for the record</summary>

W2 needs real cells in week 1 so the harness is never queued behind W1's
five-patient pilot. W4 owns both Lee cohorts (GSE132465 SMC, GSE144735 KUL3) and
does its own ingest and ambient correction in weeks 1–2. If W2 ingests Lee too,
two people run the same QC in the same fortnight.

**Proposal: split the pair rather than duplicate it.**

| | Cohort | Does |
|---|---|---|
| W4 | GSE132465 (SMC) | ingest, QC, ambient correction, labels — as planned |
| W2 | GSE144735 (KUL3) | same pipeline *shape*, then hands the cleaned object to W4 |

W4 needs both cohorts eventually and gets one for free. W2 gets real cells in
week 1 without a redundant fortnight. Pipeline shape is coordinated; code is not
shared yet — the same rule W4 already has with W1.

**Fallback if W4 objects:** W2 waits for W4's SMC object and spends week 1 on the
design spec and the data-free plumbing (`common/io.py`, `harness/results.py`,
the deconvolution adapter protocol, the NNLS baseline) — all of which are
already done, so the cost of the fallback is low.

**Note on independence:** sharing a cohort between W2 and W4 does not compromise
the harness. The harness needs cells with labels, not a dataset nobody else has
touched. What it does require is that the patients held out for a pseudobulk
sample are held out — which is enforced in `generate_pseudobulk`, not assumed.

</details>

---

## 7 · Detectable effect for cutpoint calibration — OPEN

**Raised:** W2 design spec, 2026-08-15 · **Owner:** W2 + W4 · **Needed by:** week 4,
before the attenuation sweep runs

The cutpoints are derived from the smallest mature-cell count at which the
estimator can still discriminate a pre-registered effect. The spec pre-registers
`s = 0.5` — a halving of per-cell output — because it sits mid-band in the
published attenuation range (×0.6–0.8).

If W4 expects MLH1 silencing to be closer to complete (`s` near 0), a smaller
effect is easier to detect, the cutpoints get looser and G4 gets easier to pass.
If the realistic effect is milder, the reverse. **This number must be fixed
before the sweep is looked at**, which is what makes the cutpoints derived rather
than chosen. See [harness_design_spec.md §4](harness_design_spec.md).

---

## 8 · The harness needs raw counts; `lee_io` emits CP10K — OPEN

**Raised:** W2, 2026-08-15, reviewing `w4/estimator-core` · **Owner:** W2 + W4 ·
**Needed by:** week 2, before the harness runs on real cells

`LeeCohort.expression` is **CP10K-normalised** and restricted to
genes-of-interest. That is exactly right for `decompose_cohort`, which consumes
per-patient summary means and needs a linear depth-normalised scale for the
Kitagawa identity to be additive. It is the wrong artifact for the harness,
which needs:

| | `lee_io` gives | harness needs |
|---|---|---|
| scale | CP10K floats | **raw integer counts** — binomial thinning is defined on counts |
| gene width | genes-of-interest | **the full index** — deconvolution wants 500–2000 markers |

Handing the CP10K frame to `generate_pseudobulk` used to truncate silently:
`astype(int64)` turns every value below 1.0 into 0, destroying exactly the
low-expressing cells the near-zero mature-cell edge cases are made of. **The
generator now refuses a non-integer matrix** with a message naming CP10K, so this
can no longer happen quietly — but refusing is not the same as being unblocked.

Options:

| | Approach | Cost |
|---|---|---|
| a | `lee_io` grows a `return_raw_counts=True` / second accessor on `LeeCohort` | Small, and W4 already streams the raw matrix — the counts exist before the CP10K step at `lee_io.py:322` |
| b | W2 writes its own reader against the same manifest files | Duplicates the parsing and the QC that `ingest.qc_flags` already does |
| c | `lee_io` caches a raw-count parquet next to the CP10K frame | Both arms served, one parse, costs disk |

**Recommendation: (a).** It is a keyword argument on code W4 has already
written, and the raw counts are in hand at the point of normalisation. Needs
W4's agreement since it is their module.

**Until it lands**, the harness runs on synthetic fixtures. That is not a block
on the week-5 gate — G3 is a statement about the estimator's algebra and
synthetic cells test it fine — but the attenuation curve wants real
expression distributions, so this should not drift past week 3.

---

## 9 · `doubly_robust` folds the interaction into both arms — OPEN

**Raised:** W2, 2026-08-15, reviewing `d9c08c0` · **Owner:** W4 + W2 ·
**Needed by:** before any decomposition result is written with
`weighting="doubly_robust"`

CLAUDE.md invariant 7: *the interaction term is reported separately, never
folded into either arm.* The pooled-reference split in
`kitagawa._doubly_robust_split` reports `interaction = 0.0` — but the cross term
has not vanished, it has been distributed **50/50 into the other two arms**:

```
comp_pooled = Δf·(m_n+m_t)/2 = Δf·m_n + Δf·Δm/2 = comp_normal + interaction/2
intr_pooled = (f_n+f_t)/2·Δm = f_n·Δm + Δf·Δm/2 = intr_normal + interaction/2
```

Verified numerically on `f_n=0.40, f_t=0.10, m_n=10.0, m_t=4.0`:

| weighting | compositional | intrinsic | interaction |
|---|---|---|---|
| normal | −3.0000 | −2.4000 | **+1.8000** |
| tumour | −1.2000 | −0.6000 | −1.8000 |
| doubly_robust | −2.1000 | −1.5000 | **0.0000** |

`−2.1 = −3.0 + 0.9` and `−1.5 = −2.4 + 0.9`, and `0.9` is exactly half the
interaction term. Here that shrinks the intrinsic term by **37.5%, toward
zero** — the direction of our prior hypothesis, which README calls out as the
worst way for a result to move.

This is not a claim that W4 is wrong. Kline (2011) is a real citation, the
pooled reference is a real estimator, and W4's docstring is explicit that it is
a first cut pending cell-level AIPW. The problem is narrower:

1. It is what invariant 7 was written to forbid, and invariants change by PR
   with two approvals and a written reason — not by implementation.
2. Writing `0.0` into the schema's `interaction` column reads as "no interaction
   here" when the honest statement is "interaction distributed into the arms."
   That is the invariant-1 failure mode wearing a different hat.

Options:

| | Approach | Cost |
|---|---|---|
| a | Amend invariant 7 to "never folded silently" and require the pooled split to record what it absorbed | Frozen-doc PR, 2 approvals. Keeps a citable estimator. |
| b | Keep `doubly_robust` but report the normal-weighting interaction alongside, not `0.0` | Schema already has the column; the value stops being a lie |
| c | Rename it — it is a pooled-reference split, not AIPW — and defer true doubly-robust to cell-level data | Cheapest, most honest about what it is; W4's docstring already says this |

**Recommendation: (c) plus (b).** Call it what it is, and put the real cross
term in the interaction column rather than zero. Neither needs new maths.

---

## 10 · Which interval goes in the schema's `ci_low`/`ci_high` — W2 PROPOSES, W4 TO CONFIRM

**Raised:** W4 in `bootstrap_over_patients`' docstring · **Answered by:** W2,
2026-08-16 · **Owner:** W2 + W4 · **Needed by:** week 5

W4 wrote, and was right to:

> *"which quantity belongs there (this cohort-level bootstrap band on
> `intrinsic`, a per-patient interval from within-patient cell-count
> uncertainty, or a value read off the hierarchical model instead) is an open
> call for whoever owns the week 4–5 hierarchical-model deliverable to make
> deliberately."*

**The two are different estimands and the project needs both.**

| Question | Sampling unit | Where |
|---|---|---|
| What is the cohort's intrinsic loss, and how sure are we? | **patients** (invariant 5) | W4's `bootstrap_over_patients` |
| Does *this* patient have enough mature cells for their own estimate to mean anything? | **cells**, by construction | W2's `harness/interval.py` |

The cutpoints have to be calibrated on the second. A cohort-level band
broadcast onto every patient row is *identical* for a patient with 800 mature
cells and one with 21, so coverage against `n_cells_mature` is flat and no
cutpoint exists. Measured on the new interval, CI width runs 2.583 at 5 mature
cells to 0.596 at 800 — a 4.3× range, which is the dependence a cutpoint needs.

**This is not a violation of invariant 5.** Invariant 5 exists because
resampling cells to make a *population* claim inflates n by the cells-per-patient
count. That failure mode requires the claim to be about patients. A
within-patient statement has cells as its sample by construction. The argument
is in `src/harness/interval.py`'s module docstring in full, so a reader who finds
the resampling loop finds the reasoning next to it.

**Proposal for the schema slot:** keep W4's current choice — the cohort-level
intrinsic band — because the per-patient row is read as part of a population
result, and say which one it is wherever it is presented, as
`attach_intrinsic_ci` already instructs. The within-patient interval is a
harness-side quantity for calibration and does not need a schema slot.

**W4: flag it if you disagree.** The only thing that would change is which
number lands in `ci_low`/`ci_high`; both intervals exist and are tested either
way.

---

## 11 · GSE178341 ships NO unfiltered droplets — ANSWERED, BLOCKING

> **Renumbered from #8 on merge.** W1 and W2 both reached for the next free
> number on the same day and collided. This is W1's decision; #8 is W2's
> `lee_io` raw-counts question. Anything still citing "open decision #8"
> for the ambient problem means this one.

**Raised:** W1, 2026-08-16 · **Answered:** W1, 2026-08-16 ·
**Owner:** W1 (W4 has the same exposure) · **Decision needed by:** week 1

SoupX and CellBender both learn the ambient profile from **empty droplets**.
"Raw counts" and "raw droplets" are different properties, and the repo's existing
gotcha ("verify you have raw counts, not normalised values") does not catch the
second.

**The answer is no.** GEO GSE178341's per-sample records state, in the data
processing section:

> "Cellranger count is run using cell-ranger v3.1 on firecoud. R package
> **dropletUtils is used to exclude chimeric reads and identify and exclude empty
> droplets**"

The deposit is three series-level files and nothing per-sample:

| File | Size |
|---|---|
| `GSE178341_crc10x_full_c295v4_submit.h5` | 1.1 Gb |
| `GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz` | 2.8 Mb |
| `GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz` | 2.4 Mb |

All 181 GSM records say "Supplementary data files not provided". 371,223 cells is
the post-calling count.

### The one lead — CHECKED, dead

**Resolved 2026-08-16.** `colon10x_default_dDvec_sampleID.csv.gz` has
**371,223 rows** — exactly the published post-QC cell count. The Broad-hosted
file is the same cells as GEO in a fatter HDF5 layout, not a superset. **There
are no empty droplets in any public source.** The section below is kept for the
record; the `dDvec_ensgID` and `dDvec_batchID` files are still worth having for
open decision #3 and for per-batch QC.



The authors' analysis repo (`matanhofree/crc-immune-hubs`) points at a
Broad-hosted directory that is **not linked from GEO**:

```
https://portals.broadinstitute.org/crc-immune-hubs/extra/data/colon10x_default/
  colon10x_default_dSp_rawCount.h5      2,261,717,327 bytes   (2021-08-24)
  colon10x_default_dDvec_sampleID.csv.gz        1.8 MB
  colon10x_default_dDvec_geneID.csv.gz          161 KB
  colon10x_default_dDvec_ensgID.csv.gz          188 KB
  colon10x_default_dDvec_batchID.csv.gz          18 KB
  colon10x_default_dZ_annot.csv.gz              541 KB
  colon10x_default_dZ_metatable.csv.gz          326 KB
```

That `.h5` is **1.9x the size of GEO's** (1,203,550,558 bytes) and two months
newer. Both are public and unauthenticated. The extra gigabyte could be empty
droplets — or just different HDF5 compression, since "rawCount" conventionally
means un-normalised rather than un-filtered.

**Cheap way to settle it:** the `dDvec` files are axis labels, one row per
barcode. Download `colon10x_default_dDvec_sampleID.csv.gz` (1.8 MB) and count
its rows. ~371k means filtered and the lead is dead; millions means empty
droplets survived and **CellBender is back on**. That is a 1.8 MB answer to a
2.1 GB question.

Two side benefits regardless of the outcome: `dDvec_ensgID` gives Ensembl IDs
alongside symbols, which settles open decision #3, and `dDvec_batchID` gives the
batch key that `src/reference/qc.py` needs for per-batch thresholds.

### Sources audited

| Source | Status |
|---|---|
| GEO GSE178341 supplementary | Filtered — dropletUtils excluded empty droplets |
| GSM records (all 181) | "Supplementary data files not provided" |
| Broad SCP1162 | Login-walled; reports the same 371,223 cells |
| **Broad crc-immune-hubs portal** | **2.1 GB `rawCount.h5` — unverified, see above** |
| dbGaP `phs002407.v1.p1` | Controlled access, raw FASTQ, weeks of latency |
| GitHub `matanhofree/crc-immune-hubs` | Points to the Broad portal; no cloud buckets |

### The Lee cohorts are worse — W4 should read this

- **GSE132465 (SMC):** filtered only, and *"raw data not provided due to patient
  privacy concerns"*. No SRA, no dbGaP. Dead end.
- **GSE144735 (KUL3):** filtered only on GEO, **but raw sequencing is at
  ArrayExpress `E-MTAB-8410`**, which is typically open access.

If E-MTAB-8410 is open, KUL3 becomes the **calibration cohort**: re-run
CellRanger to get genuine empty droplets, run CellBender there, and compare
against SoupX/DecontX on the same cells. That yields a measured estimate of how
much the degraded methods miss, which can then be carried onto GSE178341 as a
stated bound rather than an assumption. Worth W4's time in week 1.

### Consequences — narrower than they first look

- **CellBender cannot run at all.** Its model is a generative fit across *all*
  barcodes, learning background from the cell-free ones. There is no
  filtered-matrix mode. The GPU step goes away with it.
- **SoupX still runs, degraded.** Reading the soup off empty droplets is the
  default, not a requirement: construct with `calcSoupProfile = FALSE` and supply
  a profile estimated from the filtered cells via `setSoupProfile()`.
  `autoEstCont()` derives the contamination fraction from clusters and marker
  genes without touching empties.
- **DecontX (`bioconductor-celda`) is the drop-in second method.** Built for this
  case — models each cell as a mixture of its own cluster's distribution and
  contamination from the others, no empty droplets required. **Not currently in
  `env/w1_reference.yml`.**
- **G1 still has substrate.** It needs a pre/post retention ratio per gene and a
  check of whether that tracks abundance across tiers. Any correction that runs
  provides it. Week 2's "two methods, compared" survives as **SoupX (degraded)
  vs DecontX**.
- **W4 is likely in the same position** on GSE132465 / GSE144735. Check before
  building `correct_ambient()`.

### The real cost, and why G1 matters more now

Both surviving methods estimate contamination **from the cells themselves**
rather than from true empty droplets. DecontX defines contamination as counts
that look like they came from other clusters — and this project's hypothesis is
that mature markers sit at unexpectedly low levels in cells that are still
present. A rare mature cell's genuine low-level GUCA2A is exactly what a
cluster-based method may reassign to "contamination from elsewhere." **The
correction can absorb the intrinsic signal it exists to protect.**

That makes G1 more load-bearing, not less: it is the check that catches this.
Note the mild circularity in the methods write-up — G1 asks whether retention
tracks abundance, and DecontX's correction is itself driven by cluster and
abundance structure.

### Quality assurance that does not need empty droplets

Losing CellBender is not losing quality control. Four routes remain, and the
first two are independent of any correction method:

1. **Measure contamination with impossible genes.** You do not need empty
   droplets to *quantify* soup — only cells where a transcript is biologically
   impossible. Haemoglobin (HBB, HBA1) in epithelium, immunoglobulin in
   non-plasma cells, PTPRC in EPCAM+ cells. Residual counts there *are* the
   ambient rate. Implemented as `contamination_fraction()` and
   `contamination_by_sample()` in
   [src/reference/ambient.py](../src/reference/ambient.py), with recovery tests
   at five known contamination levels. Use it to set DecontX priors **and** to
   audit the correction afterwards.
2. **The panel's tier structure is already an ambient control.** Tier D
   (MS4A12) is colonocyte-restricted yet usually retained — expected *neither*
   compositional nor intrinsic. Ambient contamination is abundance-driven, so
   soup driving the results would move tier D too. The falsification rule
   (A, B and D agreeing means the estimator is broken) was designed in before
   any of this surfaced and doubles as an ambient check.
3. **Ambient dose-response in W2's harness.** Spike simulated soup into
   pseudobulk at known contamination levels and measure where the intrinsic term
   breaks. Converts "we could not correct properly" into a quantified
   sensitivity bound, in the same spirit as the §2.2 attenuation curve — a
   publishable object rather than a caveat. **Proposed addition to W2's scope.**
4. **The ambient-free substrates already in the plan.** The plate-based subset
   from ICBI has no soup by construction; spatial (week 13+) has no droplet
   encapsulation at all.

### Options

| | Approach | Cost |
|---|---|---|
| a | **Proceed with SoupX (degraded) + DecontX.** Week 2 and G1 both run. Add `bioconductor-celda` to `env/w1_reference.yml`. | Cheap and available today. Both methods infer the ambient profile from the cells, so the correction may absorb real low-level expression — stated as a limitation, and G1 is the check on it. |
| b | **Add the plate-based ambient-free subset** from the ICBI atlas as corroboration. Plate protocols have essentially no soup, so an intrinsic signal surviving there is strong evidence it is not contamination (§8.2). | Free. Independent of the correction method, which is exactly why it is worth having alongside (a). Power depends on the subset's size — the week-1 ICBI metadata pull tells us. |
| c | **Ask the authors** for the pre-dropletUtils matrices. | Free, uncertain, slow. Worth sending regardless; would restore CellBender if it lands. |
| d | **dbGaP `phs002407.v1.p1`** (controlled access) → re-run CellRanger without cell calling to regenerate raw droplet matrices. | Application through your signing official, typically weeks, plus CellRanger over 181 samples. Only worth it if CellBender specifically is judged necessary. |

**Recommendation: (a) + (b), with (c) sent in parallel.** Week 2 proceeds on
schedule with SoupX and DecontX, the plate-based subset gives an independent
corroboration that does not depend on either correction, and the author email
costs nothing. (d) only if the week-5 gate discussion concludes CellBender is
required — the dbGaP latency means that would move the gate, not fit inside it.

`assert_unfiltered_droplets()` in
[src/reference/ingest.py](../src/reference/ingest.py) implements the check; run it
on the downloaded `.h5` to confirm on our own bytes rather than on the metadata.

---

## 9 · Only 36 of 62 patients have matched normal — OPEN

**Raised:** W1, 2026-08-17 · **Owner:** W1 + W4 · **Needed by:** week 3

Measured from the deposit, not assumed. `patient_cohort_table()` over all
370,115 cells:

| | matched | unmatched | total |
|---|---|---|---|
| MMRd | 21 | 13 | 34 |
| MMRp | 15 | 13 | 28 |
| **All** | **36** | **26** | **62** |

The compositional term is Δ(mature fraction) against the patient's **own**
normal, so an unmatched patient contributes to neither arm — not partially, not
with a wide interval. Absent.

execution_plan.md §8.2 listed exactly this as a week-1 check that could
invalidate the plan, and §8.4 states "n≈60 supports the primary paired
analysis". **The real figure is n=36.** Consequences:

- The primary paired analysis is powered by 36, not 62.
- The pre-registered MMR subgroup contrast is **21 vs 15**. §8.4 says interaction
  contrasts need roughly 4x the primary; at n=36 that would want ~144. This
  supports the existing commitment to report MMR as an estimate with an interval
  rather than as a test — it does not rescue it.
- Matching is not differential by MMR status (62% of MMRd, 54% of MMRp), so the
  cohort is not distorted in that dimension.

**The decision: what happens to the other 26?**

| | Approach | Cost |
|---|---|---|
| a | Drop them from the decomposition entirely. | Cleanest. n=36 everywhere. Discards 42% of the cohort, including their tumour cells. |
| b | Include them against a **pooled** normal reference, reported separately and flagged. | Uses all 62 for a weaker, clearly-labelled secondary analysis. Breaks the within-patient design, so the pooled estimate is not comparable to the paired one and must never be merged with it. |
| c | Treat them as `estimability="not_estimable"` rows and emit them. | Fits the frozen schema and the project's own three-way framing: the honest answer for these patients is that the split is not identifiable. Makes the 26 visible in the output rather than absent from it. |

**Corrected recommendation (see CORRECTIONS #4): (c) is right in spirit but is
NOT free, and cannot be done by W1 alone.** Verified against the code:

> ### The arithmetic, so this is not decided on intuition — 2026-08-22
>
> Let **X** be how many of the 36 matched patients fall below G4's mature-cell
> threshold. The 26 unmatched patients would enter at `n_cells_mature = 0`, so
> they are below it by construction.
>
> | population | size | G4 fails when |
> |---|---|---|
> | matched only | 36 | X ≥ **18** |
> | all patients | 62 | X ≥ **5** |
>
> **The two definitions give opposite verdicts for any X between 5 and 17** — 13
> of the 37 possible values, and the entire plausible middle of the range. This
> is not a corner case that might bite; it is the likely outcome.
>
> G4's failure consequence is pre-committed: *"Non-identifiability finding with
> diagnostics becomes the headline result, not a caveat."* So under the
> all-patients definition, a cohort-design fact — how many patients had a normal
> sample taken — decides what the paper is about.
>
> **This has to be settled before any decomposition runs.** Settled afterwards it
> is unfalsifiable, whichever way it goes.

1. **There is nowhere to record the reason.** `ESTIMABILITY` already contains
   `not_estimable`, so no schema change is needed to *emit* these rows — but
   `COLUMNS` has no reason field. In the output, "no normal arm" and "too few
   mature cells" become the same value. Those have opposite implications: one is
   a cohort-design fact known in week 1, the other is a power finding. A reader
   could not tell 26 unmatched patients from 26 depleted ones.
2. **It can flip gate criterion G4.** `gate_g4_verdict` computes the fraction of
   patients whose `n_cells_mature` is below `cutpoints.wide`. Unmatched patients
   would enter with 0. Adding 26 guaranteed-below rows to 36 real ones moves that
   fraction toward the 50% line, and G4's pre-committed failure consequence is
   *"non-identifiability with diagnostics becomes the headline result"*. A cohort
   fact would be reported as a positivity finding.

**So, three things, in this order:**

- **Now, W1 alone:** emit a **cohort-coverage artifact** — one row per patient
  with matched status and cell counts, versioned like any other result. Full
  visibility of the 26, no schema change, no gate contamination.
  `patient_cohort_table()` already computes it.
- **W2 to confirm:** G4's population is **matched patients only**. G4 asks about
  mature-cell depletion, not about cohort matching, and mixing them makes it
  answer a question it was not pre-committed to.
- **Only then, if the team still wants them in the results frame:** a
  `not_estimable_reason` column on `src/schema.py` — frozen, PR + 2 approvals +
  written reason. Worth doing, but it is a shared change, not a W1 one.

---

## 10 · Pre-register the refined tier-B (MLH1) test — OPEN

**Raised:** W1, 2026-08-17 · **Owner:** W1 + W2 + W4 · **Needed by:** before any
expression is examined — the value is entirely in committing first

GSE178341's metadata carries per-patient `MLH1Status` and free-text `MMR_IHC`.
Together they turn tier B from "does MLH1 come out intrinsic on average" into a
directional, mechanism-specific prediction with MMR status held fixed.

`assign_mlh1_strata()` in [src/reference/ingest.py](../src/reference/ingest.py):

| Stratum | n | matched | Predicted MLH1 intrinsic loss |
|---|---|---|---|
| `mlh1_methylated` | 22 | **12** | **High** — transcriptionally silenced |
| `mlh1_intact_mmrd` | 10 | **7** | **Near zero** — MMRd via MSH2/MSH6/PMS2 |
| `mlh1_deficient_unmethylated` | 2 | 2 | Ambiguous — report separately |
| `mmr_proficient` | 28 | 15 | Near zero |

Three things make this stronger than it looks:

1. **The contrast is within MMRd.** `MLH1Meth` is strictly nested in MMRd (22 of
   22, zero MMRp), but the 12-vs-7 comparison holds MMR status constant, so it
   is not a restatement of the pre-registered MMR contrast.
2. **The negative control is mechanistic.** Those 7 are MMR-deficient by the
   same MSI-H phenotype but through other genes, leaving MLH1 transcription
   untouched. Same disease biology, MLH1 specifically spared.
3. **C115 and C132 are excluded from it.** Their IHC reads "MLH1 and PMS2
   deficient" — MLH1 protein lost *without* methylation, so probably a germline
   MLH1 variant whose transcript may or may not survive NMD. They are
   indistinguishable from the negative-control group on `MMRStatus` and
   `MLH1Status` alone, and would have diluted the arm where near-zero loss is
   predicted.

**The decision: adopt this as supporting evidence for G2, not as its primary
basis.** 12 vs 7 is the ceiling and positivity (n≥50 mature cells) will reduce
both arms — possibly to 6 vs 3, which carries no weight. G2's primary test stays
tier separation (A compositional, B intrinsic, D neither) across all matched
patients. This is a second, mechanistically sharper line of evidence that costs
nothing to commit to now and is worthless if committed to later.

G2 is a pre-committed gate criterion, so **this needs the team, not just W1.**

---

## 11 · Half the cells come from sorted samples — OPEN, AFFECTS THE COMPOSITIONAL ARM

**Raised:** W1, 2026-08-17 · **Owner:** W1 + W4 · **Needed by:** before the pilot

`PROCESSING_TYPE` in the metatables takes four values, and only one leaves
cell-type composition untouched:

| PROCESSING_TYPE | What it is | normal | tumour |
|---|---|---|---|
| `unsorted` | untouched | 67,125 | 164,612 |
| `CD45pMACS` | **CD45+ magnetic sorting — immune enrichment** | 38,726 | 84,040 |
| `mixUnsortCD45MACS` | deliberate mixture | 6,782 | 3,888 |
| `LiveMACS` | viability selection | 231 | 4,711 |

The compositional term is Δ(mature epithelial fraction). A CD45-enriched sample
has had its epithelial fraction driven down **by the protocol**, so comparing it
against an unsorted sample measures the sort, not the tumour.

**Pairing does not absorb this.** Unlike chemistry, which is constant within 61
of 62 patients, `PROCESSING_TYPE` is mixed within **45 of 62**. A patient with an
unsorted normal and a CD45-sorted tumour would show a large apparent loss of
epithelium that is entirely artefact — and in the direction of the prior
hypothesis.

### What it costs

`PROCESSING_TYPE` is a **per-sample** property and patients have several
samples, so this filters samples rather than dropping patients:

- 36 patients matched
- **32 have unsorted cells in both arms**
- **30 also clear 500 tumour / 300 normal unsorted cells**

So the compositional arm is **n≈30**, down from 36. Smaller than feared.

### It broke the first pilot selection

The original five were chosen on total counts, before this was checked:

| | unsorted normal | unsorted tumour | verdict |
|---|---|---|---|
| C114 | **0** (all mixUnsortCD45MACS) | 1,390 | unusable — **and it was the only `mlh1_methylated` patient**, i.e. tier B's positive control |
| C115 | **0** | **0** | unusable in both arms |
| C140 | 1,677 | 2,417 | fine |
| C142 | 2,038 | 2,042 | fine |
| C162 | 6,726 | 8,103 | fine |

`select_pilot.py` now computes eligibility on unsorted cells only.

### The decision

| | Approach | Cost |
|---|---|---|
| a | **Compositional estimates from `unsorted` samples only.** Sorted samples excluded from that arm. | n≈30. Clean, and the exclusion is principled rather than empirical. Recommended. |
| b | Also allow `mixUnsortCD45MACS` if the mixing ratio is documented. | Recovers a few patients, but the ratio is not in the metadata and would have to be assumed. |
| c | Model sorting as a covariate and keep everything. | Keeps all cells, but requires believing a linear adjustment can undo a physical enrichment. It cannot for a fraction. |

**Recommendation: (a).** Also worth stating explicitly that the *intrinsic* arm
is far less affected — per-cell expression among surviving epithelial cells does
not depend on the mixture the way a fraction does — so sorted samples may still
contribute there, flagged, if the extra cells are wanted. **That asymmetry needs
W4's sign-off**, since it means the two arms of the same decomposition run on
different cell sets.

---

## 12 · The 20% mitochondrial cap cuts normal harder than tumour — ANSWERED BY DATA

> **For W4 (CORRECTIONS #5):** the ask is not "adopt 50". 29.8% is the epithelial
> median measured on **GSE178341**; W4's cohort is Lee, with different chemistry
> and dissociation. What transfers is the *reasoning* — 20% is a lymphocyte
> number, colonic epithelium runs far higher, and a cap that cuts one arm harder
> than the other biases the compositional term directly. **Run the same
> per-compartment `pct_mito` breakdown on Lee and set the cap from it.** If Lee's
> epithelium genuinely sits near 20%, keeping 20 is the right answer there, and
> the two cohorts differing for a measured reason is fine. The two cohorts
> differing because one number was never checked is not.
>
> Note also that W4's cap is a single hard threshold by design
> (`src/estimator/ingest.py`), while W1's is per batch. That difference is
> defensible — but it should be a decision, not an accident.

> **Resolved 2026-08-17 from the pilot.** The `pct_mito` distribution settles it.
> Colonic epithelium runs at a **median of 29.8% in normal** and 21.1% in tumour,
> against 4–11% for every immune and stromal compartment — so a 20% cap sits
> *below the median for normal epithelium*.
>
> | cap | epithelium kept | tumour/normal gap |
> |---|---|---|
> | 20% | 40.7% | **22.7 pts** |
> | 30% | 66.1% | 21.7 pts |
> | 40% | 85.9% | 11.4 pts |
> | **50%** | **100%** | **0.0 pts** |
>
> And the deposit is **already filtered at 50%** upstream — observed max 49.976
> (normal) / 49.988 (tumour) — so a 50 cap is a no-op rather than an opinion, and
> anything lower double-filters on top of the authors' cut.
>
> **`DEFAULT_MAX_PCT_MITO` is now 50.0.** W4 uses 20 on the Lee cohorts and
> **must match or justify diverging**, or the cohorts are not comparable at the
> gate. That is the part still open.
>
> One correction to the original reasoning below: the mitochondrial content is
> **not** mostly ambient. Contamination is ~2.7% of counts and MT genes are ~18%
> of the soup, so ambient contributes roughly 0.5% of a cell's counts — nowhere
> near a 30% observed fraction. Colonocytes are simply metabolically active,
> which strengthens the case for the high cap rather than weakening it.



**Raised:** W1, 2026-08-17 (from the pilot run) · **Owner:** W1 + W4 ·
**Needed by:** before any compositional estimate

The pilot retained **61.5% of 44,794 cells**, and `pct_mito > 20%` is doing
almost all of the cutting — over half the cells in several samples
(C107_N 1,022/1,814; C122_N_1_1_1 1,078/2,155; C162_T_c1_v3 2,343/4,128), while
the MAD-based count and gene bounds cut very few.

**It is not cutting the two arms equally.** Retention by tissue:

| patient | normal | tumour | gap |
|---|---|---|---|
| C107 | 42.8% | 68.2% | **-25.4** |
| C122 | 49.4% | 65.2% | **-15.8** |
| C165 | 52.2% | 60.4% | -8.2 |
| C138 | 69.9% | 79.3% | -9.4 |
| C162 | 71.3% | 57.4% | +13.9 |

Four of five patients lose substantially more of their **normal** than their
tumour.

### Why this is not an ordinary QC quibble

The compositional term is Δ(mature fraction) between a patient's tumour and
their own normal. Mature colonocytes are metabolically active, fragile, and
carry high mitochondrial content — they are exactly what a mitochondrial cap
removes first. Cutting normal harder than tumour **understates the normal mature
fraction**, which **inflates the apparent compositional loss**.

That is a bias in the direction of the prior hypothesis, which README's opening
identifies as the worst kind of result. It would be produced by a QC parameter,
not by biology.

The soup profile makes the mechanism plain: **eight of the ten most abundant
ambient genes are mitochondrial** (MT-CO2, MT-CO3, MT-CO1, MT-ATP6, MT-ND4,
MT-RNR2, MT-CYB, MT-ND3). Much of the measured mitochondrial fraction is ambient
contamination from lysed cells, not cell-intrinsic stress — so the cap is partly
filtering on soup, and soup load differs by sample.

### Options

| | Approach | Cost |
|---|---|---|
| a | **Raise the cap for epithelium** to 30-50%. Colonic epithelium routinely exceeds 20% in published work; 20% is a lymphocyte-appropriate number. | Keeps more real cells. Needs a value chosen from the data rather than convention — plot the per-compartment distribution first. |
| b | **Per-batch MAD on pct_mito**, like the count and gene bounds, instead of a hard cap. | Adapts to differing soup load. Departs from W4's hard cap, so the two cohorts diverge unless W4 follows. |
| c | Apply the cap **after** ambient correction, so it filters on cell-intrinsic mitochondrial content rather than on soup. | Most defensible; the soup profile says the cap is partly measuring contamination. Reorders the pipeline. |
| d | Keep 20% and carry differential retention as a stated limitation. | Cheapest, and wrong — the bias points at the conclusion. |

**Recommendation: (c) then (a).** Filter mitochondria after ambient correction,
and set the threshold from the observed per-compartment distribution rather than
from convention. Whatever is chosen, `differential_retention()` in
[src/reference/qc.py](../src/reference/qc.py) must show no systematic tumour /
normal gap before any compositional number is believed. **W4 needs to match**,
or GSE178341 and the Lee cohorts are not comparable at the gate.

---

## 13 · W1 and W4 label cells differently — OPEN, BLOCKS COMPARABILITY

**Raised:** W1, 2026-08-17 · **Owner:** W1 + W4 (+ W2 consumes both) ·
**Needed by:** before either cohort's decomposition is compared

Both workstreams now have label code, which §4 sanctions ("same pipeline shape as
W1 — coordinate, do not share code prematurely"). But the *definitions* diverge,
and the mature fraction is the compositional term, so divergence is not cosmetic.

| | W4 `src/estimator/labels.py` | W1 `src/reference/labels.py` |
|---|---|---|
| Column | `mature__{axis}__{rung}`, boolean | `label_{axis}_{rung}`, categorical |
| epithelial rung | top 50% of score | **all epithelium** |
| lineage | top 35% | top 50% (median split) |
| crypt_position | top 15% | top 33% (tertiles) |
| best4 | top 5% | top 5% (marker-gated, not score-binned) |
| Scoring | mean raw expression | depth-normalised, z-scored per gene |
| Cut points computed | over the whole input | **within each patient's NORMAL arm** |
| Non-epithelial cells | caller must pre-filter | labelled `non_epithelial` explicitly |

Three of these matter:

1. **The rung schedule.** Different mature fractions per rung mean
   Δ(mature fraction) is not measuring the same thing in the two cohorts, so the
   granularity curve (§6.2) would not be comparable and the week-5 gate would be
   comparing incomparable numbers.
2. **Depth normalisation.** GSE178341 mixes v2 and v3 chemistry, which differ in
   capture efficiency, so a raw-mean score ranks deep cells as more mature for
   technical reasons. W1 normalises for this. Lee may not need it; the two should
   still agree on whether it is done.
3. **Where the cut points come from.** Three options, and they fail differently.

   *Per-sample quantiles* are the worst and were W1's bug: a within-sample
   quantile cannot express a between-sample difference, so the mature fraction
   equals the quantile in every arm and Δ is **identically zero by construction**.
   The pilot showed it — every `opposite_lineage` arm returned exactly 0.500 at
   `lineage`, 0.333 at `crypt_position`.

   *Pooled quantiles* are what W4 does today (`score.quantile(...)` over the whole
   input). **This does not zero the term** — the cut is global, so the two arms
   are free to differ. The problem is subtler: the threshold depends on the
   cohort's tumour:normal cell mix, so the same biology gives a different mature
   fraction in Pelka and in Lee. Invariant 4 has us estimate per study and then
   meta-analyse, which requires the per-study numbers to be on a comparable scale.
   Pooled cuts are not.

   *Reference-arm cuts* — W1's current approach — take the threshold from the
   patient's own normal and apply that absolute value to their tumour. Anchored to
   a biologically meaningful reference, and a within-patient contrast, which is
   what the decomposition is built on. **This is the ask of W4**, and the reason
   is comparability across studies, not a zeroed term.

W4's own docstring defers on the numbers — *"Swap them for whatever W1/W2 lands on
once real cells are in hand; this file's job is the plumbing... not the final
number."* So this is a handoff waiting to happen rather than a disagreement.

**Recommendation: adopt W1's definitions, now that they have been run on real
cells.** Specifically: all-epithelium at the coarsest rung (it is the lower bound
of the granularity curve and is *supposed* to look degenerate), **cut points
from the patient's own normal arm — not per-sample quantiles**,
depth-normalised z-scored input, and marker-gating rather than
score-binning for BEST4+ since that population is discrete. W1 keeps the
categorical bin names — they carry which bin, not merely mature-or-not, which the
granularity curve needs.

**Interop already exists**, so nothing is blocked while this is decided:
`cell_type_vector()` emits the `cell_type` array W2's `generate_pseudobulk`
consumes (renaming whichever bin is mature to `mature_colonocyte`), and
`maturity_summary()` emits every column `decompose_cohort` requires except
`gene`, `mean_normal` and `mean_tumour`, which need the target gene's expression.
Both are covered by tests that call W2's and W4's real functions rather than mocks.

---

## 14 · Neither labelling axis is a clean maturity measure — OPEN, BLOCKS COMPOSITION

> ### UPDATE 2026-08-20 — the kappa has been read
>
> Options (a) and (c) below are **both implemented** (`unresolved_depth` label,
> binomial depth matching). The run that followed says three things, and the
> second was not anticipated.
>
> **1 · Axis 1 carries real but *fair* signal.** Cohen's kappa against the
> authors' `cl295v11SubFull`, at `depth_quantile=0.10` (target 1,352 UMIs):
>
> | axis | rung | kappa | agreement | sens | spec |
> |---|---|---|---|---|---|
> | `stem_pole` | `lineage` | **0.313** | 62.1% | 0.93 | 0.44 |
> | `stem_pole` | `crypt_position` | **0.313** | 62.1% | 0.93 | 0.44 |
> | `stem_pole` | `best4` | 0.030 | 63.0% | 0.04 | 0.98 |
> | `opposite_lineage` | `lineage` | −0.271 | 32.4% | 0.43 | 0.26 |
> | `opposite_lineage` | `crypt_position` | −0.291 | 33.7% | 0.35 | 0.33 |
>
> 0.313 is "fair" on Landis–Koch, not "good". Backing the 2x2 out of the
> published margins: of 10,276 scored epithelial cells the authors call **6,327
> stem/TA-like** and 3,949 mature; we call **7,216** mature. We recover 93% of
> their mature cells and also call **56% of their stem cells mature**. The error
> is one-directional and it points the same way dropout does.
>
> **2 · `lineage` and `crypt_position` are the same partition on axis 1.**
> Identical mature counts in all ten pilot arms (885, 765, 342, 694, 149 normal;
> 271, 862, 137, 2354, 757 tumour), `crypt_middle` absent entirely, and the
> per-patient "supports only 2 of 3 bins" note fired for every patient. Reading
> `_bin_against_reference`, the mechanism is that the reference-arm quantiles
> both land on the tied score and `searchsorted(..., side="right")` puts the
> whole tied block on the mature side — which is also why the **normal** arm
> reads 87–96% mature at a rung whose nominal target is 50% (and 33%).
>
> So axis 1's mature call is, operationally, **"no stem marker detected at 1,352
> UMIs"** — a detection gate, not a median split and not a tertile. The
> granularity curve has **three points on axis 1, not four**, and two of them
> were being reported as different resolutions of the same boundary.
> `rung_degeneracy()` now reports this rather than leaving it to be noticed.
>
> **3 · The depth flag cannot clear on this axis, and should stop being the
> gate.** After thinning, every scored cell sits at the same expected depth, so
> detection probability depends only on a marker's *fraction* of the
> transcriptome. The route from depth to label is closed; the residual
> association (`depth_auc` 0.378) is then a correlation between RNA content and
> stem-marker fraction, which is what proliferating cells are expected to show.
> `label_depth_confounding` cannot tell that from the artifact — which is the
> reason `annotation_concordance` exists — so **the depth target must be chosen
> by kappa, not by the flag.** The sweep now prints kappa per target.
>
> ### RESOLVED IN CODE 2026-08-22 — the parts W1 could settle alone
>
> **Depth target set to q=0.25 (3,281 UMIs).** The full sweep: kappa 0.247 /
> 0.313 / 0.444 / 0.495 / 0.343 at q = 0.05 / 0.10 / 0.25 / 0.50 / 0.65. Kappa
> tracking detection depth is the signature of a detection-limited measure, which
> settles the SIGNAL-or-NOISE question in favour of signal. 0.25 dominates the
> previous 0.10 on every axis at once, so the change is free. **The peak at 0.50
> was deliberately not taken:** it buys +0.05 kappa for a quarter of the
> epithelium, a depth floor removes shallow *samples* rather than a random slice,
> and kappa is measured on the survivors — so a higher floor can raise it without
> the labels improving. The sweep now prints a `paired` column (patients keeping
> both arms) so that trade is visible next to the number it improves.
>
> **The measurement definition now travels with the number.** `maturity_summary`
> carries `depth_target` and `mature_definition`. A mature fraction gated at one
> depth and another gated at a different depth look comparable and are not, and
> invariant 4 has us meta-analyse across studies.
>
> **"Drop `crypt_position` from axis 1" was too blunt — corrected.** The collapse
> is a property of *this data at this depth target, per patient*, not of the rung
> definition: at a higher target the tied block shrinks and three bins may be
> supported. Hard-coding the removal would bake a data artifact into a frozen
> vocabulary. Instead `maturity_summary` carries `degenerate_with`, naming the
> coarser rung a finer one duplicated, so the granularity curve is built from
> non-degenerate points and the collapse is visible rather than silent.
>
> **Still needs the team:** whether a detection gate is quotable under an honest
> name, what happens to `best4`, and what axis 2 is for.
>
> ### CONFIRMED ON THE CLUSTER 2026-08-22 — and one new problem
>
> **q=0.25 was the right call, and by more than expected.** The `paired` column:
> 5 patients at q=0.05/0.10/0.25, **3 at q=0.50, 2 at q=0.65**. The +0.05 kappa
> at the peak costs 40% of the pilot's compositional n. Specificity also rose
> from 0.44 to 0.62 as detection improved (sens held at 0.90), which is what a
> detection-limited measure should do and is further evidence the signal is real.
>
> **The rung collapse persists at the higher target.** `stem_pole`
> `lineage` == `crypt_position`, Jaccard 1.000. Three points on axis 1, not four.
>
> **NEW — the depth floor cuts the two arms unequally, and it can flip the sign
> of the compositional term.** Unresolved fractions at q=0.25:
>
> | patient | normal | tumour | gap |
> |---|---|---|---|
> | C165 | **65.2%** | **0.6%** | **−64.6** |
> | C138 | 58.1% | 62.6% | +4.5 |
> | C107 | 23.1% | 40.1% | +17.0 |
> | C122 | 39.3% | 34.1% | −5.2 |
> | C162 | 24.2% | 13.3% | −10.9 |
>
> C165's Δ(mature fraction) on `stem_pole`/`lineage` was **+0.140 at q=0.10 and
> −0.053 at q=0.25** — a sign change. The other four held their sign across both
> targets. C165 is also the patient with the 64.6-point resolution gap and the
> deepest tumour sample in the pilot (upper QC bound 162,736 against a normal arm
> at 15,100), so the floor bites one arm and not the other.
>
> This is decision #12's problem one stage later: **the depth floor is QC by
> another name**, and QC that cuts one arm harder than the other moves the
> compositional term directly. A `paired` count cannot catch it — C165 kept both
> arms comfortably. `differential_resolution()` now reports it per patient, and
> `run_pilot.py` prints it beside the mature-cell counts.
>
> **What this does NOT change:** q=0.25 is still better than q=0.10 on every
> measured axis. The floor has to exist. **What it adds to the decision:** a
> patient with a large resolution gap should probably be excluded from the
> compositional arm, or the target set per patient rather than globally. Neither
> is obviously right and both are the team's call.
>
> **What this does not yet settle.** The kappa above was measured at the sweep's
> *worst* setting: q=0.10 gives a 70.2% tied block, against 33.3% at q=0.50,
> which also has the same usable fraction (33.4%) and the `depth_auc` closest to
> 0.5. Better detection should mean better agreement. **Rerun the pilot and read
> the kappa column of the sweep before resolving this.** If kappa rises
> materially at q=0.50, option (b) becomes defensible on a stated-limitation
> basis; if it stays near 0.3 at every target, axis 1 is detection-limited on 10x
> data and option (d) is the honest answer.
>
> **Axis 2's negative kappa is not noise and not a bug — it is a criterion
> mismatch, and it costs the project an argument.** Systematic anti-agreement
> (−0.29, not ~0) is what you get when two labels disagree *by definition*: axis
> 2 calls goblet cells immature and stem cells mature, the authors' stem-vs-rest
> criterion says the reverse. This test therefore does not measure whether axis 2
> works. But it does show that **README design decision 2's "agreement across
> structurally different axes" was never a testable claim with these two axes** —
> they cannot agree, because they do not measure the same quantity. Confirm by
> cross-tabulating `label_opposite_lineage_lineage` against `cl295v11SubFull`
> (expect axis-2-immature to be enriched for Goblet clusters); then either
> re-scope axis 2 as a secretory-composition control rather than a maturity axis,
> or replace it. **That is a decision for the team, and it is arguably larger
> than #14 itself.**
>
> **`best4` should not be quoted at any resolution.** Sensitivity 0.04 — the gate
> recovers 4% of the authors' BEST4+ cells while calling 279 cells BEST4+. It is
> the finest rung and therefore the upper bound of the granularity curve, so its
> failure truncates the curve at both ends.


**Raised:** W1, 2026-08-18 (from the pilot) · **Owner:** W1 + W2 + W4 ·
**Needed by:** before any compositional number is quoted

Measured on 16,955 QC-passing epithelial cells:

| axis | tied fraction | largest tied block | distinct scores |
|---|---|---|---|
| `stem_pole` | **44.8%** | 7,593 | 9,329 |
| `opposite_lineage` | 13.5% | 2,296 | 14,373 |

**Axis 1's mature bin IS the tied block, exactly.** 16,955 − 7,593 = 9,362, which
is precisely the `stem_like` / `crypt_bottom` count from the same run. The
threshold lands on the tie boundary, so "mature" on axis 1 means nothing more
than *no stem marker was detected in this cell*. It is a detection split, not a
graded maturity call.

That is a problem because **zero counts stay zero after depth normalisation**. A
shallower cell is likelier to have none of LGR5/ASCL2/MKI67/OLFM4/SMOC2 and so be
called mature, and per-sample count thresholds on this deposit span 5,140 to
162,736. `label_depth_confounding()` in
[src/reference/labels.py](../src/reference/labels.py) measures it directly.

**Axis 2 has the opposite problem**: it resolves well (13.5% tied) but measures
*absence of the goblet program*, and "not secretory" is not "mature". On the
pilot C162's tumour read 0.856 mature on axis 2 against 0.200 on axis 1, and
C165's tumour read exactly 0.000 — a tumour with few goblet cells scores as
uniformly mature.

So the two axes fail in opposite ways: axis 1 is well-conceived but poorly
resolved on this data; axis 2 is well-resolved but measuring the wrong quantity.
README design decision 2 chose them for structural independence, which they have
— but neither is clean, and "agreement across axes" is weakened when both are
individually suspect.

### Options

| | Approach | Cost |
|---|---|---|
| a | **Treat unresolvable cells as unscored** — give them their own label, exclude from the mature numerator, and report the fraction as partially unidentifiable. | Most honest, and matches the project's own three-way framing. Changes what the compositional denominator means, so W4 and W2 must agree. |
| b | Keep them as mature and carry the depth confound as a stated limitation. | Cheapest, and wrong if `label_depth_confounding()` flags — the bias would sit in the headline number. |
| c | Restrict the analysis to cells above a depth floor where absence is informative. | Defensible; costs cells, and the floor is another arbitrary cut needing justification. |
| d | Add axis 3 (chromatin/spatial) sooner than week 13+, since it is not transcript-based and has neither failure mode. | Expensive, but it is the only axis immune to both problems. |

**Recommendation: (a), with (c) as a sensitivity analysis.** A cell with no
detected stem markers might be genuinely differentiated or merely shallow, and
scoring it as maximally mature is inference from absence of evidence — the exact
move this project refuses elsewhere. **This is a scientific choice, not a default
W1 should pick**, which is why it is here rather than in the code.
## 12 · `build_signature()` asserts on the whole index; W3's matrix needs panel genes — OPEN

**Raised:** W3, 2026-08-17, building the provisional index · **Owner:** W1 + W3 ·
**Needed by:** week 4, before W1 builds an S matrix

Two requirements that are both correct and look incompatible:

| | Needs | Where |
|---|---|---|
| W1 | the index handed to `build_signature()` to contain **no** target gene | [`signature.py:96`](../src/reference/signature.py#L96), pinned by [`test_leakage.py:45`](../tests/test_leakage.py#L45) |
| W3 | the bulk matrix to **contain** GUCA2A and CDX2 | they are the outcome variables for the week-2 premise check and the Stage 4 variance question |

If the shared index is made target-free to satisfy the assertion, W3's own
deliverable becomes unrepresentable on the shared index. If it is not, every
`build_signature()` call raises.

Note also that [`config/gene_index/README.md`](../config/gene_index/README.md)
already says the opposite of what the code does — *"the index itself may contain
them for other purposes"* — so the doc and the assertion currently disagree.

**W3 has implemented the narrow fix** without touching W1's module: the
committed index carries every gene, and
[`target_free_index()`](../src/bulk/gene_index.py) produces the filtered view for
the call site:

```python
build_signature(..., gene_index=target_free_index(ids, index_map, targets))
```

Invariant 2 binds where it matters — the reference matrix — and the shared index
stays whole.

Options for the durable fix:

| | Approach | Cost |
|---|---|---|
| a | Keep it at the call site. W1 wraps every `build_signature` call in `target_free_index`. | Zero code change to frozen-ish W1 code. Relies on W1 remembering. |
| b | `build_signature()` filters targets out of `gene_index` itself, as it already does for `usable` at line 98, and drops the line-96 assertion. | One line in W1's module plus a test update. Makes the guard unforgettable. |
| c | Amputate panel genes from the shared index. | **Do not.** Breaks W3.2 and Stage 4. Listed only to rule it out in writing. |

**Recommendation: (b).** The assertion on the *reference pool* is the one that
protects anything; the one on the *whole shared index* protects nothing extra
and costs the shared index its completeness. (a) is what is in place today and
is a fine holding position.

### Why this is now urgent

[`src/reference/jobs/run_pilot.py`](../src/reference/jobs/run_pilot.py) passes
`read_gene_index(...)` straight into `build_signature_sparse` as `gene_index=`,
unfiltered, inside a `try/except` that prints `FAILED: {exc}` and continues to
the next rung. So the moment W1 has an index containing panel genes — and every
candidate index does, decision #2 measured all 23 present on both — **every rung
fails and no S matrix is written**, as a printed line rather than a raised
error. The S matrix is the W1 → W2 handoff the gate depends on.

It does not get that far today: `config/gene_index/` holds only a README, so the
job reports `!! no gene index` and skips the step. The failure is sitting on the
other side of decision #2, not in front of it.

### Ready to apply — drafted and tested against W1's branch, 2026-08-20

Option (b) as a patch:
[`decision_12_signature_filter.patch`](decision_12_signature_filter.patch).

```
git switch w1/ingest-gse178341
git apply docs/decision_12_signature_filter.patch
```

Three files, +53/−19. It replaces the assertion on the whole `gene_index` with a
filter, in both `build_signature` and `build_signature_sparse`, and updates the
two tests that pinned the old behaviour. The other three leakage guards are
untouched — `usable`, `markers`, and the post-condition on the emitted S matrix,
which is the one that still binds and is what invariant 2 is actually about.

Verified against `origin/w1/ingest-gse178341` at da6ba75: applies cleanly, full
suite **459 passed, 2 skipped**, `ruff` clean.

**W1 owns this call.** It is drafted rather than committed because it is W1's
module and W1's test. The patch exists so the decision is a review rather than a
rewrite — and so that "we agreed on (b)" and "(b) is in the tree" do not end up
three weeks apart.

---

## 15 · CNV-based malignancy calling may fail differentially by MMR status — OPEN, AFFECTS THE PRE-REGISTERED CONTRAST

**Raised:** W1, 2026-08-22 (from the inferCNV pilot) · **Owner:** W1 + W4 +
whoever owns G2 · **Needed by:** before the MMR contrast is computed

inferCNV separates malignant from normal epithelium by aneuploidy. **MMR-
deficient tumours are characteristically near-diploid** — chromosomal stability
is what MSI-H looks like, and it is one of the oldest results in the field. So
the method is expected to work *worse* in MMRd than in MMRp, and MMRd-vs-MMRp is
a pre-registered contrast.

Aneuploid fraction on the pilot, as enrichment over the 0.10 null (the fraction
of a patient's tumour epithelium above its own held-out normal epithelium's 90th
percentile):

| patient | stratum | enrichment |
|---|---|---|
| C162 | **MMRp** | **7.4x** |
| C138 | MMRd, MLH1-intact | 2.9x |
| C107 | **MMRp** | 2.0x |
| C122 | MMRd, methylated | 1.7x |
| C165 | MMRd, MLH1-intact | **0.65x** — below the null |

**This is directionally consistent and NOT established.** Two MMRp against three
MMRd, and C107 (MMRp, 2.0x) sits between two MMRd patients. Nothing here would
survive a significance test and none is claimed. What makes it worth recording
is that the biological prior is independent of these five points.

### CONFIRMED AT FULL SCALE — 2026-08-24

62 patients through inferCNV. Restricted to the **36 with a matched normal**, so
a patient failing for want of a comparator is not counted as failing for
biology:

| stratum | separable | rate | median enrichment |
|---|---|---|---|
| `mmr_proficient` | 15 / 15 | **100%** | **7.38** |
| `mlh1_methylated` | 10 / 12 | 83% | 2.50 |
| `mlh1_intact_mmrd` | 4 / 6 | 67% | 2.26 |
| `mlh1_deficient_unmethylated` | 1 / 2 | 50% | 1.52 |

**Every MMR-proficient patient separates. None fails.** The MMRd strata fall away
monotonically, and the median aneuploid enrichment in MMRp is **~3x** the
methylated stratum's.

**The statistics, stated honestly.** All five separability failures land outside
MMRp; under random assignment that has probability **0.048** one-sided. The
direction was pre-specified in this document before the run, so a one-sided read
is legitimate — but n=35 and the two-sided p is not significant. The stronger
evidence is not the 2x2: it is that **four strata order monotonically, on both
the rate and the enrichment, in the direction an independent literature
predicts**. A single test on a collapsed table throws that away.

**So the concern this decision was opened on is real.** It is no longer a
prediction from MSI biology; it is measured in this cohort.

### Two consequences that now need deciding, not noting

1. **The MMRd arm loses ~25% of its patients to `not_called`, and the survivors
   are weaker.** Median enrichment 2.5 against MMRp's 7.4 means the surviving
   MMRd calls sit closer to their threshold — so the MMRd arm is both smaller
   *and* noisier, in a contrast where it is compared against an arm that lost
   nobody.
2. **Any comparison of malignancy-filtered results between MMR strata is
   confounded by this**, and no downstream method removes it. The filtered MMRd
   arm is a biased subsample of MMRd; the filtered MMRp arm is all of MMRp.

The sensitivity run proposed below — decomposition with and without malignancy
filtering — stops being optional. It is the only way to see how much of any MMR
difference is this artifact.

### Why it matters more than a QC wrinkle

If malignant cells cannot be separated in MMRd tumours, the MMRd "tumour" arm
retains more non-malignant epithelium. Non-malignant epithelium is *mature*, so
its retention **inflates the apparent mature fraction in MMRd tumours** — making
them look less compositionally depleted than they are.

That is a bias running **along** the axis being tested, not across it. A bias
that differs between arms of a pre-registered contrast is the one kind that
cannot be argued away as noise.

### It is not fixable by a better tool

CopyKAT has the same limitation, because it is the same measurement: both infer
copy number from expression. A genuinely near-diploid tumour has no aneuploid
population to find. This is an **information** limit, not a method choice, and
"try another caller" is not a mitigation.

### Two things the first real run exposed — 2026-08-23

**1 · `MALIGNANT_QUANTILE = 0.99` is a judgement, not a constant, and it is
conservative.** Malignant fractions on the pilot: C107 5.4%, C122 3.1%, C138
7.1%, C162 10.8%. Internally consistent — C107 has ~20% of query cells above the
copy-neutral 90th percentile and ~5% above the 99th — but low for colorectal
tumour samples, where a large share of tumour-sample epithelium is usually
malignant.

The cost of under-calling is not symmetric with over-calling here. **The
malignant set defines the tumour arm**, so a conservative threshold shrinks the
denominator of every compositional estimate and pushes patients toward
`not_estimable` on mature-cell counts — which then interacts with G4. The
threshold should be reported as a **sensitivity across quantiles**, the way the
depth target now is, rather than fixed at 0.99 and forgotten.

**2 · The specificity check has become nearly circular.** The threshold is the
99th percentile of copy-neutral epithelium, so a specificity near 0.99 on other
copy-neutral epithelium is what arithmetic predicts, not evidence. Observed:
0.989–0.999 across four patients. It now confirms the threshold was *applied*
correctly and nothing more.

Recovering an independent check means **splitting the held-out normal epithelium
in two** — one half sets the threshold, the other validates. At
`HOLDOUT_FRACTION = 0.30` that leaves 15% each, which C165 (117 held-out cells)
cannot spare. So it is a real trade and a team decision, not a code fix.

### Options

| | Approach | Cost |
|---|---|---|
| a | **Emit malignancy calls with per-patient confidence, and mark patients with no separated population as `not_estimable` for the malignant/normal distinction.** | The project's own three-way framing, applied one level down. Costs patients from the malignant-only analysis, and the loss is not random — which is the point, and must be reported. |
| b | Run the decomposition twice — malignancy-filtered and unfiltered — and report whether conclusions differ by MMR stratum. | Cheap, and it measures the exposure directly rather than assuming it. Doubles the result rows. |
| c | Use sample-of-origin instead of malignancy calls, as today. | No differential bias, but the tumour arm keeps its non-malignant epithelium in every patient, which is the gap the whole malignancy stage exists to close. |
| d | Ignore it. | The MMR contrast then carries an unquantified bias in a known direction. |

**Recommendation: (a) plus (b).** (a) is honest about which patients can carry
the distinction; (b) turns the residual exposure into a measured sensitivity
rather than a caveat. Together they cost one extra run and no credibility.

**Confirm at full scale before acting on it.** 36 matched patients, ~21 MMRd and
15 MMRp, is enough to see whether the association is real. Look at it once, on
the full cohort, and pre-specify that here rather than after.

---

## 16 · Ambient contamination: measure, do not correct — and exclude above 10% — PRE-COMMITTED 2026-08-23

**Raised:** W1 · **Owner:** W1 + W2 · **Status:** threshold committed **before**
counting what it costs

Measured across all 62 patients, 84 unsorted samples with >=20 epithelial cells
(`results/2026-08-23_*/ambient_contamination.parquet`):

| | |
|---|---|
| median | **2.2%** |
| 75th | 4.7% |
| 90th | **10.2%** |
| 95th | 13.8% |
| max | **19.4%** (C132_N) |

The five-patient pilot sat near the median and **understated the tail**.

### Two decisions, and the second is the one that needs pre-committing

**1 · Measure and report; do not correct.** At a 2.2% median, correction removes
very little and risks more than it removes: DecontX defines contamination as
counts resembling other clusters, so it can absorb genuine low-level expression
of a marker in a rare population — precisely this project's signal. SoupX and
DecontX still run, and the **per-gene retention table and their correlation**
remain the week-2 deliverable. What changes is that their output is reported as
a diagnostic rather than applied to the counts.

**2 · Exclude samples above 10% contamination from the compositional arm.**

The threshold is **10%**, and the reason is not the cost:

- Above roughly a tenth of counts being ambient, the per-cell marker detection
  that axis 1's maturity call depends on is materially perturbed. That call is a
  **detection gate** at a matched depth of 3,281 UMIs, so ~330 ambient counts
  per cell is not a rounding error in it.
- Having chosen not to correct, a sample whose ambient share approaches the
  effects being estimated cannot be rescued by a caveat.
- 10% is where the cohort's own distribution turns: it is the 90th percentile,
  and the gap from 75th (4.7%) to 90th (10.2%) is where samples stop resembling
  each other.

**This is committed before counting the patients it removes**, for the same
reason G1's threshold is: a threshold chosen after seeing its cost is not a
threshold. Nine of 84 samples exceed it; how many *patients* that costs — and
whether it costs them their matched arm — is to be reported, not to be used to
revise the number.

### The sharper criterion, once it exists

Total contamination is a proxy. What actually threatens this project is
**target-gene soup share** — ambient GUCA2A gives an immature cell false mature
counts and inflates the intrinsic term. The pilot's soup was dominated by
MT-CO2/CO3/CO1, IGKC and MALAT1, not by panel genes, which is reassuring but
unmeasured cohort-wide. When that is measured, it should replace this threshold
rather than supplement it — and **that replacement must also be pre-committed.**

### ANSWERED 2026-08-23 — the asymmetry is real but has no direction

**Measured, 23 patients with both arms interpretable:** 6 differ by more than 5
points, and the median tumour-minus-normal is **+1.5%**. The flagged patients
split evenly — C106 +12.3, C155 −12.1, C140 +9.2, C132 −8.8, C135 −5.3, C130
+5.2.

**Retracting the inference that prompted this.** "Four of the worst five samples
are tumour" was an artifact of reading the top five rows of a sorted table; it is
not a cohort pattern. Tumour samples are not systematically dirtier.

That distinction matters more than the flag count:

- **Systematic** asymmetry would bias the cohort-level compositional estimate in
  a known direction — the serious case.
- **Random** asymmetry inflates per-patient error and widens intervals without
  biasing the mean.

This is the second. So the case for switching #16 to a gap-based rule is
**weaker than anticipated and rests on variance, not bias**. Six patients carry
a >5-point gap and their individual Δ is correspondingly noisy; excluding them
buys precision, not correctness. That is a judgement about how much n to spend,
not a correction — and it should be argued as one.

**Recommendation: keep the 10% level rule as committed, and report the gap
alongside it** so a reader can see which patients are noisy. Revisit only if the
gap turns out to correlate with something the analysis cares about.

### CORROBORATED 2026-08-23 — SoupX and the impossible-gene estimator agree; DecontX does not

First real sample, C122_N_1_1_0_c1_v2 (1,609 cells, 55 clusters):

| route | contamination |
|---|---|
| impossible genes | **0.8%** |
| SoupX (degraded mode) | **0.8%** — median retention 0.992 |
| DecontX | **8.5%** — median retention 0.915 |

**Two unrelated routes land on the same number. DecontX removes ten times more.**

And they are not finding different genes: Spearman between the two retention
vectors is **0.71**, and the hardest-stripped list is textbook ambient —
`MT-ATP6`, `MT-CO3`, `RPL13`, `RPL18`, `RPS2`, `RPS23`, `TMSB4X`, `PTMA`,
`Metazoa_SRP`. Both methods identify the same soup. They disagree on how much of
it there is.

**The likely cause is the one the docstring warned about.** DecontX defines
contamination as counts resembling *other clusters*, and this sample has 55 of
them. With that many, a large amount of genuine cell-type-specific expression is
indistinguishable from cross-cluster bleed. A mature marker expressed in one
small population looks exactly like soup to that model — which is this project's
signal.

**Three consequences:**

1. **The measure-and-report decision is vindicated on data, not on argument.**
   Correcting with DecontX would have removed ~8.5% of every gene's counts on
   the strength of a number that two independent routes put under 1%.
2. **The 2.2% cohort median that this decision's 10% threshold rests on is
   corroborated**, since SoupX agrees with the estimator that produced it.
3. **If anyone ever does want correction, SoupX is the defensible choice here**
   and DecontX is not — on this deposit, with this clustering. That is a
   finding about the data, not a general claim about the methods.

### SECOND SAMPLE — the cluster-count hypothesis is wrong, and the real pattern is sharper

C162_T_0_0_0_c1_v3, 4,128 cells, **84** clusters against C122's 55:

| | C122 normal (55 clusters) | C162 tumour (84 clusters) |
|---|---|---|
| impossible genes | 0.8% | **1.7%** |
| SoupX | 0.8% | **1.1%** |
| DecontX | 8.5% | **6.7%** |

**Clusters went up and DecontX's over-removal went down**, so it does not scale
with cluster count. That hypothesis is retracted.

The pattern across the two samples is more informative than the one that was
looked for. The independent estimate **doubled**, 0.8% to 1.7%. SoupX tracked it.
**DecontX did not move** — it returned 7-8% both times, on samples whose actual
contamination differs by a factor of two.

An estimate that does not respond to the quantity it is estimating is dominated
by something other than the data. Whatever the mechanism, the operational
conclusion is firmer than "DecontX over-corrects": **on this deposit DecontX's
contamination estimate does not track contamination**, and it should not be used
to size a correction here.

One further difference worth noting rather than explaining: SoupX's
hardest-stripped genes on the tumour sample are small RNAs — `SNORA81` (0.47),
`Metazoa_SRP`, `7SK`, `U6`, `SNORD70`, `SNORA46` — and DecontX leaves every one
of them untouched at retention 1.000. The two methods disagree about a whole
class of transcript, not only about magnitude.

**Two samples is two samples.** Both statements above are consistent across a
normal and a tumour arm from different patients, chemistries and cluster counts,
which is more than nothing and less than a result. They are recorded as
observations to check when the cohort-wide run happens, not as established.

### WHAT THE THRESHOLD COST — 2026-08-25, measured

The rule was committed before this was counted. Here is the count.

| arm lost | patients | cause |
|---|---|---|
| normal | C112, C114, C116, C155 | sorted-only (#11) or ambient (#16) |
| tumour | C106, C140 | **ambient (#16)** — tumour samples at 14.6% and 10.4% |

**Six patients leave the paired cohort: 34 usable become 28.**

Both tumour-arm losses are squarely #16's: those samples exceeded 10% and were
excluded, so the patients retain a normal arm and have no tumour to compare it
against. They are not "unmatched" and must not be reported as such.

**The paired n W1 delivers is 28** — against 32 matched-and-unsorted, 36 matched,
and the ~60 §8.4 assumed. Each reduction has a different cause and a different
remedy, and collapsing them into one number is how a cohort quietly shrinks:

| n | what it counts |
|---|---|
| 62 | patients in the deposit |
| 36 | with a matched normal sample (#9) |
| 32 | matched **and** unsorted in both arms (#11) |
| 28 | **and** both arms under 10% contamination (#16) |

Whether 28 is enough is a separate question from whether the threshold is right,
and it should be argued separately. Lowering the threshold to recover C106 and
C140 would be revising a rule to buy back two patients after seeing their
names — which is the move this decision was pre-committed to prevent.

### A number nobody had looked at

**Only 23 patients have both arms among the interpretable samples**, against 32
matched-and-unsorted. Nine lose an arm to the ≥20-epithelial-cell floor or to an
unestimable ratio. If the ambient exclusions then remove more, the compositional
n falls further — and 23 is already well below the 36 that decision #9 reports
and the ~60 §8.4 assumed. **Whatever the final rule, the paired n it leaves
should be stated in the same sentence as the rule.**

---

## 17 · G1's threshold, committed before any G1 number exists — PRE-COMMITTED 2026-08-25

> ### ⚠ SUPERSEDED — see [Amendment 2](prereg_amendment_2_g1_tier_d.md), 2026-08-25
>
> **Two independent defects, both arithmetic, neither found by looking at a
> result.** `checks.py` has never been run against expression.
>
> **1. The statistic cannot pass.** "Apparent loss" as a raw Δ per-cell mean
> carries abundance inside it — a gene averaging 100 counts can lose 30, a gene
> averaging 0.1 cannot. Simulated on 20,000 genes where *every* gene loses
> exactly 30% and the truth has no abundance dependence at all, the statistic
> below returns **ρ = −0.997**. G1 fails at |ρ| > 0.5, so as pre-registered it
> fails whatever the biology. A gate that always fails carries as much
> information as one that always passes.
>
> **2. Tier D holds one gene.** #17 justifies its 0.2 threshold with "n≈8 genes
> per tier". `config/panel.yaml` has 4 in tier A, 3 in tier B and **1 in tier D**
> — and a Spearman over one gene is undefined, not merely noisy. Tier D is the
> half of the gate carrying the falsification logic.
>
> Amendment 2 proposes the standard MA construction (M = log₂ ratio, A = mean
> log abundance) and moves the unit of analysis genome-wide with within-bin
> percentiles, since with one gene you cannot compute a correlation but you can
> compute a percentile. It commits three replacement thresholds.
>
> Until the team ratifies it, `g1_verdict()` returns `not_estimable` and G1 is
> undecided. **Everything below stands as the record of what was committed
> first**; it is not edited, because the ordering is the point.


**Raised:** W1 · **Owner:** W1 + whoever owns the gate · **Status:** committed
before `checks.py` was written, let alone run

G1 asks whether the residual signal is ambient RNA rather than biology. Ambient
counts are enriched for whatever is abundant, so **if a gene's apparent loss
tracks its abundance, the loss is a property of the soup and not of the tumour.**

### CORRECTION 2026-08-25 — this changed G1's statistic, and said so late

**execution_plan.md §4 specifies G1 as "post-correction *retention* vs total
abundance".** The version below uses abundance vs **apparent loss**. That is a
different measurement and the substitution was not flagged when it was written.
Recording it now rather than letting it stand.

They answer different questions, and both are worth having:

- **Retention vs abundance** (the plan's) asks whether *the correction* is
  abundance-driven. On this cohort it is close to tautological — soup is
  enriched for abundant genes, so a correction that removes soup will strip
  abundant genes hardest, and finding that tells you little.
- **Loss vs abundance** (below) asks whether *the project's signal* is
  abundance-driven. That is the question G1's own stated consequence is about:
  "if retention tracks abundance across all tiers, the residual signal is soup."

**Recommendation: run both, report both, and let the plan's version be the
named gate criterion.** Substituting a better statistic for a pre-registered one
is exactly the move this project refuses elsewhere; adding a second is not.
Where they disagree, that disagreement is the finding.

The thresholds below apply to whichever is being read. **The team should ratify
this before `checks.py` runs**, because after it runs the choice is
unfalsifiable.

### The statistic

Spearman correlation between **gene abundance** (mean expression across the
cohort) and **apparent loss** (Δ per-cell mean, tumour minus normal), computed
**within each panel tier separately** — A (compositional targets), B (intrinsic
targets), D (neither, the negative control).

Spearman rather than Pearson for the same reason as the retention comparison:
abundance spans orders of magnitude and a handful of very high genes would
otherwise decide the answer.

### The pre-committed thresholds

**G1 FAILS if either holds:**

1. **|ρ| > 0.5 in tier D.** Tier D genes are chosen to have no differentiation
   story. A strong abundance-loss relationship *there* has no biological reading
   left — it is the soup, measured.
2. **The three tier correlations fall within 0.2 of one another.** That is what
   "tracks abundance across all tiers" means: if A, B and D behave alike, the
   panel is measuring abundance and the tier structure — which is the whole
   falsification design — carries no information.

**G1 PASSES if** tier D is flat (|ρ| ≤ 0.5) **and** tiers A and B separate from
D by more than 0.2.

### Why these numbers

0.5 is the same rank-correlation line already used for method agreement in #16,
so the project uses one meaning of "these two things track each other" rather
than a different one per test. 0.2 is the smallest separation that survives
n≈8 genes per tier — below that, tier differences are not distinguishable from
noise at this panel size, and pretending otherwise would manufacture a pass.

**Neither number was chosen by looking at a G1 result, because none exists.**
When `checks.py` runs, its output is compared against this and the comparison is
reported whichever way it goes.

### What a failure would mean

Not that the project is wrong — that **this cohort cannot separate the signal
from the soup**, and the honest report is the non-identifiability, with the
diagnostics, as the result. That is the same consequence G4 carries, and the
same three-way framing the whole project rests on.

Recorded here rather than in code so the commitment has a date and a diff.

---

## Closed

*(none yet — move entries here with the date and the decision, do not delete them)*

---

## 13 · Bulk GUCA2A is continuous, not bimodal — CONFIRMED AFTER PURITY ADJUSTMENT

**Raised:** W3, 2026-08-17 · **Confirmed:** W3, 2026-08-18, purity-conditioned re-run · **Owner:** all four ·
**Needed by:** week 2 standing meeting, and before W2 calibrates cutpoints (#7)

The premise check ran. **Hartigan's dip test finds no multimodality in GUCA2A in
any stratum** — p = 0.851 in COAD tumours (n=458), 0.919 in READ (n=166), 0.982
pooled. The loss is large (median 2.8 vs 9.0 log2CPM, ~64-fold) and entirely
graded. Zero inflation is ruled out at 0.2%. CDX2 is unimodal *and* barely lost
at the median (8.08 vs 8.31).

BIC prefers two components in most strata. That is the skew artifact the brief
warned about, not evidence of two groups — see
[the note](../results/notes/w3.2_premise_check.md).

execution_plan.md §8.2 pre-registered the consequence: *"if continuous, the
two-type classification dissolves into a regression."* So:

1. **Does the bulk arm become explicitly a regression?** There is no natural cut
   point, so any threshold is a modelling choice and should be reported as one.
2. **Does this change what W2 calibrates against (#7)?** Cutpoints for the
   *positivity* threshold are about mature-cell counts and are unaffected. But
   anything calibrating a bulk-phenotype cut is calibrating against a continuum.
3. **It does not touch G1–G4.** This is a statement about bulk GUCA2A, not about
   whether the decomposition separates compositional from intrinsic.

**CONFIRMED 2026-08-18 — no longer provisional.** The re-run conditioned on
W3.3 purity is done and the finding is unchanged. Purity explains **1.9%** of
GUCA2A variance (ABSOLUTE, n=556). The smallest dip p-value across every
adjustment — residualised on ABSOLUTE, residualised on ESTIMATE, and within
each purity tertile — is **0.332**. Nothing approaches significance.

A further reason to prefer ABSOLUTE in the W3.6 lock: GUCA2A correlates with
ABSOLUTE purity at r=-0.139 (purer tumour, less GUCA2A — the direction the
compositional hypothesis predicts) but at r=+0.064 against ESTIMATE, which is
near zero and the wrong sign. See the addendum in
[the note](../results/notes/w3.2_premise_check.md).

---

## 14 · Plate explains more expression variance than any biological variable — OPEN

**Raised:** W3, 2026-08-18, from W3.4 · **Owner:** W3 + W2 · **Needed by:** week 4–5,
before the covariate set is locked (W3.6)

PVCA-style variance analysis over 624 TCGA tumours, permutation-nulled:

| Factor | Excess variance over null | p |
|---|---|---|
| **plate** | **0.132** | 0.010 |
| **TSS** | **0.084** | 0.010 |
| msi_status | 0.042 | 0.010 |
| site | 0.023 | 0.010 |
| stage | 0.007 | 0.010 |
| vial (negative control) | −0.001 | 0.683 |

**Plate is associated with three times more expression variance than MSI status
and twenty times more than stage.** And plate is confounded with MSI
(Cramér's V = 0.21, permutation p = 0.006), which is the project's single
pre-registered subgroup variable.

Two decisions follow, neither taken here:

1. **Does plate join the locked covariate set?** It is not on the current list
   (stage, age, sex, MMR/MSI, purity, site) and on this evidence has a stronger
   claim than stage. But 29 levels against 624 samples is a lot of degrees of
   freedom — a random effect or a coarser grouping is probably wanted rather
   than 28 fixed-effect dummies.
2. **How is the MSI subgroup contrast reported?** With 76 MSI patients, a plate
   confound and no batch correction, an unqualified MSI contrast overstates the
   evidence. Minimum: report the plate association alongside it.

**Invariant 4 is not in question.** The response to a measured confound is to
carry it into the model, not to remove it from the data.

See [the note](../results/notes/w3.4_batch_structure.md).

---

## 15 · TSS and COAD/READ are nearly the same variable — bears on #4

**Raised:** W3, 2026-08-18, from W3.4 · **Owner:** W3 · **Needed by:** week 3,
alongside decision #4

Tissue source site against project: **Cramér's V = 0.971, permutation p = 0.001.**
A hospital submits colon cases or rectal cases, essentially never both.

So **a COAD-vs-READ contrast is also an institution contrast**, and adding a
project covariate does not adjust for site because the two are nearly the same
column. This does not by itself settle [#4](#4--coad-and-read-pool-or-stratify--open),
but whichever way it goes should be decided knowing the two cannot be separated
in this cohort.

MSI coverage, also from W3.4, resolves the brief's other escalation: **98.6%
annotated** (523 MSS, 76 MSI, 16 conflicting, 9 missing). Coverage was never the
binding constraint — the 76-patient subgroup size is.

---

## 16 · The CDR calls DSS "approximated" for COAD/READ — invariant 9 needs a look

**Raised:** W3, 2026-08-18, from W3.5 · **Owner:** all four ·
**Needed by:** week 5, before W3.7 reports anything

[CLAUDE.md](../CLAUDE.md) invariant 9 makes **DSS and PFI primary, OS
secondary**. The brief also says to honour the CDR's own recommended-use flags.
The CDR's notes sheet says:

> "we recommend the use of **PFI** [...] and **OS** [...] Given the relatively
> short follow-up time, PFI is preferred over OS."
>
> "DSS is relatively accurate for CESC, PAAD, and UVM, and is **approximated for
> other tumor types**."

COAD and READ are "other tumor types". The CDR derives DSS as *dead and with
tumour*, and states that a with-tumour patient who dies of an unrelated cause is
"incorrectly considered as an event" — the same contamination invariant 9
objects to in OS, reduced but not removed.

**Invariant 9's reasoning survives**: DSS is still better than OS for COAD. Two
things argue for leading with PFI regardless:

1. Both the project and the CDR endorse PFI.
2. **DSS has 75 events; PFI has 156.** The locked covariate set expands to ~11
   degrees of freedom, so DSS gives under 7 events/df — below the conventional
   floor of 10. PFI gives ~14.

**Options:** (a) keep invariant 9 verbatim and lead with PFI in practice, which
is what the code does today; (b) amend invariant 9 to name PFI primary and DSS
co-primary-with-caveat, which needs a PR and two approvals.

**Recommendation: (a) now, decide (b) at the gate.** Nothing is blocked either
way; what must not happen is a DSS-led headline result that does not mention the
approximation.

See [the note](../results/notes/w3.5_clinical_table.md).

---

## 17 · W3.6 covariate set — LOCKED 2026-08-18

**Raised:** W3, 2026-08-18 · **Locked:** W3 owner (jeremy749), 2026-08-18, commit
53442b0 · **Owner:** W3 (owner confirms) · **Was needed by:** week 5, before W3.7 ran

[`config/covariate_set.yaml`](../config/covariate_set.yaml) was committed with
`status: proposed`; `src/bulk/covariates.py:require_locked` refused to let any
survival model run against it until it flipped to `locked` in its own commit.
That flip is 53442b0, authorised by the W3 owner and recorded in the file's
`lock_authorisation` field. [W3.7](../results/notes/w3.7_baseline_survival.md)
ran after it, not before.

**This closes #17 and nothing else.** The lock carries answers to three
decisions that are formally the team's, and locking the config did not ratify
them — it implemented them so W3.7 could run. **#4, #14 and #16 stay open until
the group says otherwise**, and overruling any of them is a config change with
its own commit and a stated reason:

- **#4 (pool or stratify)** → COAD/READ as **strata**, not a covariate. Rectal
  cancer usually gets neoadjuvant chemoradiation, so a shared baseline hazard
  with a proportional shift is not credible. Costs zero degrees of freedom, and
  since TSS is confounded with project at V = 0.971 (#15) it absorbs the
  institution effect too.
- **#14 (does plate join the set)** → **no** for the clinical baseline, **yes,
  required** for expression models. Plate affects expression *measurement*, and
  the clinical baseline contains no expression, so it is not a confounder there.
- **#16 (which endpoint leads)** → invariant 9 unchanged; PFI carries
  `lead: true`.

**The driving constraint:** the six covariates cost **ten** degrees of freedom,
not six. Applying all of them to DSS gives 59 events over 10 df. Excluding
purity from the clinical baseline — where it is not a confounder — and dropping
site from DSS gets every endpoint over the floor:

| Context | PFI | DSS | OS |
|---|---|---|---|
| clinical baseline | 15.8 ✅ | 11.2 ✅ | 12.2 ✅ |
| expression models | 12.4 ✅ | 8.1 ❌ | 9.7 ❌ |

Reasoning in full: [the note](../results/notes/w3.6_covariate_lock.md).

---

## 18 · The premise-check finding replicates independently — CLOSES the artifact worry on #13

**Raised:** W3, 2026-08-18, from W3.8 · **Owner:** all four · **Status:** informational

[Decision #13](#13--bulk-guca2a-is-continuous-not-bimodal--confirmed-after-purity-adjustment)
rests on a finding from one cohort. It has now been replicated in **GSE39582**
(Marisa et al., 566 tumours) on **Affymetrix microarray** rather than Illumina
RNA-seq:

| Gene | TCGA dip p | GSE39582 dip p | TCGA % of normal | GSE39582 % of normal |
|---|---|---|---|---|
| GUCA2A | 0.982 | 0.997 | 1.40% | 1.72% |
| CDX2 | 0.969 | 0.987 | 94.7% | 84.8% |
| MS4A12 | — | 0.992 | 0.82% | 0.35% |

Different patients, different country, different measurement technology — and
the effect sizes land within a factor of two on all three genes.

**The CIMP+ test is the substantive addition.** GSE39582 carries CpG island
methylator phenotype status, which TCGA's GDC clinical does not. Promoter
hypermethylation is the one mechanism that would plausibly produce a discrete
off-state, so CIMP+ tumours are where bimodality should hide. Dip p = 0.804
(GUCA2A), 0.515 (CDX2), 0.785 (MS4A12) across 91 patients. **Nothing there
either.**

So the classification-to-regression question in #13 is no longer contingent on
whether TCGA is representative. It is a decision about framing, on evidence from
two independent cohorts.

Cohorts are **not pooled** (invariant 4) — estimated separately, reported side
by side, using the same test code so a difference could not come from the
analysis. See [the note](../results/notes/w3.8_replication_gse39582.md).

---

## 20 · `unresolved_fraction` is measured, means "bounded not measured", and gates nothing — OPEN

**Raised:** W1, 2026-08-25 · **Owner:** W2 (`src/harness/positivity.py`) · **Needs:** the weekly

Numbered 20, not 18 or 19, to stay clear of the duplicate-numbering collision
flagged at the top of this file — W3 already holds #17 and #18.

### What was measured

From `results/2026-08-25_9e3ca1a/mature_cell_counts_full.parquet`, 928 rows,
30 patients (28 paired, plus C106 and C140 with one arm each):

```
unresolved_fraction   mean 0.310   median 0.257   max 0.923
rows with n_cells_resolved < 50      108 / 696   (excluding the epithelial rung)
C165 normal arm                      90.7% unresolved
C119 tumour                          n_cells_resolved = 0, mature_fraction NaN
```

### The gap

`src/reference/labels.py` computes `unresolved_fraction` and says what it means:
*"How much of the epithelium the fraction could not speak for. A large value
means the fraction is **bounded, not measured**."*

**Nothing anywhere consumes it.** Grepped across `src/`: the only other
`unresolved` hits are W3's unrelated panel resolution and W1's own pilot job.

`classify_estimability()` gates on `n_cells_mature` alone, which protects the
**intrinsic** arm — too few mature cells and you cannot ask about expression
within them. There is no matching gate on the **compositional** arm. A
`mature_fraction` of 0.92 computed on 9% of the epithelium is reported the same
way as one computed on 90%, and `classify_estimability`'s own docstring says
*"The compositional term is still estimable in that case — do not drop the row."*

The zero case is safe: `n_cells_resolved = 0` forces `n_cells_mature = 0`, which
falls below `wide=20` and classifies `not_estimable`. C119 is handled. **It is
the middle of the range that is not** — enough mature cells to pass the gate,
too few resolved cells for the fraction to mean much.

### Why it may not be benign

The unresolved cells are not missing at random. They are the cells whose labels
are ambiguous, and on a maturity axis that means intermediate and transitional
states — precisely the cells whose classification determines the fraction.

### What was checked and did NOT hold

**Tumour is not systematically harder to label.** `unresolved_fraction` is
constant within patient x arm, so the 168 paired rows are 28 patients counted six
times each. Tumour exceeds normal in 102/168 rows = **17 of 28 patients**,
binomial p ~ 0.34. Testing on rows rather than patients inflates that sixfold and
would have produced a spurious directional bias — CLAUDE.md invariant 5's
principle applied one level up.

### Recommendation

A second cutpoint on `n_cells_resolved` (or equivalently on
`unresolved_fraction`), pre-committed before it is applied, carrying the
compositional term the way `Cutpoints` carries the intrinsic one. **W1 is not
proposing the number** — `src/harness/positivity.py` is W2's file under
CONTRIBUTING §2, and the threshold should be chosen by whoever owns the gate.

W1 supplies `n_cells_resolved` and `unresolved_fraction` in the frozen output
already, so no W1 change is needed to act on this.
