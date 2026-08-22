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

**2 · The hg19 warning was overstated.** "A silent 8% gene loss on the join" was
asserted, not measured. Unversioned ENSG identifiers are largely stable across
GRCh37 and GRCh38 — that is *why* decision #3 chose them — so the expected loss
comes from GENCODE release differences, not the assembly. It is still worth
measuring, and `src/reference/jobs/check_gene_index.py` now does.

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

**Recommendation: (b).** The assertion at line 98 is the one that protects the
reference pool; line 96 protects nothing extra and costs the shared index its
completeness. (a) is what is in place today and is a fine holding position.

---

## Closed

*(none yet — move entries here with the date and the decision, do not delete them)*
