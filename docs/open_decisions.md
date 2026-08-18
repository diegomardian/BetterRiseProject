# Open decisions

Things the scaffold could not decide on its own, recorded so they get decided
deliberately rather than discovered in week 9. Close each with a line in the
weekly meeting and a PR that updates this file.

Open questions of a scientific kind live in
[execution_plan.md §10](../execution_plan.md#10-open-questions). This file is
for decisions that block or shape code.

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

## 17 · W3.6 covariate set — PROPOSED, awaiting confirmation

**Raised:** W3, 2026-08-18 · **Owner:** W3 (owner confirms) · **Needed by:** week 5,
before W3.7 runs

[`config/covariate_set.yaml`](../config/covariate_set.yaml) is committed with
`status: proposed`. `src/bulk/covariates.py:require_locked` refuses to let any
survival model run against it until that flips to `locked` in its own commit.

The proposal carries answers to three decisions that are formally the team's,
because leaving them open blocks W3.7 and each now has evidence:

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
