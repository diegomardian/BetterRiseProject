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

## 2 · Who owns the shared gene index — OPEN

**Raised:** scaffold, 2026-08-15 · **Owner:** W1 + W3 · **Needed by:** week 1

W1 emits S matrices on a fixed gene index; W3 emits bulk on the same index;
integration is a join. The plan does not say who *produces* the index, and both
need it in week 1. See [config/gene_index/README.md](../config/gene_index/README.md).

**Recommendation: W1 emits it in week 1** from the GSE178341 feature table,
committed as `config/gene_index/gene_index_1.0.0.txt`. W3 reindexes onto it
rather than the reverse — single-cell features are the more constrained set. If
W1's ingest slips, W3 emits a provisional index from the GDC gene model and W1
conforms. Decide in the week-1 meeting; do not let both build one.

---

## 3 · Symbol vs. Ensembl ID on the shared index — OPEN

**Raised:** scaffold, 2026-08-15 · **Owner:** W1 + W3 · **Needed by:** week 1

Follows from #2 and is the usual way a join silently loses 8% of genes. TCGA
STAR counts arrive as Ensembl IDs with versions; the panel and both labelling
axes are written as symbols; single-cell references are usually symbols.

**Recommendation: Ensembl IDs (unversioned) as the index, symbols as a mapped
column.** Panel and axis lookups resolve through the map, which is committed
alongside the index. Unmapped panel genes are a week-1 finding, not a week-9
surprise.

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

## 6 · Which cohort W2 ingests — OPEN

**Raised:** W2 planning, 2026-08-15 · **Owner:** W2 + W4 · **Needed by:** week 1,
before anyone downloads anything

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

## 8 · GSE178341 ships NO unfiltered droplets — ANSWERED, BLOCKING

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

**Recommendation: (c) as the default, with (b) as a pre-specified sensitivity
analysis if anyone wants it.** (c) costs nothing, is what the schema was designed
to express, and means a reader can see that 42% of the cohort could not be
decomposed instead of wondering where they went. Note that (c) needs a sign-off
that "no matched normal" is a legitimate `not_estimable` reason — the frozen
schema ties `estimability` to the intrinsic term and mature-cell counts, and
this is a different cause with the same consequence.

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

## Closed

*(none yet — move entries here with the date and the decision, do not delete them)*
