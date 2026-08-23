# Project Execution Plan

## Decomposing differentiation-marker loss in colorectal cancer

**Compositional vs. cell-intrinsic, with explicit non-identifiability. Four-person parallel execution.**

*Derived from project plan v4 + Executive Brief reconciliation*

**Team:** W1 Reference · W2 Method · W3 Bulk & Clinical · W4 Estimator & Replication

> Weeks 0–5 are branch-independent. The week-5 gate decides what the project becomes.

---

## Contents

1. [What this project is](#1-what-this-project-is)
2. [Reconciliation — read before writing code](#2-reconciliation--read-before-writing-code)
3. [Week 0 — shared setup](#3-week-0--shared-setup)
4. [Workstream specifications, weeks 1–5](#4-workstream-specifications-weeks-15)
5. [Week 5 gate](#5-week-5-gate)
6. [Weeks 5–13 and beyond](#6-weeks-513-and-beyond)
7. [Every branch ends in a result](#7-every-branch-ends-in-a-result)
8. [Feasibility audit](#8-feasibility-audit)
9. [Data inventory](#9-data-inventory)
10. [Open questions](#10-open-questions)
11. [Standing meetings](#11-standing-meetings)

---

## 1. What this project is

Bulk expression is a product of two quantities: how many cells make a gene, and how much each one makes. This project separates those terms per patient in colorectal cancer, at multiple annotation resolutions, and — the part nobody else does — reports explicitly when the separation is **not identifiable** from the data available.

> **THE DELIVERABLE**
>
> A three-part split, not a two-part one: **compositional / cell-intrinsic / not-estimable**. Every existing method returns a number. None flags that the intrinsic estimate is meaningless in a tumour with no mature cells left.

### 1.1 Out of scope

- No treatment-selection claim. No differentiation therapy exists in CRC.
- No confident intrinsic fraction from bulk alone. This is decided empirically at the gate ([§2.2](#22-the-bulk-question-is-empirical-not-assumed)), not assumed.
- No deployable clinical assay. No wet-lab validation.

---

## 2. Reconciliation — read before writing code

### 2.1 Errors in the Executive Brief that must not be carried forward

| # | Problem | Correction |
|---|---------|------------|
| 1 | GUCA2A appears in the reference matrix and is the target | Target genes are excluded from all label and reference construction. Enforced by assertion in `build_signature()`. Otherwise silenced mature cells are read as absent mature cells — the classifier cannot detect the phenomenon it was built to detect. |
| 2 | Wnt index built from CTNNB1 / TCF7L2 transcript levels | CTNNB1 mRNA is not Wnt activity — activation is post-translational. Use a target signature: AXIN2, NKD1, ASCL2, TCF7, RNF43, LGR5. These overlap the stem axis; see [§2.3](#23-the-wnt--stem-circularity). |
| 3 | Two-column reference (stem, mature) applied to bulk | Bulk CRC is 30–60% non-epithelial. Reference must include stromal, immune and endothelial columns or stromal signal is absorbed arbitrarily. This is the CMS4 failure mode. |
| 4 | 11-gene signature used for deconvolution | ν-SVR robustness comes from high dimensionality. Use 500–2000 genes to estimate fractions; the panel is for interpretation only. |
| 5 | Hard >70% classification cutoffs | Cutpoints come from the simulation harness (W2), not from round numbers. |
| 6 | "Assuming bulk GUCA2A is negligible" | Untested premise. W3 checks the distribution in week 2. If loss is continuous rather than bimodal, this is a regression, not a classification. |
| 7 | Overall survival as endpoint | Use DSS and PFI from TCGA-CDR (Liu et al. 2018). COAD OS is heavily contaminated by non-cancer death. |
| 8 | "Feasibility 10/10 on a standard laptop" | Realistic: 32 GB RAM minimum, GPU strongly preferred for CellBender, several multi-hour steps. See [§8](#8-feasibility-audit). |

### 2.2 The bulk question is empirical, not assumed

Published attenuation on cell-type-specific expression recovery from bulk: fold changes compressed to roughly 60–80% of true magnitude, for a gene-dependent subset. The bias is **directional** — it shrinks the intrinsic term and leaves the compositional term intact, pushing toward "compositional," which is also our prior hypothesis.

We do not decide in advance whether bulk can recover the intrinsic term. The harness measures it. Output is a **calibration curve**: given a known split, how well does bulk recover it as a function of mature-cell fraction. Expected shape — extremes separable, middle band not. That curve is itself a result.

### 2.3 The Wnt / stem circularity

Wnt-target genes and stem-identity genes are largely the same genes, so "high mature fraction and high Wnt activity" is close to a contradiction in a single expression space. Mitigations, in order of preference:

1. **Test Wnt activity within differentiation-matched cells** (Stage 3). This is the per-cell test and it dissolves the circularity.
2. **Use non-overlapping targets** where possible (NKD1, RNF43, NOTUM); drop ASCL2 and LGR5 from the Wnt score when the stem axis is in play.
3. **Report the correlation** between Wnt score and stem fraction explicitly, so readers see the confound rather than inferring it.

---

## 3. Week 0 — shared setup

Half a day, everyone present. Three things get frozen before anyone writes analysis code. Retrofitting any of them costs weeks.

### 3.1 Panel freeze

| Tier | Genes | Expectation | Role |
|------|-------|-------------|------|
| **A · compositional control** | GUCA2A, GUCA2B, OTOP2, CA7 | compositional | BEST4+ program; under 5% of epithelium even in healthy colon |
| **B · intrinsic control** | MLH1, SFRP1, SFRP2 | intrinsic | Methylation-silenced in CIMP/MSI, in cells that stay epithelial |
| **C · hypothesis** | CDX2 | mixed | Loss varies within tumours, concentrated at buds |
| **D · retained control** | MS4A12 | neither | Colonocyte-restricted yet frequently maintained |
| **E · exploratory** | AQP8, CA1, CA2, CA4, CEACAM7, SLC26A3, FABP1, PIGR, KRT20, LGALS4, VIL1, SATB2, MUC2, TFF3 | — | Unspecified |

Tier B is load-bearing. MLH1 is broadly expressed across colonic epithelium, so its compositional term should be structurally near zero and its loss almost entirely intrinsic. A compositional control alone validates only the half nobody doubted.

> **FALSIFICATION RULE — WRITTEN DOWN NOW**
>
> If tiers A, B and D all return the same answer, the estimator is broken and no biological claim may be made. This test exists only because the panel has controls; it is unavailable to a single-gene study.

### 3.2 Labelling axes freeze

| Axis | Genes / basis | Notes |
|------|---------------|-------|
| **1 · distance from stem pole** | LGR5, ASCL2, MKI67, OLFM4, SMOC2 | Overlaps the Wnt score — see [§2.3](#23-the-wnt--stem-circularity) |
| **2 · opposite lineage** | MUC2, TFF3, SPDEF, ITLN1 | Goblet program, orthogonal to absorptive |
| **3 · different measurement** | Chromatin accessibility (Becker/Chang multiome); crypt position (Visium HD / Xenium) | Not transcript-based — the strongest defence against leakage |

Fully panel-independent labels probably do not exist; differentiation state in colonic epithelium largely *is* this program. The claim is **agreement across structurally different axes**, not independence.

> **SEQUENCING CONSTRAINT**
>
> Tier E consumes genes that would otherwise be candidate labels (SPIB, BEST4, KRT20, CA1, VIL1, SATB2). Expanding the panel first leaves nothing to label with. **Labels are frozen first, in week 0.**

### 3.3 Output schema freeze

Enforced in code, not in the write-up. One row per (patient, gene, granularity rung, labelling axis, study):

```python
patient_id            str
study_id              str
gene                  str
granularity_rung      enum{epithelial, lineage, crypt_position, best4}
labeling_axis         enum{stem_pole, opposite_lineage, chromatin, spatial}
n_cells_mature        int
compositional         float | None
intrinsic             float | None   # None != 0.0 — this distinction is the project
interaction           float | None
weighting             enum{normal, tumour, doubly_robust}
estimability          enum{ok, wide_interval, not_estimable}
ci_low, ci_high       float | None
```

`assert not (estimability == "not_estimable" and intrinsic is not None)` — in the writer, not in review. Reporting zero when the truthful answer is "not estimable" is the single most likely route to a wrong conclusion here.

### 3.4 Repo contract

```
data/       raw/ interim/ processed/   (gitignored; manifest with checksums)
src/        reference/  W1     harness/   W2
            bulk/       W3     estimator/ W4
            schema.py   shared — changed only by PR + 2 approvals
results/    parquet, versioned by date + git sha
env/        one conda env per workstream, pinned
```

W1 emits `S_matrix_{rung}_{version}.parquet` with a fixed gene index. W3 emits bulk on the same index. Integration is then a join, not a negotiation. Random seeds fixed and logged; every result carries the git sha that produced it.

---

## 4. Workstream specifications, weeks 1–5

The split works because each stream owns a distinct data artifact with a frozen interface. Three of the four start on day one with no dependencies.

| Stream | Owns | Produces | Blocked by |
|--------|------|----------|------------|
| **W1 Reference** | Single-cell primary (GSE178341) | Cleaned AnnData, cell labels, versioned S matrices at 4 granularities | nothing |
| **W2 Method** | Simulation harness, deconvolution bake-off | Verdict on whether the estimator works; calibrated cutpoints | W1 pilot, wk 2 |
| **W3 Bulk & clinical** | TCGA-COAD/READ | Harmonised bulk matrix, purity, curated survival covariates | nothing |
| **W4 Estimator** | Lee cohorts + Kitagawa estimator | Per-patient compositional / intrinsic / not-estimable labels | its own cohort only |

### W1 — Reference

**Owner:** strongest scRNA-seq person. Owns GSE178341 (Pelka et al. 2021), ~371k cells, 62 patients, matched normal, MMRp and MMRd.

| Wk | Task | Done when |
|----|------|-----------|
| 1 | Ingest, QC (per-study thresholds, not global), doublet detection (scDblFinder or Scrublet) | Cell counts by patient and tissue tabulated; QC thresholds documented with rationale |
| 1–2 | Pilot subset: 5 patients through the full pipeline | Handed to W2 by end of week 2 |
| 2 | Ambient correction: **SoupX and DecontX**, both, compared | Per-gene retention table; correlation between methods reported. **Restated 2026-08-22: was SoupX and CellBender.** CellBender requires unfiltered droplets, and open decision #8 established that none exist in any public source for GSE178341 — GEO ran dropletUtils upstream. DecontX (`bioconductor-celda`) replaces it: it needs no empty droplets, modelling each cell as a mixture of its own cluster's distribution and contamination from the others. The deliverable is unchanged — two methods, compared. |
| 2–3 | Malignant vs. normal epithelium: inferCNV, CopyKAT as cross-check | Per-cell malignancy call with confidence; normal epithelium not misread as tumour |
| 3–4 | Cell labels, axes 1 and 2, all four granularity rungs | Labels stored as separate columns, never overwriting each other |
| 4–5 | Build S matrices — including stromal, immune, endothelial columns | Versioned parquet, fixed gene index, one per rung |
| 5 | Cross-gene ambient check: post-correction retention vs. total abundance, across all tiers | Plot and statistic. If retention tracks abundance across all tiers, the residual signal is soup. |

**Gotchas**

- GSE178341 supplementary structure is awkward — budget a day just for parsing. Verify you have raw counts, not normalised values.
- inferCNV on 371k cells is slow. Subsample per patient, or run per-patient in parallel.
- ~~CellBender wants a GPU. Without one, budget overnight runs.~~ **CellBender cannot run on this deposit at all** — no unfiltered droplets exist publicly (decision #8). DecontX replaces it and is CPU-only. Note also that inferCNV is CPU-only, so the single-cell arm needs no GPU whatever; the constraint that actually binds is **disk**, and the project filesystem is 55 GB.
- Use backed AnnData for anything that does not need the full matrix in memory.

### W2 — Method

**Owner:** strongest ML/stats person. The most underrated stream — the only place true ground truth exists, and the stream that adjudicates the week-5 gate.

| Wk | Task | Done when |
|----|------|-----------|
| 1 | Harness design spec, reviewed by all four | Written before code. Defines what "known ground truth" means here. |
| 2–3 | Pseudobulk generator: held-out patients mixed at known fractions with known per-cell shifts | Can generate arbitrary (composition, intrinsic) pairs, including near-zero mature-cell edge cases |
| 3–4 | Deconvolution bake-off: ν-SVR, CIBERSORTx, MuSiC, BayesPrism, plain NNLS baseline | Ranked on fraction recovery (r, RMSE) against harness ground truth |
| 4 | 11-gene vs. full-signature comparison, quantified | Settles the Executive Brief claim empirically rather than by argument |
| 4–5 | Attenuation curve: intrinsic-term recovery from pseudobulk vs. mature-cell fraction | The [§2.2](#22-the-bulk-question-is-empirical-not-assumed) calibration curve — a publishable object on its own |
| 5 | Calibrated positivity cutpoints replacing the provisional values below | Cutpoints derived from where CI width crosses a stated threshold, not chosen |

**Provisional positivity cutpoints, until recalibrated**

| Cells available | What gets reported |
|-----------------|--------------------|
| n ≥ 50 | Patient contributes to the intrinsic estimate |
| 20 ≤ n < 50 | Wide-interval flag; sensitivity analysis with and without |
| n < 20 | Intrinsic term **undefined — not zero**. Compositional term still estimable. |
| >50% of patients under 20 | Non-identifiability becomes the headline result |

W2 also owns the negative controls: housekeeping genes (should show neither term) and within-patient label permutation (should destroy both terms).

### W3 — Bulk & clinical

**Owner:** can be a strong analyst without scRNA-seq background. Fully independent for the first five weeks.

| Wk | Task | Done when |
|----|------|-----------|
| 1 | GDC ingest, COAD + READ, STAR counts; normalisation (TPM and log-CPM both retained) | Matrix on the shared gene index |
| 2 | Premise check: distribution of GUCA2A and CDX2 in TCGA-COAD | Histogram and bimodality test. Report to the team in week 2 — this can redirect the project. |
| 2–3 | Tumor purity: ESTIMATE; pull precomputed ABSOLUTE calls where available | Purity per sample, with method noted |
| 3 | Batch and technical structure: plate, TSS, sequencing batch | Documented; confounding with stage and MMR tested |
| 3–4 | Clinical table from TCGA-CDR — DSS and PFI primary, OS secondary | Curated, with censoring rules explicit |
| 4–5 | Covariate set pre-specified and locked: stage, age, sex, MMR/MSI, purity, tumor site | Written down before any survival model is run |
| 5 | Baseline survival models on clinical covariates alone | Sanity check — do known effects reproduce? If stage is not prognostic, something is wrong upstream. |

**Gotchas**

- COAD and READ have different treatment patterns. Decide now whether to pool or stratify, and write down the reason.
- MSI status in TCGA is incompletely annotated. Check coverage before committing to it as the subgroup variable.
- Normal-adjacent samples in TCGA-COAD are few and not matched to all tumors. Do not assume pairing.

### W4 — Estimator & replication

**Owner:** strongest methods implementer. Owns the Kitagawa decomposition and the Lee cohorts (GSE132465 SMC, GSE144735 KUL3). Develops on Lee so it is not queued behind W1, and delivers independent replication as a by-product.

| Wk | Task | Done when |
|----|------|-----------|
| 1–2 | Lee cohorts ingest, QC, ambient correction (same pipeline shape as W1 — coordinate, do not share code prematurely) | Cells labelled, axes 1 and 2 |
| 2–3 | Kitagawa standardisation: both weightings plus interaction term reported separately | Unit-tested on synthetic data with analytically known answers |
| 3–4 | Doubly-robust reweighted version | Agreement with the plain version quantified |
| 4 | Cross-check against cacoa and QuasiMed | Correlation reported. **Not** CoCoA-diff — it assumes cell fractions are not a mediator, which assumes away the compositional arm. |
| 4–5 | Patient-level bootstrap (over patients, not cells); hierarchical model with patient as grouping factor | CIs that reflect the real unit of inference |
| 5 | First decomposition results on the Lee cohorts | Independent of W1 timeline — a second opinion at the gate |

**Estimator definition**

```
compositional = Δ(mature fraction) × normal per-cell mean
intrinsic     = tumour mature fraction × Δ(per-cell mean)
interaction   = Δ(mature fraction) × Δ(per-cell mean)   [reported separately, never folded in]
```

The split is not unique — normal-weighted and tumour-weighted give different answers, and the difference lives in the interaction term. Report both plus doubly-robust. This is Kitagawa (1955) demographic standardisation, **not** regression Oaxaca–Blinder. Multiple testing: Benjamini–Hochberg within tier, reported separately for each term.

---

## 5. Week 5 gate

Go / no-go. All four present, one afternoon. Four criteria, each with a pre-committed consequence. Decide here, not in week 12.

| # | Criterion | Pass | Fail → |
|---|-----------|------|--------|
| **G1** | Ambient correction does not eliminate the intrinsic signal. Post-correction retention does not track total abundance across tiers. | Proceed | Paper becomes a caution about a widely-run analysis. Pivot to snRNA-seq and spatial, which carry what remains. |
| **G2** | Control tiers separate on pilot data: GUCA2A compositional, MLH1 intrinsic, MS4A12 neither. | Proceed | Methods and validation paper; no biological claim. If MLH1 specifically fails — harness passes means detection floor in the data, harness also fails means broken estimator. Both reportable. |
| **G3** | Harness shows the estimator recovers known ground truth on pseudobulk. | Proceed | The estimator is the problem, not the biology. Fix, or report the failure mode. |
| **G4** | Fewer than 50% of patients fall below the mature-cell positivity threshold. | Proceed | Non-identifiability finding with diagnostics becomes the headline result, not a caveat. This is a real paper. |

> **ALSO DECIDED AT THE GATE**
>
> Whether to attempt cell-type-specific expression recovery from bulk in Stage 4, based on W2's attenuation curve. Not decided now, and not decided by argument.

---

## 6. Weeks 5–13 and beyond

### 6.1 Weeks 5–10 — main decomposition and bulk fractions

W1 and W4 jointly run the full decomposition on GSE178341 across all tiers, four granularity rungs and three labelling axes.

> **ARCHITECTURE RULE**
>
> **Estimate per study and meta-analyse. Do not pool.** Batch-correction methods work by removing between-dataset variation, and between-dataset variation is where the compositional signal lives — running the decomposition on an integrated embedding risks correcting away the measurement. Use integration for label transfer only, then random-effects meta-analysis, which yields between-study heterogeneity as a free robustness statistic.

W2 and W3 jointly run bulk deconvolution on TCGA with the full-signature S — fractions only at this stage — producing mature-colonocyte fraction per patient with CIs. W3 runs survival models in parallel with purity as a pre-specified covariate.

**Output by week 10:** the per-patient split with three-way estimability at every rung and axis, with meta-analytic CIs; plus mature-colonocyte fraction for roughly 620 TCGA patients.

### 6.2 Weeks 10–13 — integration, which is the paper

- **Does the granularity or labelling choice change the answer?** Report the split as a curve across four resolutions plus a continuous maturation-score version. A single point estimate would present a modelling choice as a measurement. If it swings, that divergence is the contribution.
- **Does bulk-derived fraction predict the single-cell-derived label?** Validated on pseudobulk from held-out patients (W2 owns), reported as the calibration curve from [§2.2](#22-the-bulk-question-is-empirical-not-assumed).
- **Does the fraction — rather than the marker — carry the prognostic signal in TCGA?** How much variance in bulk GUCA2A and CDX2 is explained by mature-colonocyte fraction alone? If most of it, then CDX2, which already influences adjuvant decisions, is substantially a differentiation-content readout. This is the strongest honest clinical claim available and it needs only the reliable half.
- **One pre-registered subgroup contrast:** MMR status, reported as an estimate with an interval rather than a test. The primary paired analysis is adequately powered at n≈60; interaction contrasts need roughly four times that. Everything else is labelled exploratory, plus a design calculation for what a confirmatory study would need.

### 6.3 Weeks 13+ — long-runway additions

- **Chromatin axis (axis 3).** Becker/Chang multiome — a labelling axis not made of transcripts at all. If the decomposition agrees across a transcriptional and a chromatin definition of maturity, the leakage objection largely dissolves. Closes open question 5.
- **Stage 3 mechanism.** Within differentiation-matched cells, correlate per-cell expression against a Wnt-target signature. Never done in human tissue. This is what elevates the work from methods to biology. Spatial validation also tests the CDX2 tumour-bud prediction directly.

---

## 7. Every branch ends in a result

| If this happens | The paper is |
|-----------------|--------------|
| Ambient correction removes the intrinsic signal, or retention tracks abundance | A caution about a widely-run analysis; snRNA-seq and spatial carry what remains |
| Most patients fall below the positivity threshold | A non-identifiability finding with diagnostics — the headline, not a caveat |
| MLH1 fails to come out intrinsic | Broken estimator (harness fails too) or detection floor in the data (harness passes) — both reportable, decided by the harness |
| Tiers A/B/D do not separate | A methods and validation paper; no biological claim |
| The split swings with granularity or labelling axis | That divergence is the contribution |
| Bulk attenuation is severe | The calibration curve is the result: how far bulk can get you, measured |
| Everything works | The partition, tier contrasts, MMR stratification, and the Stage 4 prognostic result |
| Someone publishes first | Granularity-dependence, labelling-axis analysis and the simulation benchmark all remain open |

---

## 8. Feasibility audit

### 8.1 Compute

| Task | Realistic requirement |
|------|------------------------|
| GSE178341 full load | 32 GB RAM (backed mode makes 16 GB workable but painful) |
| ~~CellBender~~ → DecontX | CellBender is unusable here (no empty droplets, decision #8). DecontX is CPU-only and runs per sample. |
| inferCNV, 371k cells | Hours to a day; parallelise per patient |
| Bootstrap × rungs × axes × tiers | Embarrassingly parallel; the combinatorics are the cost, not any single fit |
| TCGA-COAD/READ bulk | Trivial by comparison — laptop-fine |

The Executive Brief's "10/10, standard consumer laptop" is wrong for the single-cell arm. Budget one machine with 32–64 GB and GPU access, or cloud equivalent. **Confirm this before week 1** — it is the only hard infrastructure dependency.

### 8.2 Things that could invalidate the plan, checkable in weeks 1–2

| Check | Owner | Consequence if it fails |
|-------|-------|--------------------------|
| ICBI atlas metadata: paired tumour/normal count, epithelial fraction by study, platform mix | W1 | Gives the real sample size. Also reveals whether the plate-based subset is large enough as an ambient-free validation set — plate protocols have essentially no soup, so an intrinsic signal surviving there is strong evidence it is not contamination. |
| GSE178341 raw counts availability and matched-normal completeness | W1 | If normals are sparse, the compositional term loses its reference |
| Bimodality of GUCA2A in TCGA | W3 | If continuous, the two-type classification dissolves into a regression |
| MSI/MMR annotation coverage in TCGA | W3 | If sparse, the pre-registered subgroup variable changes |
| Becker/Chang multiome accessibility and CRC overlap | W4 | Axis 3 may need to be spatial-only |

> **WEEK ONE, BEFORE ANYTHING ELSE**
>
> Pull the ICBI atlas **metadata table only** — not the 4.27M-cell object. That single table gives you the real sample size before you commit any compute.

### 8.3 Timeline honesty

Stage 1 was widened from four weeks to eight in plan v4 for a reason: harmonising annotation vocabularies, chemistries and QC conventions across atlases is six to eight weeks on its own. The five-week gate here is tight because W1 delivers a five-patient pilot at week 2 rather than the full object. **If that pilot slips past week 3, move the gate to week 7 rather than compressing W2's harness work** — the harness is what makes the gate meaningful.

### 8.4 Limitations, to be stated in the paper

- The composition-versus-intrinsic concept is established in neurodegeneration and more maturely in DNA methylation (TCA, CellDMC). Only the CRC application, the resolution and labelling dependence, and the identifiability treatment are new.
- BEST4+ depletion in CRC is already reported descriptively. Only the partition is new.
- Fully panel-independent labels may not exist; agreement across axes is a mitigation, not a proof.
- The split is annotation-relative by construction. There is no single true value — only a family indexed by resolution and axis.
- n≈60 supports the primary paired analysis, not multi-way stratification. Panel genes are co-regulated in blocks, so effective sample size stays small at any panel size.
- No wet lab. For the clinical framing: a trial-enrichment hypothesis, not a treatment-selection claim.

---

## 9. Data inventory

| Resource | Use | Owner |
|----------|-----|-------|
| Pelka GSE178341 | 371k cells, 62 patients, MMRp and MMRd, matched normal — primary | W1 |
| Lee GSE132465 / GSE144735 | SMC and KUL3 cohorts with matched normal — replication | W4 |
| ICBI integrated CRC atlas | 4.27M cells, 650 patients, 48 studies. Metadata first, object later. | W1 |
| HTAN / Vanderbilt polyp atlas | Conventional vs. serrated; crypt-top colonocytes pre-annotated | W1 |
| Joanito 2022 | iCMS subtypes with Wnt/MYC annotation | W4 |
| Becker/Chang multiome | Chromatin labels (axis 3), locus accessibility | W4, wk 13+ |
| Visium HD / Xenium | Spatial validation, crypt-position labels | W1, wk 13+ |
| TCGA-COAD/READ + TCGA-CDR | Bulk plus survival | W3 |

---

## 10. Open questions

1. Should positivity cutpoints come entirely from the harness rather than being set a priori? *(Leaning yes — W2 answers this at the gate.)*
2. Do SFRP1 and SFRP2 survive their stromal-expression problem, or should tier B rest on MLH1 alone?
3. Is one pre-registered subgroup contrast too conservative if iCMS is more biologically informative than MMR?
4. Does the granularity and labelling analysis stand alone as a methods note, or does splitting it weaken the main paper?
5. Is the chromatin axis worth the effort given noise levels? *(With the long runway: yes — see [§6.3](#63-weeks-13--long-runway-additions).)*

---

## 11. Standing meetings

| Cadence | Who | Purpose |
|---------|-----|---------|
| Weekly, 30 min | All four | Blockers only. Not status theatre. |
| Week 2 | W3 → all | Premise check result (GUCA2A distribution). Can redirect the project. |
| Week 2 | W1 → W2 | Pilot handoff |
| Week 5 | All four, half day | The gate. G1–G4 decided against pre-committed criteria. |
| Week 10 | All four | Integration planning |

Anything that changes `schema.py`, the frozen panel or the labelling axes needs **two approvals and a written reason**. The freeze is the point.
