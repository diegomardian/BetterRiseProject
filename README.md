# When a gene disappears from a tumor, did the cells stop making it — or did the cells leave?

A computational decomposition of differentiation-marker loss in colorectal cancer, scoped to what public data can actually support.

```
┌──────────────────────────────┬──────────────────────────┬─────────────────┐
│  COMPOSITIONAL — cells gone  │ CELL-INTRINSIC — gene off │  NOT ESTIMABLE  │
└──────────────────────────────┴──────────────────────────┴─────────────────┘
```

**Diego · Bode · Jeremy · Diego**

The bar above is the deliverable. Not a p-value, not a classifier — a per-patient split with an honest third segment for the part the data cannot resolve. Most methods in this space return two numbers and no third segment. **That omission is the opening.**

> Week-by-week execution, workstream ownership and the week-5 gate live in **[execution_plan.md](execution_plan.md)**.

---

## The problem, in one figure

Bulk RNA sequencing grinds tissue up and averages everything in it. What comes back is a single number per sample, and that number is a product:

```
bulk signal = (how many cells make the gene) × (how much each one makes)
```

A product can fall two ways. In normal colon, mature absorptive cells near the crypt top produce markers like GUCA2A. Tumors regress — they fill with stem-like cells and lose the mature ones. So when the signal vanishes, either the producers are gone, or the producers are present and silent.

| | NORMAL COLON | TUMOUR A | TUMOUR B |
|---|---|---|---|
| **Mature cells** | present | gone | present but silenced |
| **Bulk reading** | high | ~zero | ~zero |

**A and B are indistinguishable to any measurement that averages.** The two tumours have different biology and identical readings.

This is not a hypothetical failure mode. CMS4 — the poor-prognosis "mesenchymal" colorectal subtype — was defined from bulk expression and later shown to derive substantially from stromal cells rather than tumour cells. The classifier was already in wide use. The same confound is suspected in ovarian cancer subtyping.

---

## What the data can actually resolve

Single-cell sequencing separates the two factors because it measures cells individually: you can count the mature ones and read each one's output. Deconvolution attempts the same inversion on bulk data, using single-cell profiles as a reference. **It works for one half and not the other.**

| Quantity | Single-cell data | Bulk data |
|---|---|---|
| **How many cells** | Counted directly · reliable | Deconvolution · r ≈ 0.92 |
| **How much each makes** | Measured directly · reliable | Residual · fold changes ×0.6–0.8 |

Cell-fraction estimates from bulk correlate ≈0.92 with ground truth. Cell-type-specific expression estimates come back with fold changes compressed to roughly 60–80% of true size, for a subset of genes only.

> **The bias points the wrong way**
>
> Attenuation shrinks the intrinsic term while leaving the compositional term intact. A bulk-trained model would drift toward reporting "compositional" regardless of truth — which is also the prior hypothesis. A result that confirms your expectation for methodological reasons is the worst kind of result, and this is the specific mechanism by which it would happen here.

---

## What this project claims, and what it does not

**In scope**

- Per-patient compositional / intrinsic / not-estimable split from single-cell atlases
- How that split shifts with annotation granularity and labelling choice
- Compositional fraction estimated from bulk at population scale, with survival
- Whether existing bulk prognostic signals track differentiation content rather than tumour biology

**Out of scope**

- Any treatment-selection claim — no differentiation therapy exists in CRC
- Intrinsic fraction from bulk alone as a confident number
- A deployable clinical assay
- Anything requiring wet-lab validation

The actionability research was blunt: there is no 2026 clinical decision that flips on this. GUCY2C therapies currently in trials target the receptor, which tumours retain regardless of mechanism. CDX2 — a far more mature differentiation marker — is used mechanism-agnostically and isn't guideline-endorsed. The defensible destinations are **trial enrichment** and **reinterpreting existing signatures**.

---

## Prior art: the honest position

| Component | Status |
|---|---|
| Cell fractions from bulk | **Solved, crowded.** Dozens of methods. You use this, you don't invent it. |
| Cell-type expression from bulk | **Attempted, unreliable.** CIBERSORTx high-res, bMIND, swCAM, BayesPrism. Authors warn against exactly this use. |
| Composition vs. intrinsic framing | **Done elsewhere.** Alzheimer's "cell-intrinsic DEGs"; more mature still in DNA methylation (TCA, CellDMC). |
| BEST4+ depletion in CRC | **Reported.** Descriptively, in a Dec 2025 spatial atlas preprint. |
| Per-marker partition with identifiability limits, in CRC | **Open.** Not found. |

So the novelty is not "a model that decomposes bulk into two mechanisms" — that is assembly plus a borrowed framing. It is the part everyone else skips: **reporting when the answer isn't knowable.** Every existing method returns a number. None flags that the intrinsic estimate is meaningless in a tumour with no mature cells left.

---

## Design decision 1 — the panel is a control set

Not a survey of interesting genes. Tiers with written-down expectations, so the estimator can fail visibly rather than quietly.

| Tier | Genes | Expected result |
|---|---|---|
| **A · compositional control** | GUCA2A, GUCA2B, OTOP2, CA7 | **compositional** — BEST4+ program; those cells are <5% of epithelium even in healthy colon |
| **B · intrinsic control** | MLH1, SFRP1, SFRP2 | **intrinsic** — methylation-silenced in CIMP/MSI tumours, in cells that stay epithelial |
| **C · hypothesis** | CDX2 | **mixed** — loss varies within tumours, concentrated at tumour buds |
| **D · retained control** | MS4A12 | **neither** — colonocyte-restricted yet frequently maintained in tumours |
| **E · exploratory** | AQP8, CA1, CA2, CA4, CEACAM7, SLC26A3, FABP1, PIGR, KRT20, LGALS4, VIL1, SATB2, MUC2, TFF3 … | Unspecified |

Tier B is the load-bearing addition. MLH1 is expressed broadly across colonic epithelium rather than confined to mature cells, so its compositional term should be structurally near zero and its loss almost entirely intrinsic. A compositional control alone validates only the half nobody doubted.

> **Falsification rule**
>
> If tiers A, B and D all return the same answer, the estimator is broken and no biological claim may be made. This test only exists because the panel has controls — it is unavailable to a single-gene study.

---

## Design decision 2 — labels that don't leak

To ask "are the surviving mature cells silenced?", you must first decide which cells are mature — using something other than the gene you're testing. Leaving out just that one gene is not enough, because the panel is largely one co-regulated program. Hold out GUCA2A and its five partners still label the cells. CDX2 is worse: it sits upstream of MS4A12 and several tier-E genes, so using it to define maturity and then testing it is the same operation twice.

**Leaky:** panel labels panel — one co-regulated program. **Replace with three structurally different axes:**

| | Axis 1 | Axis 2 | Axis 3 |
|---|---|---|---|
| **Basis** | Distance from the stem pole | The opposite lineage | A different measurement |
| **Markers** | LGR5, ASCL2, MKI67 | MUC2, TFF3, SPDEF | chromatin · crypt position |

**Agreement = robust · divergence = the finding.**

Fully panel-independent labels may not exist — differentiation state in colonic epithelium largely *is* this program. The mitigation is agreement across structurally different axes, not a claim of independence.

> **Sequencing constraint**
>
> Labels are frozen before the panel expands. Tier E consumes genes (SPIB, BEST4, KRT20, CA1, VIL1, SATB2) that would otherwise be candidate labels. Expanding first leaves nothing to label with. **Panel freeze: end of week one.**

---

## Design decision 3 — granularity is a variable, not a setting

GUCA2A's cell of origin is genuinely disputed: single-cell atlases call it a BEST4+ marker, in-situ hybridisation puts it in goblet cells and colonocytes, and pseudotime shows it rising gradually with maturation rather than switching on at a cell-type boundary. Define the population narrowly and you are working inside a compartment that is under 5% of epithelium before any tumour depletion.

The deeper issue is that **a compositional change becomes an expression change purely by re-drawing cluster boundaries.** The split is annotation-relative and is not an objective property of the tissue. So report it as a curve across four resolutions — epithelial / lineage / crypt-position / BEST4+ — plus a continuous maturation-score version. A single point estimate would present a modelling choice as a measurement.

---

## Design decision 4 — the estimator, and the third segment

Framed as potential outcomes with cell state as mediator: compositional loss is the indirect effect through the cell-state distribution, intrinsic loss is the direct effect holding state fixed. Estimated transparently by Kitagawa standardisation (1955, demographic standardisation — **not** regression Oaxaca–Blinder, which decomposes by regression coefficients):

```
compositional = Δ(mature fraction) × normal per-cell mean
intrinsic     = tumour mature fraction × Δ(per-cell mean)
+ interaction term, reported separately
```

The split is not unique — weighting by normal or tumour values gives different answers, and the difference lives in the interaction term. Report both weightings plus a doubly robust reweighted version. Cross-check with cacoa and QuasiMed. **Not CoCoA-diff:** it explicitly assumes cell fractions are not a mediator, which would assume away the compositional arm.

### Positivity — where the grey segment comes from

If a tumour has essentially no mature cells left, "how much does each mature cell make" is not a hard question, it is an undefined one. Count the cells available, per patient, per gene, per rung, per axis, using labels that exclude the tested gene:

| Cells available | What gets reported |
|---|---|
| n ≥ 50 | Patient contributes to the intrinsic estimate |
| 20 ≤ n < 50 | Wide-interval flag; sensitivity analysis with and without |
| n < 20 | Intrinsic term **undefined — not zero.** Compositional term still estimable. |
| >50% of patients under 20 | Non-identifiability becomes the headline result |

Reporting zero when the truthful answer is "not estimable" is the single most likely route to a wrong conclusion here. The distinction is **enforced in the output schema, not left to the write-up.** Cut points are provisional and should be recalibrated from the simulation harness.

---

## Execution

| Stage | Weeks | Focus | Gate failure → |
|---|---|---|---|
| **1 · De-risk** | 1–8 | both controls gate | Ambient dominates, positivity fails, or MLH1 not intrinsic |
| **2 · Decompose** | 8–16 | all rungs, all axes | Tiers A/B/D don't separate |
| **3 · Mechanism** | 16+ | Wnt, spatial, subtype | — |
| **4 · Bulk + survival** | optional | fractions only | — |

*Both gate failures still yield a paper — see [Every branch ends in a result](#every-branch-ends-in-a-result).*

Stage 1 widened from four weeks. Harmonising annotation vocabularies, chemistries and QC conventions across atlases is six to eight weeks on its own.

### Stage 1 · weeks 1–8 · de-risk and validate both arms

Freeze panel and labelling axes by end of week one. Run SoupX and CellBender before any biological analysis — abundant epithelial transcripts leak into the droplet soup and appear in cells that never made them, which is precisely the signature that would fake an intrinsic result. Add the cross-gene check: plot post-correction retention against total abundance across the panel; **if retention tracks abundance across all tiers, the residual signal is soup.** Run inferCNV or CopyKAT so normal epithelium isn't misread as tumour. Then validate both controls — Milo/scCODA for the compositional arm, MLH1 for the intrinsic arm — plus housekeeping negatives and a within-patient label permutation.

### Stage 2 · weeks 8–16 · decomposition

Estimator across all tiers, every granularity rung, all three labelling axes. Bootstrap over patients, not cells; hierarchical model with patient as grouping factor; Benjamini–Hochberg within tier, reported separately for each term. Build the simulation harness — pseudobulk with known mixing fractions and known per-cell shifts — which is the only place true ground truth exists and the right place to calibrate the positivity cut points.

> **Architecture note**
>
> **Estimate per study and meta-analyse. Do not pool.** Batch-correction methods work by removing between-dataset variation, and between-dataset variation is where the compositional signal lives — running the decomposition on an integrated embedding risks correcting away the measurement. Use integration for label transfer only, then random-effects meta-analysis, which yields between-study heterogeneity as a free robustness statistic.

### Stage 3 · weeks 16+ · mechanism and stratification

Within differentiation-matched cells, correlate per-cell expression against a Wnt-target signature (AXIN2, NKD1, ASCL2, TCF7) — the per-cell test of the β-catenin silencing model, never done in human tissue. Note that **CTNNB1 transcript level is not Wnt activity**; β-catenin activation is post-translational. Validate spatially, which also tests the CDX2 tumour-bud prediction directly.

The unit of inference is the patient. The primary paired analysis is adequately powered at n≈60; stratification is not, since interaction contrasts need roughly four times the sample size. So: **one pre-registered subgroup contrast** — MMR status, best annotated and most balanced — reported as an estimate with an interval rather than a test. Everything else labelled exploratory, plus a design calculation for what a confirmatory study would need.

### Stage 4 · optional · bulk extension, fractions only

**No cell-type-specific expression imputation anywhere.** Deconvolve TCGA-COAD/READ for mature-colonocyte fraction, then ask how much of the variance in bulk GUCA2A and CDX2 that fraction alone explains, and whether the fraction rather than the marker carries the prognostic signal. This uses only the reliable half, needs the survival data single-cell cohorts lack, and supports the strongest honest claim available: that a marker influencing therapy decisions may be substantially a differentiation-content readout.

---

## Every branch ends in a result

| If this happens | The paper is |
|---|---|
| Ambient correction removes the intrinsic signal, or retention tracks abundance | A caution about a widely-run analysis; snRNA-seq and spatial carry what remains |
| Most patients fall below the positivity threshold | A non-identifiability finding with diagnostics — the headline, not a caveat |
| MLH1 fails to come out intrinsic | Either a broken estimator (harness fails too) or a detection floor in the data (harness passes, real data doesn't) — both reportable, decided by the harness, and written down now so it can't be argued away in week eight |
| Tiers A/B/D don't separate | A methods and validation paper; no biological claim |
| The split swings with granularity or labelling axis | That divergence is the contribution |
| Everything works | The partition, tier contrasts, MMR stratification, and the Stage 4 prognostic result |
| Someone publishes first | Granularity-dependence, labelling-axis analysis, and the simulation benchmark all remain open |

---

## Data

| Resource | Use |
|---|---|
| **Pelka GSE178341** | 371k cells, 62 patients, 28 MMRp / 34 MMRd, matched normal — primary |
| **Lee GSE132465 / GSE144735** | SMC and KUL3 cohorts with matched normal — replication |
| **HTAN / Vanderbilt polyp atlas** | Conventional vs. serrated; crypt-top colonocytes already annotated |
| **ICBI integrated CRC atlas** | 4.27M cells, 650 patients, 48 studies, cellxgene-browsable. Verify first: epithelial depth, matched-normal availability, platform mix |
| **Joanito 2022** | iCMS subtypes with Wnt/MYC annotation |
| **Becker/Chang multiome** | Chromatin labels (axis 3) and locus accessibility |
| **Visium HD / Xenium** | Spatial validation; crypt-position labels |
| **TCGA-COAD/READ** | Bulk plus survival — Stage 4 only |

> **Week one, before anything else**
>
> Pull the ICBI atlas **metadata table only.** Count patients with paired tumour and normal, tabulate epithelial fraction by study, tabulate platform. That single table gives you the real sample size, tells you whether the epithelial compartment is deep enough, and reveals whether the plate-based subset is large enough to serve as an ambient-free validation set — plate protocols have essentially no soup, so an intrinsic signal surviving there is strong evidence it isn't contamination.

---

## Open questions

1. Should the positivity cut points come entirely from the simulation harness rather than being set a priori?
2. Do SFRP1/SFRP2 survive their stromal-expression problem, or should tier B rest on MLH1 alone?
3. Is one pre-registered subgroup contrast too conservative if iCMS is more biologically informative than MMR?
4. Does the granularity-and-labelling analysis stand alone as a methods note, or does splitting it weaken the main paper?
5. Is the chromatin axis worth the effort given noise levels, or better deferred?

---

## Limitations, stated plainly

- The composition-versus-intrinsic concept is established in neurodegeneration and in methylation analysis. Only the CRC application, the resolution- and labelling-dependence, and the identifiability treatment are new.
- BEST4+ depletion in CRC is already reported descriptively. Only the partition is new.
- Fully panel-independent labels may not exist; agreement across axes is a mitigation, not a proof.
- The split is annotation-relative by construction. There is no single true value — only a family indexed by resolution and axis.
- n≈60 supports the primary paired analysis, not multi-way stratification. Gene-level covariate analysis is descriptive; panel genes are co-regulated in blocks, so effective sample size stays small at any panel size.
- No wet lab. The defensible output is a computational decomposition with stated limits — and, for the clinical framing, a trial-enrichment hypothesis rather than a treatment-selection claim.

---

## Repository

| Path | Contents |
|---|---|
| [README.md](README.md) | This document — what the project is and why it is designed this way |
| [execution_plan.md](execution_plan.md) | Four-person parallel execution: workstreams W1–W4, week-5 gate, frozen schema and repo contract |
| [CLAUDE.md](CLAUDE.md) | Invariants that must not be violated by any code change |
| [CONTRIBUTING.md](CONTRIBUTING.md) | **Start here.** Setup, branch naming, what you own, what you must not edit |
| [docs/open_decisions.md](docs/open_decisions.md) | Decisions that block code, and who closes them by when |

**Code**

| Path | Owner | Contents |
|---|---|---|
| [src/schema.py](src/schema.py) | shared · frozen | The output contract and the results writer. `None` is not `0.0`, enforced here. |
| [src/common/](src/common/) | shared | Paths, provenance stamping, loaders for the frozen panel and axes |
| [src/reference/](src/reference/README.md) | **W1** | GSE178341, QC, ambient correction, labels, S matrices |
| [src/harness/](src/harness/README.md) | **W2** | Simulation harness, deconvolution bake-off, calibrated cutpoints |
| [src/bulk/](src/bulk/README.md) | **W3** | TCGA-COAD/READ, purity, clinical table, survival |
| [src/estimator/](src/estimator/README.md) | **W4** | Kitagawa decomposition, bootstrap, Lee-cohort replication |
| [config/](config/) | shared · frozen | `panel.yaml`, `labeling_axes.yaml`, the fixed gene index |
| [tests/](tests/) | shared | The invariants, asserted. A red build means a frozen contract moved. |
| [data/](data/README.md) | — | Gitignored. Only [`manifest.csv`](data/manifest.csv) travels. |
| [results/](results/README.md) | — | Versioned parquet, one directory per date + git sha |
| [env/](env/README.md) | — | One pinned conda env per workstream |
