# W1 — Reference

**Owner:** Bode · **Env:** `env/w1_reference.yml` → `conda activate brp-w1`
**Branch prefix:** `w1/…` · **Blocked by:** nothing

> ## STATE — 2026-08-20
>
> Week 1 complete; week 2's handoff artifact exists. 525 tests, ruff clean.
> Working branch `w1/ingest-gse178341`, PR #5.
>
> **Read `docs/open_decisions.md` first.** Seven decisions are open and every one
> shapes numbers that weeks 3-5 will produce.
>
> **The one number outstanding:** the last `run_pilot.py` added
> `annotation_concordance` and its Cohen's kappa has not been read. That kappa
> decides open decision #14, which blocks every compositional number. Run the
> pilot and read that block before anything else.
>
> ### Measured, not assumed
>
> | | |
> |---|---|
> | Cells / genes | 370,115 x 43,113 (published 371,223 — 1,108 unreconciled) |
> | **Matched tumour + normal** | **36 of 62** — not the ~60 §8.4 assumes |
> | **Unsorted in both arms** | **~30** — the real compositional n |
> | Ambient contamination | median 1.6%, max 4.1% |
> | QC retention | 91.5% |
>
> ### What the deposit actually is
>
> 10x CellRanger HDF5 v2, CSC, genes x barcodes, float64-but-integral, 764M
> nonzeros. Feature ids `ENSG00000243485.5_4` on **GRCh37_liftover_v28 — hg19**,
> while TCGA is GRCh38 (tell W3). Barcodes encode patient, tissue and chemistry.
> `TA`/`TB` are two tumour regions, not a tissue type. Chemistry is mixed but
> constant within 61 of 62 patients; `PROCESSING_TYPE` is **not** — mixed within
> 45 of 62, and `CD45pMACS` is immune enrichment. Already mito-filtered at 50%
> upstream, and **no unfiltered droplets exist in any public source** (#8), so
> CellBender cannot run and the ambient arm is SoupX + DecontX.
>
> ### Still stubbed, all deliberately judgement over real data
>
> `malignancy.run_infercnv` · `qc.flag_doublets` · `ambient.run_soupx` ·
> `ambient.run_decontx`
>
> ### Before any compositional number
>
> Malignancy calls (largest gap — the "tumour" arm still holds non-malignant
> epithelium, so the contrast is sample-of-origin), ambient subtraction, doublet
> removal, the 29-cell sample pooled with samples 100x larger, and the 1,108-cell
> discrepancy.

You own GSE178341 (Pelka et al. 2021): ~371k cells, 62 patients, matched normal,
MMRp and MMRd. Everything downstream is built on what you emit.

## What you deliver

| Wk | Task | Done when |
|----|------|-----------|
| 1 | Ingest, QC (per-study thresholds, **not** global), doublets (scDblFinder or Scrublet) | Cell counts by patient and tissue tabulated; QC thresholds documented with rationale |
| 1–2 | **Pilot: 5 patients through the full pipeline** | Handed to W2 by end of week 2 |
| 2 | Ambient correction: SoupX **and** CellBender, both, compared | Per-gene retention table; correlation between methods reported |
| 2–3 | Malignant vs. normal epithelium: inferCNV, CopyKAT as cross-check | Per-cell malignancy call with confidence |
| 3–4 | Cell labels, axes 1 and 2, all four granularity rungs | Labels in **separate columns**, never overwriting each other |
| 4–5 | Build S matrices — including stromal, immune, endothelial columns | Versioned parquet, fixed gene index, one per rung |
| 5 | Cross-gene ambient check: retention vs. total abundance, all tiers | Plot + statistic. **This is gate criterion G1.** |

The week-2 pilot is the critical path for the whole project. If it slips past
week 3, the gate moves to week 7 — do not compress W2 to cover for it
(execution_plan.md §8.3).

## Week one, before anything else

Pull the **ICBI atlas metadata table only** — not the 4.27M-cell object. That
one table gives the real sample size: patients with paired tumour and normal,
epithelial fraction by study, platform mix, and whether the plate-based subset
is large enough to serve as an ambient-free validation set. Plate protocols have
essentially no soup, so an intrinsic signal surviving there is strong evidence
it is not contamination.

## What you emit, and where

```python
from src.common import s_matrix_path
s_matrix_path("lineage", "1.0.0")
# data/processed/reference/S_matrix_lineage_1.0.0.parquet
```

On the **fixed gene index** in `config/gene_index/`. W3 emits bulk on the same
index. Integration is a join, not a negotiation.

## The one thing you cannot get wrong

`build_signature()` in [signature.py](signature.py) asserts that target genes
never reach the reference matrix (CLAUDE.md invariant 2). Call
`assert_no_target_leakage()` in your **label** construction too — that function
cannot see your labels, and leaked labels break the project in exactly the same
way. If GUCA2A defines "mature" and GUCA2A is the test, a silenced mature cell
reads as an absent mature cell and the classifier cannot detect the phenomenon
it was built to detect.

`_select_markers()` is a `NotImplementedError` waiting for you. It is
deliberately not scaffolded — marker selection is a judgement call and W2's
bake-off will be interpreted against whatever you choose, so write down why in
the docstring.

## Gotchas

- GSE178341's supplementary structure is awkward — budget a day just for
  parsing. **Verify you have raw counts, not normalised values.**
- inferCNV on 371k cells is slow. Subsample per patient, or run per-patient in
  parallel.
- CellBender wants a GPU. Without one, budget overnight runs.
- Use backed AnnData for anything that does not need the full matrix in memory.
- Matched normals: if they are sparse, the compositional term loses its
  reference. Check completeness in week 1, not week 4.

## For W2 — consuming W1's labels

`src/reference/labels.py` emits eight `label_{axis}_{rung}` columns (2 axes × 4
rungs), one row per cell, categorical, none overwriting another.

```python
from src.reference.ingest import assign_compartments, read_gse178341, read_gse178341_clusters
from src.reference.labels import assign_labels, cell_type_vector, maturity_summary
from src.common.panel import tier_genes

adata    = read_gse178341(h5, patients=["C122", "C165", "C107", "C138", "C162"])
clusters = read_gse178341_clusters(cluster_csv)

labels = assign_labels(
    adata.X, adata.var["gene_symbol"],
    compartment=assign_compartments(clusters).reindex(adata.obs.index).to_numpy(),
    sample_id=adata.obs["sample_id"].to_numpy(),
    target_genes=tier_genes("A"),          # REQUIRED — see below
    index=adata.obs.index,
)

# -> the cell_type array generate_pseudobulk() expects
cell_type = cell_type_vector(labels, "stem_pole", "lineage")   # mature -> "mature_colonocyte"

# -> everything decompose_cohort() needs except gene / mean_normal / mean_tumour
summary = maturity_summary(labels, patient_id=..., tissue=..., study_id="GSE178341")
```

**`target_genes` is required and has no default.** Which genes count as targets
for a run is [open decision #1](../../docs/open_decisions.md), and a default would
bury it: passing the whole panel makes axis 2 unusable, while a permissive default
silently disables invariant 2. Consequence you will hit — **a run testing MUC2 or
TFF3 cannot use axis 2**, because those two genes are simultaneously tier-E targets
and axis-2 markers. Pass `axes=["stem_pole"]` for those runs.

**Which bin is mature** is `RUNG_SPECS[rung].mature`, and `mature_mask()` /
`cell_type_vector()` both read it from there. Do not hard-code a bin name — the
rung partitions are W1's proposal, not frozen, and they will change.

**The rungs are meant to disagree.** The coarsest calls all epithelium mature; the
finest calls ~5%. That spread *is* the granularity curve (§6.2) — a single point
estimate would present a modelling choice as a measurement.

**W1 and W4 currently define the rungs differently** — see
[open decision #13](../../docs/open_decisions.md). The interop functions above work
regardless, but the two cohorts are not comparable until that is settled.


## Lessons that cost real time

- **Per-sample quantile binning pins the mature fraction to the quantile**, making
  Δ(mature fraction) identically zero by construction. Cut points must come from the
  reference (normal) arm so the tumour arm is free to differ.
- **`np.where` evaluates both branches** — it cannot guard a division.
- **Index alignment.** `labels` is barcode-indexed; internal frames are not. Assigning a
  Series aligns on index, matches nothing, and yields a silent all-NaN column. Use
  `.to_numpy()`.
- **Never densify.** cells × genes is 8.3 GB at pilot scale and grows with the cohort.
  `build_signature_sparse` aggregates without materialising it.
- **Depth matching removes bias, not error.** It makes error unbiased with respect to
  depth; it does not make it smaller.
- **Dropout is stochastic**, so depth strata cannot separate technical from biological
  confounding. Only concordance with an independent annotation can.
- **Tests that check shapes and monotonicity pass while the measurement is meaningless.**
  Test that a known injected effect is recovered — that is the test that caught the
  quantile-pinning bug, and its absence is why the bug reached real data.
