# W1 — Reference

> **W2 needs two things from you** — see [docs/handoff_w2_to_w1_w4.md](../../docs/handoff_w2_to_w1_w4.md).

**Owner:** Bode · **Env:** `env/w1_reference.yml` → `conda activate brp-w1`
**Branch prefix:** `w1/…` · **Blocked by:** nothing

> ## STATE — 2026-08-22
>
> Weeks 1–2 complete. **Malignancy calling works end to end.** 871 tests, ruff
> clean. Working branch `w1/infercnv-mtx-loader`, PR #28 — 11 commits ahead of
> main. Work on the branch; merging per fix cost an evening.
>
> ### Do this first
>
> Five `cnv_scores.csv` files are hours of compute sitting in **gitignored**
> scratch, on a filesystem that was 84% full. The step that reads them, runs
> `call_malignancy` + `validate_normal_epithelium` and writes a versioned
> parquet under `results/` **does not exist yet**. Write it before anything
> clears `data/interim`.
>
> ### Measured
>
> | | |
> |---|---|
> | Cohort | 62 patients · **36 matched** · **32 matched + unsorted** |
> | `stem_pole` kappa | 0.444 |
> | `opposite_lineage` kappa | **0.529** against its own criterion (−0.24 against the wrong one) |
> | inferCNV | 39,516 genes joined · specificity **0.99–1.00** out-of-sample |
> | Runtime | 7k cells 14 min · 22k cells ~45 min |
>
> ### Left
>
> Malignancy results writer · ambient correction (the only stubs) · full-scale
> labels and S matrices · `checks.py` for G1.
>
> **The week-2 ambient deliverable as written is impossible.** It asks for
> "SoupX and CellBender, both, compared"; CellBender needs empty droplets that
> exist in no public source (#8). Restate as SoupX vs DecontX first.
>
> ### The finding worth reading
>
> **Open decision #15.** inferCNV separates malignant cells by aneuploidy, and
> MMR-deficient tumours are characteristically near-diploid — so the method is
> expected to fail in exactly the stratum this project compares against MMRp.
> The MMRd tumour arm would keep more non-malignant (mature) epithelium,
> inflating its apparent mature fraction. **That bias runs along the
> pre-registered contrast, not across it**, and no other caller fixes it —
> CopyKAT infers copy number from expression too.
>
> Directionally consistent on the pilot, not established: n=2 MMRp vs n=3 MMRd
> with overlap. The prior is what makes it worth pre-specifying.
>
> ### Cluster
>
> **Disk, not CPU, is the constraint** — 55 GB for the whole project. Submit
> inferCNV with **`-tc 2`**. inferCNV is **CPU-only**, despite CLAUDE.md listing
> a GPU beside it.


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
