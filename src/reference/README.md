# W1 — Reference

**Owner:** strongest scRNA-seq person · **Env:** `env/w1_reference.yml` → `conda activate brp-w1`
**Branch prefix:** `w1/…` · **Blocked by:** nothing — you start day one

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
