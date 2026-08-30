# W1 — Reference

> **W2 needs two things from you** — see [docs/handoff_w2_to_w1_w4.md](../../docs/handoff_w2_to_w1_w4.md).

**Owner:** Bode · **Env:** `env/w1_reference.yml` → `conda activate brp-w1`
**Branch prefix:** `w1/…` · **Blocked by:** nothing

> ## STATE — 2026-08-29
>
> **W1 is complete. Nothing is blocked on W1.** Everything outstanding is a
> response to something W1 raised, or another workstream's decision.
>
> ### The result, and what it may and may not be used for
>
> The decomposition ran on GSE178341 — all tiers A–D, both axes, four rungs,
> 32 patients, patient-level bootstrap. **On the depth-matched read exactly one
> rung clears every check: `lineage`.**
>
> | GUCA2A at `lineage` | compositional | intrinsic |
> |---|---|---|
> | unmatched | −2.412 | −13.821 |
> | **depth-matched** | −2.827 | **−14.741** |
>
> Tier A is the panel's designated *compositional* control, and it comes out
> predominantly **intrinsic** — the pre-registered G2 pattern failing. **Depth
> matching did not collapse it**, which is the test it had to survive: had the
> intrinsic term been the arm depth imbalance, matching would have removed it.
> W2's independent Lee result (#49) points the same way from data sharing no
> cells, no platform and no labelling step.
>
> **What must travel with that number, every time:** one study, one rung, 45.4%
> of cells after matching, and a cohort-level interval broadcast onto each
> patient row rather than a per-patient one (#10). `best4` — the cleanest rung
> on depth — cannot corroborate it, being unestimable on every patient.
>
> ### Why the other three rungs are out, each for a different reason
>
> | rung | why not |
> |---|---|
> | `epithelial` | Degenerate by construction — every scored cell is mature, so the compositional term is **exactly 0.000** and cannot move. Matching does not touch this. |
> | `crypt_position` | Depth fixed by matching, but still a two-bin split on ~90% of patients (#42), so not an independent point on the curve. |
> | `best4` | 1,134 rows unestimable, **0** `ok`, median 1 mature cell against a cutpoint of 50. This is G4 (#48) on W1's own cohort. |
>
> ### Five things a new session must NOT undo
>
> 1. **`checks.py` returns `not_estimable` and must keep doing so.** G1's
>    specification is unresolved (#46, #48). It is not stale code.
> 2. **Threshold 2 (MS4A12 ≥ 0.50) is frozen.** Committed 2026-08-25, since seen
>    to fail 31/31. Any change now is indistinguishable from tuning.
> 3. **Tier A's G1 inputs are deliberately unmeasured.** Thresholds 1 and 3 are
>    the last unseen G1 inputs; measuring them is the team's call.
> 4. **No third G1 repair was proposed.** Two of W1's own failed; a third chosen
>    after watching them fail is a rule iterated until it passes. A repair
>    *exists* — match the comparison set on cell-type restriction, measurable in
>    the normal arm alone — which is why "G1 is unrepairable" is not established
>    (see the #48 review).
> 5. **The measurement that would settle #42 is undone on purpose** — it is one
>    step from an axis-change proposal, and the axes are frozen.
>
> ### Open, none of it W1's to move
>
> | | |
> |---|---|
> | **#48** | Withdraw G1 or re-specify it. W1's review showed the withdrawal premise does not hold. No cluster dependency — G1 has never run. |
> | **#54** | Stage 4's primary prediction is satisfiable by a floor effect. Needs a low-abundance control **before the pre-spec locks** — it is `status: proposed`, and after the lock a control chosen with the result in view is not a control. |
> | **#55 #56 #57** | Need reviewers. #55 is approved by W1 and waiting on W4. |
> | **#42 #43 #8 #38** | With W3 and W4. |
>
> ### The single cheapest thing that would improve the result
>
> **Run `match_arm_depth` on Lee's `lineage`.** #49 used the unmatched read for
> the same reason W1's first decomposition did — both predated #24.1. On
> GSE178341 matching cost 55% of cells but only ~7% of estimable rows and lost no
> patient. Lee's arms are 2.36× apart, so it is the same operation.
>
> That is **one run**, and it is the difference between k = 1 and k = 2 — i.e.
> between "there is nothing to meta-analyse" (decision #25) and §6.1's stated
> output existing at all. It is W2's cohort, so W1 can ask but not do it.
>
> ### A pattern worth knowing before touching anything
>
> **Symbols and Ensembl ids meet in more places than anyone has enumerated, and
> every collision fails silently.** Seven instances this week across three
> workstreams, including two inside the guard written to end it and one in a job
> whose only purpose was auditing someone else's code. The general fix adopted:
> **a check that cannot fire must refuse, not report success.** If you add a
> comparison between two gene-name sets, assert the spaces match.
>
> The same shape recurs beyond identifiers — `nan` compared against a threshold
> returns `False` and reads as a pass; a filter that removes nothing reads as a
> filter. Both happened in merged code this week.
>
> ### Measured at full scale
>
> | | |
> |---|---|
> | Cohort | 62 patients · 36 matched · **23 with both arms ambient-interpretable** |
> | `stem_pole` | kappa 0.444 — mature means *no stem marker detected at 3,281 UMIs* |
> | `opposite_lineage` | kappa **0.529 against its own criterion** — a goblet axis, **not** maturity |
> | `best4` | kappa 0.045 — **do not quote** |
> | inferCNV | 39,516 genes · specificity 0.99–1.00 · **30/62 separable** |
> | Ambient | median **2.2%**, 9 of 84 above 10% |
>
> ### The finding that shapes everything downstream
>
> **Open decision #15, confirmed.** MMR-proficient tumours separate **15/15**;
> MMRd **15/20** — and within *callable* patients MMRd yields **3.4× fewer**
> cells called malignant. Both follow from MMRd tumours being near-diploid, and
> no other caller escapes it. The bias runs **along** the pre-registered MMR
> contrast rather than across it.
>
> `docs/prereg_amendment_1_mmr_tumour_arm.md` responds: report the contrast under
> **both** tumour-arm definitions, and treat disagreement as *not identifiable*
> rather than choosing. Written before any expression was examined — **that
> timing is the whole of its credibility.**
>
> ### Cluster
>
> Disk is the constraint, not CPU — 55 GB for the project. inferCNV `-tc 1` at
> ~25 GB free. inferCNV is **CPU-only** despite CLAUDE.md listing a GPU.



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
