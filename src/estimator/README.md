# W4 — Estimator & replication

> **W2 needs two things from you** — see [docs/handoff_w2_to_w1_w4.md](../../docs/handoff_w2_to_w1_w4.md).

**Owner:** strongest methods implementer · **Env:** `env/w4_estimator.yml` → `conda activate brp-w4`
**Branch prefix:** `w4/…` · **Blocked by:** your own cohort only

You own the Kitagawa decomposition and the Lee cohorts (GSE132465 SMC,
GSE144735 KUL3). Developing on Lee means you are not queued behind W1 — and
independent replication falls out as a by-product, so the gate gets a second
opinion.

## What you deliver

| Wk | Task | Done when |
|----|------|-----------|
| 1–2 | Lee ingest, QC, ambient correction (same pipeline **shape** as W1 — coordinate, do not share code prematurely) | Cells labelled, axes 1 and 2 |
| 2–3 | Kitagawa standardisation: both weightings + interaction reported separately | **Unit-tested on synthetic data with analytically known answers** |
| 3–4 | Doubly-robust reweighted version | Agreement with the plain version quantified |
| 4 | Cross-check against cacoa and QuasiMed | Correlation reported |
| 4–5 | Patient-level bootstrap; hierarchical model with patient as grouping factor | CIs that reflect the real unit of inference |
| 5 | First decomposition results on Lee | Independent of W1's timeline — a second opinion at the gate |

## The estimator

```
compositional = Δ(mature fraction) × normal per-cell mean
intrinsic     = tumour mature fraction × Δ(per-cell mean)
interaction   = Δ(mature fraction) × Δ(per-cell mean)   [separately, never folded in]
```

Kitagawa (1955) demographic standardisation — **not** regression
Oaxaca–Blinder, which decomposes by regression coefficients and answers a
different question. The scalar identity is implemented and unit-tested in
[kitagawa.py](kitagawa.py); the per-patient version over AnnData, the
doubly-robust reweighting and the bootstrap are yours.

The split is **not unique**. Normal-weighted and tumour-weighted give different
answers and the difference lives in the interaction term. Report both plus
doubly-robust — three rows per (patient, gene, rung, axis), one per weighting.

## Three things that will silently ruin the result

1. **`None` is not `0.0`.** Take the verdict from W2:
   ```python
   from src.harness import classify_estimability
   est = classify_estimability(n_cells_mature)
   intrinsic = None if est == "not_estimable" else estimate
   ```
   Do not reimplement the thresholds — W2 recalibrates them at week 5 and you
   want to inherit that. `src.schema.write_results()` will reject the frame if a
   `not_estimable` row carries a number, but the honest value starts in your code.
2. **Bootstrap over patients, not cells** (CLAUDE.md invariant 5). Resampling
   cells inflates n by roughly the cells-per-patient count and produces
   intervals wrong by an order of magnitude, in the flattering direction.
3. **Estimate per study, then meta-analyse. Never pool** (invariant 4).
   Batch correction removes between-dataset variation, which is where the
   compositional signal lives. Integration is for label transfer only.

## Cross-checks

cacoa and QuasiMed, correlation reported. **Not CoCoA-diff** — it explicitly
assumes cell fractions are not a mediator, which assumes away the compositional
arm entirely.

Multiple testing: Benjamini–Hochberg **within tier**, reported separately for
each term.

## Week 13+

Becker/Chang multiome for the chromatin axis (axis 3) — a labelling axis not
made of transcripts at all. Check CRC overlap and accessibility early; if it is
thin, axis 3 becomes spatial-only. That check is cheap and it closes open
question 5.
