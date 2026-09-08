# Pre-registration — Zheng three-stage descriptive gradient

**Written:** 2026-09-08 · **Status:** RAN 2026-09-08; all three fixed patients
and stages remained measurable after QC. **Scope:** Item 1c of
`docs/NEXT_AVENUES.md`: a three-patient normal → polyp → carcinoma timing
companion. It is not a replication, premise test, decomposition, or per-cell
silencing analysis.

## 1 · Fixed cohort and population

The input is the checksum-pinned ICBI atlas counts layer and its row-aligned
`icbi_obs.parquet` cache. The study is fixed as
`Zheng_2022_Signal_Transduct_Target_Ther`. The three stages are fixed from
`sample_type`: `adjacent normal` or `healthy normal` → normal, `polyp` → polyp,
and `primary tumor` → carcinoma.

Only the atlas's coarse **epithelial** annotation enters. This is deliberately
not a transcript-derived mature-cell label: the table describes the epithelial
population present at each stage, including any compositional change. The
population label uses no panel genes; the claim measurement is raw transcript
counts for the frozen six-gene panel. The provenance guard runs before any
matrix read.

A patient enters only when all three stages contain at least 100 naive,
coarse-epithelial cells before QC. The verified object supplies exactly three
such patients. Any other count stops the job; it does not silently add patients
or discard one. QC follows the repository's per-`sample_id` rules over each
patient's full cell block, before epithelial cells are summarized.

## 2 · Descriptive output

For every retained patient, stage, and frozen-panel gene, emit the number of
epithelial cells before and after QC, number detected at >=1 UMI, detection
fraction, CP10K mean, median total UMI, median genes detected, and number of
samples. A separate availability table retains each stage's before/after-QC
counts, including stages with zero post-QC epithelial cells.

The job reports rows and patient trajectories only. It computes **no** pooled
mean, stage contrast, p-value, bootstrap, confidence interval, premise verdict,
or model. It must not be combined with Chen or any other study: n=3 is below
the project’s inferential floor and a percentile interval at n=3 is known to be
anti-conservative.

## 3 · Fixed reading

The output can show whether the six measured transcripts move before, at, or
after the polyp stage within each of the three patients. It is timing context
only. It cannot establish a common trajectory, a stage effect, cell-intrinsic
loss, or a causal order. A zero after QC means unmeasurable at that stage, never
zero expression.

## 4 · Result

The run used commit `f6b22f6` and wrote clean artifacts in
`results/2026-09-08_f6b22f6/`. It produced 54 trajectory rows (three patients
× three stages × six genes) and nine availability rows; every stage retained
measurable post-QC epithelium.

GUCA2A detection fell from normal to polyp in each individual trajectory:
P1 0.579 → 0.040, P2 0.703 → 0.027, and P3 0.642 → 0.040. Its CP10K mean also
fell in all three. MS4A12 followed that normal-to-polyp direction in all three.
CDX2 did not: it rose in P1 and P2 polyp and fell in P3. Carcinoma GUCA2A
remained below normal in every trajectory but partially rebounded relative to
polyp in P1 and P2.

These rows remain descriptive timing context only. P1 and P2 normal samples
were substantially shallower than their polyp samples, and P3 retains only 150
post-QC polyp epithelial cells. No pooled trend, interval, stage contrast,
identity conclusion, or intrinsic-loss conclusion is licensed.
