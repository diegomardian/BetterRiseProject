# D2: GUCA2A variation feasibility gate

## Purpose

D2 is a future marker-to-survival analysis, not the locked Stage 4 variance
arm. Before a survival pre-specification exists, this gate asks only whether the
committed TCGA expression matrix contains usable GUCA2A variation. It must run
without reading, joining, or loading clinical outcomes.

## Fixed inputs and unit

- Expression: `data/processed/bulk/tcga_log2cpm_1.0.0.parquet` on the cluster.
- Gene mapping: `config/gene_index/gene_index_1.0.0.map.tsv`; GUCA2A must map
  one-to-one to a matrix column. An ambiguous or absent mapping is a stop.
- Unit: one primary-tumour RNA sample per TCGA participant, selected using the
  committed `sample_manifest.tsv`. Technical aliquots do not increase `n`.
- Excluded data: all TCGA-CDR fields, survival outcomes, event indicators, and
  covariates. They must not be opened for this gate.

## Pre-committed pass rule

The gate passes only if all conditions hold:

1. At least 100 unique primary-tumour participants have a mapped GUCA2A value.
2. The participant-level GUCA2A distribution has IQR at least 0.5 log2-CPM.
3. At least 10% and at least 25 participants lie strictly below the median, and
   at least 10% and at least 25 lie strictly above it.

The output records `n`, minimum, quartiles, median, maximum, IQR, exact-zero
count, and the counts on either side of the median. It records only the gate
verdict and these distributional quantities; no association, Kaplan–Meier
curve, hazard ratio, or outcome-derived threshold is calculated.

## Interpretation

`PASS` licenses drafting a separate D2 survival pre-specification; it does not
license fitting a model. `STOP` means there is insufficient variation for the
proposed continuous-marker question in this matrix and no D2 survival model is
written from it. The threshold is fixed before the cluster read and is not
retuned after seeing the distribution.

## Execution boundary

The local checkout has no committed TCGA processed matrix, so execution belongs
to W3 on the cluster. Run, from a clean checkout with `BRP_DATA_DIR` exported:

```bash
python -m src.bulk.d2_feasibility_gate
```

It writes `d2_guca2a_feasibility.parquet` through the standard results writer
and cites this document as the gate contract.
