# Pre-registration — Becker lesion-level Wnt detection gate

**Written:** 2026-09-08 · **Status:** proposed; no GSE201348 count matrix has
been read by this gate. **Scope:** the outcome-blind feasibility step for Item
4 of `docs/NEXT_AVENUES.md`, not B1 and not the completed within-cell D1 test.

## 1 · Question and fixed input

Can the five-gene Wnt signature (`AXIN2`, `NKD1`, `RNF43`, `NOTUM`, `TCF7`) be
measured in the Becker GSE201348 snRNA-seq polyp nuclei well enough to justify a
separate lesion-within-donor Wnt specification?

The only inputs are `GSE201348_RAW.tar` and
`GSE201348_series_matrix.txt.gz`, both checksum-pinned in `data/manifest.csv`.
The metadata-defined paired donors are fixed as A001, A002, A014, and A015.
Only their `Polyp` / `tumour` samples enter. CRC, unpaired polyps, normal tissue,
the differentiation panel, and every outcome are excluded from this gate.

One physical lesion is one `sample_id`; technical replicate sequencing rows are
concatenated within that lesion and their row count is retained. The donor is
the gate unit. Lesions are not promoted to independent patients.

## 2 · Measurement and pass rule

Detection is count >= 1 UMI (`DETECTION_MIN_UMI = 1`) in the raw 10x matrix.
Before a detection can be called zero, each signature symbol must be located
exactly once in the CellRanger feature-symbol column and that feature order must
match across every included sample.

For **every signature gene in every one of the four paired donors**, both must
hold after pooling that donor's polyp nuclei:

1. at least **30 detected nuclei**, matching the repository's minimum cell
   count for a Wnt correlation; and
2. detection in at least **1%** of that donor's polyp nuclei.

All 20 donor-by-gene rows must pass. These are measurement floors, not effect
sizes and not tuned after inspection.

## 3 · Fixed consequences

| gate outcome | consequence |
|---|---|
| all 20 rows pass | Draft a separate lesion-within-donor Wnt pre-specification. Do not fit it yet. |
| a feature is absent, duplicated, or feature order changes | **NO SUBSTRATE — identifier/assay failure.** Do not call the gene biologically absent. |
| any donor-by-gene row fails detection | **NO SUBSTRATE AT THIS snRNA RESOLUTION.** Report the table; do not fit a reduced signature or replace the gene. |

This gate cannot test Wnt biology, a differentiation association, or a
between-lesion effect. It decides only whether those questions can be specified
without asking a sparse nuclear assay to support a score it cannot measure.
