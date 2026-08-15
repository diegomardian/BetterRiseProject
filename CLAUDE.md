# CLAUDE.md

Project context and invariants. Read [README.md](README.md) for what the project is and [execution_plan.md](execution_plan.md) for who does what when.

**One-line summary:** decompose differentiation-marker loss in colorectal cancer into *compositional* (mature cells gone), *cell-intrinsic* (mature cells present but silenced), and *not estimable* (too few mature cells to ask) — per patient, at several annotation resolutions.

## Invariants — do not violate without an explicit written decision

1. **`None` is not `0.0`.** An unestimable intrinsic term is `None` with `estimability="not_estimable"`. Writing zero there is the single most likely route to a wrong conclusion in this project. Enforced by assertion in the results writer, not by review.
2. **Target genes never appear in labels or the reference matrix.** `build_signature()` asserts this. A silenced mature cell must not be readable as an absent mature cell.
3. **The panel and the labelling axes are frozen.** Changing either — or `src/schema.py` — needs a PR, two approvals and a written reason.
4. **Estimate per study, then meta-analyse. Never pool.** Batch correction removes between-dataset variation, which is where the compositional signal lives. Integration is for label transfer only.
5. **Bootstrap over patients, not cells.** The patient is the unit of inference.
6. **No cell-type-specific expression imputation from bulk.** Bulk gives fractions (r ≈ 0.92); intrinsic estimates from bulk are attenuated ×0.6–0.8 in the direction of the prior hypothesis. Fractions only.
7. **The interaction term is reported separately, never folded into either arm.**
8. **CTNNB1 / TCF7L2 transcript level is not Wnt activity.** Use a target signature (AXIN2, NKD1, RNF43, NOTUM, TCF7); drop ASCL2/LGR5 when the stem axis is in play.
9. **Survival endpoints are DSS and PFI from TCGA-CDR.** OS is secondary — COAD OS is contaminated by non-cancer death.
10. **Every result carries the git sha and a fixed random seed.** Results are versioned parquet under `results/`.

## Repo layout

```
data/       raw/ interim/ processed/   gitignored; manifest with checksums
src/        reference/  W1     harness/   W2
            bulk/       W3     estimator/ W4
            schema.py   shared — PR + 2 approvals to change
results/    parquet, versioned by date + git sha
env/        one conda env per workstream, pinned
```

W1 emits `S_matrix_{rung}_{version}.parquet` on a fixed gene index; W3 emits bulk on the same index. Integration is a join, not a negotiation.

## Two timelines, on purpose

README.md describes the sequential staging (Stage 1 weeks 1–8, gate at the end of Stage 1). execution_plan.md describes the four-person parallel version, where the same gate lands at week 5 because W1 hands W2 a five-patient pilot at week 2 instead of the full object. If that pilot slips past week 3, the gate moves to week 7 — the harness is what makes the gate meaningful, so it is not the thing to compress.

## Compute

32–64 GB RAM and GPU access for the single-cell arm (CellBender, inferCNV on 371k cells). The bulk/TCGA arm is laptop-fine. "Standard consumer laptop" was wrong and is corrected in both documents.
