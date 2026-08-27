# RETRACTED — the four `0.1.0-pilot` S matrices in this directory

**Retracted:** 2026-08-26 · **By:** W1 (Bode) · **Found by:** W3, [issue #35](https://github.com/diegomardian/BetterRiseProject/issues/35)

The four `S_matrix_*_0.1.0-pilot.parquet` files here **violate CLAUDE.md
invariant 2** and must not be used. Everything else in this directory stands.

| file | tier-A target genes present |
|---|---|
| `S_matrix_best4_0.1.0-pilot.parquet` | **GUCA2A, GUCA2B, CA7, OTOP2** — all of tier A |
| `S_matrix_crypt_position_0.1.0-pilot.parquet` | **GUCA2A** |
| `S_matrix_epithelial_0.1.0-pilot.parquet` | **GUCA2A** |
| `S_matrix_lineage_0.1.0-pilot.parquet` | **GUCA2A** |

## Why this is the serious kind of wrong

Invariant 2 exists because **a silenced mature cell must not be readable as an
absent mature cell** — Executive-Brief error #1. With GUCA2A in the reference
matrix, a tumour cell that has switched GUCA2A off looks to the deconvolution
like a cell that was never there. That is the compositional/intrinsic confusion
the entire project is built to separate, reintroduced inside the instrument.

These files are the W1 → W2 handoff, and **G2 asks whether the control tiers
separate on exactly these genes.**

## How it happened

`assert_no_target_leakage()` intersected target *symbols* (`GUCA2A`) with the
genes it was handed. When those genes are unversioned Ensembl ids the
intersection is empty **whatever the data**, so all four guard call sites passed
unconditionally while reading as enforced. A guard that cannot fire reports
success.

Two things kept it invisible:

1. The S-matrix build caught bare `Exception` and printed `{rung} skipped`, so a
   real leak and an unfirable guard looked identical to a skipped rung.
2. `tests/test_reference_gene_index.py` asserted the guard **passed**, with the
   comment *"Ensembl-keyed, so symbol-named panel genes cannot collide with it."*
   The defect was written down as intended behaviour and pinned there by a test.

Fixed in `59ae14f`: the guard now detects the identifier space of both sides and
**refuses** rather than passing vacuously, leakage errors are never swallowed,
and the test asserts the refusal.

## Replacements, and what they are NOT valid for

`S_matrix_{rung}_1.0.0.parquet`, built 2026-08-26 at `59ae14f`, on the
39,236-gene shared index. Verified: **0 of 4 tier-A genes** in any of them. They
live under `data/processed/reference/`, which is gitignored, so they travel by
copy — not by git.

> **They exclude tier A only.** `run_full_reference.py` passes
> `targets = tier_genes("A")`, and open decision #1's narrow reading filters the
> target set *for the run in question*. So `SFRP1`, `SFRP2` (tier B), `MS4A12`
> (tier D) and `CDX2` (tier C) **are present** in these matrices.
>
> **A run testing tier B, C or D against them violates invariant 2**, and the
> guard will not catch it, because those genes were never passed as targets.
> G2 spans tier A, B and D — so G2 cannot use one of these matrices for all
> three arms as things stand. Raised for the weekly.
