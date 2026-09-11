# Replogle sensitivity — terms, pinning and provenance gate

**Status, 2026-09-10:** pre-lock gates **1, 2 and 4 discharged**; **gate 3 is
open** and its blocker is specific. `docs/prereg_replogle_sensitivity.md` §7.
No expression value has been parsed and the protocol is **not yet locked**.

## 1 · Terms and licence — DISCHARGED

The Figshare+ deposit `10.25452/figshare.plus.20029387.v1` (Replogle, Weissman,
published 2022-06-09) is **CC BY 4.0**, stated by the repository API as
`https://creativecommons.org/licenses/by/4.0/`.

CC BY permits use, modification and derived work, including commercially, with
**attribution**. That is materially more permissive than MHIST's research-use
agreement, which prohibited derivative works. What is owed is citation of
Replogle et al., *Cell* 2022 (`10.1016/j.cell.2022.05.013`) and of the dataset
DOI in any output.

## 2 · Snapshot pinned — DISCHARGED for what is held

The deposit is **12 files**. The repository supplies an md5 per file, so files
are pinned **without downloading them**, which is what makes this gate
metadata-only.

| file | size | publisher md5 |
|---|---|---|
| `K562_essential_raw_bulk_01.h5ad` | 79.8 MB | `8321d5d3…9735` |
| `K562_essential_raw_singlecell_01.h5ad` | **10.66 GB** | `4f1122ce…266d3` |
| `rpe1_raw_singlecell_01.h5ad` | 8.70 GB | `6a2a9d0d…b268c` |
| `K562_gwps_raw_singlecell_01.h5ad` | 65.83 GB | `887e3e6a…c9546` |

**Held and verified:** `K562_essential_raw_bulk_01.h5ad` only. Downloaded
2026-09-10, its md5 reproduced the publisher's value exactly, and it is in
`data/manifest.csv` by sha256. The genome-scale files are **not** needed: §4
fixes the essential-scale dataset.

**A correction to §7's wording.** Gate 2 said "sha256 in `data/manifest.csv`",
but sha256 requires the bytes, which is not metadata-only. The pre-lock pin is
therefore the **publisher's md5**, which is stronger provenance anyway — it
comes from the depositor rather than from our copy — and sha256 is recorded when
a file lands and must agree with the pinned md5.

## 3 · Guide-assignment field at source — OPEN, and the blocker is exact

The condition requires confirming the per-cell guide field from the deposit's
own schema rather than inferring it from a column name. **It is not yet met.**

The deposit's description documents the file naming and AnnData layout but names
no `.obs` field. The pseudobulk file's schema was read directly — columns only,
no values — and it is **aggregated by target, not by cell**: indexed on
`gene_transcript`, with `control_expr`, `fold_expr`, `pct_expr`,
`num_cells_filtered`, `num_cells_unfiltered`, `core_control` and the authors'
own `anderson_darling_counts`, `mann_whitney_counts`, `energy_test_p_value`.
**There is no guide column, because a pseudobulk row is not a cell.**

**So gate 3 needs the `.obs` of `K562_essential_raw_singlecell_01.h5ad`, and
that requires the 10.66 GB download.** Downloading bytes and reading column
names parses no count and cannot shape the design, so it is inside the pre-lock
remit — it is a resource decision, not a methodological one.

### What the pseudobulk file changed about the plan

`control_expr` and `fold_expr` are the authors' **own measured baseline
abundance and on-target knockdown per target**. The prereg treats measured
knockdown as an axis (§4a) and this supplies it at 80 MB rather than 10.66 GB.
The single-cell file is still required for the curve itself — thinning, cell
counts and per-cell detection cannot come from pseudobulk — and for gate 3.

**Not read, deliberately.** The distribution of `control_expr` would say which
abundance bins are populated. It was **not** inspected, because §3 anchors the
bins to this project's own regimes and §4a already handles a sparse cell with a
not-estimable rule. Reading it pre-lock would be a value read that could shape
the grid.

## 4 · Label family declared — DISCHARGED

`src/reference/replogle_labels.py`. One family: the perturbed arm is
**`genotype` — per-cell CRISPRi guide identity from guide capture** — against a
**`transcript`** claim. Guard passes with an empty overlap.

`tests/test_replogle_labels.py` pins the failure it exists to prevent: a label
defined by the targeted transcript against a claim on the same transcript is
refused by `check_no_circular_claim`. That is the configuration §2 forbids, and
the one that would manufacture near-perfect power.

## 5 · What remains before lock

Gate 3 only. On it passing, the protocol locks immediately, and the §4
positive-control gate runs as the first post-lock artifact.
