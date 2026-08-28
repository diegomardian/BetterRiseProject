"""Open decision #2, executed: the shared gene index at 1.0.0, as the intersection.

WHY THIS EXISTS
---------------
Both arms measured the overlap independently and reached the same answer, in
writing:

- W3, `cc06981`: 39,236 genes on both, all 23 panel genes present.
- W1, [issue #7](https://github.com/diegomardian/BetterRiseProject/issues/7)
  and its comment: *"We independently reached the same conclusion: 1.0.0 is the
  intersection, each arm keeping its full matrix for its own work. So #2 and #3
  are arithmetic now rather than opinion."*

What never happened is the commit. `config/gene_index/` holds only `0.9.0`,
W3's provisional index built from the GDC gene model when the documented
fallback fired. Meanwhile `src/reference/jobs/run_full_reference.py` hard-codes
``GENE_INDEX_VERSION = "1.0.0"`` and, finding no such file, prints

    !! no gene index 1.0.0: ...
       S matrices skipped — they must sit on the shared index or integration is
       a negotiation, not a join.

and skips the step. So the agreed decision not being executed is what is
currently costing the project W1's S matrices, which are the W1 → W2 handoff the
gate reads.

WHY THE INTERSECTION AND NOT EITHER ARM'S NATIVE SET
----------------------------------------------------
Deconvolution needs a gene present in **both** matrices. A row of the S matrix
for a gene the single-cell data never measured is undefined; a bulk gene with no
single-cell counterpart contributes nothing. So the operative index is the
intersection whatever file we bless — the only question is whether that is
written down or rediscovered later by whoever hits it.

Adopting either arm wholesale discards what the other measured: W3 conforming to
W1 loses 21,380 bulk-only genes, W1 conforming to W3 loses 3,877. Each arm keeps
its own full matrix for its own work; the shared index is what integration joins
on, and 39,236 is far more than the 500–2,000 markers deconvolution needs
(execution_plan.md §2.1 error #4).

WHAT IS CARRIED OVER, AND WHAT IS NOT
--------------------------------------
The 1.0.0 map is W3's 0.9.0 map **filtered**, not rebuilt. Symbols, gene types,
`symbol_ambiguous` and `on_panel` are the same columns decided in W3.1 and
already reviewed; regenerating them here would create a second source of truth
for identifier decisions that decision #3 settled once.

Ordering is W3's 0.9.0 ordering, preserved. `build_signature` reindexes onto the
shared index and `test_rows_follow_the_shared_gene_index_order` pins that the
row order is the index's — so the order is part of the contract, not incidental.

PANEL GENES ARE IN THIS INDEX ON PURPOSE
-----------------------------------------
All 23 of them, and the emitter refuses to write an index that has lost one.
They are W3's **outcome variables** for the premise check and the Stage 4
variance question; an index without them makes W3's own deliverable
unrepresentable. That is invariant 2's boundary as decision #12 draws it —
invariant 2 governs what may enter the *reference matrix*, not what may sit on
the shared index.

**Emitting 1.0.0 does NOT unblock the S matrices on its own, and the reason is
worse than decision #12.** Measured, not assumed:

1. `run_full_reference.py:208` passes ``adata.var["gene_symbol"]`` as
   ``gene_names`` while ``gene_index`` is Ensembl-keyed. Symbols never intersect
   ENSG ids, so `build_signature_sparse` raises *"no gene in the matrix appears
   on the shared gene index"* — the exact error decision #3 anticipated. The
   pilot job got this right (`run_pilot.py` used ``adata.var["ensembl_id"]``);
   the full-scale job regressed to symbols.
2. **Every invariant-2 guard is inert in Ensembl space.**
   ``assert_no_target_leakage`` intersects target *symbols* with a *gene id*
   list, and `set(["GUCA2A"]) & set(["ENSG00000197273"])` is empty. All four
   call sites — the index, the reference pool, the marker set, and the emitted
   S matrix — pass unconditionally.

That is not hypothetical. **GUCA2A (`ENSG00000197273`) is in all four committed
pilot S matrices**, and `S_matrix_best4_0.1.0-pilot.parquet` carries all four
tier-A targets: GUCA2A, GUCA2B, CA7, OTOP2. That is CLAUDE.md invariant 2 and
Executive-Brief error #1 — *a silenced mature cell readable as an absent mature
cell* — in a committed artifact.

Neither is fixed here: `src/reference/` is W1's (CONTRIBUTING §2). Both are
raised as an issue. Note the fix is bigger than
`docs/decision_12_signature_filter.patch`, which filters `gene_index` by
`targets` and is therefore *also* inert in Ensembl space — the guard has to
resolve panel symbols into the index's key space before comparing, which is what
`src.bulk.gene_index.resolve_symbols` and the 1.0.0 map exist to do.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from src.bulk.gene_index import (
    MAP_COLUMNS,
    GeneIndexError,
    load_gene_index_map,
    write_gene_index,
)
from src.common.panel import panel_genes

#: W3's provisional index, built from the GDC gene model (GENCODE v36, GRCh38).
BULK_VERSION = "0.9.0"

#: The shared index both arms join on. Decision #2.
SHARED_VERSION = "1.0.0"


def intersect(bulk_map: pd.DataFrame, reference_ids: list[str]) -> pd.DataFrame:
    """W3's map, restricted to identifiers the single-cell deposit also carries.

    Joined on the **unversioned** Ensembl id, which is decision #3's whole
    point: the deposit is GRCh37_liftover_v28 and TCGA is GRCh38, and the
    unversioned id is the only form both builds share. Joining on the versioned
    string drops ~9% of genes silently.
    """
    if "ensembl_id" not in bulk_map.columns:
        raise GeneIndexError("the bulk map has no `ensembl_id` column")
    shared = set(reference_ids)
    if not shared:
        raise GeneIndexError("the reference identifier list is empty")

    keep = bulk_map.loc[bulk_map["ensembl_id"].isin(shared)].copy()
    if keep.empty:
        raise GeneIndexError(
            "no identifier is on both sides. Check the identifier FORM — the "
            "index is unversioned Ensembl ids (decision #3), and a versioned "
            "or symbol-keyed list will intersect to nothing rather than fail."
        )
    return keep.loc[:, list(MAP_COLUMNS)].reset_index(drop=True)


def assert_panel_survives(index: pd.DataFrame) -> None:
    """Refuse to emit an index that has lost a frozen panel gene.

    A panel gene absent from the join is untestable in the integrated analysis,
    and the panel is frozen (invariant 3) precisely so that is not renegotiated
    after the fact. If this ever fires it is a reason to revisit the index, not
    a reason to drop the gene.
    """
    panel = set(panel_genes())
    present = set(index.loc[index["gene_symbol"].isin(panel), "gene_symbol"])
    missing = sorted(panel - present)
    if missing:
        raise GeneIndexError(
            f"{len(missing)} frozen panel gene(s) absent from the intersection: "
            f"{missing}. Do not emit this index. A panel gene lost in the join is "
            f"untestable in the integrated analysis (invariant 3)."
        )


def reconciliation(
    bulk_map: pd.DataFrame, reference_ids: list[str], shared: pd.DataFrame
) -> pd.DataFrame:
    """One row per side, so the losses are counted rather than described."""
    bulk = set(bulk_map["ensembl_id"])
    reference = set(reference_ids)
    both = bulk & reference
    rows = [
        {"side": "bulk (W3, GENCODE v36)", "n_genes": len(bulk),
         "on_both": len(both), "lost": len(bulk - reference)},
        {"side": "reference (W1, GSE178341)", "n_genes": len(reference),
         "on_both": len(both), "lost": len(reference - bulk)},
        {"side": "shared index 1.0.0", "n_genes": len(shared),
         "on_both": len(both), "lost": 0},
    ]
    table = pd.DataFrame(rows)
    table["lost_fraction"] = (table["lost"] / table["n_genes"]).round(4)
    return table


def build_shared_index(
    reference_ids: list[str],
    *,
    bulk_version: str = BULK_VERSION,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """The 1.0.0 map and its reconciliation. Writes nothing."""
    bulk_map = load_gene_index_map(bulk_version)
    shared = intersect(bulk_map, reference_ids)
    assert_panel_survives(shared)
    return shared, reconciliation(bulk_map, reference_ids, shared)


def emit(
    reference_ids: list[str],
    *,
    version: str = SHARED_VERSION,
    bulk_version: str = BULK_VERSION,
    config_dir: Path | None = None,
) -> tuple[Path, Path, pd.DataFrame]:
    """Build and write the shared index. Refuses to overwrite an existing version."""
    shared, report = build_shared_index(reference_ids, bulk_version=bulk_version)
    index_path, map_path = write_gene_index(
        shared, version=version, config_dir=config_dir
    )
    return index_path, map_path, report
