"""The shared gene index at 1.0.0 — open decision #2, executed.

The committed file is checked here as well as the code that builds it, because
the file is the contract: `build_signature` reindexes onto it, W3 emits bulk on
it, and integration is a join only if both arms sit on the same one.
"""

from __future__ import annotations

import pandas as pd
import pytest

from src.bulk.gene_index import (
    MAP_COLUMNS,
    GeneIndexError,
    load_gene_index,
    load_gene_index_map,
)
from src.bulk.shared_index import (
    BULK_VERSION,
    SHARED_VERSION,
    assert_panel_survives,
    build_shared_index,
    emit,
    intersect,
    reconciliation,
)
from src.common.panel import panel_genes


def _map(ids, symbols=None):
    symbols = symbols or [f"SYM{i}" for i in range(len(ids))]
    return pd.DataFrame(
        {
            "ensembl_id": ids,
            "ensembl_version": ["1"] * len(ids),
            "gene_symbol": symbols,
            "gene_type": ["protein_coding"] * len(ids),
            "symbol_ambiguous": [False] * len(ids),
            "on_panel": [False] * len(ids),
        }
    )


# ---------------------------------------------------------------------------
# The build
# ---------------------------------------------------------------------------


def test_intersect_keeps_only_identifiers_on_both_sides():
    bulk = _map(["ENSG1", "ENSG2", "ENSG3"])
    assert list(intersect(bulk, ["ENSG2", "ENSG3", "ENSG9"])["ensembl_id"]) == [
        "ENSG2",
        "ENSG3",
    ]


def test_intersect_preserves_the_bulk_index_order():
    """`test_rows_follow_the_shared_gene_index_order` pins that an S matrix's rows
    follow the index, so the index's own order is part of the contract."""
    bulk = _map(["ENSG3", "ENSG1", "ENSG2"])
    got = list(intersect(bulk, ["ENSG1", "ENSG2", "ENSG3"])["ensembl_id"])
    assert got == ["ENSG3", "ENSG1", "ENSG2"]


def test_intersect_keeps_the_committed_map_columns():
    bulk = _map(["ENSG1", "ENSG2"])
    assert tuple(intersect(bulk, ["ENSG1"]).columns) == MAP_COLUMNS


def test_an_empty_overlap_fails_loudly_rather_than_emitting_nothing():
    """The likely cause is a versioned or symbol-keyed list, which intersects to
    nothing rather than erroring — decision #3's silent-loss failure mode."""
    with pytest.raises(GeneIndexError, match="identifier FORM"):
        intersect(_map(["ENSG1"]), ["ENSG00000141510.16"])


def test_an_empty_reference_list_is_refused():
    with pytest.raises(GeneIndexError, match="empty"):
        intersect(_map(["ENSG1"]), [])


# ---------------------------------------------------------------------------
# The panel guard — the one that must not be waved through
# ---------------------------------------------------------------------------


def test_a_lost_panel_gene_refuses_the_emit():
    """A frozen panel gene absent from the join is untestable in the integrated
    analysis. That is a reason to revisit the index, not to drop the gene."""
    frame = _map(["ENSG1"], symbols=["NOT_A_PANEL_GENE"])
    with pytest.raises(GeneIndexError, match="panel gene"):
        assert_panel_survives(frame)


def test_the_panel_guard_passes_when_every_gene_is_present():
    panel = sorted(panel_genes())
    frame = _map([f"ENSG{i}" for i in range(len(panel))], symbols=panel)
    assert_panel_survives(frame)  # does not raise


# ---------------------------------------------------------------------------
# Reconciliation and writing
# ---------------------------------------------------------------------------


def test_reconciliation_counts_both_sides_losses():
    bulk = _map(["ENSG1", "ENSG2", "ENSG3"])
    shared = intersect(bulk, ["ENSG2", "ENSG3", "ENSG9"])
    table = reconciliation(bulk, ["ENSG2", "ENSG3", "ENSG9"], shared)
    by_side = dict(zip(table["side"], table["lost"], strict=True))
    assert by_side["bulk (W3, GENCODE v36)"] == 1  # ENSG1
    assert by_side["reference (W1, GSE178341)"] == 1  # ENSG9
    assert by_side["shared index 1.0.0"] == 0


def test_emit_refuses_to_overwrite_an_existing_version(tmp_path, monkeypatch):
    """An in-place edit makes every earlier result unreproducible without saying
    so — config/gene_index/README.md."""
    panel = sorted(panel_genes())
    ids = [f"ENSG{i}" for i in range(len(panel))]
    monkeypatch.setattr(
        "src.bulk.shared_index.load_gene_index_map", lambda _v: _map(ids, panel)
    )
    emit(ids, version="9.9.9", config_dir=tmp_path)
    with pytest.raises(GeneIndexError, match="already exists"):
        emit(ids, version="9.9.9", config_dir=tmp_path)


def test_build_returns_the_map_and_the_reconciliation_without_writing(monkeypatch):
    panel = sorted(panel_genes())
    ids = [f"ENSG{i}" for i in range(len(panel))]
    monkeypatch.setattr(
        "src.bulk.shared_index.load_gene_index_map", lambda _v: _map(ids, panel)
    )
    shared, report = build_shared_index(ids)
    assert len(shared) == len(ids)
    assert set(report["side"]) == {
        "bulk (W3, GENCODE v36)",
        "reference (W1, GSE178341)",
        "shared index 1.0.0",
    }


# ---------------------------------------------------------------------------
# The committed artifact
# ---------------------------------------------------------------------------


def test_the_committed_shared_index_is_a_subset_of_the_bulk_index():
    """1.0.0 is 0.9.0 filtered, not rebuilt. If it ever stops being a subset,
    a second source of truth for identifier decisions has appeared."""
    shared = set(load_gene_index(SHARED_VERSION))
    bulk = set(load_gene_index(BULK_VERSION))
    assert shared <= bulk
    assert shared  # and it is not empty


def test_the_committed_shared_index_carries_every_panel_gene():
    """They are W3's outcome variables. An index without them makes the premise
    check and the Stage 4 variance question unrepresentable (decision #12)."""
    assert_panel_survives(load_gene_index_map(SHARED_VERSION))


def test_the_committed_index_and_its_map_agree():
    index = load_gene_index(SHARED_VERSION)
    mapping = load_gene_index_map(SHARED_VERSION)
    assert list(mapping["ensembl_id"]) == index
    assert mapping["ensembl_id"].is_unique


def test_the_committed_index_is_read_identically_by_both_arms():
    """W1's reader and W3's reader must return the same list, in the same order.
    This is the join that integration rests on."""
    from src.reference.gene_index import read_gene_index

    assert read_gene_index(SHARED_VERSION) == load_gene_index(SHARED_VERSION)
