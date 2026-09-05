"""The producer for the two tables that had none.

`purity_conditioned_check` and `purity_association` were tested from the day
they were written; what was missing was anything that CALLED them, so their
output reached results/ from an uncommitted script and stayed there as the only
copy. These tests cover the driver, and pin the shape of what it emits against
the committed tables it supersedes.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.bulk.run_purity_conditioned import (
    GENES,
    newest_purity_table,
    purity_series,
    resolve_genes,
)
from src.common.paths import RESULTS_DIR

#: The cluster run that replaced the 2026-08-18 originals.
PRODUCED = RESULTS_DIR / "2026-09-05_40b823b"


def test_the_two_genes_resolve_unambiguously():
    resolved = resolve_genes(GENES)
    assert set(resolved) == set(GENES)
    assert all(v.startswith("ENSG") for v in resolved.values())


def test_an_unresolvable_symbol_stops_the_run():
    with pytest.raises(SystemExit, match="cannot resolve"):
        resolve_genes(["NOT_A_GENE_SYMBOL_XYZ"])


def test_purity_series_selects_one_method_and_drops_missing_calls():
    purity = pd.DataFrame({
        "barcode": ["a", "b", "c", "d"],
        "method": ["absolute", "absolute", "aran_cpe", "absolute"],
        "purity": [0.5, np.nan, 0.7, 0.9],
    })
    got = purity_series(purity, "absolute")
    assert list(got.index) == ["a", "d"]
    assert got.tolist() == [0.5, 0.9]


def test_it_finds_a_committed_purity_table_to_condition_on():
    path = newest_purity_table()
    if path is None:
        pytest.skip("no committed purity table")
    assert path.name == "tcga_purity.parquet"


# ---------------------------------------------------------------------------
# Against the tables this supersedes


def test_the_orphan_tables_are_gone_and_their_replacements_are_clean():
    """The housekeeping this driver existed to make safe, now done.

    Sixteen of the eighteen dirty tables had a clean twin and were deleted on
    2026-09-05. These two did not, and they were the same two with no producer,
    so sweeping all eighteen would have destroyed the only copies of results
    nothing could regenerate. The driver ran on the cluster, its output
    reproduced the originals (all four association rows within 0.0003 of their
    r-squared, all thirty-two conditioned rows at 100% verdict agreement), and
    the originals are now deleted.

    This asserts the end state: no dirty copy survives, and every remaining copy
    is clean-stamped.
    """
    import json

    for name in ("tcga_premise_purity_conditioned", "tcga_purity_expression_association"):
        copies = sorted(RESULTS_DIR.glob(f"*/{name}.parquet"))
        assert copies, f"{name} has no copy at all"
        for path in copies:
            meta = json.loads(
                path.with_suffix("").with_suffix(".meta.json").read_text()
            )
            assert not meta["git_dirty"], (
                f"{path.parent.name}/{name} is dirty. Every copy of these two "
                f"should now be clean-stamped -- see 6ad0bbb."
            )


def test_the_driver_emits_the_columns_the_committed_tables_carry():
    """Shape parity, so a re-run is comparable to what it replaces."""
    path = PRODUCED / "tcga_premise_purity_conditioned.parquet"
    if not path.exists():
        pytest.skip("produced table not present")
    committed = pd.read_parquet(path)
    assert {"gene", "stratum", "purity_method", "purity_adjusted"} <= set(committed.columns)
    assert set(committed["gene"]) == set(GENES)
    assert committed["purity_adjusted"].all()
    # The driver runs every method with >= 30 calls, so the 2026-09-05 table
    # also carries aran_cpe, which the August run skipped. Require the two the
    # comparison depends on rather than pinning the set -- a third method is
    # more coverage, not a regression.
    assert {"absolute", "estimate_affy_extrapolated"} <= set(committed["purity_method"])


def test_expression_derived_purity_explains_more_variance_than_absolute():
    """The circularity Stage 4's gate excludes, visible in the committed table.

    Purity explains 0.104 of bulk CDX2 against the EXPRESSION-derived ESTIMATE
    score and 0.042 against copy-number ABSOLUTE -- two and a half times as
    much, same tumours, same gene. Part of that gap is the shared derivation,
    not confounding, which is why `src/bulk/instrument.py` filters on
    `expression_derived` rather than on a method name.
    """
    path = PRODUCED / "tcga_purity_expression_association.parquet"
    if not path.exists():
        pytest.skip("produced table not present")
    table = path and pd.read_parquet(path)
    cdx2 = table[table["gene"] == "CDX2"].set_index("purity_method")["r_squared"]
    assert cdx2["estimate_affy_extrapolated"] > 2 * cdx2["absolute"], (
        "the expression-derived purity call no longer explains materially more "
        "variance than the copy-number one -- re-derive the instrument gate's "
        "justification before relying on it"
    )


def test_the_driver_refuses_missing_inputs_rather_than_writing_a_partial_table(tmp_path):
    from src.bulk.run_purity_conditioned import main

    with pytest.raises(SystemExit, match="not found"):
        main(["--expression", str(tmp_path / "nope.parquet"),
              "--results-dir", str(tmp_path), "--allow-dirty"])
    assert not list(tmp_path.glob("*/*.parquet"))
