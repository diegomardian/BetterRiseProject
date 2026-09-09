"""D2's feasibility gate is outcome-blind and patient-level by construction."""

import json

import numpy as np
import pandas as pd
import pytest

from src.bulk.d2_feasibility_gate import (
    D2FeasibilityError,
    evaluate,
    main,
    primary_tumour_values,
    resolve_gene,
)

GENE_ID = "ENSG00000112345"


def _manifest(barcodes, *, sample_types=None, patients=None):
    n = len(barcodes)
    return pd.DataFrame({
        "barcode": barcodes,
        "patient_id": patients or [f"P{i:03d}" for i in range(n)],
        "sample_type": sample_types or ["01"] * n,
    })


def _expression(values, *, barcodes=None):
    barcodes = barcodes or [f"S{i:03d}" for i in range(len(values))]
    return pd.DataFrame({GENE_ID: values}, index=barcodes)


def _map():
    return pd.DataFrame({"gene_symbol": ["GUCA2A"], "ensembl_id": [GENE_ID]})


def test_pass_requires_patient_level_spread_on_both_sides_of_median():
    values = pd.Series([0.0] * 50 + [1.0] * 50, index=[f"P{i}" for i in range(100)])

    got = evaluate(values, gene_id=GENE_ID).iloc[0]

    assert got["verdict"] == "PASS"
    assert got["n_participants"] == 100
    assert got["iqr"] == pytest.approx(1.0)
    assert got["n_below_median"] == got["n_above_median"] == 50


def test_iqr_failure_is_a_stop_even_when_both_median_tails_exist():
    values = pd.Series(np.linspace(0.0, 0.4, 100), index=[f"P{i}" for i in range(100)])

    got = evaluate(values, gene_id=GENE_ID).iloc[0]

    assert got["verdict"] == "STOP"
    assert not got["condition_iqr"]
    assert got["condition_below_median"] and got["condition_above_median"]


def test_duplicate_primary_participant_is_refused_not_selected():
    expression = _expression([0.0, 1.0], barcodes=["S1", "S2"])
    manifest = _manifest(["S1", "S2"], patients=["P1", "P1"])

    with pytest.raises(D2FeasibilityError, match="more than one primary-tumour"):
        primary_tumour_values(expression, manifest, GENE_ID)


def test_primary_manifest_rows_missing_from_expression_are_refused():
    expression = _expression([0.0], barcodes=["S1"])
    manifest = _manifest(["S1", "S2"])

    with pytest.raises(D2FeasibilityError, match="absent from expression"):
        primary_tumour_values(expression, manifest, GENE_ID)


def test_mapping_must_be_one_to_one():
    gene_map = pd.DataFrame({
        "gene_symbol": ["GUCA2A", "GUCA2A"],
        "ensembl_id": ["ENSG1", "ENSG2"],
    })

    with pytest.raises(D2FeasibilityError, match="one-to-one"):
        resolve_gene(gene_map)


def test_cli_refuses_missing_expression_before_writing(tmp_path):
    with pytest.raises(SystemExit, match="not found"):
        main([
            "--expression", str(tmp_path / "missing.parquet"),
            "--results-dir", str(tmp_path / "results"),
            "--allow-dirty",
        ])
    assert not list((tmp_path / "results").glob("*/*.parquet"))


def test_cli_writes_only_outcome_blind_gate_table(tmp_path):
    values = [0.0] * 50 + [1.0] * 50
    barcodes = [f"S{i:03d}" for i in range(100)]
    expression = _expression(values, barcodes=barcodes)
    manifest = _manifest(barcodes)
    expression_path = tmp_path / "expression.parquet"
    manifest_path = tmp_path / "manifest.tsv"
    map_path = tmp_path / "map.tsv"
    expression.to_parquet(expression_path)
    manifest.to_csv(manifest_path, sep="\t", index=False)
    _map().to_csv(map_path, sep="\t", index=False)

    assert main([
        "--expression", str(expression_path),
        "--manifest", str(manifest_path),
        "--gene-index-map", str(map_path),
        "--results-dir", str(tmp_path / "results"),
        "--allow-dirty",
    ]) == 0

    result = next((tmp_path / "results").glob("*/d2_guca2a_feasibility.parquet"))
    table = pd.read_parquet(result)
    meta = json.loads(result.with_suffix(".meta.json").read_text())
    assert table.loc[0, "verdict"] == "PASS"
    assert meta["outcomes_read"] is False
    assert meta["contract"] == "docs/d2_feasibility_gate.md"


# ---------------------------------------------------------------------------
# The two conditions that had no forcing input
# ---------------------------------------------------------------------------


def test_too_few_participants_is_a_stop_even_with_ample_spread():
    """`condition_n` had no failing test. A gate condition that is never shown
    to fail is a condition nobody has checked — this repository's own rule
    (`tests/test_checks_can_fail.py`)."""
    values = pd.Series([0.0] * 40 + [5.0] * 40, index=[f"P{i}" for i in range(80)])

    got = evaluate(values, gene_id=GENE_ID).iloc[0]

    assert got["verdict"] == "STOP"
    assert not got["condition_n"], "80 participants is below the 100 floor"
    # the spread conditions are satisfied, so only `n` is carrying the STOP
    assert got["condition_iqr"]
    assert got["condition_below_median"] and got["condition_above_median"]


def test_a_gene_at_floor_in_most_tumours_stops_on_the_median_tail():
    """`condition_below_median` had no failing test, and this is the failure
    mode D2 actually risks.

    If GUCA2A is exactly zero in more than half the tumours the median IS zero,
    nothing lies strictly below it, and there is no usable low arm for a
    continuous-marker model. That is precisely what the gate exists to catch —
    `docs/d2_feasibility_gate.md` was written because "GUCA2A loss is
    near-universal" was named as a risk and never checked.
    """
    values = pd.Series([0.0] * 70 + list(np.linspace(1.0, 6.0, 50)),
                       index=[f"P{i}" for i in range(120)])

    got = evaluate(values, gene_id=GENE_ID).iloc[0]

    assert got["median"] == 0.0
    assert got["n_below_median"] == 0
    assert got["verdict"] == "STOP"
    assert not got["condition_below_median"]
    assert got["condition_n"], "120 participants clears the floor"
    # and the exact-zero count is on the record, which is how a reader sees why
    assert got["n_exact_zero"] == 70


def test_the_gate_module_imports_nothing_that_carries_an_outcome():
    """`outcomes_read: False` in the sidecar is a DECLARATION, not a check.

    The existing CLI test asserts the sidecar says false, which it would say
    whatever the module imported. This asserts the thing itself: no module that
    reads TCGA-CDR, clinical fields or survival endpoints is reachable from
    this job's imports.
    """
    import ast

    from src.common.paths import REPO_ROOT

    source = (REPO_ROOT / "src" / "bulk" / "d2_feasibility_gate.py").read_text()
    imported: set[str] = set()
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module)
        elif isinstance(node, ast.Import):
            imported.update(alias.name for alias in node.names)

    forbidden = {"src.bulk.cdr", "src.bulk.clinical", "src.bulk.survival",
                 "src.bulk.covariates", "src.bulk.run_survival",
                 "src.bulk.run_clinical", "src.bulk.prespec"}
    assert not (imported & forbidden), (
        f"the outcome-blind gate imports {sorted(imported & forbidden)}")
