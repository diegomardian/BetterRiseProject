"""The output contract. CLAUDE.md invariants 1 and 10.

If a test in this file goes red, do not adjust the test. The contract is frozen
and the assertion is the point.
"""

from __future__ import annotations

import json

import pandas as pd
import pytest

from src.schema import (
    COLUMNS,
    ESTIMABILITY,
    GRANULARITY_RUNGS,
    LABELING_AXES,
    WEIGHTINGS,
    SchemaViolation,
    coerce_results,
    empty_results_frame,
    read_results,
    validate_results,
    write_results,
)


def _row(**overrides):
    base = {
        "patient_id": "P001",
        "study_id": "GSE178341",
        "gene": "GUCA2A",
        "granularity_rung": "lineage",
        "labeling_axis": "stem_pole",
        "n_cells_mature": 120,
        "compositional": -0.8,
        "intrinsic": -0.1,
        "interaction": 0.02,
        "weighting": "normal",
        "estimability": "ok",
        "ci_low": -0.3,
        "ci_high": 0.1,
    }
    base.update(overrides)
    return base


def _frame(*rows):
    return coerce_results(pd.DataFrame(list(rows) or [_row()]))


# ---------------------------------------------------------------------------
# Frozen vocabulary — execution_plan.md §3.3
# ---------------------------------------------------------------------------


def test_schema_columns_are_exactly_the_frozen_set():
    assert list(COLUMNS) == [
        "patient_id",
        "study_id",
        "gene",
        "granularity_rung",
        "labeling_axis",
        "n_cells_mature",
        "compositional",
        "intrinsic",
        "interaction",
        "weighting",
        "estimability",
        "ci_low",
        "ci_high",
    ]


def test_frozen_enums():
    assert GRANULARITY_RUNGS == ("epithelial", "lineage", "crypt_position", "best4")
    assert LABELING_AXES == ("stem_pole", "opposite_lineage", "chromatin", "spatial")
    assert WEIGHTINGS == ("normal", "tumour", "doubly_robust")
    assert ESTIMABILITY == ("ok", "wide_interval", "not_estimable")


def test_empty_frame_round_trips():
    df = empty_results_frame()
    assert list(df.columns) == list(COLUMNS)
    validate_results(df)  # an empty frame is valid, just uninformative


# ---------------------------------------------------------------------------
# INVARIANT 1 — None is not 0.0. The one that matters.
# ---------------------------------------------------------------------------


def test_not_estimable_with_none_intrinsic_is_valid():
    df = _frame(_row(n_cells_mature=3, estimability="not_estimable", intrinsic=None))
    validate_results(df)  # compositional term is still reported — the row stays


def test_not_estimable_with_zero_intrinsic_is_rejected():
    """The single most likely route to a wrong conclusion in this project."""
    df = _frame(_row(n_cells_mature=3, estimability="not_estimable", intrinsic=0.0))
    with pytest.raises(SchemaViolation, match="None is not 0.0"):
        validate_results(df)


def test_not_estimable_with_any_intrinsic_value_is_rejected():
    df = _frame(_row(n_cells_mature=3, estimability="not_estimable", intrinsic=-0.42))
    with pytest.raises(SchemaViolation, match="not_estimable"):
        validate_results(df)


def test_ok_row_missing_its_intrinsic_estimate_is_rejected():
    """The converse leak: an estimate that vanished but kept its 'ok' label."""
    df = _frame(_row(estimability="ok", intrinsic=None))
    with pytest.raises(SchemaViolation, match="null intrinsic"):
        validate_results(df)


# ---------------------------------------------------------------------------
# Structural checks
# ---------------------------------------------------------------------------


def test_unknown_enum_value_is_rejected():
    df = _frame(_row(granularity_rung="single_cell"))
    with pytest.raises(SchemaViolation, match="frozen set"):
        validate_results(df)


def test_unexpected_column_is_rejected():
    df = pd.DataFrame([_row() | {"pvalue": 0.01}])
    with pytest.raises(SchemaViolation, match="schema is frozen"):
        coerce_results(df)


def test_missing_column_is_rejected():
    df = pd.DataFrame([{k: v for k, v in _row().items() if k != "intrinsic"}])
    with pytest.raises(SchemaViolation, match="missing required columns"):
        coerce_results(df)


def test_duplicate_keys_are_rejected():
    df = _frame(_row(), _row())
    with pytest.raises(SchemaViolation, match="duplicate key"):
        validate_results(df)


def test_three_weightings_are_three_distinct_rows():
    """Both weightings plus doubly-robust — not folded together."""
    df = _frame(*(_row(weighting=w) for w in WEIGHTINGS))
    validate_results(df)


def test_inverted_interval_is_rejected():
    df = _frame(_row(ci_low=0.5, ci_high=-0.5))
    with pytest.raises(SchemaViolation, match="ci_low > ci_high"):
        validate_results(df)


def test_half_an_interval_is_rejected():
    df = _frame(_row(ci_high=None))
    with pytest.raises(SchemaViolation, match="only one CI bound"):
        validate_results(df)


# ---------------------------------------------------------------------------
# INVARIANT 10 — every result carries a git sha and a fixed seed
# ---------------------------------------------------------------------------


def test_writer_stamps_provenance_and_round_trips(tmp_path):
    df = _frame(
        _row(),
        _row(patient_id="P002", n_cells_mature=4, estimability="not_estimable", intrinsic=None),
    )
    path = write_results(
        df, name="unit_test", seed=20260815, results_dir=tmp_path, allow_dirty=True
    )

    assert path.exists()
    meta = json.loads(path.with_suffix("").with_suffix(".meta.json").read_text())
    assert meta["seed"] == 20260815
    assert meta["git_sha"]
    assert meta["n_rows"] == 2
    assert meta["n_not_estimable"] == 1
    assert path.parent.name.endswith(meta["git_sha_short"])

    back = read_results(path)
    assert len(back) == 2
    assert back.loc[back["patient_id"] == "P002", "intrinsic"].isna().all()


def test_writer_refuses_an_invalid_frame(tmp_path):
    """Bad rows never reach disk — the assertion is in the writer, not in review."""
    df = _frame(_row(n_cells_mature=0, estimability="not_estimable", intrinsic=0.0))
    with pytest.raises(SchemaViolation):
        write_results(df, name="bad", seed=1, results_dir=tmp_path, allow_dirty=True)
    assert not list(tmp_path.rglob("*.parquet"))


def test_writer_requires_an_explicit_seed(tmp_path):
    with pytest.raises(TypeError):
        write_results(_frame(), name="no_seed", results_dir=tmp_path)  # type: ignore[call-arg]
