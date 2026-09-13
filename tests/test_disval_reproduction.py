"""Forcing inputs for the DIS/VAL reproduction comparison.

A reproduction is only a check if a changed cell flips it to False and a missing
cell raises. These tests force a changed mean, a changed sample size, a changed
interval verdict, and a cell present on only one side.
"""

from __future__ import annotations

import pandas as pd
import pytest

from src.reference.disval_reproduction import (
    DisvalReproductionError,
    compare_disval,
)


def _table(overrides: dict | None = None) -> pd.DataFrame:
    """One row per (half, contrast); ``overrides`` keyed by that pair.

    Keyed by a tuple rather than keyword arguments: ``**{("a","b"): ...}`` is a
    TypeError, so the obvious signature cannot express the case under test.
    """
    overrides = overrides or {}
    rows = []
    for half, mean in (("VUMC_HTAN_discovery", 0.881),
                       ("VUMC_HTAN_validation", 0.952)):
        for contrast in ("GUCA2A - ACTB", "GUCA2A - KRT8"):
            row = {
                "half": half, "statistic": "log_ratio", "contrast": contrast,
                "summary": "mean", "n_patients": 15,
                "mean": mean, "ci_low": mean - 0.2, "ci_high": mean + 0.2,
                "excludes_zero": True,
            }
            row.update(overrides.get((half, contrast), {}))
            rows.append(row)
    return pd.DataFrame(rows)


def test_an_exact_reproduction_is_flagged_true():
    out = compare_disval(_table(), _table())
    assert out["reproduces"].all()
    assert out["mean_diff"].max() == 0.0


def test_a_changed_mean_fails():
    reproduced = _table({("VUMC_HTAN_discovery", "GUCA2A - ACTB"): {"mean": 0.5}})
    out = compare_disval(reproduced, _table()).set_index(["half", "contrast"])
    assert not out.loc[("VUMC_HTAN_discovery", "GUCA2A - ACTB"), "reproduces"]
    assert out.loc[("VUMC_HTAN_validation", "GUCA2A - ACTB"), "reproduces"]


def test_a_changed_sample_size_fails():
    reproduced = _table(
        {("VUMC_HTAN_discovery", "GUCA2A - ACTB"): {"n_patients": 14}}
    )
    out = compare_disval(reproduced, _table())
    assert not out["reproduces"].all()


def test_a_changed_interval_verdict_fails():
    reproduced = _table(
        {("VUMC_HTAN_discovery", "GUCA2A - ACTB"): {"excludes_zero": False}}
    )
    out = compare_disval(reproduced, _table()).set_index(["half", "contrast"])
    assert not out.loc[("VUMC_HTAN_discovery", "GUCA2A - ACTB"), "reproduces"]


def test_a_cell_on_only_one_side_raises():
    reproduced = _table().iloc[:1]
    with pytest.raises(DisvalReproductionError, match="only one table"):
        compare_disval(reproduced, _table())


def test_missing_columns_are_refused():
    bad = _table().drop(columns=["mean"])
    with pytest.raises(DisvalReproductionError, match="missing"):
        compare_disval(bad, _table())
