"""The locked §5 contrast: the floor bites, None is not 0.0, branches are fixed."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.reference.jobs.chen_subtype_contrast import (
    MIN_ESTIMABLE_PER_ARM,
    ContrastError,
    falsifier_branch,
    per_gene,
    welch_difference,
)


def _values(ad: list[float | None], ser: list[float | None], gene: str = "CDX2") -> pd.DataFrame:
    rows = []
    for arm, vals in (("AD", ad), ("SER", ser)):
        for i, v in enumerate(vals):
            rows.append({
                "patient_id": f"{arm}{i}", "arm": arm, "gene": gene,
                "weighting": "doubly_robust", "intrinsic": v,
                # Every term needs variation: a zero-variance arm has no
                # standard error, and per_gene walks all three.
                "compositional": 0.5 + 0.01 * i, "interaction": 0.25 - 0.01 * i,
            })
    return pd.DataFrame(rows)


def test_welch_uses_unequal_variances_and_widens_when_an_arm_is_noisy():
    tight = welch_difference(np.array([1.0, 1.1, 0.9, 1.0]), np.array([0.0, 0.1, -0.1, 0.0]))
    noisy = welch_difference(np.array([1.0, 1.1, 0.9, 1.0]), np.array([-5.0, 5.0, -4.0, 4.0]))
    assert noisy["ci_high"] - noisy["ci_low"] > tight["ci_high"] - tight["ci_low"]
    assert tight["interval_method"] == "welch_t"


def test_an_arm_below_the_floor_reports_no_interval_on_any_statistic():
    """§7's last row, fixed before the estimability counts were seen."""
    ad = [1.0] * MIN_ESTIMABLE_PER_ARM
    ser = [1.0] * (MIN_ESTIMABLE_PER_ARM - 1)
    out = per_gene(_values(ad, ser))
    row = out[(out.term == "intrinsic")].iloc[0]
    assert row["estimability"] == "not_estimable"
    assert row["mean_difference"] is None
    assert row["ci_low"] is None
    assert "floor" in row["reason"]


def test_a_not_estimable_patient_is_dropped_and_counted_not_zeroed():
    """Invariant 1. Zeroing would drag the arm mean toward zero."""
    # Values must vary: a zero-variance arm has no standard error and is refused.
    ad = [2.0 + 0.01 * i for i in range(MIN_ESTIMABLE_PER_ARM)] + [None]
    ser = [1.0 + 0.01 * i for i in range(MIN_ESTIMABLE_PER_ARM)]
    out = per_gene(_values(ad, ser))
    row = out[out.term == "intrinsic"].iloc[0]
    assert row["n_AD"] == MIN_ESTIMABLE_PER_ARM + 1
    assert row["n_AD_estimable"] == MIN_ESTIMABLE_PER_ARM
    assert row["n_AD_not_estimable"] == 1
    assert row["mean_AD"] == pytest.approx(2.0 + 0.01 * (MIN_ESTIMABLE_PER_ARM - 1) / 2)


def test_the_compositional_term_never_carries_a_test():
    """§4: it is confounded with the selection criterion by construction."""
    out = per_gene(_values([1.0 + 0.01 * i for i in range(10)],
                           [0.01 * i for i in range(10)]))
    assert set(out.loc[out.term == "compositional", "standing"]) == {"descriptive_only"}
    assert not out.loc[out.term == "compositional", "tested"].any()
    assert out.loc[out.term == "intrinsic", "tested"].all()


@pytest.mark.parametrize(
    ("cdx2", "cdx2_diff", "guca2a", "expected"),
    [
        (True, +1.0, False, "PREDICTION REPLICATES"),
        (False, 0.0, False, "NOT SEPARABLE AT THIS RESOLUTION"),
        (False, 0.0, True, "CONTRADICTS THE PUBLISHED PREDICTION"),
        (True, +1.0, True, "BOTH EXCLUDE ZERO"),
        (True, -1.0, False, "CDX2 SEPARATES IN THE OPPOSITE DIRECTION"),
        (True, -1.0, True, "CDX2 SEPARATES IN THE OPPOSITE DIRECTION"),
    ],
)
def test_every_falsifier_branch_is_pre_committed(
    cdx2: bool, cdx2_diff: float, guca2a: bool, expected: str
):
    table = pd.DataFrame([
        {"gene": "CDX2", "term": "intrinsic", "weighting": "w",
         "excludes_zero": cdx2, "mean_difference": cdx2_diff, "estimability": "estimated"},
        {"gene": "GUCA2A", "term": "intrinsic", "weighting": "w",
         "excludes_zero": guca2a, "mean_difference": 0.0, "estimability": "estimated"},
    ])
    assert falsifier_branch(table, "w")["branch"] == expected


def test_a_downward_cdx2_exclusion_is_never_read_as_replication():
    """§6 predicts CDX2 falls in SER, so AD minus SER is positive."""
    table = pd.DataFrame([
        {"gene": "CDX2", "term": "intrinsic", "weighting": "w",
         "excludes_zero": True, "mean_difference": -2.0, "estimability": "estimated"},
        {"gene": "GUCA2A", "term": "intrinsic", "weighting": "w",
         "excludes_zero": False, "mean_difference": 0.0, "estimability": "estimated"},
    ])
    out = falsifier_branch(table, "w")
    assert out["branch"] != "PREDICTION REPLICATES"
    assert out["cdx2_in_predicted_direction"] is False


def test_a_floor_failure_short_circuits_the_branch():
    table = pd.DataFrame([
        {"gene": "CDX2", "term": "intrinsic", "weighting": "w",
         "excludes_zero": None, "mean_difference": None, "estimability": "not_estimable"},
        {"gene": "GUCA2A", "term": "intrinsic", "weighting": "w",
         "excludes_zero": True, "mean_difference": 1.0, "estimability": "estimated"},
    ])
    assert falsifier_branch(table, "w")["branch"] == "NOT ESTIMABLE"


def test_welch_refuses_a_single_observation_arm():
    with pytest.raises(ContrastError):
        welch_difference(np.array([1.0]), np.array([0.0, 1.0]))
