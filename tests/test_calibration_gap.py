"""The calibration sweep, and the two properties it exists to expose.

Two claims in the paper are properties of *this* module rather than of the data,
and both are asserted here rather than read off a run:

1. The committed grid cannot reach any mature count between 40 and 100, so a
   cutpoint returned from it is the first grid point past the criteria and not
   where the criteria first hold.
2. The oracle arm's recovery ratio carries no information about the estimator,
   because the estimate reproduces the sweep's own realised summary statistics
   exactly. ``recovery_summary`` persists that as a column so the property is
   visible in the result table instead of having to be re-derived.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.harness.calibration_gap import (
    COMMITTED_FRACTIONS,
    EXTENDED_FRACTIONS,
    GAP_FILLING_FRACTIONS,
    N_CELLS,
    _pool_mask,
    counts_reachable,
    recovery_summary,
)

# ---------------------------------------------------------------------------
# The grid reports its own resolution
# ---------------------------------------------------------------------------


def test_the_committed_grid_cannot_reach_the_interval_the_answer_lives_in():
    """The defect, stated as a property rather than observed as an outcome."""
    reachable = counts_reachable(COMMITTED_FRACTIONS)
    assert reachable == [0, 20, 40, 100, 200, 400, 800]
    assert not [n for n in reachable if 40 < n < 100], (
        "the committed grid must have no point strictly between 40 and 100 — "
        "that gap is why a returned cutpoint of 100 was indistinguishable from "
        "any value in it"
    )


def test_the_extended_grid_fills_that_interval_and_changes_nothing_else():
    extended = counts_reachable(EXTENDED_FRACTIONS)
    committed = counts_reachable(COMMITTED_FRACTIONS)
    assert set(committed) < set(extended), "the extended grid must be a superset"
    added = sorted(set(extended) - set(committed))
    assert added == [50, 60, 70, 80, 90, 110, 120, 130, 150]
    assert len([n for n in added if 40 < n < 100]) == 5


def test_the_added_fractions_are_the_only_difference():
    assert set(EXTENDED_FRACTIONS) == set(COMMITTED_FRACTIONS) | set(
        GAP_FILLING_FRACTIONS
    )


def test_counts_are_a_property_of_the_fixed_cell_budget():
    """A grid over fractions is only a grid over counts because ``n_cells`` is
    held fixed. Change the budget and the same fractions reach other counts —
    which is the whole reason the resolution was invisible."""
    assert counts_reachable((0.05,), n_cells=N_CELLS) == [100]
    assert counts_reachable((0.05,), n_cells=4000) == [200]


# ---------------------------------------------------------------------------
# The draw pool
# ---------------------------------------------------------------------------


def test_the_pooled_reading_draws_from_both_tissues():
    tissue = pd.Series(["normal", "tumour", "normal", "border"])
    assert _pool_mask(tissue, "pooled").tolist() == [True, True, True, False]


def test_the_reference_reading_draws_from_one():
    tissue = pd.Series(["normal", "tumour", "normal", "border"])
    assert _pool_mask(tissue, "reference").tolist() == [True, False, True, False]


def test_an_unknown_pool_is_refused_rather_than_silently_defaulted():
    with pytest.raises(ValueError, match="unknown pool"):
        _pool_mask(pd.Series(["normal"]), "whatever")


# ---------------------------------------------------------------------------
# The recovery curve cannot see the estimator
# ---------------------------------------------------------------------------


def _sweep(*, residual: float) -> pd.DataFrame:
    """A sweep whose estimate misses realised truth by exactly ``residual``."""
    n = 8
    realised = np.linspace(-5.0, -1.0, n)
    return pd.DataFrame(
        {
            "arm": ["oracle"] * n,
            "replicate": range(n),
            "frac_mature_tumour": [0.05] * n,
            "shift": [0.5] * n,
            "n_cells_mature": [100] * n,
            "intrinsic_true_realised": realised,
            "intrinsic_true_parametric": realised - 0.1,
            "intrinsic_hat": realised + residual,
            "attenuation_ratio": (realised + residual) / (realised - 0.1),
        }
    )


def test_a_zero_residual_is_reported_as_zero_not_hidden():
    """THE check. An estimator reproducing the generator's own sufficient
    statistics gives a recovery curve that contains no estimator, and the only
    visible signature is that this column is identically zero."""
    out = recovery_summary(_sweep(residual=0.0))
    assert out["max_abs_residual_vs_realised"].max() == 0.0


def test_an_estimator_that_does_not_reproduce_its_input_is_distinguished():
    """The column has to be capable of being non-zero, or it is not a check."""
    out = recovery_summary(_sweep(residual=0.25))
    assert out["max_abs_residual_vs_realised"].max() == pytest.approx(0.25)


def test_the_null_arm_ratio_is_undefined_and_counted_rather_than_dropped():
    """At ``shift = 1.0`` the parametric intrinsic term is exactly zero, so the
    ratio divides by zero. The design's own null control is the one point where
    this validation statistic cannot be evaluated — it is counted, not silently
    dropped, and the surviving median is not computed over an inf."""
    sweep = _sweep(residual=0.0)
    sweep["shift"] = 1.0
    sweep["intrinsic_true_parametric"] = 0.0
    sweep["attenuation_ratio"] = np.inf

    out = recovery_summary(sweep)
    assert out["n_ratio_undefined"].sum() == len(sweep)
    assert out["n_replicates"].sum() == len(sweep)
    assert out["ratio_median"].isna().all()


def test_undefined_ratios_do_not_contaminate_a_defined_median():
    defined = _sweep(residual=0.0)
    undefined = _sweep(residual=0.0)
    undefined["frac_mature_tumour"] = 0.0
    undefined["attenuation_ratio"] = np.inf

    out = recovery_summary(pd.concat([defined, undefined], ignore_index=True))
    row = out[out["frac_mature_tumour"] == 0.05].iloc[0]
    assert np.isfinite(row["ratio_median"])
    assert row["n_ratio_undefined"] == 0


def test_only_the_oracle_arm_is_summarised():
    """The bulk arm has no per-patient interval and is not what the cutpoints
    are calibrated on. Folding it in would average two different questions."""
    sweep = pd.concat(
        [_sweep(residual=0.0), _sweep(residual=9.0).assign(arm="bulk")],
        ignore_index=True,
    )
    out = recovery_summary(sweep)
    assert out["max_abs_residual_vs_realised"].max() == 0.0
