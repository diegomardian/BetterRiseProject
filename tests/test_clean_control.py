from __future__ import annotations

import numpy as np

from src.harness.clean_control import REFERENCE_TOLERANCE, run, summarise


def _table():
    return summarise(run(seed=42, n_replicates=800, n_patients=200))


def test_equal_residual_class_contains_good_and_bad_intervals():
    table = _table().set_index(["requested_effect", "estimator"])
    for effect in (0.0, 0.5):
        calibrated = table.loc[(effect, "empirical-mean-calibrated")]
        narrow = table.loc[(effect, "same-point-narrow-interval")]
        assert calibrated["reference_audit_outcome"] == "numerical_equality"
        assert narrow["reference_audit_outcome"] == "numerical_equality"
        assert calibrated["interval_coverage"] > 0.90
        assert narrow["interval_coverage"] < 0.50


def test_departure_class_contains_valid_and_biased_estimators():
    table = _table().set_index(["requested_effect", "estimator"])
    for effect in (0.0, 0.5):
        valid = table.loc[(effect, "independent-half-mean")]
        biased = table.loc[(effect, "known-bias-plus-0.5")]
        assert valid["reference_audit_outcome"] == "numerical_departure"
        assert biased["reference_audit_outcome"] == "numerical_departure"
        assert valid["interval_coverage"] > 0.90
        assert abs(biased["bias"] - 0.5) < 0.02
        assert biased["interval_coverage"] < 0.10


def test_null_effect_keeps_difference_metrics_and_no_ratio_is_needed():
    table = _table()
    null = table[table["requested_effect"] == 0.0]
    assert np.isfinite(null["bias"]).all()
    assert np.isfinite(null["rmse"]).all()
    assert not any("ratio" in column for column in table.columns)


def test_invalid_outputs_are_counted_and_cannot_certify_equality():
    runs = run(seed=7, n_replicates=20, n_patients=50, effects=(0.0,))
    target = runs["estimator"] == "empirical-mean-calibrated"
    first = runs[target].index[0]
    runs.loc[first, ["estimate", "ci_low", "ci_high", "residual_vs_reference"]] = np.nan
    runs.loc[first, "valid_output"] = False
    summary = summarise(runs, tolerance=REFERENCE_TOLERANCE)
    row = summary[summary["estimator"] == "empirical-mean-calibrated"].iloc[0]
    assert row["n_invalid"] == 1
    assert row["reference_audit_outcome"] == "insufficient_valid_pairs"
