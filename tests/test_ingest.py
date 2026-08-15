"""Per-study QC flagging. execution_plan.md W1 task table: "per-study
thresholds, not global" -- the same shape W4 mirrors for the Lee cohorts.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.estimator.ingest import correct_ambient, flag_doublets, qc_flags


def _metrics(seed=0):
    """Two studies with different library-size scales -- a global cutoff
    would either over-filter the low-depth study or under-filter the
    high-depth one. A few planted outliers in each."""
    rng = np.random.default_rng(seed)
    smc = pd.DataFrame(
        {
            "study_id": "GSE132465",
            "n_genes": rng.normal(1500, 200, 200),
            "n_counts": rng.normal(5000, 800, 200),
            "pct_mito": rng.normal(5, 2, 200).clip(0),
        }
    )
    kul3 = pd.DataFrame(
        {
            "study_id": "GSE144735",
            "n_genes": rng.normal(3000, 300, 200),
            "n_counts": rng.normal(9000, 1000, 200),
            "pct_mito": rng.normal(5, 2, 200).clip(0),
        }
    )
    metrics = pd.concat([smc, kul3], ignore_index=True)
    # Plant unambiguous outliers in each study.
    metrics.loc[0, "n_genes"] = 5  # SMC: near-empty droplet
    metrics.loc[200, "n_counts"] = 200000  # KUL3: doublet-scale library
    metrics.loc[5, "pct_mito"] = 80  # dying cell, either study
    return metrics


def test_qc_flags_catches_planted_outliers_in_both_studies():
    flags = qc_flags(_metrics())
    assert flags.loc[0]  # near-empty SMC droplet
    assert flags.loc[200]  # doublet-scale KUL3 library
    assert flags.loc[5]  # high-mito cell


def test_qc_flags_thresholds_are_computed_per_study_not_globally():
    """A cell with GSE132465-typical depth must not be flagged just because
    GSE144735 runs deeper -- that would be a global, not per-study, cutoff."""
    metrics = _metrics()
    smc_typical = metrics[(metrics["study_id"] == "GSE132465") & ~metrics.index.isin([0, 5])]
    flags = qc_flags(metrics)
    assert not flags.loc[smc_typical.index].any()


def test_qc_flags_requires_the_expected_columns():
    bad = _metrics().drop(columns=["pct_mito"])
    with pytest.raises(ValueError, match="missing column"):
        qc_flags(bad)


def test_qc_flags_mito_cap_is_a_single_hard_threshold():
    metrics = _metrics()
    lenient = qc_flags(metrics, max_pct_mito=100)
    strict = qc_flags(metrics, max_pct_mito=1)
    assert strict.sum() > lenient.sum()


def test_doublet_and_ambient_correction_are_explicit_stubs_pending_real_data():
    """These need actual GSE132465/GSE144735 matrices to fit -- not fakeable
    without them. Confirm they fail loudly, not silently."""
    with pytest.raises(NotImplementedError, match="real Lee count matrices"):
        flag_doublets()
    with pytest.raises(NotImplementedError, match="real Lee count matrices"):
        correct_ambient()
