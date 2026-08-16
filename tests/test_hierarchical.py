"""Hierarchical (mixed-effects) patient-clustered CIs -- the model-based
counterpart to kitagawa.bootstrap_over_patients.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.estimator.hierarchical import hierarchical_intrinsic_ci


def _synthetic_cells(
    n_patients=15,
    n_cells_per_patient=120,
    frac_normal=0.4,
    frac_tumour=0.4,
    mean_normal=10.0,
    mean_tumour=6.0,
    patient_sd=0.5,
    noise_sd=1.0,
    seed=0,
):
    """Cell-level data with a genuine patient random effect on expression, so
    a plain per-cell model (ignoring clustering) would understate uncertainty
    relative to a model that accounts for it."""
    rng = np.random.default_rng(seed)
    frames = []
    tissue_specs = (
        ("normal", frac_normal, mean_normal),
        ("tumour", frac_tumour, mean_tumour),
    )
    for p in range(n_patients):
        patient_effect = rng.normal(0, patient_sd)
        for tissue, frac, mean in tissue_specs:
            mature = rng.random(n_cells_per_patient) < frac
            expression = np.where(
                mature,
                mean + patient_effect + rng.normal(0, noise_sd, n_cells_per_patient),
                rng.normal(0.5, 0.2, n_cells_per_patient),
            )
            frames.append(
                pd.DataFrame(
                    {
                        "patient_id": f"P{p}",
                        "tissue": tissue,
                        "mature": mature.astype(int),
                        "expression": expression,
                    }
                )
            )
    return pd.concat(frames, ignore_index=True)


def test_recovers_a_pure_intrinsic_ground_truth_within_tolerance():
    """frac unchanged (0.4 -> 0.4), mean drops 10 -> 6: true intrinsic =
    0.4*(6-10) = -1.6, true compositional = 0."""
    cells = _synthetic_cells(frac_normal=0.4, frac_tumour=0.4, mean_normal=10.0, mean_tumour=6.0)
    result = hierarchical_intrinsic_ci(cells, seed=1, n_draws=1000)
    assert result["intrinsic"] == pytest.approx(-1.6, abs=0.5)
    assert result["compositional"] == pytest.approx(0.0, abs=0.1)
    assert result["ci_low"] < result["ci_high"]


def test_recovers_a_pure_compositional_ground_truth_within_tolerance():
    """frac drops 0.4 -> 0.1, mean unchanged: true compositional =
    (0.1-0.4)*10 = -3.0, true intrinsic = 0."""
    cells = _synthetic_cells(frac_normal=0.4, frac_tumour=0.1, mean_normal=10.0, mean_tumour=10.0)
    result = hierarchical_intrinsic_ci(cells, seed=1, n_draws=1000)
    assert result["compositional"] == pytest.approx(-3.0, abs=0.5)
    assert result["intrinsic"] == pytest.approx(0.0, abs=0.3)


def test_wider_patient_variance_widens_the_interval():
    """More between-patient noise in the mean should widen the CI on
    intrinsic -- the whole reason to model patient as a grouping factor."""
    tight = _synthetic_cells(patient_sd=0.05, seed=2)
    loose = _synthetic_cells(patient_sd=3.0, seed=2)
    tight_result = hierarchical_intrinsic_ci(tight, seed=1, n_draws=1000)
    loose_result = hierarchical_intrinsic_ci(loose, seed=1, n_draws=1000)
    tight_width = tight_result["ci_high"] - tight_result["ci_low"]
    loose_width = loose_result["ci_high"] - loose_result["ci_low"]
    assert loose_width > tight_width


def test_is_reproducible_given_the_same_seed():
    cells = _synthetic_cells()
    a = hierarchical_intrinsic_ci(cells, seed=7, n_draws=200)
    b = hierarchical_intrinsic_ci(cells, seed=7, n_draws=200)
    assert a == b


def test_rejects_cells_missing_required_columns():
    cells = _synthetic_cells().drop(columns=["mature"])
    with pytest.raises(ValueError, match="missing column"):
        hierarchical_intrinsic_ci(cells, seed=1)


def test_rejects_a_tissue_value_that_is_not_normal_or_tumour():
    cells = _synthetic_cells()
    cells.loc[0, "tissue"] = "metastasis"
    with pytest.raises(ValueError, match="normal.*tumour"):
        hierarchical_intrinsic_ci(cells, seed=1)


def test_rejects_non_positive_n_draws():
    with pytest.raises(ValueError, match="n_draws"):
        hierarchical_intrinsic_ci(_synthetic_cells(), seed=1, n_draws=0)


def test_requires_an_explicit_seed():
    with pytest.raises(TypeError):
        hierarchical_intrinsic_ci(_synthetic_cells())  # type: ignore[call-arg]
