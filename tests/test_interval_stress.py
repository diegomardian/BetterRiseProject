"""Forcing inputs for the interval stress calibration.

Each test makes one claim fail if it stops holding: that the effect shapes are
standardised (so tau means the same thing under every regime), that the
extraction of the thinning did not change ``simulate_deltas``, that the Wilson
interval does not put mass below zero, and that the audit's rule needs both the
tolerance and the interval before it fires.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.reference.interval_calibration import (
    CalibrationError,
    simulate_deltas,
    simulate_deltas_from_log_fc,
)
from src.reference.interval_stress import (
    REGIMES,
    STRESS_TAU,
    StressRegime,
    audit_flags,
    draw_log_fc,
    input_validation_problems,
    original_check_flags,
    standardised_effects,
    stress_calibration_table,
    summarise_stress,
    wilson_interval,
)


def _regime(name: str) -> StressRegime:
    return next(r for r in REGIMES if r.name == name)


@pytest.mark.parametrize("name", [r.name for r in REGIMES])
def test_effect_shapes_are_standardised(name):
    """Mean 0 and variance 1, whatever the shape.

    If a regime had a different variance then ``tau`` would mean a different
    spread in each row, and the stress comparison would confound shape with
    scale.
    """
    rng = np.random.default_rng(0)
    draws = standardised_effects(_regime(name), 400_000, rng)
    assert abs(float(draws.mean())) < 0.01
    assert abs(float(draws.var()) - 1.0) < 0.02


def test_skewed_is_actually_skewed():
    rng = np.random.default_rng(1)
    skewed = standardised_effects(_regime("skewed"), 200_000, rng)
    gaussian = standardised_effects(_regime("gaussian"), 200_000, rng)
    from scipy import stats

    assert float(stats.skew(skewed)) > 0.5
    assert abs(float(stats.skew(gaussian))) < 0.05


def test_heavy_tailed_has_heavier_tails_than_gaussian():
    rng = np.random.default_rng(2)
    heavy = standardised_effects(_regime("heavy_tailed"), 200_000, rng)
    gaussian = standardised_effects(_regime("gaussian"), 200_000, rng)
    from scipy import stats

    # Excess kurtosis measures the tail directly; a single high quantile is
    # noisy and, after standardising to unit variance, does not separate them
    # reliably. Student-t(4) has very heavy tails; the normal has none.
    assert float(stats.kurtosis(heavy)) > 3.0
    assert abs(float(stats.kurtosis(gaussian))) < 0.2


def test_spike_slab_is_half_zeros():
    rng = np.random.default_rng(3)
    draws = standardised_effects(_regime("spike_slab"), 200_000, rng)
    assert 0.49 < float(np.mean(draws == 0.0)) < 0.51


def test_draw_log_fc_has_sd_tau_for_every_regime():
    for name in ("gaussian", "skewed", "heavy_tailed", "spike_slab"):
        rng = np.random.default_rng(4)
        draws = draw_log_fc(_regime(name), STRESS_TAU, 400_000, rng)
        assert abs(float(draws.std()) - STRESS_TAU) < 0.02 * STRESS_TAU


def test_regime_parameters_are_checked():
    with pytest.raises(CalibrationError, match="heavy-tailed"):
        StressRegime("bad", "heavy_tailed", "", cp10k=3.0, df=1.5)
    with pytest.raises(CalibrationError, match="sigma"):
        StressRegime("bad", "skewed", "", cp10k=3.0, sigma=0.0)
    with pytest.raises(CalibrationError, match="spike_prob"):
        StressRegime("bad", "spike_slab", "", cp10k=3.0, spike_prob=1.0)
    with pytest.raises(CalibrationError, match="unknown kind"):
        StressRegime("bad", "wobbly", "", cp10k=3.0)


def test_thinning_extraction_preserves_simulate_deltas():
    """The refactor that let stress generators reuse the thinning is exact.

    ``simulate_deltas`` now draws the log fold change and delegates. If the
    split changed the order of the RNG calls or dropped the boundary rule, every
    committed calibration would move. Reconstruct the same draw by hand and
    check the two paths agree to the last bit.
    """
    n = np.array([12, 40, 7, 133])
    depth = np.array([900.0, 2100.0, 400.0, 3300.0])

    r_direct = np.random.default_rng(7)
    direct = simulate_deltas(
        n_cells=n, depth=depth, cp10k=3.0, fold_change=0.5, tau=0.2,
        rng=r_direct,
    )

    r_split = np.random.default_rng(7)
    log_fc = np.full(n.size, np.log(0.5)) + r_split.normal(0.0, 0.2, n.size)
    split = simulate_deltas_from_log_fc(
        n_cells=n, depth=depth, cp10k=3.0, log_fc=log_fc, rng=r_split,
    )
    assert np.array_equal(direct, split)


def test_thinning_refuses_a_broadcast_effect():
    with pytest.raises(CalibrationError, match="per patient"):
        simulate_deltas_from_log_fc(
            n_cells=np.array([10, 20]), depth=np.array([1.0, 1.0]), cp10k=3.0,
            log_fc=np.array([0.0]), rng=np.random.default_rng(0),
        )


def test_wilson_interval_matches_a_known_value_and_stays_in_range():
    lo, hi = wilson_interval(5, 100)
    assert lo == pytest.approx(0.0216, abs=1e-3)
    assert hi == pytest.approx(0.1117, abs=1e-3)
    # A zero count must not produce a negative lower bound: the Wald interval
    # does, and a "worse than impossible" rate is how noise becomes a finding.
    lo0, hi0 = wilson_interval(0, 100)
    assert lo0 >= 0.0
    assert lo0 == pytest.approx(0.0, abs=1e-12)
    assert hi0 > 0.0


def test_input_validation_catches_degenerate_samples():
    assert input_validation_problems(np.array([1.0, 2.0, 3.0])) == []
    assert "no spread" in input_validation_problems(np.array([2.0, 2.0, 2.0]))
    assert "non-finite value" in input_validation_problems(
        np.array([1.0, np.nan, 3.0])
    )
    assert "fewer than two patients" in input_validation_problems(np.array([1.0]))


def test_original_check_is_undefined_where_the_closed_form_is_not():
    assert original_check_flags(0.10, float("nan")) is None
    assert original_check_flags(0.10, 0.05) is True
    assert original_check_flags(0.06, 0.05) is False


def test_audit_needs_both_the_tolerance_and_the_interval():
    """A rate above tolerance whose interval still contains nominal is noise.

    This is the forcing input for the one rule that keeps Monte Carlo spread
    from being recorded as a miscalibration finding.
    """
    assert audit_flags(0.08, 0.06) is True
    assert audit_flags(0.08, 0.049) is False   # interval includes nominal
    assert audit_flags(0.06, 0.06) is False    # inside tolerance


def test_stress_cell_is_deterministic_and_in_range():
    cohort = (np.array([12, 30, 45, 80]), np.array([1200.0, 2000.0, 900.0, 1500.0]))
    frame = stress_calibration_table(
        cohorts={"synthetic": cohort},
        regimes=(_regime("gaussian"),),
        methods=("percentile", "student_t"),
        seeds=(11, 11), n_trials=50,
    )
    # Two seeds are two cells per method; the same seed twice must be
    # byte-identical within a method, and methods differ by construction.
    assert len(frame) == 4
    for _, block in frame.groupby("method"):
        assert block["false_positive_rate"].nunique() == 1
    assert frame["false_positive_rate"].between(0, 1).all()


def test_summary_counts_the_families_separately():
    rows = []
    for seed, fpr, verdict, inp, orig, audit in (
        (1, 0.12, "MISCALIBRATED", 0, True, True),
        (2, 0.11, "MISCALIBRATED", 0, False, True),
    ):
        rows.append({
            "cohort": "c", "n_patients": 10, "regime": "skewed",
            "regime_kind": "skewed", "method": "percentile", "seed": seed,
            "n_trials": 100, "false_positive_rate": fpr,
            "input_validation_failures": inp,
            "original_check_flags": orig, "audit_flags": audit,
            "verdict": verdict,
        })
    # A non-percentile row where the original check is undefined (None).
    rows.append({
        "cohort": "c", "n_patients": 10, "regime": "skewed",
        "regime_kind": "skewed", "method": "student_t", "seed": 1,
        "n_trials": 100, "false_positive_rate": 0.05,
        "input_validation_failures": 0,
        "original_check_flags": None, "audit_flags": False,
        "verdict": "CALIBRATED",
    })
    summary = summarise_stress(pd.DataFrame(rows))
    pct = summary[summary.method == "percentile"].iloc[0]
    assert pct["detection_rate_input_validation"] == 0.0
    assert pct["detection_rate_original_check"] == 0.5
    assert pct["detection_rate_audit_check"] == 1.0
    stu = summary[summary.method == "student_t"].iloc[0]
    # None is "undefined", not a fire: the original check is percentile-only.
    assert stu["detection_rate_original_check"] == 0.0
