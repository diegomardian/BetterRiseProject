"""The interval calibration, against things whose answer is known independently.

The load-bearing tests here are the ones where two routes must agree:

* the simulated false-positive rate against the closed form, which contains no
  simulation;
* the estimated between-patient heterogeneity against a ``tau`` the generator
  was given.

A calibration measurement checked only against itself is the thing this
repository keeps finding. These check it against arithmetic and against known
truth.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.reference.interval_calibration import (
    CALIBRATED_METHOD,
    INTERVAL_METHODS,
    CalibrationError,
    bca_interval,
    calibration_verdict,
    excludes_zero,
    expected_false_positive_rate,
    heterogeneity_tau,
    percentile_interval,
    power_curve,
    rejection_rate,
    simulate_deltas,
    student_t_interval,
    width_ratio,
)

SEED = 20260906


def _cohort(n_patients: int, cells: int = 262, depth: float = 8000.0):
    return np.full(n_patients, cells), np.full(n_patients, depth)


# ---------------------------------------------------------------------------
# The closed form, and the simulation that must agree with it
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "n_patients, expected", [(4, 0.188), (10, 0.096), (20, 0.071), (44, 0.059)]
)
def test_the_closed_form_matches_the_published_rates(n_patients, expected):
    """The rates quoted in the prereg and the handoff, re-derived here.

    These are the numbers the MLH1 design and the best4 caveat are written
    against. If this test moves, those documents are wrong.
    """
    assert expected_false_positive_rate(n_patients) == pytest.approx(
        expected, abs=0.002
    )


def test_the_percentile_interval_is_narrower_than_the_calibrated_one():
    """The mechanism, measured: a normal quantile and a divide-by-n variance.

    The ratio is ``z*sqrt((n-1)/n)/t(n-1)`` and it is a function of n alone --
    no data, no gene, no scale. Measuring it on real draws is what says the
    arithmetic describes the implementation rather than a different bootstrap.
    """
    rng = np.random.default_rng(SEED)
    n_cells, depth = _cohort(10)
    widths = {"percentile": [], "student_t": []}
    for _ in range(200):
        values = simulate_deltas(
            n_cells=n_cells, depth=depth, cp10k=0.039,
            fold_change=1.0, tau=0.2, rng=rng,
        )
        for name in widths:
            lo, hi = INTERVAL_METHODS[name](values, rng=rng)
            widths[name].append(hi - lo)
    measured = np.mean(widths["percentile"]) / np.mean(widths["student_t"])
    assert measured == pytest.approx(width_ratio(10), abs=0.03)


def test_the_simulated_rate_agrees_with_the_closed_form():
    """Two independent routes to one number.

    The closed form is a floor -- it assumes the per-patient values are normal
    enough for the bootstrap mean to be, and a rare transcript's delta is also
    skewed. So the simulation may sit above it, and must not sit below it by
    more than Monte Carlo error.
    """
    rng = np.random.default_rng(SEED)
    n_cells, depth = _cohort(10)

    def null(r):
        return simulate_deltas(
            n_cells=n_cells, depth=depth, cp10k=3.0,
            fold_change=1.0, tau=0.0, rng=r,
        )

    measured = rejection_rate(null, method="percentile", rng=rng, n_trials=800)
    closed = expected_false_positive_rate(10)
    assert measured > 0.05, "the defect must reproduce, not merely be asserted"
    assert measured == pytest.approx(closed, abs=0.035)


def test_the_calibrated_interval_is_calibrated_where_the_other_is_not():
    """The comparison the choice of interval rests on, run at MLH1's n."""
    rng = np.random.default_rng(SEED)
    n_cells, depth = _cohort(10)

    def null(r):
        return simulate_deltas(
            n_cells=n_cells, depth=depth, cp10k=0.039,
            fold_change=1.0, tau=0.2, rng=r,
        )

    percentile = rejection_rate(null, method="percentile", rng=rng, n_trials=800)
    student = rejection_rate(null, method=CALIBRATED_METHOD, rng=rng, n_trials=800)
    assert calibration_verdict(percentile) == "MISCALIBRATED"
    assert calibration_verdict(student) == "CALIBRATED"


# ---------------------------------------------------------------------------
# Heterogeneity, against a tau the generator was given
# ---------------------------------------------------------------------------


def _synthetic_deltas(tau: float, *, n_patients=40, cells=250, p=0.4, seed=1):
    rng = np.random.default_rng(seed)
    log_fc = rng.normal(0.0, tau, n_patients)
    mu = -np.log1p(-p)
    k_normal = rng.binomial(cells, p, n_patients)
    k_tumour = rng.binomial(cells, 1 - np.exp(-mu * np.exp(log_fc)), n_patients)
    return pd.DataFrame({
        "gene": "SYNTH", "n_normal": cells, "n_tumour": cells,
        "detect_normal": k_normal / cells, "detect_tumour": k_tumour / cells,
    })


@pytest.mark.parametrize("tau", [0.0, 0.2, 0.5])
def test_heterogeneity_recovers_the_tau_it_was_given(tau):
    """Known truth in, the same number out -- the only check that means much.

    A tau estimator that always returned zero would make every power figure in
    this project optimistic and nothing would raise; one that always returned a
    large number would make them all pessimistic. Both are excluded by running
    the same code over three generators whose answer is known.
    """
    out = heterogeneity_tau(_synthetic_deltas(tau), seed=SEED)
    assert out["tau"].iloc[0] == pytest.approx(tau, abs=0.08)


def test_heterogeneity_reports_the_raw_excess_so_a_negative_one_is_visible():
    """``tau`` is floored at zero; ``tau_squared_raw`` is not.

    A systematically negative excess would mean the sampling model is wrong,
    not that patients agree. Flooring without reporting would hide that.
    """
    out = heterogeneity_tau(_synthetic_deltas(0.0), seed=SEED)
    assert "tau_squared_raw" in out.columns
    assert out["tau"].iloc[0] >= 0.0


# ---------------------------------------------------------------------------
# Power
# ---------------------------------------------------------------------------


def test_power_rises_with_the_effect_and_falls_with_heterogeneity():
    n_cells, depth = _cohort(10)
    curve = power_curve(
        n_cells=n_cells, depth=depth, cp10k=0.039,
        fold_changes=(1.0, 0.5, 0.25), taus=(0.0, 0.4),
        seed=SEED, n_trials=300, cohort="test",
    )
    at = {(r.tau, r.fold_change): r.power for r in curve.itertuples()}
    assert at[(0.0, 0.25)] > at[(0.0, 0.5)] > at[(0.0, 1.0)]
    assert at[(0.4, 0.5)] <= at[(0.0, 0.5)]


def test_every_power_row_carries_its_own_false_positive_rate():
    """The structural guard: there is no way to get power without calibration."""
    n_cells, depth = _cohort(10)
    curve = power_curve(
        n_cells=n_cells, depth=depth, cp10k=0.039,
        fold_changes=(0.5,), taus=(0.0,), seed=SEED, n_trials=200, cohort="test",
    )
    assert "false_positive_rate" in curve.columns
    assert curve["false_positive_rate"].notna().all()
    assert (curve["method"] == CALIBRATED_METHOD).all()


# ---------------------------------------------------------------------------
# Boundaries
# ---------------------------------------------------------------------------


def test_a_t_interval_refuses_a_single_observation():
    with pytest.raises(CalibrationError, match="at least 2"):
        student_t_interval(np.array([0.3]), rng=np.random.default_rng(0))


def test_mismatched_cohort_vectors_are_refused():
    with pytest.raises(CalibrationError, match="per-patient vectors"):
        simulate_deltas(
            n_cells=np.array([100, 100, 100]), depth=np.array([5000.0, 5000.0]),
            cp10k=0.039, fold_change=1.0, tau=0.0, rng=np.random.default_rng(0),
        )


def test_a_negative_tau_is_refused():
    with pytest.raises(CalibrationError, match="standard deviation"):
        simulate_deltas(
            n_cells=np.array([100]), depth=np.array([5000.0]), cp10k=0.039,
            fold_change=1.0, tau=-0.1, rng=np.random.default_rng(0),
        )


def test_an_unknown_interval_method_is_named_rather_than_defaulted():
    with pytest.raises(CalibrationError, match="unknown interval method"):
        rejection_rate(
            lambda r: np.array([0.1, 0.2, 0.3]), method="jackknife",
            rng=np.random.default_rng(0), n_trials=1,
        )


def test_all_three_methods_share_one_call_shape():
    """A comparison between methods must not also be a comparison of call sites."""
    rng = np.random.default_rng(SEED)
    values = np.array([0.1, -0.2, 0.35, 0.02, -0.11, 0.4, -0.3, 0.15, 0.0, 0.22])
    for name, fn in INTERVAL_METHODS.items():
        lo, hi = fn(values, rng=rng, alpha=0.05)
        assert lo < hi, name
        assert isinstance(excludes_zero((lo, hi)), bool)


def test_bca_does_not_return_an_infinite_bound_when_every_draw_is_below():
    """``norm.ppf(0)`` is ``-inf`` and would silently become an order statistic."""
    rng = np.random.default_rng(SEED)
    lo, hi = bca_interval(np.full(10, 0.5), rng=rng)
    assert np.isfinite(lo) and np.isfinite(hi)


def test_the_percentile_interval_is_the_one_the_repository_actually_uses():
    """Pins the implementation this calibration is a statement about.

    If ``coexpression_silencing`` ever changes how it forms its interval, this
    module's measurements stop describing it, and the MLH1 design's reason for
    departing from it stops being a reason.
    """
    from src.reference.jobs.coexpression_silencing import N_BOOTSTRAP

    rng = np.random.default_rng(11)
    values = rng.normal(0.0, 1.0, 30)
    draws = np.random.default_rng(3).choice(
        values, size=(4000, values.size), replace=True
    ).mean(axis=1)
    reference = np.percentile(draws, [2.5, 97.5])
    lo, hi = percentile_interval(values, rng=np.random.default_rng(3), n_boot=4000)
    assert N_BOOTSTRAP == 10_000
    assert (lo, hi) == pytest.approx(tuple(reference), abs=1e-9)
