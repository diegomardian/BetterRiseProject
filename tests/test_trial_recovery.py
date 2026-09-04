"""The blind spot, asserted rather than described.

Every claim §2 makes about the trial simulator is a property of these functions,
so each one is tested here against an input that could falsify it. The central
claim -- that an estimator sharing sufficient statistics with the generator
returns the realised effect *identically*, not approximately -- is worth nothing
as prose and everything as an assertion.
"""

from __future__ import annotations

import numpy as np
import pytest

from src.harness.trial_recovery import (
    ESTIMATORS,
    IS_DEGENERATE,
    gcomp_from_generator,
    ipw_cross_fitted,
    ipw_saturated,
    ols_stratum_dummies,
    run,
    simulate_trial,
    summarise,
    unadjusted,
)


def _trial(n=2000, theta=3.0, seed=0):
    return simulate_trial(n, theta, rng=np.random.default_rng(seed))


# ---------------------------------------------------------------------------
# The degenerate estimators
# ---------------------------------------------------------------------------


def test_gcomp_returns_the_realised_effect_exactly_not_approximately():
    """THE claim. Not ``approx`` — exactly, or the argument is weaker than
    stated and the recovery curve retains some information."""
    for seed in range(20):
        trial = _trial(seed=seed)
        assert gcomp_from_generator(trial) == trial.theta_realised


def test_saturated_ipw_is_the_same_functional_wearing_different_code():
    """The finding we did not expect. Different literature, different code, same
    sufficient statistics — so it is just as blind, and nothing about reading it
    would tell you."""
    for seed in range(20):
        trial = _trial(seed=seed)
        assert ipw_saturated(trial) == pytest.approx(trial.theta_realised, abs=1e-9)


def test_the_two_degenerate_estimators_agree_with_each_other_to_machine_precision():
    for seed in range(20):
        trial = _trial(seed=seed)
        assert gcomp_from_generator(trial) == pytest.approx(
            ipw_saturated(trial), abs=1e-9
        )


# ---------------------------------------------------------------------------
# The non-degenerate ones. These must NOT reproduce the realised effect, or the
# demonstration has no contrast in it.
# ---------------------------------------------------------------------------


def test_cross_fitting_breaks_the_degeneracy():
    """Weights applied to a patient come from data that patient is not in, which
    makes it a different functional. If this ever starts returning the realised
    effect, the arm has stopped demonstrating anything."""
    residuals = [
        abs(ipw_cross_fitted(t) - t.theta_realised)
        for t in (_trial(seed=s) for s in range(20))
    ]
    assert min(residuals) > 1e-6, "cross-fitted IPW must not reproduce realised truth"


def test_the_unadjusted_estimator_does_not_reproduce_realised_truth_either():
    residuals = [
        abs(unadjusted(t) - t.theta_realised)
        for t in (_trial(seed=s) for s in range(20))
    ]
    assert min(residuals) > 1e-6


def test_the_unadjusted_estimator_is_biased_because_assignment_depends_on_stratum():
    """The curve has to be able to catch something, or its blindness is not
    interesting. This is the arm it catches."""
    trial = _trial(n=20000, seed=1)
    assert unadjusted(trial) > trial.theta_realised + 1.0


def test_cross_fitted_ipw_is_consistent_so_its_curve_still_sits_near_one():
    """Non-degenerate must not mean wrong — otherwise the contrast is between a
    blind curve and a broken estimator rather than a blind curve and a working
    one."""
    trial = _trial(n=50000, seed=2)
    assert ipw_cross_fitted(trial) == pytest.approx(trial.theta_requested, abs=0.6)


def test_ols_is_non_degenerate_without_splitting_the_sample():
    """The specificity claim: the residual check detects a shared *functional*,
    not a shared sample.

    ``ipw_cross_fitted`` breaks the degeneracy by holding data out, which leaves
    the reading that the check merely notices sample splitting. This estimator
    sees every record and splits nothing, and is still non-degenerate — because
    OLS weights strata by within-stratum treatment variance where
    standardisation weights them by prevalence, and the propensities differ.
    """
    for seed in range(10):
        trial = _trial(n=2000, seed=seed)
        residual = abs(ols_stratum_dummies(trial) - trial.theta_realised)
        assert residual > 1e-6, (
            f"seed {seed}: OLS reproduced the realised effect to {residual:.3g}, "
            f"which would make it degenerate and this test's premise wrong"
        )


def test_ols_is_consistent_so_its_curve_also_sits_near_one():
    """Non-degenerate must not mean wrong here either. If this estimator were
    biased it would be a second `unadjusted`, and the contrast would collapse
    back to blind-versus-broken."""
    trial = _trial(n=50000, seed=3)
    assert ols_stratum_dummies(trial) == pytest.approx(trial.theta_requested, abs=0.3)


def test_the_conventional_estimator_is_the_non_degenerate_one():
    """The inversion worth stating: the two estimators that look like careful
    causal inference are blind, and the one everybody actually writes is not."""
    trial = _trial(n=5000, seed=4)
    blind = abs(gcomp_from_generator(trial) - trial.theta_realised)
    conventional = abs(ols_stratum_dummies(trial) - trial.theta_realised)
    assert blind == 0.0
    assert conventional > 1e-6


# ---------------------------------------------------------------------------
# The declared table has to match what the estimators actually do
# ---------------------------------------------------------------------------


def test_the_degeneracy_table_is_measured_not_asserted():
    """``IS_DEGENERATE`` is a claim about each estimator. Check every entry
    against the residual rather than trusting the dict."""
    trials = [_trial(seed=s) for s in range(15)]
    for name, estimator in ESTIMATORS.items():
        worst = max(abs(estimator(t) - t.theta_realised) for t in trials)
        degenerate = worst < 1e-9
        assert degenerate == IS_DEGENERATE[name], (
            f"{name}: measured degenerate={degenerate}, table says "
            f"{IS_DEGENERATE[name]} (worst residual {worst:.3g})"
        )


# ---------------------------------------------------------------------------
# The generator
# ---------------------------------------------------------------------------


def test_the_realised_effect_differs_from_the_requested_one_by_sampling():
    """If these were equal the whole distinction would be vacuous."""
    trial = _trial(n=200, seed=3)
    assert trial.theta_realised != trial.theta_requested
    assert abs(trial.theta_realised - trial.theta_requested) < 3.0


def test_a_stratum_missing_an_arm_gives_no_number(monkeypatch):
    """Positivity fails and the estimand has no referent, so it returns NaN
    rather than a number — the same rule the rest of this project follows."""
    trial = _trial(n=40, seed=4)
    trial.records.loc[trial.records["stratum"] == 0, "treated"] = 1
    assert np.isnan(gcomp_from_generator(trial))
    assert np.isnan(ipw_saturated(trial))


def test_the_sweep_is_reproducible_under_a_fixed_seed():
    a = run(seed=7, cohort_sizes=(200,), n_replicates=5)
    b = run(seed=7, cohort_sizes=(200,), n_replicates=5)
    assert a.equals(b)


def test_different_seeds_give_different_draws():
    a = run(seed=7, cohort_sizes=(200,), n_replicates=5)
    b = run(seed=8, cohort_sizes=(200,), n_replicates=5)
    assert not a["estimate"].equals(b["estimate"])


def test_the_summary_carries_the_residual_so_the_table_shows_the_blindness():
    out = summarise(run(seed=9, cohort_sizes=(500,), n_replicates=20))
    degenerate = out[out["shares_sufficient_statistics"]]
    other = out[~out["shares_sufficient_statistics"]]
    assert (degenerate["max_residual_vs_realised"] < 1e-9).all()
    assert (other["max_residual_vs_realised"] > 1e-6).all()


def test_the_blind_curves_are_indistinguishable_from_each_other():
    """The reason the residual column has to exist: on the curve itself, the two
    degenerate estimators are the same line."""
    out = summarise(run(seed=11, cohort_sizes=(500, 2000), n_replicates=40))
    a = out[out["estimator"] == "gcomp-from-generator"]["ratio_median"].to_numpy()
    b = out[out["estimator"] == "ipw-saturated"]["ratio_median"].to_numpy()
    assert np.allclose(a, b, atol=1e-9)
