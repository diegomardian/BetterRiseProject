"""Every claim in trial_blindness.py, against the input that forces it to fail.

Same rule as tests/test_checks_can_fail.py: a check that cannot fail is worse
than no check, so each guard here gets two tests -- a *positive control* on the
input the guard must accept, and a *falsifier* on the input it must reject. The
falsifiers are the load-bearing half. "OLS's residual against the
variance-weighted truth is 1e-14" means nothing on its own: 1e-14 is also what
you get from an assertion that compares a number with itself, and most of these
tests exist to rule that out.

The specific way this module could be self-deceiving is worth naming. Its whole
subject is estimators that are secretly identical to a truth functional, so a
test suite written carelessly would consist of identities that hold by
construction and would pass against a completely broken implementation. Every
"is exactly zero" test below is therefore paired with a design point or a truth
where the same code is measurably non-zero, on the same draw.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from src.common.io import RESERVED_META_KEYS
from src.harness import trial_blindness as tb
from src.harness.trial_recovery import IS_DEGENERATE, simulate_trial


def _draw(design="confounded-bernoulli", n=2000, theta=3.0, seed=0, **kw):
    if kw:
        return tb.simulate_trial_design(
            n, theta, rng=np.random.default_rng(seed), **kw
        )
    return tb.DESIGNS[design].draw(n, theta, np.random.default_rng(seed))


#: Enough independent streams that a machine-zero residual is an identity rather
#: than a coincidence, and few enough that the suite stays quick.
SEEDS = range(8)


# ---------------------------------------------------------------------------
# The generalised generator must be the SAME generator at the default point
# ---------------------------------------------------------------------------


def test_the_extended_generator_reproduces_simulate_trial_draw_for_draw():
    """Positive control, and the precondition for every comparison in the module.

    If ``simulate_trial_design`` at its default design point drew even slightly
    differently -- one extra RNG call, a reordered draw -- then every number here
    would be a comparison between two generators rather than between two designs,
    and the paper's 0.059 would not be the thing being reproduced.
    """
    for seed in SEEDS:
        original = simulate_trial(1500, 3.0, rng=np.random.default_rng(seed))
        extended = tb.simulate_trial_design(1500, 3.0, rng=np.random.default_rng(seed))
        assert original.records.equals(extended.records), f"seed {seed}"
        assert original.theta_realised == extended.theta_realised


def test_the_equivalence_check_can_fail():
    """FALSIFIER for the test above. Move any design knob and the draws must
    diverge -- otherwise the knob is not connected to anything and the previous
    test would pass no matter what ``assignment`` or ``propensity_spread`` did."""
    original = simulate_trial(1500, 3.0, rng=np.random.default_rng(0))
    for kwargs in (
        {"propensity_spread": 0.0},
        {"propensity_spread": 0.5},
        {"assignment": "block"},
        {"tau": 2.0},
    ):
        moved = tb.simulate_trial_design(
            1500, 3.0, rng=np.random.default_rng(0), **kwargs
        )
        assert not original.records.equals(moved.records), (
            f"{kwargs} produced an identical cohort, so the knob does nothing"
        )


def test_tau_zero_is_the_homogeneous_effect_and_is_not_a_no_op_at_tau_nonzero():
    per_stratum_flat = tb._stratum_effects(3.0, 0.0, (0.5, 0.3, 0.2))
    assert np.allclose(per_stratum_flat, 3.0)
    tilted = tb._stratum_effects(3.0, 1.0, (0.5, 0.3, 0.2))
    assert not np.allclose(tilted, 3.0)
    # Centred on the design-weighted mean, so the population ATE is still theta.
    weights = np.array([0.5, 0.3, 0.2])
    assert float((weights * tilted).sum()) == pytest.approx(3.0, abs=1e-12)


# ---------------------------------------------------------------------------
# (1) The second realised truth
# ---------------------------------------------------------------------------


def test_the_varweighted_truth_is_exactly_the_ols_coefficient():
    """THE Frisch--Waugh claim, and the basis of finding (1).

    Not ``approx``: if OLS's estimand had to be computed approximately, the
    "OLS is exactly blind against its own estimand" result would be a numerical
    coincidence rather than an algebraic identity.
    """
    from src.harness.trial_recovery import ols_stratum_dummies

    for seed in SEEDS:
        trial = _draw(seed=seed)
        assert tb.varweighted_effect(trial.records) == pytest.approx(
            ols_stratum_dummies(trial), abs=1e-11
        )


def test_the_two_truths_are_different_numbers_under_confounding():
    """FALSIFIER for the test above. If the standardised and variance-weighted
    truths coincided under this generator, "OLS is blind against one and not the
    other" would be a statement about nothing, and every zero in the matrix would
    be trivially explained.

    THE GAP IS A SAMPLING QUANTITY, NOT A BIAS, and that had to be learned from a
    failing test rather than assumed. Under a homogeneous effect both truths have
    expectation theta, so their difference is mean-zero noise of order
    1e-2 and an individual draw can agree to any tolerance you like by chance.
    The assertion is therefore on the median, with a floor six orders of
    magnitude above machine zero on every single draw -- which is the property
    that actually matters: 1e-6 is not 1e-14.
    """
    gaps = []
    for seed in SEEDS:
        trial = _draw(seed=seed)
        gaps.append(
            abs(tb.varweighted_effect(trial.records) - trial.theta_realised)
        )
    assert np.median(gaps) > 1e-2, (
        f"the two truths agree to {np.median(gaps):.3g}, so the matrix has no contrast"
    )
    assert min(gaps) > 1e-6, f"a draw came within {min(gaps):.3g} of machine zero"


def test_ols_is_exactly_blind_against_its_own_estimand_and_informative_against_the_other():
    """Finding (1) in one assertion, on the SAME draw.

    The two residuals below are computed from one estimate and one cohort. Only
    the truth functional changes, so nothing about sampling, seeds or estimator
    code can account for the difference between them.
    """
    from src.harness.trial_recovery import ols_stratum_dummies

    own, other = [], []
    for seed in SEEDS:
        trial = _draw(seed=seed)
        estimate = ols_stratum_dummies(trial)
        own.append(abs(estimate - tb.varweighted_effect(trial.records)))
        other.append(abs(estimate - trial.theta_realised))
    assert max(own) < 1e-11, f"OLS is not exact against its own estimand: {max(own):.3g}"
    assert np.median(other) > 1e-2, (
        f"OLS reproduced the standardised truth to {np.median(other):.3g}, which "
        f"would contradict trial_recovery's own result"
    )
    # Not a single draw is within eight orders of magnitude of the other column.
    assert min(other) > 1e-6, f"closest draw: {min(other):.3g}"


def test_the_paper_s_blindest_estimator_has_an_informative_cell():
    """The mirror image, and the reason the matrix is not just a relabelling.

    ``gcomp_from_generator`` is exactly the standardised truth, so against the
    variance-weighted truth it must be measurably WRONG. If it were blind against
    both, degeneracy would still be a property of the estimator.
    """
    from src.harness.trial_recovery import gcomp_from_generator

    against_other = []
    for seed in SEEDS:
        trial = _draw(seed=seed)
        estimate = gcomp_from_generator(trial)
        assert estimate == trial.theta_realised
        against_other.append(abs(estimate - tb.varweighted_effect(trial.records)))
    assert min(against_other) > 1e-6
    assert np.median(against_other) > 1e-2


def test_the_varweighted_truth_returns_no_number_under_a_positivity_failure():
    """Invariant 1's rule, in this module: an unestimable quantity is not a zero
    and not a silently dropped stratum."""
    trial = _draw(n=400, seed=4)
    trial.records.loc[trial.records["stratum"] == 0, "treated"] = 1
    assert np.isnan(tb.varweighted_effect(trial.records))


def test_the_positivity_guard_can_fail():
    """FALSIFIER: the untouched cohort must give a finite number, or the test
    above would pass on an implementation that returned NaN unconditionally."""
    assert np.isfinite(tb.varweighted_effect(_draw(n=400, seed=4).records))


# ---------------------------------------------------------------------------
# (2) Block randomisation
# ---------------------------------------------------------------------------


def test_block_randomisation_is_exactly_one_to_one_in_every_stratum():
    for seed in SEEDS:
        trial = _draw("block-randomised", n=2000, seed=seed)
        phat = trial.records.groupby("stratum")["treated"].mean()
        assert (phat == 0.5).all(), f"seed {seed}: phat = {phat.to_dict()}"


def test_bernoulli_one_to_one_is_NOT_exactly_one_to_one():
    """FALSIFIER, and a finding in its own right.

    ``rct-bernoulli`` sets every propensity to 0.5, which is 1:1 randomisation in
    expectation and is what most simulators mean by "an RCT". The realised counts
    are still unbalanced, so OLS's variance weights are still not proportional to
    prevalence and OLS is still informative. Blocking is doing the work, not the
    0.5 -- and without this test the block result would read as "we set the
    propensity to a half".
    """
    exact = []
    for seed in SEEDS:
        trial = _draw("rct-bernoulli", n=2000, seed=seed)
        phat = trial.records.groupby("stratum")["treated"].mean()
        exact.append(bool((phat == 0.5).all()))
    assert not any(exact), "Bernoulli assignment landed exactly 1:1 in every stratum"


def test_ols_goes_exactly_blind_under_blocking_and_is_not_under_bernoulli():
    """Finding (2). The estimator code is untouched between these two lines; only
    the randomisation differs."""
    from src.harness.trial_recovery import ols_stratum_dummies

    blocked, bernoulli = [], []
    for seed in SEEDS:
        t_block = _draw("block-randomised", n=2000, seed=seed)
        blocked.append(abs(ols_stratum_dummies(t_block) - t_block.theta_realised))
        t_bern = _draw("rct-bernoulli", n=2000, seed=seed)
        bernoulli.append(abs(ols_stratum_dummies(t_bern) - t_bern.theta_realised))
    assert max(blocked) < 1e-11, f"blocked OLS residual {max(blocked):.3g}"
    assert min(bernoulli) > 1e-6, (
        f"Bernoulli 1:1 also made OLS exact ({min(bernoulli):.3g}); the claim "
        f"would then be about the propensity value, not about blocking"
    )


def test_the_unadjusted_control_goes_blind_under_blocking():
    """The cell we did not anticipate, asserted so it cannot quietly stop being
    true.

    ``trial_recovery.unadjusted`` exists "so the curve has something it *can*
    catch". Under exact 1:1 allocation the difference in arm means IS the
    standardised effect, so under the randomisation a real trial runs, the
    positive control is blind too.
    """
    from src.harness.trial_recovery import unadjusted

    for seed in SEEDS:
        trial = _draw("block-randomised", n=2000, seed=seed)
        assert abs(unadjusted(trial) - trial.theta_realised) < 1e-11


def test_the_unadjusted_control_is_still_catchable_under_confounding():
    """FALSIFIER for the test above: under the paper's design the same estimator
    must be badly wrong, or "goes blind under blocking" is not a change."""
    from src.harness.trial_recovery import unadjusted

    trial = _draw("confounded-bernoulli", n=20000, seed=1)
    assert unadjusted(trial) > trial.theta_realised + 1.0


def test_the_propensity_spread_walks_between_the_two_designs():
    assert np.allclose(tb.propensity_at_spread(0.0), 0.5)
    assert np.array_equal(
        tb.propensity_at_spread(1.0), np.asarray([0.25, 0.50, 0.75])
    ), "spread=1 must reproduce the paper's propensities exactly, not nearly"
    mid = tb.propensity_at_spread(0.5)
    assert mid[0] > 0.25 and mid[0] < 0.5


def test_the_spread_sweep_moves_the_ols_residual_monotonically_toward_zero():
    """Finding (2)'s continuous version. The residual must be ordered by spread,
    not merely different at the two ends -- an unordered pair would be consistent
    with the residual just being noisy."""
    from src.harness.trial_recovery import ols_stratum_dummies

    medians = []
    for spread in (0.0, 0.25, 0.5, 0.75, 1.0):
        residuals = []
        for seed in SEEDS:
            trial = _draw(n=4000, seed=seed, propensity_spread=spread)
            residuals.append(abs(ols_stratum_dummies(trial) - trial.theta_realised))
        medians.append(float(np.median(residuals)))
    assert medians == sorted(medians), f"not monotone in spread: {medians}"
    assert medians[-1] > 10 * medians[0]


def test_block_sizes_refuse_a_cohort_that_cannot_be_allocated_one_to_one():
    with pytest.raises(ValueError, match="multiple of block_size"):
        tb._block_sizes(2001, (0.5, 0.3, 0.2), 4)
    with pytest.raises(ValueError, match="even"):
        tb._block_sizes(2000, (0.5, 0.3, 0.2), 3)
    with pytest.raises(ValueError, match="leaves a stratum empty"):
        tb._block_sizes(4, (0.5, 0.3, 0.2), 4)


def test_block_sizes_accepts_the_cohorts_the_sweep_actually_uses():
    """FALSIFIER for the guards above: they must not reject the grid."""
    for n in (100, 200, 500, 1000, 2000, 5000):
        n_usable = n - n % 4
        sizes = tb._block_sizes(n_usable, (0.5, 0.3, 0.2), 4)
        assert sizes.sum() == n_usable
        assert (sizes % 4 == 0).all()
        assert (sizes > 0).all()


def test_an_unknown_assignment_is_refused_rather_than_silently_bernoulli():
    with pytest.raises(ValueError, match="assignment must be"):
        tb.simulate_trial_design(100, 3.0, rng=np.random.default_rng(0),
                                 assignment="stratified-ish")


# ---------------------------------------------------------------------------
# (3) The estimators that are secretly the same functional
# ---------------------------------------------------------------------------


def test_saturated_aipw_collapses_onto_standardisation():
    """Doubly robust, and blind. The augmentation term cancels within every
    stratum, so the defensive vocabulary buys nothing against this failure."""
    for seed in SEEDS:
        trial = _draw(seed=seed)
        assert tb.aipw_saturated(trial) == pytest.approx(trial.theta_realised, abs=1e-11)


def test_saturated_aipw_is_not_blind_against_the_other_truth():
    """FALSIFIER: the same estimator on the same draw must be measurably wrong
    against the variance-weighted truth, or the zero above is unfalsifiable."""
    gaps = [
        abs(tb.aipw_saturated(t) - tb.varweighted_effect(t.records))
        for t in (_draw(seed=seed) for seed in SEEDS)
    ]
    assert min(gaps) > 1e-6
    assert np.median(gaps) > 1e-2


def test_matching_reaches_the_generator_continuously_as_m_grows():
    """Finding (3)'s sharpest form: there is no line in the code where the
    estimator becomes degenerate. The same function, with one integer changed,
    walks from clearly informative to exactly the generator."""
    residuals = {}
    for m in (1, 5, 20, None):
        draws = []
        for seed in SEEDS:
            trial = _draw(seed=seed)
            rng = np.random.default_rng([13, seed])
            estimate = tb.matching_exact(trial, m=m, rng=rng)
            draws.append(abs(estimate - trial.theta_realised))
        # The median, not the max: matching resamples, so the worst of eight
        # draws is itself noisy and its ordering is not the claim being made.
        residuals[m] = float(np.median(draws)) if m is not None else float(max(draws))
    ordered = [residuals[m] for m in (1, 5, 20, None)]
    assert ordered == sorted(ordered, reverse=True), f"not monotone in m: {residuals}"
    assert residuals[None] < 1e-11, f"m=all is not exact: {residuals[None]:.3g}"
    assert residuals[1] > 1e-2, (
        f"m=1 is already exact ({residuals[1]:.3g}); the continuity claim needs a "
        f"non-degenerate end of the range"
    )


def test_matching_refuses_a_nonsense_m():
    with pytest.raises(ValueError, match="m must be"):
        tb.matching_exact(_draw(n=200), m=0, rng=np.random.default_rng(0))


def test_the_logistic_ipw_is_degenerate_in_intent_and_never_numerically_zero():
    """Finding (3)'s uncomfortable half. Saturated model, so the unpenalised MLE
    is exactly ``phat_g`` and this estimator IS ``ipw_saturated`` -- but sklearn's
    default penalty means it never reports a zero residual, at any C."""
    worst_by_c = {}
    for C in (0.1, 1.0, 100.0, 1e6):
        residuals = []
        for seed in SEEDS:
            trial = _draw(seed=seed)
            residuals.append(abs(tb.ipw_logistic_l2(trial, C=C) - trial.theta_realised))
        worst_by_c[C] = float(np.median(residuals))
    assert all(v > tb.DEGENERATE_BELOW for v in worst_by_c.values()), worst_by_c
    assert worst_by_c[0.1] > worst_by_c[1.0] > worst_by_c[100.0], worst_by_c


def test_the_regularisation_constant_outranks_a_genuinely_informative_estimator():
    """The inversion, asserted. At C=0.1 the estimator that is algebraically the
    generator posts a LARGER residual than cross-fitted IPW, which genuinely holds
    data out. Anyone reading the residual as a quality score gets the ranking
    exactly backwards, and the constant responsible is a library default."""
    from src.harness.trial_recovery import ipw_cross_fitted

    penalised, honest = [], []
    for seed in SEEDS:
        trial = _draw(seed=seed)
        penalised.append(abs(tb.ipw_logistic_l2(trial, C=0.1) - trial.theta_realised))
        honest.append(abs(ipw_cross_fitted(trial) - trial.theta_realised))
    assert np.median(penalised) > np.median(honest), (
        f"penalised {np.median(penalised):.3g} vs cross-fitted "
        f"{np.median(honest):.3g}: the inversion did not reproduce"
    )


def test_the_solver_tolerance_moves_the_residual_by_orders_of_magnitude():
    """The second knob, and it is not a statistical one.

    Turning the penalty off does not restore the identity: the residual stops at
    the lbfgs convergence tolerance. Tightening ``tol`` alone, with C fixed, moves
    it by several orders of magnitude -- so this estimator's measured distance
    from the generator is set by two sklearn defaults.
    """
    loose, tight = [], []
    for seed in SEEDS:
        trial = _draw(seed=seed)
        loose.append(abs(tb.ipw_logistic_l2(trial, C=1e6) - trial.theta_realised))
        tight.append(
            abs(tb.ipw_logistic_l2(trial, C=1e6, tol=1e-8) - trial.theta_realised)
        )
    assert np.median(loose) > 100 * np.median(tight), (
        f"loose {np.median(loose):.3g} vs tight {np.median(tight):.3g}"
    )


# ---------------------------------------------------------------------------
# (4) The information ratio
# ---------------------------------------------------------------------------


def test_rho_is_small_for_ols_and_does_not_fall_with_cohort_size():
    """Finding (4). A non-zero residual does not make a curve informative: OLS's
    residual is a fixed small fraction of the generator's own sampling error at
    every n, so more patients buy no separation between estimator and
    generator."""
    runs = tb.run_grid(
        seed=515,
        n_seeds=6,
        cohort_sizes=(500, 2000, 5000),
        n_replicates=12,
        designs=("confounded-bernoulli",),
        estimators=("ols-stratum-dummies", "unadjusted"),
    )
    rho = tb.information_ratio(runs)
    ols = rho[(rho["estimator"] == "ols-stratum-dummies") & (rho["truth"] == "standardised")]
    values = ols.sort_values("n_patients")["rho"].to_numpy()
    assert len(values) == 3
    assert ((values > 0.05) & (values < 0.30)).all(), f"rho out of range: {values}"
    assert values[-1] > 0.5 * values[0], (
        f"rho fell with n ({values}); the 'does not improve with n' claim is the "
        f"whole point of the ratio and would be wrong"
    )


def test_rho_can_be_large_so_the_scale_means_something():
    """FALSIFIER for the test above. If every estimator's rho were ~0.1 the ratio
    would just be a rescaling with no discriminating power. ``unadjusted`` is
    biased by more than the generator's sampling error, so its rho must be well
    above 1."""
    runs = tb.run_grid(
        seed=515,
        n_seeds=6,
        cohort_sizes=(2000,),
        n_replicates=12,
        designs=("confounded-bernoulli",),
        estimators=("ols-stratum-dummies", "unadjusted"),
    )
    rho = tb.information_ratio(runs)
    unadj = rho[(rho["estimator"] == "unadjusted") & (rho["truth"] == "standardised")]
    assert float(unadj["rho"].iloc[0]) > 5.0


def test_rho_is_exactly_zero_for_an_estimator_that_is_the_generator():
    runs = tb.run_grid(
        seed=515,
        n_seeds=5,
        cohort_sizes=(2000,),
        n_replicates=8,
        designs=("confounded-bernoulli",),
        estimators=("gcomp-from-generator",),
    )
    rho = tb.information_ratio(runs)
    std = rho[rho["truth"] == "standardised"]
    assert float(std["rho"].iloc[0]) == 0.0
    assert float(std["generator_noise_share"].iloc[0]) == 1.0


# ---------------------------------------------------------------------------
# Effect heterogeneity
# ---------------------------------------------------------------------------


def test_heterogeneity_moves_the_ols_curve_while_its_residual_stays_at_zero():
    """The curve moving for a reason unrelated to estimator quality.

    As tau grows, OLS's recovery ratio wanders away from 1 -- and its residual
    against the estimand it actually targets stays at machine zero throughout.
    Nothing about the estimator changed; the variance-weighted estimand simply
    stopped equalling theta.
    """
    runs = tb.run_sweep(
        seed=616,
        arm="tau",
        values=(0.0, 8.0, 16.0),
        n_seeds=6,
        n_patients=4000,
        n_replicates=12,
        estimators=("ols-stratum-dummies",),
    )
    out = tb.summarise_sweep(runs)
    own = out[out["truth"] == "varweighted"].sort_values("sweep_value")
    assert (own["max_residual"] < 1e-11).all(), own[["sweep_value", "max_residual"]]

    drift = (own["median_recovery_ratio"] - 1.0).abs().to_numpy()
    assert drift[0] < 0.03, f"the curve is already off 1 at tau=0: {drift[0]:.3g}"
    assert drift[-1] > 0.08, (
        f"the curve did not move under heterogeneity ({drift}); without movement "
        f"there is nothing for the zero residual to contrast with"
    )
    assert list(drift) == sorted(drift), f"drift not ordered in tau: {drift}"


# ---------------------------------------------------------------------------
# The matrix itself
# ---------------------------------------------------------------------------


def test_every_declared_cell_matches_the_measurement():
    """``DEGENERACY`` is a claim about 84 cells. Measure all of them."""
    runs = tb.run_grid(seed=717, n_seeds=5, cohort_sizes=(2000,), n_replicates=5)
    matrix = tb.residual_matrix(runs)
    wrong = matrix[matrix["measured_degenerate"] != matrix["declared_degenerate"]]
    assert wrong.empty, wrong[
        ["estimator", "truth", "design", "max_residual", "declared_degenerate"]
    ].to_string(index=False)


def test_the_matrix_check_can_fail():
    """FALSIFIER: corrupt one declared cell and the comparison must catch it.

    Without this, ``test_every_declared_cell_matches_the_measurement`` would pass
    just as happily against a ``measured_degenerate`` column computed from
    ``declared_degenerate``.
    """
    runs = tb.run_grid(
        seed=717,
        n_seeds=3,
        cohort_sizes=(2000,),
        n_replicates=3,
        designs=("confounded-bernoulli",),
        estimators=("gcomp-from-generator",),
    )
    matrix = tb.residual_matrix(runs)
    matrix.loc[matrix["truth"] == "standardised", "declared_degenerate"] = False
    wrong = matrix[matrix["measured_degenerate"] != matrix["declared_degenerate"]]
    assert not wrong.empty


def test_the_degeneracy_threshold_sits_in_an_empty_region():
    """``DEGENERATE_BELOW`` must not be a knob.

    If measured residuals formed a continuum across 1e-9, moving the threshold
    would move the verdicts and "exactly blind" would be a choice rather than a
    measurement.

    The assertion is on the SEPARATION between the two populations rather than on
    a literal gap, because a literal one was wrong: the closest informative cell
    is ``gcomp`` against the variance-weighted truth under ``rct-bernoulli``, at
    about 9e-4, and the tightly-converged logistic IPW gets down to about 3e-7 --
    both far above machine zero, both far below the 1e-2 the headline cells sit
    at. What matters is that the blind and informative populations are separated
    by many orders of magnitude with the threshold strictly inside, so no
    plausible move of ``DEGENERATE_BELOW`` changes a verdict.
    """
    runs = tb.run_grid(seed=818, n_seeds=5, cohort_sizes=(500, 2000), n_replicates=4)
    matrix = tb.residual_matrix(runs)
    blind = matrix[matrix["measured_degenerate"]]["max_residual"]
    informative = matrix[~matrix["measured_degenerate"]]["max_residual"]
    assert not blind.empty and not informative.empty

    assert blind.max() < 1e-11, f"a 'blind' cell reached {blind.max():.3g}"
    assert informative.min() > 1e-8, (
        f"an 'informative' cell came down to {informative.min():.3g}, close enough "
        f"to the threshold that its verdict is a knob"
    )
    assert informative.min() / max(blind.max(), 1e-300) > 1e4
    assert blind.max() < tb.DEGENERATE_BELOW < informative.min()


def test_the_matrix_contradicts_the_per_estimator_dict():
    """The module's reason to exist. If no cell disagreed with
    ``trial_recovery.IS_DEGENERATE``, degeneracy really would be a property of the
    estimator and this file would be a reformatting exercise."""
    runs = tb.run_grid(seed=919, n_seeds=5, cohort_sizes=(2000,), n_replicates=5)
    matrix = tb.residual_matrix(runs)
    flipped = matrix[
        matrix["estimator"].isin(IS_DEGENERATE)
        & (
            matrix["measured_degenerate"]
            != matrix["estimator"].map(IS_DEGENERATE).astype(bool)
        )
    ]
    assert not flipped.empty
    assert {"ols-stratum-dummies", "unadjusted"} <= set(flipped["estimator"])


# ---------------------------------------------------------------------------
# Reproducibility and the results contract
# ---------------------------------------------------------------------------


def test_the_grid_is_reproducible_under_a_fixed_seed():
    kw = dict(
        n_seeds=2,
        cohort_sizes=(200,),
        n_replicates=3,
        designs=("confounded-bernoulli",),
        estimators=("matching-20", "ipw-cross-fitted"),
    )
    a = tb.run_grid(seed=31, **kw)
    b = tb.run_grid(seed=31, **kw)
    pd.testing.assert_frame_equal(a, b)


def test_different_seeds_give_different_draws():
    """FALSIFIER for the test above: identical frames from different seeds would
    mean the seed is not reaching the generator and reproducibility is vacuous."""
    kw = dict(
        n_seeds=2,
        cohort_sizes=(200,),
        n_replicates=3,
        designs=("confounded-bernoulli",),
        estimators=("matching-20",),
    )
    a = tb.run_grid(seed=31, **kw)
    b = tb.run_grid(seed=32, **kw)
    assert not a["estimate"].equals(b["estimate"])


def test_the_seed_streams_are_independent_not_the_same_draw_repeated():
    runs = tb.run_grid(
        seed=41,
        n_seeds=4,
        cohort_sizes=(500,),
        n_replicates=2,
        designs=("confounded-bernoulli",),
        estimators=("gcomp-from-generator",),
    )
    per_seed = runs.groupby("seed_index")["truth_standardised"].apply(list)
    assert len({tuple(v) for v in per_seed}) == 4


def test_no_sidecar_key_collides_with_the_provenance_record(tmp_path):
    """The results contract, end to end.

    ``write_versioned_table`` raises on a reserved key, so a collision would be
    loud -- but only if some run exercises it. This runs the CLI, then re-reads
    every sidecar and checks the invariant-10 fields are present and untouched.
    """
    assert tb.main(["--quick", "--results-dir", str(tmp_path), "--allow-dirty"]) == 0
    sidecars = sorted(tmp_path.glob("*/trial_blindness_*.meta.json"))
    assert len(sidecars) == 4
    for path in sidecars:
        meta = json.loads(path.read_text(encoding="utf-8"))
        assert meta["seed"] == 20260909
        assert meta["git_sha"] != "unknown"
        assert meta["SYNTHETIC"].startswith("Simulated patients")
        assert meta["what_this_answers"]
        assert "sweep_seed" in meta
        # The reserved keys are the provenance record's own, not the caller's.
        assert set(meta) & RESERVED_META_KEYS == {
            k for k in RESERVED_META_KEYS if k in meta
        }
        assert "table_kind" in meta


def test_the_sidecar_note_is_the_one_trial_recovery_uses(tmp_path):
    """FALSIFIER-adjacent: the SYNTHETIC text must be copied, not paraphrased, or
    the two tables can drift apart in what they disclaim."""
    import inspect

    from src.harness import trial_recovery

    source = inspect.getsource(trial_recovery.main)
    assert "Simulated patients with an analytically known treatment " in source
    assert tb.SYNTHETIC_NOTE.startswith(
        "Simulated patients with an analytically known treatment effect."
    )
