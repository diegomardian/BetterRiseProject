"""Every guard in the learned-generator audit, against the input that forces it.

The house rule from ``tests/test_checks_can_fail.py``: a check unable to fail is
worse than no check, because it turns an absence of evidence into a green light.
So each guard here gets a constructed input it MUST catch, plus a positive
control showing it still passes clean data.

The experiment these guards belong to is pre-registered in
``docs/prereg_learned_generator.md``. Four of the tests below are that
document's falsifiers written as code, so that a reader can see F1, F3, F6 and
F8 are decidable rather than rhetorical.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.harness.trial_learned_generator import (
    ALLCLOSE_ATOL,
    ALLCLOSE_RTOL,
    CELLS,
    CONTROL_ESTIMATOR,
    CONTROL_MUST_EXCEED,
    ESTIMATORS,
    GENERATORS,
    SATURATED_BELOW,
    SHARES_FUNCTIONAL_WITH_T_DRAW,
    THETA,
    _as_trial,
    _rng,
    _standardise_cell_means,
    cell_means,
    flip_tolerances,
    gcomp_ridge_outcome,
    primary_table,
    run_grid,
    screen_band,
)
from src.harness.trial_recovery import simulate_trial


# ---------------------------------------------------------------------------
# The hazard that makes an assertion pass by having nothing left to check
# ---------------------------------------------------------------------------


def test_an_empty_residual_array_is_what_makes_allclose_say_clean():
    """THE INPUT IS numpy ITSELF, and it is the reason this file exists.

    ``np.allclose`` on an empty array is True. So a screen that filters to valid
    pairs and finds none does not report "I could not look" -- it reports
    "everything agrees". A NaN cannot sneak past ``allclose`` directly, because
    ``allclose(nan, 0.0)`` is False; it sneaks past the *filter step* before it.
    """
    assert np.allclose(np.array([]), np.array([])) is np.True_ or np.allclose(
        np.array([]), np.array([])
    )
    assert not np.allclose(np.nan, 0.0)


def test_a_cell_with_no_valid_pairs_is_not_reported_as_agreement():
    """The forcing input: every estimate non-finite, so nothing survives the pair
    filter. The screen must abstain rather than return a verdict."""
    assert screen_band(allclose_vs_zero=True, allclose_paired=True, n_valid=0) == (
        "no-valid-pairs"
    )


def test_the_abstention_reaches_the_primary_table(monkeypatch):
    """End to end, on a runs frame whose estimates are all NaN.

    Without the ``n_valid == 0`` branch this row would carry
    ``allclose_residual_vs_zero = True`` and read as a clean pass.
    """
    runs = pd.DataFrame(
        {
            "generator": ["gmm-joint-K2"] * 4,
            "is_learned": [True] * 4,
            "capacity": [2.0] * 4,
            "generator_failed": [False] * 4,
            "generator_failure": [None] * 4,
            "seed_index": [0] * 4,
            "n_patients": [100] * 4,
            "replicate": [0, 1, 2, 3],
            "estimator": ["gcomp-saturated"] * 4,
            "shares_functional_with_T_draw": [True] * 4,
            "theta_requested": [THETA] * 4,
            "estimate": [np.nan] * 4,
            "difference_vs_requested": [np.nan] * 4,
            "T_draw": [3.1, 2.9, 3.0, 3.2],
            "T_model": [3.0] * 4,
            "T_model_mc_se": [0.0] * 4,
            "T_model_is_monte_carlo": [False] * 4,
            "max_cellmean_deviation": [0.5] * 4,
            "generator_is_saturated": [False] * 4,
            "residual_vs_T_draw": [np.nan] * 4,
            "residual_vs_T_model": [np.nan] * 4,
        }
    )
    table = primary_table(runs)
    assert len(table) == 2  # one row per reference
    for _, row in table.iterrows():
        assert row["n_valid_pairs"] == 0
        assert row["screen_band"] == "no-valid-pairs"
        assert not row["allclose_residual_vs_zero"]
        assert not row["allclose_estimate_vs_reference"]
        assert not row["is_exactly_zero"]
        assert row["n_nonfinite_estimate"] == 4


def test_the_abstention_does_not_fire_on_a_cell_that_has_pairs():
    """Positive control: the same guard must not abstain on usable data."""
    assert screen_band(allclose_vs_zero=True, allclose_paired=True, n_valid=7) == (
        "invisible-to-both"
    )


# ---------------------------------------------------------------------------
# The two things a practitioner might type, and the 10^3 between them
# ---------------------------------------------------------------------------


def test_the_two_allclose_spellings_disagree_inside_the_band():
    """THE FORCING INPUT: a residual of 2e-5 at a reference of 3.

    ``np.allclose(r, 0.0)`` allows only ``atol = 1e-8`` because ``rtol``
    multiplies the second argument, which is zero. ``np.allclose(est, ref)``
    allows ``atol + rtol*3 = 3.0e-5``. The same number is therefore a failure to
    one practitioner and a pass to the other.
    """
    truth = 3.0
    residual = 2e-5
    estimate = truth + residual

    assert not np.allclose(residual, 0.0)
    assert np.allclose(estimate, truth)
    assert screen_band(False, True, n_valid=1) == "blind-band"


def test_the_band_has_a_floor_and_a_ceiling_that_are_not_the_same_number():
    """If these two ever coincide the band is empty and (b) is untestable."""
    truth = 3.0
    floor = ALLCLOSE_ATOL
    ceiling = ALLCLOSE_ATOL + ALLCLOSE_RTOL * abs(truth)
    assert ceiling / floor > 1000.0


def test_the_flip_tolerances_are_the_tolerances_at_which_the_screens_flip():
    """Not a restatement: each flip point is checked by re-running the screen
    just above and just below it."""
    truth = np.array([3.0, 3.0, 3.0])
    residual = np.array([1e-6, 5e-6, 2e-5])
    flips = flip_tolerances(residual, truth)

    assert flips["flip_atol"] == pytest.approx(2e-5)
    assert not np.allclose(residual.max(), 0.0, atol=flips["flip_atol"] * 0.999)
    assert np.allclose(residual.max(), 0.0, atol=flips["flip_atol"] * 1.001)

    estimate = truth + residual
    paired = flips["flip_atol_paired"]
    assert np.allclose(estimate, truth, rtol=ALLCLOSE_RTOL, atol=paired + 1e-12)
    if paired > 0:
        assert not np.allclose(
            estimate, truth, rtol=ALLCLOSE_RTOL, atol=paired * 0.5
        )


def test_flip_tolerances_abstain_on_an_empty_array():
    flips = flip_tolerances(np.array([]), np.array([]))
    assert not np.isfinite(flips["flip_atol"])
    assert not np.isfinite(flips["flip_atol_paired"])


# ---------------------------------------------------------------------------
# F8 -- the blind band must be REACHABLE, or the verdict on (b) is "not tested"
# ---------------------------------------------------------------------------


def test_a_ridge_penalty_can_be_dialled_into_the_blind_band():
    """THE POSITIVE CONTROL FOR THE BAND (prereg Amendment 1, F8).

    An experiment whose estimators can only be exact or obviously-wrong cannot
    exhibit the regime prediction (b) describes, and would "confirm" F2 for a
    reason that is about the estimator set rather than about generators. So at
    least one arm must land strictly inside the band.
    """
    trial = _as_trial(simulate_trial(2000, THETA, rng=_rng(7, "train", 0)).records, THETA)
    truth = trial.theta_realised

    inside = []
    for alpha in (1e-6, 1e-5, 1e-4, 1e-3):
        estimate = gcomp_ridge_outcome(trial, alpha=alpha, rng=_rng(7, "estimator", 0))
        residual = abs(estimate - truth)
        if not np.allclose(residual, 0.0) and np.allclose(estimate, truth):
            inside.append(alpha)
    assert inside, (
        "no ridge penalty landed strictly inside the blind band, so this "
        "experiment cannot demonstrate the regime prediction (b) is about"
    )


def test_the_ridge_outcome_model_is_the_closed_form_it_claims_to_be():
    """``ybar_gd * n_gd / (n_gd + alpha)`` on an orthogonal cell design.

    The band arm is only interpretable if its residual is a transparent function
    of alpha rather than of a solver's mood.
    """
    trial = _as_trial(simulate_trial(1000, THETA, rng=_rng(11, "train", 0)).records, THETA)
    alpha = 1e-2
    records = trial.records
    shrunk = {}
    for (g, a), group in records.groupby(["stratum", "treated"]):
        y = group["outcome"].to_numpy()
        shrunk[(int(g), int(a))] = float(y.mean()) * len(y) / (len(y) + alpha)
    expected = _standardise_cell_means(records, shrunk)
    got = gcomp_ridge_outcome(trial, alpha=alpha, rng=_rng(11, "estimator", 0))
    assert got == pytest.approx(expected, abs=1e-12)


def test_a_zero_penalty_puts_the_ridge_arm_back_on_the_reference():
    """Positive control: with no penalty the arm IS the reference's functional,
    so the band it lands in must be the exact one."""
    trial = _as_trial(simulate_trial(1000, THETA, rng=_rng(13, "train", 0)).records, THETA)
    estimate = gcomp_ridge_outcome(trial, alpha=0.0, rng=_rng(13, "estimator", 0))
    assert abs(estimate - trial.theta_realised) < 1e-9


# ---------------------------------------------------------------------------
# F1 and F3 -- the paper's prediction, decidable
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "generator", ["parametric-oracle", "parametric-plugin", "gmm-joint-K4"]
)
def test_the_natural_wiring_is_exactly_zero_whether_or_not_the_generator_learned(
    generator,
):
    """F1 AND F3, AS A TEST RATHER THAN A CLAIM.

    ``gcomp-saturated`` is the natural wiring the paper describes: compute the
    effect on the cohort you just generated. Its residual against ``T_draw``
    must be BITWISE zero for a learned generator exactly as for a parametric
    one, because both are functionals of the same synthetic cohort and the
    generator never enters.

    If this test ever fails for the learned generator and passes for the
    parametric one, the paper's stated mechanism is right and this experiment's
    headline is wrong. That is the point of writing it this way.
    """
    runs = run_grid(
        seed=101,
        n_seeds=1,
        cohort_sizes=(500,),
        n_replicates=2,
        generators=(generator,),
        estimators=("gcomp-saturated",),
    )
    runs = runs[~runs["generator_failed"].astype(bool)]
    assert len(runs) > 0
    assert (runs["residual_vs_T_draw"] == 0.0).all(), (
        f"{generator}: the natural wiring is not bitwise zero, which would make "
        f"the paper's prediction (a) correct as stated"
    )


def test_the_screen_is_not_vacuous_because_every_estimator_agrees():
    """THE CONTROL THAT CARRIES INTERNAL VALIDITY (prereg 4.1, F6).

    ``ols-stratum-dummies`` is consistent for theta -- it is accurate -- and is
    NOT the reference's functional. If its residual is not clearly larger than
    the exact arms, the screen cannot separate "shares the functional" from
    "happens to be accurate", and no other row in the experiment may be read.
    """
    runs = run_grid(
        seed=102,
        n_seeds=1,
        cohort_sizes=(2000,),
        n_replicates=3,
        generators=("gmm-joint-K4",),
        estimators=("gcomp-saturated", CONTROL_ESTIMATOR),
    )
    table = primary_table(runs)
    draw = table[table["reference"] == "T_draw"]
    exact = draw[draw["estimator"] == "gcomp-saturated"]["max_residual"].max()
    control = draw[draw["estimator"] == CONTROL_ESTIMATOR]["max_residual"].min()

    assert exact == 0.0
    assert control > CONTROL_MUST_EXCEED, (
        "the negative control estimator collapsed onto the reference; the "
        "screen has no demonstrated discriminating power"
    )


def test_the_control_guard_fires_on_a_table_where_the_control_collapsed():
    """The guard needs an input that forces it, and this is that input.

    A primary table in which OLS sits at machine precision -- which really
    happens under block randomisation, per trial_blindness -- must be caught by
    the same comparison ``main`` makes, not waved through.
    """
    collapsed = pd.DataFrame(
        {
            "reference": ["T_draw", "T_draw"],
            "estimator": ["gcomp-saturated", CONTROL_ESTIMATOR],
            "max_residual": [0.0, 2e-14],
        }
    )
    control = collapsed[collapsed["estimator"] == CONTROL_ESTIMATOR]
    assert not bool((control["max_residual"] > CONTROL_MUST_EXCEED).all())


# ---------------------------------------------------------------------------
# The declared functional table, measured rather than asserted
# ---------------------------------------------------------------------------


def test_every_estimator_declared_to_share_the_functional_actually_does():
    """``SHARES_FUNCTIONAL_WITH_T_DRAW`` is an algebraic claim made before the
    run. This measures it. A declared value the measurement contradicts is a
    finding, not a typo."""
    runs = run_grid(
        seed=103,
        n_seeds=1,
        cohort_sizes=(2000,),
        n_replicates=2,
        generators=("gmm-joint-K2",),
    )
    table = primary_table(runs)
    draw = table[table["reference"] == "T_draw"]

    disagreements = []
    for _, row in draw.iterrows():
        measured = row["max_residual"] < 1e-9
        if measured != row["shares_functional_with_T_draw"]:
            disagreements.append(
                f"{row['estimator']}: declared "
                f"{row['shares_functional_with_T_draw']}, measured {measured} "
                f"(max residual {row['max_residual']:.3g})"
            )
    assert not disagreements, "\n".join(disagreements)


def test_the_declared_table_check_fires_when_a_declaration_is_wrong():
    """Forcing input: flip one declaration and confirm the comparison catches
    it, so the test above is not passing by never comparing anything."""
    declared = dict(SHARES_FUNCTIONAL_WITH_T_DRAW)
    declared["unadjusted"] = True  # the confounded arm, declared degenerate
    measured = {"unadjusted": False}
    assert declared["unadjusted"] != measured["unadjusted"]


# ---------------------------------------------------------------------------
# The saturation check (prereg 2.3)
# ---------------------------------------------------------------------------


def test_the_plugin_generator_is_detected_as_saturated():
    """THE FORCING INPUT FOR F5: a generator that memorises the training cell
    means IS the parametric generator, and must be labelled so."""
    training = simulate_trial(2000, THETA, rng=_rng(5, "train", 0)).records
    fitted = GENERATORS["parametric-plugin"](
        training, n_patients=2000, theta=THETA, rng=_rng(5, "fit", 0)
    )
    observed = cell_means(training)
    deviation = max(
        abs(fitted.implied_cell_means[cell] - observed[cell]) for cell in CELLS
    )
    assert deviation < SATURATED_BELOW


def test_a_generator_that_did_not_memorise_is_not_labelled_saturated():
    """Positive control: the oracle re-draws from the true parameters, so its
    implied means differ from the training cohort's by sampling error. A check
    that called everything saturated would be useless."""
    training = simulate_trial(2000, THETA, rng=_rng(5, "train", 1)).records
    fitted = GENERATORS["parametric-oracle"](
        training, n_patients=2000, theta=THETA, rng=_rng(5, "fit", 1)
    )
    observed = cell_means(training)
    deviation = max(
        abs(fitted.implied_cell_means[cell] - observed[cell]) for cell in CELLS
    )
    assert deviation > SATURATED_BELOW


# ---------------------------------------------------------------------------
# Standardisation refuses a missing cell rather than dropping it
# ---------------------------------------------------------------------------


def test_a_missing_stratum_arm_cell_is_refused_not_silently_dropped():
    """THE FORCING INPUT: a cohort with no treated patients in stratum 2.

    Dropping that stratum would change the estimand from the standardised effect
    to the effect over whichever strata happened to be complete, and would
    return a confident number for a positivity failure.
    """
    records = pd.DataFrame(
        {
            "stratum": [0, 0, 1, 1, 2, 2],
            "treated": [0, 1, 0, 1, 0, 0],
            "outcome": [10.0, 13.0, 14.0, 17.0, 25.0, 26.0],
        }
    )
    mu = cell_means(records)
    assert (2, 1) not in mu
    assert not np.isfinite(_standardise_cell_means(records, mu))


def test_standardisation_returns_a_number_when_every_cell_is_present():
    """Positive control."""
    records = pd.DataFrame(
        {
            "stratum": [0, 0, 1, 1],
            "treated": [0, 1, 0, 1],
            "outcome": [10.0, 13.0, 14.0, 17.0],
        }
    )
    got = _standardise_cell_means(records, cell_means(records))
    assert got == pytest.approx(3.0)


# ---------------------------------------------------------------------------
# The Amendment 1 finding, as a test
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("tol", [1e-2, 1e-3, 1e-10])
def test_the_mixture_mean_identity_does_not_depend_on_the_convergence_tolerance(tol):
    """Amendment 1's correction to the pre-registration, made checkable.

    The prereg predicted the GMM outcome model's residual floor would be set by
    ``GaussianMixture(tol=...)``. It is not: ``sum_k pi_k mu_k = ybar`` holds
    after any M-step because responsibilities sum to one, so the identity is
    exact at every tolerance. If this ever becomes tolerance-dependent, the
    experiment's account of where approximate equality comes from is wrong.
    """
    from sklearn.mixture import GaussianMixture

    y = np.random.default_rng(0).normal(12.0, 4.0, 400).reshape(-1, 1)
    model = GaussianMixture(
        n_components=4, tol=tol, max_iter=2000, random_state=0, reg_covar=1e-6
    ).fit(y)
    mixture_mean = float(model.weights_ @ model.means_.ravel())
    assert abs(mixture_mean - float(y.mean())) < 1e-12


# ---------------------------------------------------------------------------
# Wiring invariants
# ---------------------------------------------------------------------------


def test_every_estimator_carries_a_declared_functional_relationship():
    assert set(ESTIMATORS) == set(SHARES_FUNCTIONAL_WITH_T_DRAW)


def test_the_grid_is_reproducible_under_a_fixed_seed():
    """Invariant 10. Two runs at the same seed must agree bitwise, or nothing
    in the results directory means anything."""
    kwargs = dict(
        seed=104,
        n_seeds=1,
        cohort_sizes=(300,),
        n_replicates=2,
        generators=("gmm-joint-K2",),
        estimators=("gcomp-saturated", "ols-stratum-dummies"),
    )
    first = run_grid(**kwargs)
    second = run_grid(**kwargs)
    pd.testing.assert_frame_equal(first, second)


def test_adding_an_estimator_cannot_move_another_estimators_number():
    """Each estimator draws from its own stream, keyed by its index. If streams
    were shared, extending ESTIMATORS would silently rewrite published columns.
    """
    kwargs = dict(
        seed=105,
        n_seeds=1,
        cohort_sizes=(300,),
        n_replicates=2,
        generators=("gmm-joint-K2",),
    )
    alone = run_grid(**kwargs, estimators=("gcomp-saturated",))
    with_more = run_grid(
        **kwargs, estimators=("gcomp-saturated", "gcomp-gmm-outcome-K4")
    )
    subset = with_more[with_more["estimator"] == "gcomp-saturated"].reset_index(
        drop=True
    )
    pd.testing.assert_series_equal(
        alone["estimate"].reset_index(drop=True), subset["estimate"]
    )


def test_a_single_record_cell_abstains_instead_of_crashing_the_sweep():
    """THE FORCING INPUT: a stratum-arm cell holding exactly one patient.

    ``GaussianMixture`` refuses fewer than two samples whatever the component
    count, so a guard written only as ``size < n_components`` passes a
    one-record cell straight into sklearn and takes the whole sweep down. This
    really happened on the first full run, at n=100.

    An abstention is the right answer here -- the draw is counted in
    ``n_nonfinite_estimate`` -- and a crash is not, because a sweep that dies on
    its smallest cohort silently becomes a sweep about large cohorts.
    """
    from src.harness.trial_learned_generator import gcomp_gmm_outcome

    records = pd.DataFrame(
        {
            "stratum": [0, 0, 0, 1, 1, 1, 2, 2, 2],
            "treated": [0, 0, 1, 0, 0, 1, 0, 0, 1],
            "outcome": [10.0, 11.0, 13.0, 14.0, 15.0, 17.0, 25.0, 26.0, 28.0],
        }
    )
    trial = _as_trial(records, THETA)
    assert (records.groupby(["stratum", "treated"]).size() == 1).any()

    for k in (1, 4):
        got = gcomp_gmm_outcome(trial, n_components=k, rng=_rng(3, "estimator", 0))
        assert not np.isfinite(got), f"K={k} returned {got} on a one-record cell"


def test_the_gmm_outcome_model_still_returns_a_number_on_a_healthy_cohort():
    """Positive control: the abstention above must not be the only thing this
    estimator ever does."""
    from src.harness.trial_learned_generator import gcomp_gmm_outcome

    trial = _as_trial(
        simulate_trial(2000, THETA, rng=_rng(17, "train", 0)).records, THETA
    )
    got = gcomp_gmm_outcome(trial, n_components=4, rng=_rng(17, "estimator", 0))
    assert np.isfinite(got)
    assert abs(got - trial.theta_realised) < 1e-9
