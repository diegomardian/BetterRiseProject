"""The censored degeneracy claim, against the inputs that would falsify it.

Same rule as ``tests/test_checks_can_fail.py``: a check unable to fail is worse
than no check. Every claim :mod:`src.harness.trial_survival` makes is paired
here with a positive control that must pass and a constructed input that must
make it fail. The claims are

  (i)   at 0% censoring the residual check catches the exponential MLE exactly;
  (ii)  censoring silences that check against the LATENT truth while the
        estimator stays exactly the generator's functional on the observed data;
  (iii) by roughly half censoring the ordering inverts, so the paper's own
        column ranks the blind estimator below the informative one;

and the negative controls are, respectively, the stratified Cox model (which
must NOT be caught), the observed-data truth (against which the residual must
stay identically zero), and the 0% row (where the ordering must be correct).

The Cox implementation is cross-checked against ``lifelines`` in one dev-only
test. ``lifelines`` is deliberately absent from ``env/w2_harness.yml`` -- the
degeneracy claims need exact control of the functional -- and is used here only
as an independent implementation to check the hand-rolled one against.
"""

from __future__ import annotations

import math
import warnings

import numpy as np
import pandas as pd
import pytest

from src.common.io import RESERVED_META_KEYS
from src.harness.trial_survival import (
    BETA_REQUESTED,
    CENSORING_TARGETS,
    DEGENERATE_ON_OBSERVED,
    ESTIMATORS,
    NO_CENSORING_TIME,
    RMST_TAU,
    _cox,
    _extra_meta,
    _partial_likelihood,
    censor,
    censoring_time_for,
    empirical_rmst_standardised,
    exponential_mle_standardised,
    exponential_rmst_standardised,
    headline,
    km_rmst_standardised,
    marginal_censored_fraction,
    run,
    simulate_latent_cohort,
    stratified_cox_breslow,
    summarise,
    unadjusted_cox,
)

LOG_HR = "log_hazard_ratio"
RMST = "rmst_difference"


def _cohort(n=3000, seed=0, beta=BETA_REQUESTED):
    return simulate_latent_cohort(n, beta, rng=np.random.default_rng(seed))


def _trial(target, n=3000, seed=0):
    return censor(_cohort(n=n, seed=seed), censoring_time_for(target))


# ---------------------------------------------------------------------------
# The censoring grid. It is the sweep's x-axis, so a horizon that silently
# misses its target would move every claim below.
# ---------------------------------------------------------------------------


def test_the_grid_hits_the_censored_fractions_it_claims():
    for target in CENSORING_TARGETS:
        achieved = marginal_censored_fraction(censoring_time_for(target))
        assert achieved == pytest.approx(target, abs=1e-6)


def test_a_horizon_that_misses_its_target_is_caught_rather_than_reported():
    """The failing input: the 12% horizon used where 49% was wanted.

    Without this the grid could be mislabelled and every row of the sweep would
    still look plausible -- the residuals grow smoothly, so a wrong x-axis does
    not announce itself.
    """
    wrong = censoring_time_for(0.12)
    assert marginal_censored_fraction(wrong) != pytest.approx(0.49, abs=0.01)


def test_the_simulated_censored_fraction_matches_the_closed_form():
    """Closed form and generator must agree, or the grid describes a different
    generator from the one that ran."""
    for target in CENSORING_TARGETS:
        trial = _trial(target, n=20000, seed=3)
        assert trial.censored_fraction == pytest.approx(target, abs=0.01)


def test_the_no_censoring_horizon_censors_nothing_at_all():
    """Not 'almost nothing'. Claim (i) is an exact-zero claim, and one censored
    draw anywhere in the 0% arm would make it approximate."""
    for seed in range(5):
        trial = censor(_cohort(n=5000, seed=seed), NO_CENSORING_TIME)
        assert trial.censored_fraction == 0.0


def test_an_impossible_censored_fraction_is_refused():
    with pytest.raises(ValueError):
        censoring_time_for(1.0)


# ---------------------------------------------------------------------------
# Claim (i) and the degeneracy itself
# ---------------------------------------------------------------------------


def test_the_exponential_mle_is_the_observed_truth_exactly_at_every_censoring():
    """THE claim. Not ``approx`` -- identically, at every censoring level, or
    the residual check against the observed truth retains some information and
    the rule this module licenses is weaker than stated."""
    for target in CENSORING_TARGETS:
        for seed in range(10):
            trial = _trial(target, seed=seed)
            estimate = exponential_mle_standardised(trial.records)
            assert estimate == trial.theta_realised_observed[LOG_HR]


def test_at_zero_censoring_the_check_fires_against_the_latent_truth_too():
    """Claim (i). With nothing censored the two truths coincide, so the check
    works exactly as ``trial_recovery`` describes."""
    for seed in range(10):
        trial = _trial(0.0, seed=seed)
        estimate = exponential_mle_standardised(trial.records)
        assert estimate == trial.theta_realised_latent[LOG_HR]


def test_the_negative_control_is_not_caught_by_the_check():
    """The failing input for the degeneracy test: the stratified Cox model must
    NOT reproduce either truth. If it did, the check would be flagging
    everything and 'flagged' would carry no information."""
    for target in CENSORING_TARGETS:
        for seed in range(5):
            trial = _trial(target, seed=seed)
            estimate = stratified_cox_breslow(trial.records)
            assert abs(estimate - trial.theta_realised_observed[LOG_HR]) > 1e-6
            assert abs(estimate - trial.theta_realised_latent[LOG_HR]) > 1e-6


def test_the_degeneracy_table_is_measured_not_asserted():
    """``DEGENERATE_ON_OBSERVED`` is a claim about each estimator. Check every
    entry against the residual rather than trusting the dict."""
    trials = [_trial(t, n=2000, seed=s) for t in CENSORING_TARGETS for s in range(3)]
    for name, arm in ESTIMATORS.items():
        worst = max(
            abs(arm.fn(t.records) - t.theta_realised_observed[arm.scale]) for t in trials
        )
        assert (worst < 1e-9) == DEGENERATE_ON_OBSERVED[name], (
            f"{name}: measured degenerate={worst < 1e-9}, table says "
            f"{DEGENERATE_ON_OBSERVED[name]} (worst residual {worst:.3g})"
        )


# ---------------------------------------------------------------------------
# Claim (ii): censoring silences the check against the latent truth
# ---------------------------------------------------------------------------


def _median_latent_residual(target, *, estimator, n=3000, n_seeds=40):
    seeds = range(n_seeds)
    values = [
        abs(estimator(t.records) - t.theta_realised_latent[LOG_HR])
        for t in (_trial(target, n=n, seed=s) for s in seeds)
    ]
    return float(np.median(values))


def test_censoring_silences_the_check_against_the_latent_truth():
    """Claim (ii). The estimator does not change -- it is still exactly the
    generator's functional on the observed data, asserted above -- but the
    residual against the latent truth stops being zero, so the check reads it as
    independent of the generator."""
    silenced = _median_latent_residual(0.78, estimator=exponential_mle_standardised)
    assert silenced > 1e-3


def test_the_silencing_grows_monotonically_with_the_censored_fraction():
    """The dose-response. A single non-zero residual could be a bug; a residual
    that tracks the censored fraction is the mechanism."""
    medians = [
        _median_latent_residual(t, estimator=exponential_mle_standardised)
        for t in CENSORING_TARGETS
    ]
    assert medians[0] == 0.0
    assert all(b > a for a, b in zip(medians, medians[1:], strict=False)), medians


def test_the_residual_against_the_observed_truth_does_not_grow_at_all():
    """The other half of claim (ii), and the control that makes it a statement
    about the truth column rather than about the estimator: the same estimator,
    the same draws, the same censoring, evaluated against the observed truth,
    stays at exactly zero."""
    for target in CENSORING_TARGETS:
        worst = max(
            abs(
                exponential_mle_standardised(t.records)
                - t.theta_realised_observed[LOG_HR]
            )
            for t in (_trial(target, seed=s) for s in range(20))
        )
        assert worst == 0.0


# ---------------------------------------------------------------------------
# Claim (iii): the ordering inversion
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_the_ordering_inverts_by_half_censoring_and_not_before():
    """Claim (iii), with the 0% row as its own failing input.

    ``inversion_rate`` is the paired share of replicates on which the blind
    estimator's latent-truth residual is the larger. It is used instead of the
    ratio of maxima because a maximum over more replicates is mechanically
    larger, so a max-based verdict moves with the replicate count and this one
    does not.
    """
    runs = run(seed=20260909, n_seeds=1, cohort_sizes=(3000,), n_replicates=120)
    head = headline(runs).set_index("censoring_target")["inversion_rate"]
    assert head[0.0] == 0.0, "with nothing censored the blind estimator must win"
    assert head[0.12] < 0.5, "at 12% the informative estimator must still win"
    assert head[0.49] > 0.5, "by 49% the ranking must have inverted"
    assert head[0.78] > 0.5, "and it must not come back at 78%"
    # NOT asserted: that the rate keeps rising from 49% to 78%. At this
    # replicate count the two are equal to the digit, and writing a strict
    # inequality here would be tuning the test to one seed. Whether the
    # inversion deepens is measured in the sweep and reported from there.


def test_the_informative_estimator_is_not_merely_worse():
    """The inversion is only interesting if the estimator it demotes is a good
    one. If Cox were biased this would be blind-versus-broken rather than
    blind-versus-informative."""
    for target in CENSORING_TARGETS:
        trial = _trial(target, n=40000, seed=7)
        assert stratified_cox_breslow(trial.records) == pytest.approx(
            BETA_REQUESTED, abs=0.06
        )


def test_the_confounded_control_is_biased_so_the_check_has_something_to_catch():
    """``unadjusted-cox`` ignores strata that drive both assignment and hazard.
    Its score test at zero is the log-rank test, so this is the log-rank arm."""
    for target in CENSORING_TARGETS:
        trial = _trial(target, n=40000, seed=8)
        assert unadjusted_cox(trial.records) > BETA_REQUESTED + 0.2


# ---------------------------------------------------------------------------
# The positive control: censoring does NOT silence everything
# ---------------------------------------------------------------------------


def test_the_horizon_bounded_functional_stays_caught_at_every_censoring_level():
    """Administrative censoring at ``c >= tau`` never bites inside ``[0, tau]``,
    so the Kaplan-Meier RMST equals the empirical truncated mean of the LATENT
    times exactly. The check keeps working there.

    This is what stops claim (ii) from being read as 'censoring breaks the
    check', which is false. It breaks it for functionals that read the
    unobserved tail."""
    for target in CENSORING_TARGETS:
        for seed in range(5):
            trial = _trial(target, seed=seed)
            estimate = km_rmst_standardised(trial.records)
            assert estimate == pytest.approx(
                trial.theta_realised_latent[RMST], abs=1e-12
            )


def test_the_tail_reading_functional_on_the_same_scale_is_silenced():
    """The failing input for the control above: same scale, same horizon, same
    data, but it reaches ``tau`` through a fitted rate. If this were also exact,
    the RMST result would be about the scale rather than about the tail."""
    for target in CENSORING_TARGETS:
        trial = _trial(target, seed=1)
        assert (
            abs(
                exponential_rmst_standardised(trial.records)
                - trial.theta_realised_latent[RMST]
            )
            > 1e-6
        )


def test_a_horizon_past_the_last_observation_is_not_estimable():
    """RMST beyond the data is an extrapolation. It returns NaN, not a number --
    the same rule invariant 1 states for the intrinsic term."""
    trial = _trial(0.78, seed=2)
    beyond = trial.censoring_time * 2
    assert np.isnan(km_rmst_standardised(trial.records, tau=beyond))
    assert not np.isnan(km_rmst_standardised(trial.records, tau=RMST_TAU))


# ---------------------------------------------------------------------------
# Estimability. None is not 0.0 here either.
# ---------------------------------------------------------------------------


def test_a_cell_with_no_events_has_no_rate_rather_than_a_rate_of_zero():
    """The failing input: censor a stratum-arm cell so hard that it has no
    events. ``d / sum(t)`` is then 0, whose log is -inf, and a standardised
    effect built from it would be a large finite-looking number."""
    trial = _trial(0.78, n=400, seed=4)
    records = trial.records.copy()
    cell = (records["stratum"] == 0) & (records["treated"] == 1)
    records.loc[cell, "event"] = 0
    assert np.isnan(exponential_mle_standardised(records))
    assert np.isnan(exponential_rmst_standardised(records))
    assert not np.isnan(exponential_mle_standardised(trial.records))


def test_a_stratum_missing_an_arm_gives_no_number():
    """Positivity fails and the estimand has no referent in that stratum."""
    trial = _trial(0.49, n=400, seed=5)
    records = trial.records.copy()
    records.loc[records["stratum"] == 0, "treated"] = 1
    assert np.isnan(exponential_mle_standardised(records))
    assert np.isnan(km_rmst_standardised(records))


def test_a_cohort_with_no_events_at_all_gives_no_cox_estimate():
    trial = _trial(0.49, n=200, seed=6)
    records = trial.records.copy()
    records["event"] = 0
    assert np.isnan(stratified_cox_breslow(records))
    assert not np.isnan(stratified_cox_breslow(trial.records))


# ---------------------------------------------------------------------------
# The estimators must not be able to see the latent times
# ---------------------------------------------------------------------------


def test_the_observed_records_do_not_carry_the_latent_column():
    """Structural, and the reason it is worth a test: if ``latent_time`` were
    ever added to ``records`` for convenience, an estimator could read it, the
    residual against the latent truth would collapse to zero, and the sweep
    would report that censoring does not silence the check."""
    trial = _trial(0.49, n=100, seed=0)
    assert set(trial.records.columns) == {"stratum", "treated", "time", "event"}


def test_two_cohorts_with_the_same_observed_data_give_the_same_estimates():
    """The behavioural version of the test above. Two draws that differ only
    beyond the horizon must produce identical estimates and different latent
    truths -- that gap is exactly what the sweep measures."""
    cohort = _cohort(n=2000, seed=9)
    horizon = censoring_time_for(0.78)
    stretched = type(cohort)(
        stratum=cohort.stratum,
        treated=cohort.treated,
        latent_time=np.where(
            cohort.latent_time > horizon, cohort.latent_time * 3.0, cohort.latent_time
        ),
        beta_requested=cohort.beta_requested,
    )
    a, b = censor(cohort, horizon), censor(stretched, horizon)
    assert a.records.equals(b.records)
    for arm in ESTIMATORS.values():
        assert arm.fn(a.records) == arm.fn(b.records)
    assert a.theta_realised_latent[LOG_HR] != b.theta_realised_latent[LOG_HR]
    assert a.theta_realised_observed[LOG_HR] == b.theta_realised_observed[LOG_HR]


# ---------------------------------------------------------------------------
# The Cox implementation, against an independent one
# ---------------------------------------------------------------------------


def test_the_cox_implementation_agrees_with_lifelines():
    """Dev-only cross-check. ``lifelines`` is not in ``env/w2_harness.yml`` and
    must not be added there: these claims need exact control of the functional.
    It is in the ``[dev]`` extra, so it can be an independent implementation to
    check against.

    Event times here are continuous, so ties have probability zero and Breslow
    and lifelines' Efron coincide. Tie handling is checked separately below,
    against a brute-force Breslow likelihood, because this test cannot see it.
    """
    lifelines = pytest.importorskip("lifelines")

    for target in (0.0, 0.12, 0.49, 0.78):
        trial = _trial(target, n=1500, seed=5)
        frame = trial.records.rename(columns={"time": "T", "event": "E"}).astype(
            {"treated": float}
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            stratified = lifelines.CoxPHFitter().fit(
                frame[["T", "E", "treated", "stratum"]],
                duration_col="T",
                event_col="E",
                strata=["stratum"],
            )
            pooled = lifelines.CoxPHFitter().fit(
                frame[["T", "E", "treated"]], duration_col="T", event_col="E"
            )
        assert stratified_cox_breslow(trial.records) == pytest.approx(
            float(stratified.params_["treated"]), abs=1e-5
        )
        assert unadjusted_cox(trial.records) == pytest.approx(
            float(pooled.params_["treated"]), abs=1e-5
        )


def _brute_force_breslow(beta, x, time, event):
    """The Breslow partial log-likelihood written the slow, obvious way.

    One term per event, each dividing by the sum over everyone still at risk.
    No cumulative sums, so it shares no code with the implementation it checks.
    """
    total = 0.0
    for i in range(len(time)):
        if event[i]:
            at_risk = time >= time[i]
            total += beta * x[i] - math.log(np.exp(beta * x[at_risk]).sum())
    return total


def test_tie_handling_matches_a_brute_force_breslow_likelihood():
    """The case lifelines cannot check, because it uses Efron.

    Times are rounded to force heavy ties, then the Newton fit is compared to a
    grid maximum of an independently written likelihood.
    """
    trial = _trial(0.49, n=600, seed=11)
    records = trial.records.copy()
    records["time"] = np.round(records["time"], 1) + 0.1  # heavy ties, no zeros
    records = records.sort_values("time").reset_index(drop=True)
    assert records["time"].duplicated().sum() > 100, "the tie test needs ties"

    stratum, treated, time, event = (
        records["stratum"].to_numpy(),
        records["treated"].to_numpy(),
        records["time"].to_numpy(float),
        records["event"].to_numpy(float),
    )
    fitted = _cox(stratum, treated, time, event)

    grid = np.linspace(fitted - 0.3, fitted + 0.3, 601)
    loglik = [
        sum(
            _brute_force_breslow(
                b,
                treated[stratum == g].astype(float),
                time[stratum == g],
                event[stratum == g],
            )
            for g in np.unique(stratum)
        )
        for b in grid
    ]
    assert fitted == pytest.approx(float(grid[int(np.argmax(loglik))]), abs=2e-3)


def test_the_partial_likelihood_would_be_wrong_on_unsorted_input_so_cox_sorts():
    """The failing input for the defensive sort. ``_partial_likelihood`` builds
    risk sets from reverse cumulative sums and returns a wrong number rather
    than an error when its input is not ascending -- so this checks both that
    the low-level function really is order-dependent (or the guard is
    pointless) and that ``_cox`` is not."""
    trial = _trial(0.49, n=800, seed=12)
    stratum, treated, time, event = (
        trial.records["stratum"].to_numpy(),
        trial.records["treated"].to_numpy(),
        trial.records["time"].to_numpy(float),
        trial.records["event"].to_numpy(float),
    )
    shuffle = np.random.default_rng(0).permutation(len(time))

    sorted_score = _partial_likelihood(0.1, treated.astype(float), time, event)[1]
    shuffled_score = _partial_likelihood(
        0.1, treated[shuffle].astype(float), time[shuffle], event[shuffle]
    )[1]
    assert sorted_score != pytest.approx(shuffled_score, abs=1e-6)

    assert _cox(stratum, treated, time, event) == pytest.approx(
        _cox(stratum[shuffle], treated[shuffle], time[shuffle], event[shuffle]),
        abs=1e-9,
    )


# ---------------------------------------------------------------------------
# The sweep, the summary, and the sidecar
# ---------------------------------------------------------------------------


def test_the_sweep_is_reproducible_under_a_fixed_seed():
    a = run(seed=7, n_seeds=1, cohort_sizes=(200,), n_replicates=4)
    b = run(seed=7, n_seeds=1, cohort_sizes=(200,), n_replicates=4)
    assert a.equals(b)


def test_different_seeds_give_different_draws():
    a = run(seed=7, n_seeds=1, cohort_sizes=(200,), n_replicates=4)
    b = run(seed=8, n_seeds=1, cohort_sizes=(200,), n_replicates=4)
    assert not a["estimate"].equals(b["estimate"])


def test_the_censoring_sweep_is_paired_within_a_replicate():
    """The four horizons must be four views of ONE draw, or a difference across
    censoring levels could be a difference in what was drawn."""
    runs = run(seed=7, n_seeds=1, cohort_sizes=(500,), n_replicates=3)
    latent = runs[runs["scale"] == LOG_HR].groupby(["replicate"])[
        "theta_realised_latent"
    ].nunique()
    assert (latent == 1).all()


def test_the_summary_carries_both_verdicts_so_the_table_shows_the_reversal():
    runs = run(seed=9, n_seeds=1, cohort_sizes=(1000,), n_replicates=25)
    table = summarise(runs).set_index(["estimator", "censoring_target"])
    blind = "exponential-mle-standardised"
    assert table.loc[(blind, 0.0), "verdict_vs_latent"].startswith("flagged")
    assert table.loc[(blind, 0.78), "verdict_vs_latent"].startswith("passed")
    for target in CENSORING_TARGETS:
        assert table.loc[(blind, target), "verdict_vs_observed"].startswith("flagged")


def test_the_sidecar_does_not_overwrite_the_provenance_record():
    """Invariant 10's guard, checked before a long sweep discovers it. The
    writer raises on a clash, so a table would be lost at the last step."""
    runs = run(seed=9, n_seeds=1, cohort_sizes=(500,), n_replicates=3)
    meta = _extra_meta(seeds=[9], n_replicates=3, runs=runs, head=headline(runs))
    assert not (RESERVED_META_KEYS & set(meta))
    assert "SYNTHETIC" in meta and "what_this_answers" in meta


def test_the_sidecar_guard_would_catch_a_clash():
    """The failing input for the test above: without the guard, a sidecar key
    named ``seed`` would silently replace the seed invariant 10 records."""
    assert RESERVED_META_KEYS & {"seed", "git_sha"}


def test_every_estimator_reports_a_scale_and_the_truths_are_scale_matched():
    """A residual is only meaningful against a truth on its own scale. RMST is
    in time units and the hazard-ratio arms are on the log scale; comparing them
    would put a 0.09 beside a 0.0006 and read it as an improvement."""
    runs = run(seed=13, n_seeds=1, cohort_sizes=(500,), n_replicates=2)
    for scale, group in runs.groupby("scale"):
        assert group["theta_realised_latent"].nunique() <= group["replicate"].nunique()
        assert scale in {LOG_HR, RMST}
    empirical = empirical_rmst_standardised(_trial(0.0, n=500, seed=0).records)
    assert 0.0 < empirical < RMST_TAU


def test_the_headline_table_names_the_rule_it_licenses():
    runs = run(seed=14, n_seeds=1, cohort_sizes=(500,), n_replicates=5)
    head = headline(runs)
    assert head["rule"].str.contains("OBSERVED").all()
    assert set(head["censoring_target"]) == set(CENSORING_TARGETS)


def test_the_runs_table_has_no_silently_missing_estimates():
    """A NaN estimate is legitimate (positivity, no events) but must be rare at
    these cohort sizes. A column that is quietly half NaN would make every
    maximum below a statement about the surviving half."""
    runs = run(seed=15, n_seeds=1, cohort_sizes=(1000,), n_replicates=10)
    assert runs["estimate"].isna().mean() < 0.01
    assert not isinstance(runs["estimate"].iloc[0], pd.Series)
