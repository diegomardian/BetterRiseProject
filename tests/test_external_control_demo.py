"""The external-control harm demonstration, against the inputs that force it to fail.

House rule, from ``tests/test_checks_can_fail.py``: *a check unable to fail is
worse than no check, because it turns an absence of evidence into a green
light.* This demonstration is itself a check -- it claims an observed-data
reference can mislead an analyst into deploying a defective estimator -- so it
owes the same debt. Every guard and every falsifier in
``docs/prereg_external_control_demo.md`` §8 gets an input here that MAKES IT
FIRE.

The falsifier tests are the load-bearing ones. ``evaluate_falsifiers`` reads a
result table and returns a verdict per falsifier; a version of it that returned
``fired = False`` on everything would let this experiment claim success on any
numbers at all. So each falsifier is handed a frame constructed to trip it.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.harness.external_control_demo import (
    EQUALITY_TOL,
    MCID,
    POOL_MIX,
    SMD_THRESHOLD,
    Registry,
    ate_ipw,
    ate_ipw_logistic,
    ate_standardisation,
    att_ipw,
    att_matching_all,
    att_standardisation,
    balance_table,
    build_registry,
    evaluate_falsifiers,
    run,
    run_shift_sweep,
    shifted_trial_mix,
    simulate_study,
    standardised_contrast,
    summarise,
)
from src.harness.trial_recovery import _standardised_effect


@pytest.fixture(scope="module")
def registry() -> Registry:
    return build_registry()


def _small_draw(registry: Registry, *, n_trial: int = 400, shift: float = 1.0, seed: int = 7):
    return simulate_study(
        n_trial, registry=registry, rng=np.random.default_rng(seed), shift=shift
    )


# ---------------------------------------------------------------------------
# The registry is a FIXED pool. That is the plasmode mechanic, not a detail.
# ---------------------------------------------------------------------------


def test_the_pool_is_the_same_object_every_time():
    """If the pool were re-drawn per replicate this would not be a plasmode design.

    The whole premise is a finite set of records that enrolment samples FROM.
    Redrawing it each replicate would restore a parametric generator and with it
    the reviewer's ``theta``, which is the assumption under test.
    """
    a = build_registry().records
    b = build_registry().records
    pd.testing.assert_frame_equal(a, b)


def test_a_pool_too_small_for_both_arms_is_refused(registry: Registry):
    """FORCING INPUT for the size guard: ask for more records than exist."""
    tiny = Registry(records=registry.records.iloc[:100].reset_index(drop=True))
    with pytest.raises(ValueError, match="pool of 100"):
        simulate_study(50, registry=tiny, rng=np.random.default_rng(0))


def test_the_arms_are_disjoint_records(registry: Registry):
    """Contemporaneous but distinct patients. Overlap would make both arms read
    the same record's outcome under two different arms."""
    draw = _small_draw(registry)
    assert len(np.unique(draw.trial_index)) == len(draw.trial_index)


# ---------------------------------------------------------------------------
# F3b, asserted rather than asserted-about: obs-pooled IS the paper's functional
# ---------------------------------------------------------------------------


def test_obs_pooled_is_the_functional_the_paper_already_analyses(registry: Registry):
    """FALSIFIER F3b. If `obs-pooled` were some bespoke definition invented to
    make this demonstration work, the contrivance charge would land. It is not:
    it is bit-identical to ``trial_recovery._standardised_effect``, the reference
    ``blind.tex`` analyses and the one the repository already shipped.
    """
    draw = _small_draw(registry)
    assert standardised_contrast(draw.records, weights_from="pooled") == (
        _standardised_effect(draw.records)
    )


def test_obs_trial_is_not_that_functional(registry: Registry):
    """The paired half. If the two weightings agreed there would be no defect to
    demonstrate, and the previous test would be vacuous."""
    draw = _small_draw(registry)
    assert abs(
        standardised_contrast(draw.records, weights_from="trial")
        - standardised_contrast(draw.records, weights_from="pooled")
    ) > 0.1


# ---------------------------------------------------------------------------
# Positivity. An unestimable contrast is NaN, never a number.
# ---------------------------------------------------------------------------


def test_a_stratum_missing_an_arm_gives_nan_not_a_number():
    """FORCING INPUT for the positivity guard: stratum 1 has no control.

    CLAUDE.md invariant 1 in a different vocabulary -- an unestimable quantity is
    not zero and is not the estimate over the strata that happened to work.
    """
    records = pd.DataFrame(
        {
            "stratum": [0, 0, 1, 1],
            "treated": [1, 0, 1, 1],
            "outcome": [5.0, 3.0, 9.0, 8.0],
        }
    )
    assert np.isnan(standardised_contrast(records, weights_from="pooled"))
    assert np.isnan(standardised_contrast(records, weights_from="trial"))


@pytest.mark.parametrize(
    "estimator",
    [att_standardisation, att_ipw, att_matching_all, ate_standardisation,
     ate_ipw, ate_ipw_logistic],
)
def test_every_estimator_abstains_on_a_degenerate_stratum(estimator):
    """FORCING INPUT: the same frame, through every estimator. A single arm that
    returns a finite number here would put a positivity failure into the sweep as
    an estimate."""
    from src.harness.external_control_demo import PlasmodeDraw

    records = pd.DataFrame(
        {
            "stratum": [0, 0, 1, 1],
            "treated": [1, 0, 1, 1],
            "outcome": [5.0, 3.0, 9.0, 8.0],
        }
    )
    draw = PlasmodeDraw(records=records, trial_index=np.array([0, 2, 3]), references={})
    assert np.isnan(estimator(draw))


# ---------------------------------------------------------------------------
# The identities the demonstration rests on
# ---------------------------------------------------------------------------


def test_the_three_att_estimators_are_one_functional(registry: Registry):
    """Three literatures, three vocabularies, the same arithmetic -- the paper's
    own point, made at the correct target population."""
    draw = _small_draw(registry)
    a, b, c = att_standardisation(draw), att_ipw(draw), att_matching_all(draw)
    assert abs(a - b) < 1e-12 and abs(a - c) < 1e-12


def test_the_defective_estimator_reproduces_the_mistaken_reference(registry: Registry):
    """The equality that makes the validation unable to fail."""
    draw = _small_draw(registry)
    assert abs(ate_standardisation(draw) - draw.references["obs-pooled"]) == 0.0
    assert abs(ate_ipw(draw) - draw.references["obs-pooled"]) < EQUALITY_TOL


def test_the_defective_estimator_does_not_reproduce_the_honest_reference(registry: Registry):
    """FORCING INPUT for the demonstration's central claim: if the defective
    estimator also reproduced ``po-trial``, the honest reference would buy
    nothing and there would be no counterfactual to show."""
    draw = _small_draw(registry)
    assert abs(ate_ipw(draw) - draw.references["po-trial"]) > 0.5


def test_no_estimator_can_reproduce_a_potential_outcome_reference(registry: Registry):
    """``po-trial`` reads both arms of the same record. Nothing that sees the
    analysis file can be that functional, and this is why it is the honest one."""
    draw = _small_draw(registry, n_trial=800)
    for fn in (att_standardisation, att_ipw, att_matching_all,
               ate_standardisation, ate_ipw):
        assert abs(fn(draw) - draw.references["po-trial"]) > EQUALITY_TOL


# ---------------------------------------------------------------------------
# The demonstration's own null. This is the check that lets it fail.
# ---------------------------------------------------------------------------


def test_at_zero_case_mix_shift_the_harm_vanishes(registry: Registry):
    """FALSIFIER F4, as a unit test.

    At ``shift = 0`` the trial enrols the registry's own case mix, the pooled and
    enrolled populations coincide, and the defective estimator is no longer
    defective in consequence -- though it is still the same functional as the
    mistaken reference. **Equality is not itself harm**, and an experiment whose
    alleged harm survived its own null would be measuring something other than
    case-mix shift.
    """
    assert shifted_trial_mix(0.0) == POOL_MIX
    draw = _small_draw(registry, n_trial=1600, shift=0.0)
    assert abs(ate_ipw(draw) - draw.references["po-trial"]) < 0.25
    # ...and the equality is still there, which is the point of the sentence above.
    assert abs(ate_ipw(draw) - draw.references["obs-pooled"]) < EQUALITY_TOL


def test_at_full_shift_the_harm_is_present(registry: Registry):
    """The paired half, so the previous test is not passing for want of signal."""
    draw = _small_draw(registry, n_trial=1600, shift=1.0)
    assert ate_ipw(draw) - draw.references["po-trial"] > 0.5


# ---------------------------------------------------------------------------
# The balance table an analyst actually runs
# ---------------------------------------------------------------------------


def test_the_routine_balance_table_is_clean_under_the_defective_weights(registry: Registry):
    """FALSIFIER F3e, as a unit test. The diagnostic checks balance BETWEEN ARMS;
    the defect is in balance against the ENROLLED population, which no protocol
    tabulates. Both columns are computed here so the contrast is visible."""
    table = balance_table(_small_draw(registry, n_trial=1600))
    ate = table[table["weighting"] == "ate"]
    assert ate["smd_between_arms"].abs().max() < SMD_THRESHOLD
    assert ate["smd_vs_enrolled"].abs().max() > 0.5


def test_the_balance_table_can_report_imbalance(registry: Registry):
    """FORCING INPUT for the balance computation itself: unweighted, where the
    case-mix shift is plainly visible. A balance function that returned ~0 on
    everything would make the previous test meaningless."""
    table = balance_table(_small_draw(registry, n_trial=1600))
    unweighted = table[table["weighting"] == "unweighted"]
    assert unweighted["smd_between_arms"].abs().max() > SMD_THRESHOLD


# ---------------------------------------------------------------------------
# THE FALSIFIERS. Each one, against a frame built to make it fire.
# ---------------------------------------------------------------------------


def _frames(n_sizes=(50, 100, 200, 400, 800, 1600)):
    """A minimal primary/shift/balance triple on which NO falsifier fires.

    The baseline the forcing tests mutate. Numbers are the shape the real run
    produces, not the real run's values.
    """
    from src.harness.external_control_demo import PRIMARY_ESTIMATORS, REFERENCES, TARGETS

    rows = []
    for est in list(PRIMARY_ESTIMATORS) + ["ate-ipw-logistic"]:
        defective = est.startswith("ate-")
        for ref in REFERENCES:
            for n in n_sizes:
                if ref == "obs-pooled":
                    resid = 0.0 if defective else 1.2
                    bias = 0.0 if defective else -0.87
                    ratio = 1.0 if defective else 0.53
                elif ref == "po-trial":
                    resid = 1.1
                    bias = 0.88 if defective else 0.01
                    ratio = 1.88 if defective else 1.0
                else:
                    resid, bias, ratio = 0.9, 0.0, 1.0
                rows.append({
                    "estimator": est, "target_population": TARGETS[est],
                    "reference": ref, "n_trial": n, "n_valid": 200, "n_excluded": 0,
                    "max_residual": resid, "equality_flagged": resid <= EQUALITY_TOL,
                    "bias": bias, "rmse": abs(bias) + (0.0 if defective else 0.01),
                    "ratio_median": ratio, "ratio_q25": ratio, "ratio_q75": ratio,
                    "estimate_median": 1.89 if defective else 1.03,
                    "reference_median": 0.99 if ref == "po-trial" else 1.02,
                })
    primary = pd.DataFrame(rows)

    shift = pd.DataFrame([
        {"shift": s, "estimator": e, "target_population": TARGETS[e], "n_trial": 1600,
         "trial_mix": "x", "estimate_median": 1.89, "po_trial_median": t,
         "gap_vs_po_trial": 1.89 - t,
         "bias_vs_po_trial": (0.02 if s == 0.0 else 0.88) if e.startswith("ate-") else 0.0,
         "max_residual_vs_obs_pooled": 0.0, "equality_flagged_vs_obs_pooled": True,
         "mcid": MCID, "estimate_clears_mcid": True, "truth_clears_mcid": t > MCID,
         "decision_flips": e.startswith("ate-") and s >= 0.75}
        for s, t in [(0.0, 2.13), (0.25, 1.85), (0.5, 1.57), (0.75, 1.29), (1.0, 1.0)]
        for e in ["ate-ipw", "att-standardisation"]
    ])

    balance = pd.DataFrame([
        {"weighting": w, "stratum": g, "stratum_name": str(g), "n_replicates": 50,
         "mean_smd_between_arms": 0.0, "max_abs_smd_between_arms": 1e-16,
         "mean_smd_vs_enrolled": 0.8, "max_abs_smd_vs_enrolled": 0.96,
         "mean_share_pseudo_population": 0.4, "mean_share_enrolled": 0.15,
         "smd_threshold": SMD_THRESHOLD, "between_arms_passes": True,
         "vs_enrolled_passes": False}
        for w in ("att", "ate", "unweighted") for g in range(3)
    ])
    return primary, shift, balance


def _fired(frames) -> set[str]:
    out = evaluate_falsifiers(*frames)
    return set(out.loc[out["fired"], "falsifier"])


def test_the_baseline_fires_no_falsifier():
    """If this frame tripped something, every forcing test below would be
    ambiguous about what it had actually forced."""
    assert _fired(_frames()) == set()


def test_f1_fires_when_the_validation_does_not_pass():
    """FORCING INPUT for F1: give the defective arm a visible residual against
    ``obs-pooled``. The analyst would then have been warned, and this would not
    be a harm case."""
    primary, shift, balance = _frames()
    mask = (primary["reference"] == "obs-pooled") & (primary["estimator"] == "ate-ipw")
    primary.loc[mask, ["max_residual", "rmse"]] = [0.4, 0.4]
    primary.loc[mask, "equality_flagged"] = False
    assert "F1" in _fired((primary, shift, balance))


def test_f1_fires_when_the_correct_estimator_wins_the_bakeoff():
    """FORCING INPUT for F1(c) specifically: the defect is still invisible in the
    residual, but the validation no longer RANKS the defective arm first. That is
    a materially weaker claim and must be caught separately."""
    primary, shift, balance = _frames()
    mask = primary["reference"] == "obs-pooled"
    primary.loc[mask & (primary["estimator"] == "att-ipw"), "rmse"] = -1.0
    assert "F1" in _fired((primary, shift, balance))


def test_f2_fires_when_the_honest_reference_sees_nothing():
    """FORCING INPUT for F2: shrink the defective arm's bias against ``po-trial``
    below the pre-registered 0.2 months. The audit would then buy nothing."""
    primary, shift, balance = _frames()
    mask = (primary["reference"] == "po-trial") & primary["estimator"].str.startswith("ate-")
    primary.loc[mask, "bias"] = 0.05
    assert "F2" in _fired((primary, shift, balance))


def test_f2_fires_when_the_bias_shrinks_like_sampling_error():
    """FORCING INPUT for F2's second clause: a bias that halves with cohort size
    is an ordinary convergence story, not a defect."""
    primary, shift, balance = _frames()
    mask = (primary["reference"] == "po-trial") & primary["estimator"].str.startswith("ate-")
    primary.loc[mask, "bias"] = 4.0 / np.sqrt(primary.loc[mask, "n_trial"])
    assert "F2" in _fired((primary, shift, balance))


def test_f3a_fires_when_the_defect_needs_the_ipw_weight_slip():
    """FORCING INPUT for F3a: make the plain-standardisation route stop
    exhibiting the defect. The failure would then be a weight-formula curiosity
    rather than a way of thinking, and the contrivance charge would land."""
    primary, shift, balance = _frames()
    mask = (primary["reference"] == "po-trial") & (primary["estimator"] == "ate-standardisation")
    primary.loc[mask, "bias"] = 0.01
    assert "F3a" in _fired((primary, shift, balance))


def test_f3c_fires_when_only_hand_counted_weights_are_equal():
    """FORCING INPUT for F3c: push the fitted-propensity arm's residual above any
    tolerance an analyst would choose. The equality would then be an artefact of
    the saturated-discrete implementation and the harm case would narrow."""
    primary, shift, balance = _frames()
    mask = (primary["reference"] == "obs-pooled") & (primary["estimator"] == "ate-ipw-logistic")
    primary.loc[mask, "max_residual"] = 1e-3
    assert "F3c" in _fired((primary, shift, balance))


def test_f3d_fires_when_only_the_most_extreme_shift_flips_the_decision():
    """FORCING INPUT for F3d: confine the decision flip to ``shift = 1.0``. A
    flip that needs the most extreme mix we wrote down and nothing less is a
    designed coincidence."""
    primary, shift, balance = _frames()
    shift.loc[(shift["shift"] < 1.0), "decision_flips"] = False
    assert "F3d" in _fired((primary, shift, balance))


def test_f3d_also_fires_when_the_decision_never_flips():
    """The other end of the same charge: no flip at all is no actionable harm."""
    primary, shift, balance = _frames()
    shift["decision_flips"] = False
    assert "F3d" in _fired((primary, shift, balance))


def test_f3e_fires_when_the_routine_balance_table_would_have_caught_it():
    """FORCING INPUT for F3e: make the post-ATE-weighting between-arm SMD exceed
    the conventional 0.1. The analyst would have had an independent warning and
    the audit would be redundant."""
    primary, shift, balance = _frames()
    balance.loc[balance["weighting"] == "ate", "max_abs_smd_between_arms"] = 0.4
    assert "F3e" in _fired((primary, shift, balance))


def test_f4_fires_when_the_harm_survives_the_null():
    """FORCING INPUT for F4, and the most important one in this file. If the
    defective arm is still biased at ``shift = 0``, where the pooled and enrolled
    populations coincide, then whatever is being measured is not case-mix shift
    and the whole demonstration is mis-attributed."""
    primary, shift, balance = _frames()
    shift.loc[(shift["shift"] == 0.0) & (shift["estimator"] == "ate-ipw"),
              "bias_vs_po_trial"] = 0.9
    assert "F4" in _fired((primary, shift, balance))


def test_f5_fires_when_the_realised_truth_is_the_parametric_one():
    """FORCING INPUT for F5: collapse ``po-trial`` onto ``induced-trial``. The
    reviewer's assumption that ``theta`` is available would then survive largely
    intact, and §1.1 of the pre-registration would be unsupported."""
    primary, shift, balance = _frames()
    primary.loc[primary["reference"] == "induced-trial", "reference_median"] = 0.99
    primary.loc[primary["reference"] == "po-trial", "reference_median"] = 0.99
    assert "F5" in _fired((primary, shift, balance))


def test_an_unevaluable_falsifier_counts_as_fired():
    """FORCING INPUT for the evaluator itself: hand it empty frames. A falsifier
    that cannot be evaluated is not a passed one, and a report that silently
    treated missing evidence as a clean bill of health would be the exact failure
    mode this project is about."""
    empty = pd.DataFrame(columns=_frames()[0].columns)
    out = evaluate_falsifiers(empty, pd.DataFrame(columns=_frames()[1].columns),
                              pd.DataFrame(columns=_frames()[2].columns))
    # F3b is asserted in this file rather than computed, so it is the one
    # exception; everything else must fire on no evidence.
    assert set(out.loc[~out["fired"], "falsifier"]) == {"F3b"}


# ---------------------------------------------------------------------------
# End to end, small
# ---------------------------------------------------------------------------


def test_the_sweep_produces_the_pre_registered_shape(registry: Registry):
    runs = run(seed=1, registry=registry, cohort_sizes=(200, 400), n_replicates=3)
    table = summarise(runs)
    assert set(table["reference"]) == {
        "obs-pooled", "obs-trial", "po-trial", "po-pooled", "induced-trial"
    }
    assert {"max_residual", "equality_flagged", "bias", "rmse",
            "ratio_median", "n_valid", "n_excluded"} <= set(table.columns)
    assert (table["n_valid"] + table["n_excluded"] == 3).all()


def test_the_shift_sweep_carries_its_own_null(registry: Registry):
    table = run_shift_sweep(
        seed=1, registry=registry, shifts=(0.0, 1.0), n_trial=400, n_replicates=3
    )
    assert set(table["shift"]) == {0.0, 1.0}
    assert "decision_flips" in table.columns
