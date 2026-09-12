"""The challenge harness, and the inputs that force each of its guards to fail.

Same rule as `tests/test_claim_check_inventory.py`: a validator with no failing
input is untested. The harness exists to say whether a check behaved as sealed,
so the thing that must not silently pass is a case whose expectation was wrong.
"""

from __future__ import annotations

import numpy as np
import pytest

from src.reference.guard_challenge import (
    CASES,
    ChallengeCase,
    ChallengeError,
    check_branch_is_directional,
    check_counts_three_outcomes,
    check_interval_is_calibrated,
    check_recovery_curve_is_informative,
    run_challenges,
    stop_rule_summary,
)


def _case(**kwargs) -> ChallengeCase:
    base = dict(
        case_id="X01",
        target_check="stub",
        description="stub",
        expected="passes",
        consequence="decision_relevant",
        run=lambda: False,
    )
    return ChallengeCase(**{**base, **kwargs})


# ---------------------------------------------------------------------------
# The committed set
# ---------------------------------------------------------------------------


def test_every_committed_case_behaves_as_sealed():
    frame = run_challenges(CASES)
    unexpected = frame.loc[~frame["as_expected"], "case_id"].tolist()
    assert not unexpected, f"cases whose sealed expectation was falsified: {unexpected}"


def test_the_set_carries_clean_controls_not_only_violations():
    """A challenge set of violations alone cannot show the audit discriminates."""
    frame = run_challenges(CASES)
    assert frame["is_clean_control"].sum() >= 4
    assert (~frame["is_clean_control"]).sum() >= 4


def test_the_set_is_not_a_held_out_evaluation_and_says_so():
    summary = stop_rule_summary(run_challenges(CASES))
    assert summary["any_case_blind_to_the_ledger"] is False


def test_running_twice_gives_the_same_answers():
    """Fixed seeds. A challenge set that moves between runs pins nothing."""
    first, second = run_challenges(CASES), run_challenges(CASES)
    assert first["observed"].tolist() == second["observed"].tolist()


def test_a_revision_after_falsification_is_recorded_in_the_table():
    """C05's construction error must stay visible, not be tidied away."""
    frame = run_challenges(CASES)
    revised = frame[frame["revised_after_falsification"]]
    assert len(revised) == 1
    assert revised.iloc[0]["case_id"] == "C05"
    assert "FALSIFIED" in revised.iloc[0]["revision_note"]


# ---------------------------------------------------------------------------
# The checks themselves, at their boundaries
# ---------------------------------------------------------------------------


def test_a_calibrated_rate_does_not_fire_and_a_miscalibrated_one_does():
    assert not check_interval_is_calibrated(0.050)
    assert check_interval_is_calibrated(0.107)


def test_a_zero_residual_curve_fires_and_a_real_one_does_not():
    assert check_recovery_curve_is_informative(np.zeros(50))
    assert not check_recovery_curve_is_informative(np.array([0.1, -0.2, 0.05]))


def test_the_three_outcome_check_needs_an_abstention_to_fire():
    """It must not fire on well-formed input: that would be flagging everything."""
    assert check_counts_three_outcomes(np.array([0.9, np.nan, 0.2]), 0.5)
    assert not check_counts_three_outcomes(np.array([0.9, 0.4, 0.2]), 0.5)


def test_the_directional_check_fires_only_on_a_wrong_sign():
    assert check_branch_is_directional(True, effect=-0.4, predicted_sign=+1)
    assert not check_branch_is_directional(True, effect=+0.4, predicted_sign=+1)
    assert not check_branch_is_directional(False, effect=-0.4, predicted_sign=+1)


# ---------------------------------------------------------------------------
# Forcing inputs for the harness's own validation
# ---------------------------------------------------------------------------


def test_an_unsealed_outcome_forces_a_refusal():
    with pytest.raises(ChallengeError, match="expected must be one of"):
        _case(expected="probably fine")


def test_an_unknown_consequence_forces_a_refusal():
    with pytest.raises(ChallengeError, match="consequence must be one of"):
        _case(consequence="interesting")


def test_a_silent_revision_forces_a_refusal():
    with pytest.raises(ChallengeError, match="unrecorded revision"):
        _case(revised_after_falsification=True, revision_note=None)


def test_duplicate_case_ids_force_a_refusal():
    with pytest.raises(ChallengeError, match="duplicate case ids"):
        run_challenges((_case(case_id="X01"), _case(case_id="X01")))


# ---------------------------------------------------------------------------
# The stop-rule arithmetic, which must exclude what it says it excludes
# ---------------------------------------------------------------------------


def test_a_decision_relevant_surprise_counts_toward_the_stop_rule():
    frame = run_challenges((_case(case_id="X01", expected="fires", run=lambda: False),))
    assert bool(frame.iloc[0]["counts_toward_stop_rule"])


def test_an_algebraic_identity_surprise_does_not_count():
    frame = run_challenges(
        (
            _case(
                case_id="X02",
                expected="passes",
                consequence="algebraic_identity",
                run=lambda: True,
            ),
        )
    )
    assert not bool(frame.iloc[0]["as_expected"])
    assert not bool(frame.iloc[0]["counts_toward_stop_rule"])


def test_an_input_validation_surprise_does_not_count():
    frame = run_challenges(
        (
            _case(
                case_id="X03",
                expected="passes",
                consequence="input_validation",
                run=lambda: True,
            ),
        )
    )
    assert not bool(frame.iloc[0]["as_expected"])
    assert not bool(frame.iloc[0]["counts_toward_stop_rule"])
