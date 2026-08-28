"""The diagnostic that would have caught issue #44 without anyone looking.

Both directions matter and are tested separately: it must fire on the Lee
pattern, and it must stay quiet on a labeller that does not read dropout.
A diagnostic that always says "confounded" is as useless as the guards this
project keeps withdrawing.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.harness.depth_confound import (
    DEPTH_RATIO_TOLERANCE,
    MATURITY_DEPTH_RHO_TOLERANCE,
    depth_confound_report,
    match_arm_depth,
    mature_share_by_depth,
    max_attainable_rho,
)


def _lee_shaped(seed: int = 0, n: int = 4000):
    """The SMC pattern: shallow normal, deep tumour, maturity = marker dropout.

    Depths and the 4.3x arm imbalance are taken from the measured cohort, and
    the maturity call is generated the way the labeller generates it — a cell is
    "mature" when it happens to sample zero of its markers.
    """
    rng = np.random.default_rng(seed)
    arm = np.where(rng.random(n) < 0.15, "normal", "tumour")
    depth = np.where(
        arm == "normal",
        rng.lognormal(np.log(4519), 0.8, n),
        rng.lognormal(np.log(19244), 0.8, n),
    )
    # P(all five markers sampled zero) falls with depth — dropout, nothing else.
    p_zero = np.exp(-depth / 9000.0)
    is_mature = rng.random(n) < p_zero
    return depth, is_mature, arm


def _clean(seed: int = 0, n: int = 4000):
    """Depth-matched arms and a maturity call independent of depth."""
    rng = np.random.default_rng(seed)
    arm = np.where(rng.random(n) < 0.5, "normal", "tumour")
    depth = rng.lognormal(np.log(15000), 0.8, n)
    base = np.where(arm == "normal", 0.60, 0.35)  # a real biological difference
    is_mature = rng.random(n) < base
    return depth, is_mature, arm


# ---------------------------------------------------------------------------
# It fires on the real pattern
# ---------------------------------------------------------------------------


def test_it_fires_on_the_lee_pattern():
    report = depth_confound_report(*_lee_shaped())
    assert report["maturity_tracks_depth"]
    assert not report["arms_are_depth_matched"]
    assert report["confounded"]
    assert "CONFOUNDED" in report["reading"]


def test_the_lee_pattern_reproduces_the_direction_that_matters():
    """Shallow normal + dropout-as-maturity => normal looks MORE mature.

    That is the sign that makes it dangerous: it is the project's hypothesis.
    """
    depth, is_mature, arm = _lee_shaped()
    report = depth_confound_report(depth, is_mature, arm)
    assert (
        report["per_arm"]["normal"]["mature_share"]
        > report["per_arm"]["tumour"]["mature_share"]
    )
    assert report["depth_ratio_between_arms"] > 3.0


def test_the_dose_response_is_monotone_when_dropout_drives_the_call():
    depth, is_mature, _ = _lee_shaped()
    table = mature_share_by_depth(depth, is_mature)
    shares = table["mature_share"].to_numpy()
    assert shares[0] > shares[-1]
    # The claim is MONOTONICITY, so it is a rank correlation. Pearson on the
    # raw axis reads -0.73 here only because depth is lognormal and the top
    # decile sits far out — that is a fact about the x-scale, not about the
    # trend. Measured on the real SMC epithelium the rank correlation is -0.98.
    import pandas as pd

    rank_depth = pd.Series(table["median_depth"]).rank()
    rank_share = pd.Series(shares).rank()
    assert float(np.corrcoef(rank_depth, rank_share)[0, 1]) < -0.9


# ---------------------------------------------------------------------------
# It stays quiet when it should
# ---------------------------------------------------------------------------


def test_it_does_not_fire_on_a_clean_labeller():
    """A real between-arm maturity difference, depth-matched, must read clean."""
    report = depth_confound_report(*_clean())
    assert not report["maturity_tracks_depth"]
    assert report["arms_are_depth_matched"]
    assert not report["confounded"]
    assert "Clean on both counts" in report["reading"]


def test_depth_sensitivity_alone_is_not_called_confounded():
    """Reading dropout is only fatal once the arms differ in depth.

    Kept separate because the fixes differ: one is the labeller's, the other is
    the cohort's.
    """
    rng = np.random.default_rng(3)
    n = 4000
    arm = np.where(rng.random(n) < 0.5, "normal", "tumour")
    depth = rng.lognormal(np.log(12000), 0.8, n)  # matched arms
    is_mature = rng.random(n) < np.exp(-depth / 9000.0)  # but depth-driven
    report = depth_confound_report(depth, is_mature, arm)
    assert report["maturity_tracks_depth"]
    assert report["arms_are_depth_matched"]
    assert not report["confounded"]
    assert "does not convert" in report["reading"]


def test_depth_imbalance_alone_is_not_called_confounded():
    rng = np.random.default_rng(4)
    n = 4000
    arm = np.where(rng.random(n) < 0.3, "normal", "tumour")
    depth = np.where(arm == "normal", 4000.0, 20000.0) * rng.lognormal(0, 0.3, n)
    is_mature = rng.random(n) < 0.4  # independent of everything
    report = depth_confound_report(depth, is_mature, arm)
    assert not report["maturity_tracks_depth"]
    assert not report["arms_are_depth_matched"]
    assert not report["confounded"]
    assert "route from it" in report["reading"]


# ---------------------------------------------------------------------------
# The thresholds and the arithmetic
# ---------------------------------------------------------------------------


def test_the_tolerances_are_what_the_docstring_says():
    assert DEPTH_RATIO_TOLERANCE == 1.5
    assert MATURITY_DEPTH_RHO_TOLERANCE == 0.20


def test_rho_is_computed_within_arm_not_pooled():
    """Pooling lets a real between-arm difference look like a depth effect.

    Construct exactly that: depth is perfectly separated by arm, maturity is
    perfectly separated by arm, and *within* either arm nothing tracks anything.
    Pooled, the correlation is near 1; within-arm it is nan/0.
    """
    n = 2000
    rng = np.random.default_rng(5)
    arm = np.where(np.arange(n) < n // 2, "normal", "tumour")
    depth = np.where(arm == "normal", 4000.0, 20000.0) * rng.lognormal(0, 0.05, n)
    is_mature = arm == "normal"  # entirely determined by arm
    report = depth_confound_report(depth, is_mature, arm)
    # Within each arm is_mature is constant, so rho is undefined, not 1.0.
    assert not report["maturity_tracks_depth"]
    assert not report["confounded"]


@pytest.mark.parametrize(
    ("depth", "mature", "arm"),
    [
        (np.array([1.0, 2.0]), np.array([True]), np.array(["a", "b"])),
        (np.array([]), np.array([]), np.array([])),
        (np.array([-1.0, 2.0]), np.array([True, False]), np.array(["a", "b"])),
    ],
)
def test_it_refuses_malformed_input(depth, mature, arm):
    with pytest.raises(ValueError):
        depth_confound_report(depth, mature, arm)


def test_share_by_depth_refuses_mismatched_lengths():
    with pytest.raises(ValueError):
        mature_share_by_depth(np.array([1.0, 2.0]), np.array([True]))


# ---------------------------------------------------------------------------
# Matching the arms, not just flooring them
# ---------------------------------------------------------------------------


def test_matching_equalises_the_arms_depth_distributions():
    """The point of matching: afterwards there is no depth difference to blame."""
    depth, is_mature, arm = _lee_shaped()
    before = depth_confound_report(depth, is_mature, arm)
    keep = match_arm_depth(depth, arm, seed=1)
    after = depth_confound_report(depth[keep], is_mature[keep], arm[keep])

    assert before["depth_ratio_between_arms"] > 3.0
    assert after["depth_ratio_between_arms"] < 1.2
    assert after["arms_are_depth_matched"]
    assert not after["confounded"]


def test_matching_keeps_equal_numbers_from_each_arm():
    depth, _, arm = _lee_shaped()
    keep = match_arm_depth(depth, arm, seed=1)
    kept = pd.Series(arm[keep]).value_counts()
    assert kept["normal"] == kept["tumour"]


def test_matching_discards_most_of_the_larger_arm_and_that_is_the_cost():
    """A matched subsample is a smaller cohort. Callers must report the size."""
    depth, _, arm = _lee_shaped()
    keep = match_arm_depth(depth, arm, seed=1)
    assert keep.sum() < len(depth) / 2


def test_matching_removes_a_gap_that_was_only_depth():
    """A world where maturity is dropout and NOTHING else: the gap must go."""
    depth, is_mature, arm = _lee_shaped()
    keep = match_arm_depth(depth, arm, seed=1)
    per_arm = pd.DataFrame({"arm": arm, "m": is_mature, "k": keep})
    before = per_arm.groupby("arm", observed=True)["m"].mean()
    after = per_arm[per_arm.k].groupby("arm", observed=True)["m"].mean()
    gap_before = abs(before["normal"] - before["tumour"])
    gap_after = abs(after["normal"] - after["tumour"])
    assert gap_before > 0.15
    assert gap_after < gap_before / 2, (gap_before, gap_after)


def test_matching_preserves_a_gap_that_is_not_depth():
    """The other direction, and the one the Lee result rests on.

    Depth-matched arms with a genuine biological difference: matching must NOT
    erase it, or the method would launder real signal away.
    """
    rng = np.random.default_rng(11)
    n = 6000
    arm = np.where(rng.random(n) < 0.3, "normal", "tumour")
    depth = np.where(
        arm == "normal",
        rng.lognormal(np.log(4519), 0.8, n),
        rng.lognormal(np.log(19244), 0.8, n),
    )
    truth = np.where(arm == "normal", 0.80, 0.55)  # real, depth-independent
    is_mature = rng.random(n) < truth
    keep = match_arm_depth(depth, arm, seed=2)
    frame = pd.DataFrame({"arm": arm, "m": is_mature})[keep]
    after = frame.groupby("arm", observed=True)["m"].mean()
    assert (after["normal"] - after["tumour"]) == pytest.approx(0.25, abs=0.08)


def test_matching_refuses_anything_other_than_two_arms():
    depth = np.array([1.0, 2.0, 3.0])
    with pytest.raises(ValueError, match="two arms"):
        match_arm_depth(depth, np.array(["a", "b", "c"]), seed=0)


# ---------------------------------------------------------------------------
# The prevalence ceiling — a rare label cannot reach the tolerance
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("p", "expected"),
    [(0.50, 0.866), (0.25, 0.750), (0.05, 0.3775), (0.0086, 0.1599)],
)
def test_max_attainable_rho_matches_the_closed_form(p, expected):
    assert max_attainable_rho(p) == pytest.approx(expected, abs=1e-3)


def test_the_ceiling_is_actually_attained_under_perfect_separation():
    """Not a formula asserted — the bound is reached by a constructed extreme."""
    n = 20000
    for p in (0.5, 0.05, 0.0086):
        y = np.zeros(n)
        y[: max(1, round(p * n))] = 1.0
        x = np.arange(n, dtype=float)
        rho = abs(np.corrcoef(pd.Series(x).rank(), pd.Series(y).rank())[0, 1])
        assert rho == pytest.approx(max_attainable_rho(p), abs=2e-3)


def test_a_rare_label_is_reported_as_not_testable_rather_than_clean():
    """The Lee best4 case: 8 mature cells in 5,564 caps |rho| at 0.066.

    Reporting that as "clean" was reading a ceiling as evidence. The diagnostic
    now says the test could not have fired.
    """
    rng = np.random.default_rng(7)
    n = 6000
    arm = np.where(rng.random(n) < 0.5, "normal", "tumour")
    depth = rng.lognormal(np.log(15000), 0.8, n)
    is_mature = np.zeros(n, dtype=bool)
    # ~0.9% prevalence, and entirely determined by depth
    is_mature[np.argsort(depth)[: int(0.009 * n)]] = True

    report = depth_confound_report(depth, is_mature, arm)
    assert not report["tolerance_is_reachable"]
    assert report["max_attainable_rho"] < MATURITY_DEPTH_RHO_TOLERANCE
    assert not report["maturity_tracks_depth"]  # it cannot fire
    assert "NOT TESTABLE" in report["reading"]


def test_rho_vs_ceiling_makes_rungs_comparable():
    """Raw |rho| across rungs whose prevalence differs by orders of magnitude is
    not one criterion. The normalised value is what can be compared."""
    rng = np.random.default_rng(3)
    n = 8000
    arm = np.where(rng.random(n) < 0.5, "normal", "tumour")
    depth = rng.lognormal(np.log(15000), 0.8, n)

    common = np.zeros(n, dtype=bool)
    common[np.argsort(depth)[: int(0.5 * n)]] = True
    rare = np.zeros(n, dtype=bool)
    rare[np.argsort(depth)[: int(0.02 * n)]] = True

    r_common = depth_confound_report(depth, common, arm)
    r_rare = depth_confound_report(depth, rare, arm)
    # Both are perfectly depth-determined, so both sit at their own ceiling...
    assert r_common["rho_vs_ceiling"] == pytest.approx(1.0, abs=0.05)
    assert r_rare["rho_vs_ceiling"] == pytest.approx(1.0, abs=0.05)
    # ...even though the raw values differ by a lot.
    assert r_common["worst_within_arm_rho"] > 3 * r_rare["worst_within_arm_rho"]
