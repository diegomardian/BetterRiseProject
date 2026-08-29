"""W2's ratification evidence for prereg Amendment 2. Issue #37.

These tests ARE the adjudication. If one fails, the answer W2 gave the team on
#37 no longer holds and the gate memo's §0 needs rewriting — that is the point of
putting it in the suite rather than in a scratch script.

Deterministic by fixed seed (CLAUDE.md invariant 10), so the pass-rate bands are
exact facts about a stated construction, not flaky thresholds.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.harness.g1_amendment import (
    G1_WORLDS,
    MIN_SEPARATION,
    TIER_A_MAX_MEDIAN,
    TIER_D_MIN_PERCENTILE,
    amendment2_verdict,
    ma_transform,
    null_pass_rate,
    run_world,
    simulate_arms,
    simulate_compartment_world,
    within_bin_percentile,
)

SEED = 20260827
N_REPLICATES = 60


@pytest.fixture(scope="module")
def worlds() -> pd.DataFrame:
    """Every world, replicated. One fixture so the suite pays for it once."""
    return pd.DataFrame(
        [
            run_world(world, seed=SEED, replicate=r)
            for world in G1_WORLDS
            for r in range(N_REPLICATES)
        ]
    )


# ---------------------------------------------------------------------------
# The ratification: can the proposed G1 both pass and fail, for the right reasons
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("world", G1_WORLDS, ids=lambda w: w.name)
def test_every_world_lands_in_its_pass_rate_band(world, worlds):
    """The whole answer to #37, world by world.

    Three worlds must almost never pass, one must almost always pass, and one —
    ``isolated_tier_a_loss`` — is a coin flip, which is the finding rather than
    an accident. See ``AMENDMENT_2_POWER_CAVEAT``.
    """
    rows = worlds[worlds["world"] == world.name]
    assert len(rows) == N_REPLICATES
    rate = float((rows["verdict"] == "PASS").mean())
    low, high = world.expected_pass_rate
    assert low <= rate <= high, f"{world.name}: P(PASS)={rate:.3f} outside [{low}, {high}]"


def test_the_superseded_statistic_fails_a_null_and_the_proposed_one_does_not(worlds):
    """Decision #17's ρ ≈ −1 on uniform proportional loss, checked not accepted.

    W1's claim, reproduced independently through W2's own construction: on data
    with no abundance dependence anywhere, the pre-registered statistic returns
    a near-perfect correlation, so G1 could never have passed.
    """
    null = worlds[worlds["world"] == "uniform_loss"]
    assert null["rho_decision_17"].mean() < -0.95
    assert abs(null["rho_g1a"].mean()) < 0.10


def test_the_proposed_statistic_detects_abundance_dependence_when_it_is_real(worlds):
    """The check the amendment did NOT do — that the replacement can fire.

    Showing the old statistic fails a null establishes the old one is broken. It
    says nothing about the new one. On a world where loss is a function of
    abundance and nothing else, G1a must see it.
    """
    soup = worlds[worlds["world"] == "pure_soup"]
    assert abs(soup["rho_g1a"].mean()) > 0.50
    # ...and it must not fire on the world where the signal is genuine biology.
    real = worlds[worlds["world"] == "broad_loss_tier_d_retained"]
    assert abs(real["rho_g1a"].mean()) < 0.50


def test_the_retained_control_separates_only_against_a_lost_background(worlds):
    """The power caveat, as a number rather than a warning.

    MS4A12 sits at ~0.89 when there is broad loss to be retained against, and at
    ~0.50 — a coin flip on threshold 2 — when nothing else moves.
    """
    broad = worlds[worlds["world"] == "broad_loss_tier_d_retained"]["tier_d_pct"]
    isolated = worlds[worlds["world"] == "isolated_tier_a_loss"]["tier_d_pct"]
    assert broad.mean() > 0.80
    assert 0.40 < isolated.mean() < 0.60


# ---------------------------------------------------------------------------
# The thresholds themselves
# ---------------------------------------------------------------------------


def test_the_thresholds_are_w1s_numbers_and_w2_has_not_tuned_them():
    """Testing whether a pre-committed threshold discriminates is legitimate.

    Adjusting it to make the ratification come out well is the move the whole
    amendment exists to avoid, so it fails a test rather than passing review.
    """
    assert (TIER_A_MAX_MEDIAN, TIER_D_MIN_PERCENTILE, MIN_SEPARATION) == (0.20, 0.50, 0.30)


def test_null_pass_rate_shows_threshold_two_is_a_coin_flip():
    """Under a null the gate passes 1 time in 40 — but tier D alone is 50/50.

    The joint rate is respectable; the per-threshold rates are where the
    falsification logic actually lives, and threshold 2 carries none of it on
    its own.
    """
    rates = null_pass_rate(seed=7, n_trials=100_000)
    assert rates["p_threshold_2_tier_d"] == pytest.approx(0.50, abs=0.01)
    assert rates["p_threshold_1_tier_a"] == pytest.approx(0.051, abs=0.01)
    assert rates["p_g1_passes_on_a_null"] == pytest.approx(0.026, abs=0.01)


def test_a_nan_percentile_is_refused_rather_than_compared():
    """The defect that would have let decision #17 pass tier D silently.

    ``scipy.stats.spearmanr`` on one observation returns nan without raising,
    and ``abs(nan) > 0.5`` is False — so "fail if |ρ| > 0.5" would have reported
    PASS for the one tier carrying the falsification logic. Refuse, don't
    compare.
    """
    with pytest.raises(ValueError, match="nan"):
        amendment2_verdict(np.array([0.1, 0.1, 0.1, 0.1]), float("nan"))
    with pytest.raises(ValueError, match="nan"):
        amendment2_verdict(np.array([0.1, np.nan, 0.1, 0.1]), 0.9)


def test_verdict_refuses_an_empty_tier():
    with pytest.raises(ValueError, match="not computable"):
        amendment2_verdict(np.array([]), 0.9)


def test_verdict_names_every_threshold_it_failed():
    result = amendment2_verdict(np.array([0.9, 0.9, 0.9, 0.9]), 0.1)
    assert result["verdict"] == "FAIL"
    assert len(result["failures"]) == 3


# ---------------------------------------------------------------------------
# The construction
# ---------------------------------------------------------------------------


def test_within_bin_percentile_is_flat_under_uniform_proportional_loss():
    """Amendment 2 §4's own claim: the median within-bin percentile is 0.500."""
    rng = np.random.default_rng(11)
    abundance = 10 ** rng.uniform(-2, 2, 20_000)
    normal, tumour = simulate_arms(
        abundance, np.full(20_000, 0.7), n_cells=3_000, rng=rng
    )
    a, m, _ = ma_transform(normal, tumour)
    assert float(np.median(within_bin_percentile(a, m))) == pytest.approx(0.5, abs=0.01)


def test_percentiles_survive_library_size_normalisation():
    """A strength of the rank construction, worth stating because it is not obvious.

    A global rescale of one arm shifts every M by the same constant, and the
    percentile is a within-bin rank, so it does not move. G1 therefore does not
    depend on whether the arms were CP10K-normalised — which matters, because a
    broad loss is partly re-centred by normalisation and a mean-based rule would
    have been sensitive to it.
    """
    rng = np.random.default_rng(3)
    abundance = 10 ** rng.uniform(-2, 2, 20_000)
    fold = np.exp(rng.normal(-0.55, 0.45, 20_000))
    normal, tumour = simulate_arms(abundance, fold, n_cells=3_000, rng=rng)

    a_raw, m_raw, keep_raw = ma_transform(normal, tumour)
    a_cp, m_cp, keep_cp = ma_transform(
        normal / normal.sum() * 1e4, tumour / tumour.sum() * 1e4
    )
    assert np.array_equal(keep_raw, keep_cp)
    moved = np.abs(within_bin_percentile(a_raw, m_raw) - within_bin_percentile(a_cp, m_cp))
    assert moved.max() < 0.05


def test_ma_transform_drops_rather_than_floors_a_zero_arm():
    """A log ratio against zero is not a number, and a pseudocount invents one."""
    normal = np.array([1.0, 0.0, 4.0])
    tumour = np.array([2.0, 3.0, 0.0])
    a, m, keep = ma_transform(normal, tumour)
    assert list(keep) == [True, False, False]
    assert len(m) == 1


def test_within_bin_percentile_refuses_a_panel_sized_input():
    """The rule is genome-wide. Eleven genes cannot populate twenty bins."""
    with pytest.raises(ValueError, match="genome-wide"):
        within_bin_percentile(np.arange(11.0), np.arange(11.0))


def test_simulate_arms_rejects_a_cohort_with_no_cells():
    rng = np.random.default_rng(0)
    with pytest.raises(ValueError, match="n_cells"):
        simulate_arms(np.array([1.0]), np.array([1.0]), n_cells=0, rng=rng)


# ---------------------------------------------------------------------------
# The gap the ratification had — issue #46
# ---------------------------------------------------------------------------


def test_the_ratification_worlds_have_no_compartment_structure():
    """Why W2's ratification passed a criterion that fails when we are right.

    Every world in ``G1_WORLDS`` varies per-gene fold change against a flat
    background, so a gene's mean over all epithelium is its mean everywhere.
    The failure mode in #46 lives entirely in the difference between those two,
    and a harness that cannot represent a confound cannot rule it out.

    This test does not assert a fix. It pins the limitation so it cannot be
    rediscovered as news.
    """
    for world in G1_WORLDS:
        rng = np.random.default_rng(0)
        abundance = 10 ** rng.uniform(-2, 2, 500)
        panel = np.array([10, 20, 30, 40, 50])
        fold = world.fold_change(abundance, panel, 500, rng)
        # A compartment world would need a per-cell-type structure; a plain
        # per-gene fold vector cannot express one.
        assert np.asarray(fold).shape == (500,)


def test_g1_fails_the_world_where_the_hypothesis_is_true():
    """#46, reproduced through W2's own harness rather than accepted from W1.

    Mature cells deplete and tier A is silenced inside the survivors — the
    project's claim, exactly. G1 is owed a PASS and returns FAIL, because
    MS4A12 is mature-restricted and sinks with its compartment while its
    abundance-matched comparison set does not.
    """
    result = simulate_compartment_world(seed=7, tier_a_silencing=0.15)
    assert result["hypothesis_is_true"]
    assert result["verdict"] == "FAIL"
    assert result["tier_d_pct"] < TIER_D_MIN_PERCENTILE
    assert any("MS4A12" in f for f in result["failures"])


def test_g1_also_fails_the_composition_only_world_so_it_cannot_discriminate():
    """The two worlds owe different verdicts and G1 returns the same one.

    That is the whole argument for withdrawal in docs/g1_withdrawal_case.md: a
    criterion returning FAIL either way carries no information, and its
    pre-committed consequence must not fire on it.
    """
    true_world = simulate_compartment_world(seed=7, tier_a_silencing=0.15)
    comp_only = simulate_compartment_world(seed=7, tier_a_silencing=1.0)
    assert comp_only["verdict"] == true_world["verdict"] == "FAIL"
    # Tier A does separate; it is threshold 2 that sinks both.
    assert true_world["tier_a_median_pct"] < comp_only["tier_a_median_pct"]
    assert true_world["tier_d_pct"] == pytest.approx(comp_only["tier_d_pct"], abs=0.05)
