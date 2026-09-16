"""Pool substitution really removes what it says it removes.

The whole diagnosis rests on ``build_pool`` producing a population that differs
from the empirical one in exactly one named respect. If a family called
``zeros_removed_mean_matched`` shifted the mean, or one called
``skew_matched_no_zeros`` did not match the skew, the percentage points
attributed to zero inflation and to skew would be attributing them to something
else. So the families are checked against their own names.
"""

from __future__ import annotations

import numpy as np
import pytest

from src.harness.interval_diagnosis import (
    DESIGNS,
    POOL_FAMILIES,
    DiagnosisError,
    arm_sizes,
    attribute,
    build_pool,
    check_design_effect_vanishes_where_the_arms_coincide,
    check_gaussian_floor_matches_its_arithmetic,
    closure_verdict,
    diagnosis_closes,
    gaussian_floor,
    normal_quantile_shortfall,
    null_rejection,
)


def _zero_inflated(rng, n=4000):
    """A pool shaped like a target gene's counts: mostly zero, long tail."""
    detected = rng.binomial(1, 0.25, size=n)
    return (detected * rng.negative_binomial(0.4, 0.05, size=n)).astype(float)


# --------------------------------------------------------------------------
# each family is what its name says
# --------------------------------------------------------------------------


def test_zeros_removed_keeps_the_mean_and_drops_the_zeros():
    rng = np.random.default_rng(1)
    values = _zero_inflated(rng)
    assert (values == 0).mean() > 0.5, "fixture must actually be zero inflated"
    sampler = build_pool(values, "zeros_removed_mean_matched")
    drawn = sampler.draw(200_000, np.random.default_rng(2))
    assert (drawn > 0).all(), "the zero atom survived the family that removes it"
    assert sampler.moments.mean == pytest.approx(values.mean(), rel=1e-12)
    assert drawn.mean() == pytest.approx(values.mean(), rel=0.03)


def test_skew_matched_family_matches_all_three_moments_or_says_it_did_not():
    rng = np.random.default_rng(3)
    values = _zero_inflated(rng)
    sampler = build_pool(values, "skew_matched_no_zeros")
    drawn = sampler.draw(400_000, np.random.default_rng(4))
    assert (drawn > 0).all(), "a strictly positive family produced non-positive draws"
    assert drawn.mean() == pytest.approx(values.mean(), rel=0.05)
    assert drawn.var() == pytest.approx(values.var(), rel=0.20)
    if sampler.moments.skew_matched:
        target = float(
            np.mean((values - values.mean()) ** 3) / values.var() ** 1.5
        )
        realised = float(np.mean((drawn - drawn.mean()) ** 3) / drawn.var() ** 1.5)
        assert realised == pytest.approx(target, rel=0.25)
    else:
        # the honest fallback: it must SAY it only matched mean and variance
        assert sampler.moments.form == "gamma_mean_var_only"


def test_gaussian_family_has_no_skew_and_no_atom():
    rng = np.random.default_rng(5)
    values = _zero_inflated(rng)
    sampler = build_pool(values, "gaussian_matched")
    drawn = sampler.draw(100_000, np.random.default_rng(6))
    skew = float(np.mean((drawn - drawn.mean()) ** 3) / drawn.var() ** 1.5)
    assert abs(skew) < 0.05
    assert drawn.mean() == pytest.approx(values.mean(), abs=0.05 * values.std())
    assert (drawn == 0).sum() == 0


def test_empirical_family_draws_only_values_the_pool_contains():
    rng = np.random.default_rng(7)
    values = _zero_inflated(rng)
    drawn = build_pool(values, "empirical").draw(5000, np.random.default_rng(8))
    assert np.isin(drawn, values).all()


def test_undefined_families_return_none_rather_than_a_substitute():
    """A family that cannot be built has no rate. It must not invent one."""
    constant = np.zeros(50)
    assert build_pool(constant, "zeros_removed_mean_matched") is None
    assert build_pool(constant, "skew_matched_no_zeros") is None
    assert build_pool(constant, "gaussian_matched") is None
    assert build_pool(np.array([1.0]), "empirical") is None


def test_unknown_family_and_design_are_refused():
    with pytest.raises(ValueError, match="unknown family"):
        build_pool(np.arange(10.0), "no_zeros_please")
    with pytest.raises(ValueError, match="unknown design"):
        arm_sizes(50, "whatever")


def test_designs_coincide_at_the_reference_count():
    assert arm_sizes(800, "fixed_fraction") == arm_sizes(800, "balanced") == (800, 800)
    assert arm_sizes(50, "fixed_fraction") == (800, 50)
    assert arm_sizes(50, "balanced") == (50, 50)


# --------------------------------------------------------------------------
# the measurement itself
# --------------------------------------------------------------------------


def test_two_sample_floor_reduces_to_the_one_sample_form_on_equal_arms():
    """At n_n = n_t the two-sample shrink must be sqrt((n-1)/n) at df 2(n-1)."""
    from scipy import stats as st

    for n in (10, 50, 800):
        z = st.norm.ppf(0.975)
        expected = 2 * st.t.sf(z * np.sqrt((n - 1) / n), 2 * (n - 1))
        assert gaussian_floor(n, n) == pytest.approx(expected, rel=1e-12)


def test_two_sample_floor_is_large_at_five_cells_and_negligible_at_800():
    """Where the closed form IS most of the story, and where it is none of it.

    This is the number that decides whether the reviewer's objection stands. At
    (800, 5) the plug-in-variance and z-vs-t term alone gives 15.3% against a
    nominal 5%; at (800, 50) it gives 5.7% and at (800, 800) 5.0%. So the
    arithmetic covers a good share of the five-cell rate and essentially none of
    the 12-30% at 50 or the 9.5% at 800 -- which is exactly why the paper cannot
    answer the reviewer by pointing at it.
    """
    assert gaussian_floor(800, 5) == pytest.approx(0.1533, abs=0.001)
    assert gaussian_floor(800, 50) == pytest.approx(0.0573, abs=0.0005)
    assert gaussian_floor(800, 800) == pytest.approx(0.0503, abs=0.0005)


def test_closure_verdict_separates_over_explaining_from_leaving_it_open():
    """Summing to more than the excess is not success."""
    assert closure_verdict({"excess_over_nominal": 0.25,
                            "share_unexplained": 0.02}) == "closes"
    assert closure_verdict({"excess_over_nominal": 0.25,
                            "share_unexplained": 0.8}) == "leaves_most_unexplained"
    assert closure_verdict({"excess_over_nominal": 0.25,
                            "share_unexplained": -1.3}) == "over_explained"
    assert closure_verdict({"excess_over_nominal": 0.001,
                            "share_unexplained": 5.0}) == "no_excess_to_explain"
    assert diagnosis_closes({"excess_over_nominal": 0.25,
                             "share_unexplained": -1.3}) is False


def test_gaussian_floor_lands_near_its_closed_form():
    """The apparatus, checked against arithmetic that contains no data."""
    rng = np.random.default_rng(20260915)
    pools = [rng.normal(4.0, 2.0, size=500) for _ in range(400)]
    row = null_rejection(
        pools, count=200, family="gaussian_matched", design="balanced",
        n_boot=800, seeds=list(range(400)),
    )
    assert row["n_scored"] == 400
    assert row["null_rejection"] < 0.10, row["null_rejection"]
    check_gaussian_floor_matches_its_arithmetic([row], tolerance_mcse=6.0)


def test_null_rejection_refuses_mismatched_pools_and_seeds():
    with pytest.raises(ValueError, match="one pool per"):
        null_rejection(
            [np.arange(10.0)], count=5, family="empirical", design="balanced",
            n_boot=10, seeds=[1, 2],
        )


def test_normal_quantile_shortfall_is_a_floor_not_an_explanation():
    """The term interval_calibration quantifies is small at these cell counts."""
    assert normal_quantile_shortfall(50) == pytest.approx(0.966, abs=0.002)
    assert normal_quantile_shortfall(800) == pytest.approx(0.99786, abs=0.0002)


# --------------------------------------------------------------------------
# attribution arithmetic
# --------------------------------------------------------------------------


def _rows(p0, p1, p2, p3, balanced):
    base = dict(n_cells_mature=50, n_replicates=2000, n_scored=2000,
                n_pool_undefined=0, median_ci_width=1.0, null_rejection_mcse=0.005)
    return [
        base | dict(family="empirical", design="fixed_fraction", null_rejection=p0,
                    n_normal=800, n_tumour=50),
        base | dict(family="zeros_removed_mean_matched", design="fixed_fraction",
                    null_rejection=p1, n_normal=800, n_tumour=50),
        base | dict(family="skew_matched_no_zeros", design="fixed_fraction",
                    null_rejection=p2, n_normal=800, n_tumour=50),
        base | dict(family="gaussian_matched", design="fixed_fraction",
                    null_rejection=p3, n_normal=800, n_tumour=50),
        base | dict(family="empirical", design="balanced", null_rejection=balanced,
                    n_normal=50, n_tumour=50),
    ]


def test_attribution_sums_to_the_excess_when_the_causes_are_additive():
    a = attribute(_rows(0.20, 0.12, 0.11, 0.058, 0.17), count=50)
    assert a["excess_over_nominal"] == pytest.approx(0.15)
    assert a["pp_zero_inflation"] == pytest.approx(0.08)
    assert a["pp_skew_given_no_zeros"] == pytest.approx(0.052)
    assert a["pp_fixed_fraction_design"] == pytest.approx(0.03)
    assert a["pp_floor_z_vs_t"] == pytest.approx(0.008)
    assert a["pp_residual_unexplained"] == pytest.approx(
        0.15 - (0.08 + 0.052 + 0.03 + 0.008)
    )


def test_attribution_reports_an_open_residual_rather_than_hiding_it():
    """A cause nobody named must land in the residual, not in a named bucket.

    Here removing the zero atom buys 1pp and removing the skew buys 4.2pp, but
    replacing the discrete pool with a continuous one buys 19pp. That 19pp is
    discreteness and the moments above the third -- not one of the reviewer's
    three candidates -- so it must show up as unexplained rather than be
    absorbed into the skew number.
    """
    a = attribute(_rows(0.30, 0.29, 0.10, 0.058, 0.30), count=50)
    assert a["pp_zero_inflation"] == pytest.approx(0.01)
    assert a["pp_skew_given_no_zeros"] == pytest.approx(0.042)
    assert a["pp_discreteness_and_higher_moments"] == pytest.approx(0.19)
    assert a["share_unexplained"] > 0.5
    assert diagnosis_closes(a) is False


def test_attribution_refuses_a_missing_family():
    rows = [r for r in _rows(0.2, 0.12, 0.11, 0.058, 0.17)
            if r["family"] != "skew_matched_no_zeros"]
    with pytest.raises(DiagnosisError, match="no rate for"):
        attribute(rows, count=50)


def test_design_guard_passes_where_the_designs_genuinely_coincide():
    rows = [
        dict(family="empirical", design=d, n_cells_mature=800, n_normal=800,
             n_tumour=800, null_rejection=r, null_rejection_mcse=0.005)
        for d, r in (("fixed_fraction", 0.071), ("balanced", 0.069))
    ]
    check_design_effect_vanishes_where_the_arms_coincide(rows)


def test_every_family_and_design_name_is_covered_by_a_test():
    """If the frozen lists grow, this file must grow with them."""
    assert POOL_FAMILIES == (
        "empirical", "zeros_removed_mean_matched",
        "skew_matched_no_zeros", "gaussian_matched",
    )
    assert DESIGNS == ("fixed_fraction", "balanced")
