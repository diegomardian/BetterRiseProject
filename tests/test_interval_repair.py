"""The candidate intervals, and the shortcuts they are allowed to take.

The load-bearing test here is ``test_affine_shortcut_equals_the_real_estimator``.
``interval_repair`` computes ``f_n * (mean_t - mean_n)`` inline instead of
calling ``kitagawa.decompose`` once per bootstrap draw, because at B=2000 by
R=2000 by 32 settings that call is the entire runtime. A shortcut around the
real estimator is exactly the thing that goes stale without anyone noticing, so
it is pinned by assertion rather than by comment.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy import stats

from src.estimator.kitagawa import decompose
from src.harness.interval import within_patient_intrinsic_ci
from src.harness.interval_repair import (
    ABSTAIN,
    CANDIDATES,
    ArmResamples,
    _jackknife_centred,
    all_intervals,
    bca_interval,
    covers,
    delta_method_se,
    excludes_zero,
    intrinsic_point,
    mcse,
    patient_cluster_interval,
    percentile_interval,
    pseudobulk_patient_t_interval,
    studentised_interval,
    welch_t_interval,
    width,
)


# --------------------------------------------------------------------------
# the shortcut is pinned to the real estimator
# --------------------------------------------------------------------------


def test_affine_shortcut_equals_the_real_estimator():
    """``intrinsic_point`` must be ``decompose(...).intrinsic``, exactly.

    If W4 changes the normal-weighted intrinsic term, this fails here rather
    than the repair study silently ceasing to describe the estimator it claims
    to calibrate.
    """
    rng = np.random.default_rng(20260915)
    for _ in range(200):
        f_n, f_t = rng.uniform(0.01, 0.9, size=2)
        m_n, m_t = rng.normal(0, 10, size=2)
        expected = decompose(
            f_n, f_t, m_n, m_t, n_cells_mature=7, weighting="normal"
        ).intrinsic
        got = intrinsic_point(m_n, m_t, frac_mature_normal=f_n)
        assert got == pytest.approx(expected, rel=0, abs=0), (
            "the affine shortcut has drifted from decompose(). Either W4 "
            "changed the normal-weighted intrinsic term or the shortcut was "
            "edited; either way interval_repair is no longer measuring the "
            "estimator."
        )


def test_percentile_candidate_reproduces_the_committed_interval():
    """R0 must be the status quo, not a reimplementation that resembles it.

    ``ArmResamples`` draws the normal arm then the tumour arm from one
    generator, in that order, which is what ``within_patient_intrinsic_ci``
    does. Same seed, same n_boot, same numbers -- otherwise the baseline the
    candidates are compared against is not the interval the paper reports.
    """
    rng = np.random.default_rng(11)
    normal = rng.poisson(3.0, size=120).astype(float)
    tumour = rng.poisson(3.0, size=45).astype(float)
    f_n, f_t, seed, n_boot = 0.4, 0.0225, 4242, 500

    expected = within_patient_intrinsic_ci(
        normal, tumour,
        frac_mature_normal=f_n, frac_mature_tumour=f_t,
        n_boot=n_boot, seed=seed, weighting="normal",
    )
    shared = np.random.default_rng(seed)
    got = percentile_interval(
        ArmResamples(normal, n_boot=n_boot, rng=shared),
        ArmResamples(tumour, n_boot=n_boot, rng=shared),
        frac_mature_normal=f_n,
    )
    np.testing.assert_allclose(got, expected, rtol=0, atol=1e-12)


def test_jackknife_closed_form_matches_the_naive_loop():
    """The O(n) form must equal ``interval_calibration``'s O(n^2) one."""
    rng = np.random.default_rng(7)
    values = rng.gamma(0.3, 5.0, size=60)
    jack = np.array([np.delete(values, i).mean() for i in range(values.size)])
    naive = jack.mean() - jack
    np.testing.assert_allclose(_jackknife_centred(values), naive, rtol=1e-12)


# --------------------------------------------------------------------------
# each candidate does what its name says
# --------------------------------------------------------------------------


def test_studentised_is_not_the_percentile_reflected():
    """The bootstrap-t sign convention: the UPPER pivot sets the LOWER bound.

    Getting this backwards yields an interval reflected about the estimate --
    plausible-looking, wrong, and invisible on symmetric data. So the check runs
    on a deliberately skewed arm, where the two conventions disagree, and pins
    the direction against a hand-computed bound.
    """
    rng = np.random.default_rng(3)
    normal = rng.gamma(0.2, 20.0, size=400)
    tumour = rng.gamma(0.2, 20.0, size=40)
    f_n, n_boot, alpha = 0.4, 4000, 0.05

    shared = np.random.default_rng(99)
    boot_n = ArmResamples(normal, n_boot=n_boot, rng=shared)
    boot_t = ArmResamples(tumour, n_boot=n_boot, rng=shared)
    lo, hi = studentised_interval(
        normal, tumour, boot_n, boot_t, frac_mature_normal=f_n, alpha=alpha
    )

    theta = intrinsic_point(normal.mean(), tumour.mean(), frac_mature_normal=f_n)
    se = delta_method_se(normal, tumour, frac_mature_normal=f_n)
    draws = f_n * (boot_t.means - boot_n.means)
    se_star = f_n * np.sqrt(boot_t.variances / 40 + boot_n.variances / 400)
    pivots = (draws - theta) / se_star
    q_lo, q_hi = np.percentile(pivots, [2.5, 97.5])

    assert lo == pytest.approx(theta - q_hi * se)
    assert hi == pytest.approx(theta - q_lo * se)
    # and the pivot distribution really is asymmetric here, so the test has
    # something to detect
    assert abs(q_hi + q_lo) > 0.05


def test_welch_t_matches_scipy_on_the_unscaled_difference():
    """R3 is Welch's t, scaled. Checked against scipy's own degrees of freedom."""
    rng = np.random.default_rng(5)
    normal = rng.normal(10.0, 3.0, size=200)
    tumour = rng.normal(10.0, 8.0, size=25)
    lo, hi = welch_t_interval(normal, tumour, frac_mature_normal=1.0)

    res = stats.ttest_ind(tumour, normal, equal_var=False)
    ref_lo, ref_hi = res.confidence_interval(0.95)
    assert lo == pytest.approx(ref_lo, rel=1e-9)
    assert hi == pytest.approx(ref_hi, rel=1e-9)


def test_bca_reduces_to_percentile_when_there_is_nothing_to_correct():
    """With z0 ~ 0 and acceleration ~ 0 BCa must not move the bounds much."""
    rng = np.random.default_rng(13)
    normal = rng.normal(0.0, 1.0, size=500)
    tumour = rng.normal(0.0, 1.0, size=500)
    shared = np.random.default_rng(21)
    boot_n = ArmResamples(normal, n_boot=4000, rng=shared)
    boot_t = ArmResamples(tumour, n_boot=4000, rng=shared)
    pct = percentile_interval(boot_n, boot_t, frac_mature_normal=0.4)
    bca = bca_interval(normal, tumour, boot_n, boot_t, frac_mature_normal=0.4)
    np.testing.assert_allclose(bca, pct, rtol=0.12)


def test_pseudobulk_is_paired_over_the_shared_patients():
    """R5 pairs on patient; an unpaired t would not see a constant offset."""
    patients = np.repeat(["p1", "p2", "p3"], 10)
    level = np.repeat([0.0, 100.0, 200.0], 10)
    normal = level + 1.0
    tumour = level + 1.5
    lo, hi = pseudobulk_patient_t_interval(
        normal, tumour, patients, patients, frac_mature_normal=1.0
    )
    # every patient's difference is exactly 0.5, so the paired interval is a
    # point at 0.5. An unpaired comparison would drown in the 0..200 level.
    assert lo == pytest.approx(0.5) and hi == pytest.approx(0.5) or (lo, hi) == ABSTAIN
    assert lo is None or (lo <= 0.5 <= hi and hi - lo < 1.0)


def test_cluster_and_pseudobulk_abstain_on_one_patient():
    """Fewer than two clusters is not a small sample; it is no sample."""
    cells = np.arange(20.0)
    one = np.repeat(["p1"], 20)
    assert patient_cluster_interval(
        cells, cells, one, one,
        frac_mature_normal=0.4, rng=np.random.default_rng(0), n_boot=50,
    ) == ABSTAIN
    assert pseudobulk_patient_t_interval(
        cells, cells, one, one, frac_mature_normal=0.4
    ) == ABSTAIN


def test_abstention_is_none_and_never_zero():
    """CLAUDE.md invariant 1, at the level of a single interval."""
    out = all_intervals(
        np.empty(0), np.arange(10.0), frac_mature_normal=0.4, seed=1, n_boot=50
    )
    assert set(out) == set(CANDIDATES)
    for name, interval in out.items():
        assert interval == ABSTAIN, name
        assert excludes_zero(interval) is None, name
        assert covers(interval, 0.0) is None, name
        assert np.isnan(width(interval)), name


def test_cluster_candidates_abstain_rather_than_raise_without_labels():
    rng = np.random.default_rng(2)
    out = all_intervals(
        rng.normal(size=40), rng.normal(size=30),
        frac_mature_normal=0.4, seed=3, n_boot=100,
    )
    assert out["patient_cluster"] == ABSTAIN
    assert out["pseudobulk_patient_t"] == ABSTAIN
    assert out["percentile"][0] is not None


def test_all_intervals_shares_one_set_of_resamples():
    """The standalone call and the bundled call must see identical draws."""
    rng = np.random.default_rng(31)
    normal, tumour = rng.poisson(2.0, 300).astype(float), rng.poisson(2.0, 60).astype(float)
    bundled = all_intervals(
        normal, tumour, frac_mature_normal=0.4, seed=808, n_boot=300,
        candidates=("percentile", "studentised", "bca"),
    )
    shared = np.random.default_rng(808)
    boot_n = ArmResamples(normal, n_boot=300, rng=shared)
    boot_t = ArmResamples(tumour, n_boot=300, rng=shared)
    np.testing.assert_allclose(
        bundled["percentile"],
        percentile_interval(boot_n, boot_t, frac_mature_normal=0.4),
        rtol=0, atol=1e-12,
    )
    np.testing.assert_allclose(
        bundled["bca"],
        bca_interval(normal, tumour, boot_n, boot_t, frac_mature_normal=0.4),
        rtol=0, atol=1e-12,
    )


def test_unknown_candidate_is_refused():
    with pytest.raises(ValueError, match="unknown candidate"):
        all_intervals(
            np.arange(5.0), np.arange(5.0),
            frac_mature_normal=0.4, seed=1, candidates=("percentile", "jackknife"),
        )


# --------------------------------------------------------------------------
# calibration sanity: on Gaussian data the t-interval must hit nominal
# --------------------------------------------------------------------------


@pytest.mark.parametrize("n_tumour", [50, 400])
def test_welch_t_is_calibrated_on_gaussian_arms(n_tumour):
    """The floor case. If this is not ~5% the apparatus is wrong, not the data.

    400 trials, so the binomial MCSE is 1.1pp; the assertion band is +/- 3.5pp
    and is deliberately loose -- this is a smoke test for the harness, and the
    real precision comes from the 2000-replicate study.
    """
    rng = np.random.default_rng(20260915)
    rejects = 0
    trials = 400
    for _ in range(trials):
        normal = rng.normal(5.0, 2.0, size=800)
        tumour = rng.normal(5.0, 2.0, size=n_tumour)
        interval = welch_t_interval(normal, tumour, frac_mature_normal=0.4)
        rejects += bool(excludes_zero(interval))
    rate = rejects / trials
    assert abs(rate - 0.05) < 0.035, f"Welch-t rejected {rate:.1%} on Gaussian arms"


def test_mcse_matches_the_binomial_formula():
    assert mcse(0.05, 2000) == pytest.approx(np.sqrt(0.05 * 0.95 / 2000))
    assert mcse(0.095, 200) == pytest.approx(np.sqrt(0.095 * 0.905 / 200))
    assert np.isnan(mcse(0.05, 0))
