"""The §2.2 sweep, bulk recovery and cutpoint calibration.

The sweep is what G3 and the week-5 cutpoints are read off, so the checks here
are about the sweep being trustworthy, not about the curve having any particular
shape. The shape is a result.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.harness.attenuation import (
    DEFAULT_SHIFTS,
    SweepConfig,
    SweepGrid,
    bulk_overconfidence,
    null_arm_noise_ratio,
    null_arm_recovers_zero,
    run_sweep,
    summarise_sweep,
)
from src.harness.bulk_recovery import (
    attenuated_mature_mean,
    attenuation_ratio,
    reference_profiles,
)
from src.harness.calibration import (
    CalibrationCriteria,
    calibrate_cutpoints,
    coverage_and_discrimination,
)
from src.harness.deconvolve import NuSVRDeconvolver
from src.harness.results import validate_harness_table

MATURE = "mature_colonocyte"
TARGET = "GUCA2A"
GENES = [TARGET, "MLH1", "ACTB", "F1", "F2", "F3"]
TYPES = [MATURE, "stem", "stromal", "immune"]

#: Per-cell-type mean of each gene. GUCA2A is mature-restricted, which is the
#: whole reason it is a compositional control.
PROFILE = {
    MATURE: [60.0, 8.0, 50.0, 5.0, 5.0, 5.0],
    "stem": [2.0, 8.0, 50.0, 5.0, 5.0, 5.0],
    "stromal": [1.0, 6.0, 55.0, 6.0, 4.0, 5.0],
    "immune": [1.0, 7.0, 45.0, 4.0, 6.0, 5.0],
}


@pytest.fixture(scope="module")
def cohort():
    """10 patients x 4 types x 60 cells. Enough patients to hold two out."""
    rng = np.random.default_rng(11)
    rows, ctypes, patients = [], [], []
    for p in range(10):
        for t in TYPES:
            for _ in range(60):
                rows.append(rng.poisson(PROFILE[t]))
                ctypes.append(t)
                patients.append(f"P{p:02d}")
    return np.array(rows), np.array(ctypes), np.array(patients)


@pytest.fixture(scope="module")
def small_sweep(cohort):
    counts, ctypes, patients = cohort
    grid = SweepGrid(
        mature_fractions=(0.40, 0.10, 0.02, 0.0),
        shifts=(1.0, 0.5),
        n_replicates=4,
        n_cells=600,
    )
    return run_sweep(
        SweepConfig(counts, ctypes, patients, GENES, TARGET), grid, seed=2026
    )


# ---------------------------------------------------------------------------
# The grid
# ---------------------------------------------------------------------------


def test_grid_must_contain_the_null():
    with pytest.raises(ValueError, match="null arm"):
        SweepGrid(shifts=(0.5, 0.25))


def test_default_shifts_are_the_preregistered_set():
    """docs/harness_design_spec.md §4. Changing these changes the cutpoints."""
    assert DEFAULT_SHIFTS == (1.0, 0.8, 0.5, 0.25)


def test_grid_points_cover_the_cross_product():
    grid = SweepGrid(mature_fractions=(0.4, 0.1), shifts=(1.0, 0.5), n_replicates=1)
    assert len(grid.points()) == 4
    assert len({gid for gid, _, _ in grid.points()}) == 4


# ---------------------------------------------------------------------------
# The sweep
# ---------------------------------------------------------------------------


def test_sweep_is_shaped_like_the_harness_table(small_sweep):
    validate_harness_table(small_sweep, "attenuation")


def test_sweep_has_both_arms_for_every_row(small_sweep):
    assert set(small_sweep["arm"]) == {"oracle", "bulk"}
    counts = small_sweep.groupby(["grid_id", "replicate"])["arm"].nunique()
    assert (counts == 2).all()


def test_null_arm_parametric_truth_is_exactly_zero(small_sweep):
    """If this fails the harness is broken and no other number here means
    anything. It is the first thing to read off a sweep."""
    assert null_arm_recovers_zero(small_sweep)
    null = small_sweep[small_sweep["shift"] == 1.0]
    assert (null["intrinsic_true_parametric"] == 0.0).all()


def test_null_arm_realised_truth_is_noisy_and_that_is_correct(small_sweep):
    """The realised null is NOT zero, and asserting it were would be wrong.

    Normal and tumour samples are different draws, so their empirical means
    differ even when nothing was silenced. What must hold is that the noise is
    small next to a real effect.
    """
    null = small_sweep[
        (small_sweep["shift"] == 1.0) & (small_sweep["arm"] == "oracle")
    ]
    assert (null["intrinsic_true_realised"] != 0.0).any()
    assert null_arm_noise_ratio(small_sweep) < 0.20


def test_sweep_is_reproducible(cohort):
    counts, ctypes, patients = cohort
    grid = SweepGrid(
        mature_fractions=(0.2,), shifts=(1.0, 0.5), n_replicates=2, n_cells=400
    )
    cfg = SweepConfig(counts, ctypes, patients, GENES, TARGET)
    a = run_sweep(cfg, grid, seed=7)
    b = run_sweep(cfg, grid, seed=7)
    pd.testing.assert_frame_equal(a, b)


def test_mature_count_and_estimability_track_the_swept_fraction(small_sweep):
    oracle = small_sweep[small_sweep["arm"] == "oracle"]
    medians = oracle.groupby("frac_mature_tumour")["n_cells_mature"].median()
    assert medians.is_monotonic_increasing
    assert (oracle.loc[oracle["frac_mature_tumour"] == 0.0, "estimability"]
            == "not_estimable").all()


def test_bulk_returns_confident_numbers_where_the_truth_is_undefined(small_sweep):
    """THE finding, measured rather than asserted rhetorically.

    At zero mature cells the intrinsic term is undefined. The oracle arm knows
    that — it counts cells. Bulk cannot count cells: deconvolution assigns a
    non-zero mature fraction anyway, the division goes through, and out comes a
    confident wrong answer. This is the structural reason the third segment
    needs single-cell data, and the reason invariant 6 exists.
    """
    overconfident = bulk_overconfidence(small_sweep)
    assert len(overconfident) > 0, (
        "bulk did not over-report on this sweep — either the grid no longer "
        "reaches an empty mature compartment, or the bulk arm gained a "
        "positivity check it is not supposed to have"
    )
    # The numbers it invents are not small, either.
    assert overconfident["intrinsic_hat"].abs().median() > 1.0


def test_oracle_arm_does_know_the_term_is_undefined(small_sweep):
    """The contrast that makes the finding a finding."""
    gone = small_sweep[
        (small_sweep["frac_mature_tumour"] == 0.0) & (small_sweep["arm"] == "oracle")
    ]
    assert (gone["n_cells_mature"] == 0).all()
    assert (gone["estimability"] == "not_estimable").all()


def test_oracle_arm_tracks_truth_where_the_estimate_is_defined(small_sweep):
    """The oracle arm is the reliable half — this is G3 in miniature."""
    ok = small_sweep[
        (small_sweep["arm"] == "oracle") & (small_sweep["estimability"] == "ok")
    ]
    assert len(ok) > 0
    np.testing.assert_allclose(
        ok["intrinsic_hat"], ok["intrinsic_true_realised"], atol=1e-9
    )


def test_summarise_collapses_replicates(small_sweep):
    summary = summarise_sweep(small_sweep)
    assert set(summary["arm"]) == {"oracle", "bulk"}
    assert (summary["n_replicates"] == 4).all()
    assert "attenuation_median" in summary.columns


def test_sweep_rejects_a_target_not_in_the_gene_list(cohort):
    counts, ctypes, patients = cohort
    grid = SweepGrid(mature_fractions=(0.2,), shifts=(1.0,), n_replicates=1, n_cells=200)
    with pytest.raises(KeyError, match="not in the gene list"):
        run_sweep(SweepConfig(counts, ctypes, patients, GENES, "NOTAGENE"), grid, seed=1)


# ---------------------------------------------------------------------------
# Bulk recovery — the part invariant 6 forbids using
# ---------------------------------------------------------------------------


def test_reference_profiles_can_exclude_the_target(cohort):
    counts, ctypes, _ = cohort
    sig = reference_profiles(counts, ctypes, GENES, exclude_genes=[TARGET])
    assert TARGET not in sig.index
    assert set(sig.columns) == set(TYPES)
    # GUCA2A is mature-restricted, so it must dominate there when included.
    full = reference_profiles(counts, ctypes, GENES)
    assert full.loc[TARGET, MATURE] > 10 * full.loc[TARGET, "stem"]


def test_attenuated_mean_recovers_the_truth_when_everything_is_known():
    """With exact fractions and an exact reference, the algebra is exact. Any
    attenuation in the sweep therefore comes from estimation error, not here."""
    n_cells = 1000
    fractions = pd.Series({MATURE: 0.25, "stem": 0.25, "stromal": 0.25, "immune": 0.25})
    profile = pd.Series({t: PROFILE[t][0] for t in TYPES})
    bulk = np.zeros(len(GENES))
    bulk[0] = n_cells * sum(fractions[t] * profile[t] for t in TYPES)
    got = attenuated_mature_mean(
        bulk, GENES, gene=TARGET, fractions=fractions,
        n_cells=n_cells, target_profile=profile,
    )
    assert got == pytest.approx(PROFILE[MATURE][0])


def test_attenuated_mean_is_undefined_when_the_compartment_is_empty():
    fractions = pd.Series({MATURE: 0.0, "stem": 1.0, "stromal": 0.0, "immune": 0.0})
    got = attenuated_mature_mean(
        np.ones(len(GENES)), GENES, gene=TARGET, fractions=fractions,
        n_cells=100, target_profile=pd.Series({t: 1.0 for t in TYPES}),
    )
    assert np.isnan(got)  # not 0.0


def test_attenuation_ratio_is_undefined_against_a_zero_truth():
    assert np.isnan(attenuation_ratio(0.5, 0.0))
    assert attenuation_ratio(0.6, 1.2) == pytest.approx(0.5)


# ---------------------------------------------------------------------------
# Calibration
# ---------------------------------------------------------------------------


def _sweep_with_cis(n_per_bin=40, seed=0):
    """Synthetic sweep where interval width grows as mature cells run out —
    the behaviour the cutpoints are supposed to detect."""
    rng = np.random.default_rng(seed)
    rows = []
    for n_mature in (5, 15, 30, 60, 150, 400):
        width = 6.0 / np.sqrt(n_mature)
        for rep in range(n_per_bin):
            truth = -1.0
            hat = truth + rng.normal(0, width / 4)
            rows.append(
                {
                    "grid_id": 0, "replicate": rep, "arm": "oracle", "gene": TARGET,
                    "weighting": "normal", "frac_mature_tumour": 0.1, "shift": 0.5,
                    "n_cells_mature": n_mature, "estimability": "ok",
                    "compositional_true_parametric": 0.0,
                    "intrinsic_true_parametric": truth,
                    "compositional_true_realised": 0.0,
                    "intrinsic_true_realised": truth,
                    "compositional_hat": 0.0, "intrinsic_hat": hat,
                    "interaction_hat": 0.0, "attenuation_ratio": hat / truth,
                    "ci_low": hat - width, "ci_high": hat + width, "seed": seed,
                }
            )
    return pd.DataFrame(rows)


def test_calibration_needs_confidence_intervals(small_sweep):
    """run_sweep leaves CIs null — attaching them is what turns a sweep into a
    calibration, and doing it silently would be worse than refusing."""
    with pytest.raises(ValueError, match="no confidence intervals"):
        calibrate_cutpoints(small_sweep)


def test_calibration_refuses_a_shift_other_than_the_registered_one():
    sweep = _sweep_with_cis()
    criteria = CalibrationCriteria(detectable_shift=0.25)
    with pytest.raises(ValueError, match="pre-registered"):
        calibrate_cutpoints(sweep, criteria)


def test_coverage_and_discrimination_are_computed_per_bin():
    table = coverage_and_discrimination(_sweep_with_cis())
    assert {"coverage", "discrimination", "median_ci_width", "verdict"} <= set(table.columns)
    assert table["n_cells_mature"].is_monotonic_increasing
    # Intervals narrow as cells accumulate, which is the point.
    assert table["median_ci_width"].iloc[0] > table["median_ci_width"].iloc[-1]


def test_calibrated_cutpoints_are_ordered_and_traceable():
    report = calibrate_cutpoints(_sweep_with_cis(), source="unit-test sweep")
    assert report.cutpoints.ok >= report.cutpoints.wide > 0
    assert "unit-test" in report.cutpoints.source
    assert set(report.comparison()["cutpoint"]) == {"ok", "wide"}


def test_calibration_reports_g4_rather_than_inventing_a_cutpoint():
    """No n meets both targets -> that IS the finding, not an error to smooth."""
    hopeless = _sweep_with_cis()
    hopeless["ci_low"] = hopeless["intrinsic_hat"] - 100.0
    hopeless["ci_high"] = hopeless["intrinsic_hat"] + 100.0
    with pytest.raises(ValueError, match="headline result"):
        calibrate_cutpoints(hopeless)


def test_criteria_are_recorded_with_the_result():
    report = calibrate_cutpoints(_sweep_with_cis())
    assert report.criteria.detectable_shift == 0.5
    assert report.criteria.coverage_target == 0.90
    assert report.criteria.discrimination_target == 0.80


# ---------------------------------------------------------------------------
# nu-SVR adapter
# ---------------------------------------------------------------------------


def test_nusvr_recovers_fractions_from_a_clean_mixture():
    rng = np.random.default_rng(3)
    sig = pd.DataFrame(
        rng.gamma(2.0, 5.0, size=(800, 4)),
        columns=TYPES,
        index=[f"G{i:04d}" for i in range(800)],
    )
    true = np.array([0.4, 0.3, 0.2, 0.1])
    bulk = sig.to_numpy() @ true
    hat = NuSVRDeconvolver().fit_predict(bulk, sig)
    assert hat.sum() == pytest.approx(1.0)
    np.testing.assert_allclose(hat.loc[TYPES].to_numpy(), true, atol=0.08)


def test_nusvr_reports_availability():
    ok, reason = NuSVRDeconvolver().is_available()
    assert ok is True and reason == "available"
