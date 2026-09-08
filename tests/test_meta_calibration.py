"""The inverse-variance weight, its moments, and the null the ceiling never had.

Pre-registration: ``docs/prereg_meta_weight_calibration.md``.

The pass paths are here; the guards' failing inputs are in
``tests/test_checks_can_fail.py``, including the one for the defect this job's
own first draft shipped — a floor comparison that compared verdict LABELS when
``premise_verdict`` reaches UNRESOLVED down two different routes.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.harness.meta import HETEROGENEITY_ALPHA
from src.reference.jobs.meta_floor_sensitivity import (
    NULL_ALPHA,
    floor_curve,
    influence,
    read_verdict,
    verdict_cause,
)
from src.reference.meta_calibration import (
    FLOORS_REPORTED,
    PATIENT_FLOOR,
    STATUS_QUO_PATIENT_FLOOR,
    WEAK_PATIENT_FLOOR,
    MetaCalibrationError,
    calibrated_p,
    min_patients_for_finite_weight_mean,
    min_patients_for_finite_weight_variance,
    null_heterogeneity_table,
    null_i_squared,
    weight_inflation,
)

#: The eleven patient counts of the committed ICBI meta.
ICBI_N = (20, 3, 11, 29, 12, 11, 15, 5, 3, 4, 5)


# --------------------------------------------------------------------------
# The closed form
# --------------------------------------------------------------------------

@pytest.mark.parametrize("n,expected", [(4, 3.0), (5, 2.0), (6, 5 / 3),
                                        (10, 9 / 7), (20, 19 / 17), (29, 28 / 26)])
def test_weight_inflation_matches_the_closed_form(n, expected):
    assert weight_inflation(n) == pytest.approx(expected)


@pytest.mark.parametrize("n", [2, 3])
def test_weight_has_no_finite_mean_at_or_below_three_patients(n):
    """``E[1/chi2(v)]`` diverges for ``v <= 2``. Not "large" -- absent.

    This is the whole pre-registration in one assertion: the committed floor is
    3, and 3 is inside the region where the quantity the estimator averages has
    no average.
    """
    assert weight_inflation(n) == float("inf")
    assert n < min_patients_for_finite_weight_mean()


def test_a_standard_error_below_two_patients_is_refused_not_infinite():
    with pytest.raises(MetaCalibrationError, match="at least 2"):
        weight_inflation(1)


def test_the_two_floors_are_the_chi_square_df_conditions():
    assert min_patients_for_finite_weight_mean() == 4      # df > 2
    assert min_patients_for_finite_weight_variance() == 6  # df > 4
    assert WEAK_PATIENT_FLOOR == min_patients_for_finite_weight_mean()
    assert PATIENT_FLOOR == min_patients_for_finite_weight_variance()
    assert STATUS_QUO_PATIENT_FLOOR == 3
    assert FLOORS_REPORTED == (3, 4, 6)


def test_closed_form_agrees_with_simulation_where_the_mean_exists():
    """A closed form nobody checked against draws is an assertion, not a result."""
    rng = np.random.default_rng(4)
    for n in (5, 8, 15):
        draws = 1.0 / (rng.chisquare(n - 1, 400_000) / (n - 1))
        assert draws.mean() == pytest.approx(weight_inflation(n), rel=0.02)


# --------------------------------------------------------------------------
# The null the ceiling was compared against without anyone computing it
# --------------------------------------------------------------------------

def test_null_i_squared_is_far_from_zero_at_the_committed_patient_counts():
    """Under EXACT homogeneity, I^2 at the ICBI n's has a median near 0.27.

    The generative model has zero between-study variance by construction, so
    every bit of this is manufactured by estimating the weights.
    """
    null = null_i_squared(ICBI_N, n_trials=40_000)
    assert null.k == 11
    assert 0.22 < null.median < 0.32
    assert null.q95 > 0.75                   # homogeneity alone clears the old ceiling
    assert 0.04 < null.p_exceeds_ceiling < 0.10


def test_cochrans_q_rejects_far_above_nominal_at_these_patient_counts():
    null = null_i_squared(ICBI_N, n_trials=40_000)
    assert null.q_rejection_rate > 0.25      # nominal is 0.05


def test_the_heterogeneity_is_manufactured_by_the_small_studies():
    """Give all eleven studies the largest study's n and the null collapses."""
    inflated = null_i_squared(ICBI_N, n_trials=40_000)
    even = null_i_squared([29] * 11, n_trials=40_000)
    assert even.median == pytest.approx(0.0, abs=1e-9)
    assert even.p_exceeds_ceiling < 0.01
    assert inflated.median > 0.2


def test_null_shrinks_monotonically_as_the_floor_rises():
    table = null_heterogeneity_table(ICBI_N, n_trials=20_000)
    medians = table.sort_values("patient_floor")["null_i_squared_median"].tolist()
    assert medians == sorted(medians, reverse=True)
    assert table.loc[table.patient_floor == 3, "weight_mean_finite"].item() is False
    assert table.loc[table.patient_floor == 6, "weight_variance_finite"].item() is True


def test_null_is_deterministic_under_its_seed():
    """Invariant 10. Two calls at one seed are one number."""
    a = null_i_squared(ICBI_N, n_trials=5_000, seed=99)
    b = null_i_squared(ICBI_N, n_trials=5_000, seed=99)
    assert a == b


def test_null_refuses_fewer_studies_than_meta_analyse_would_accept():
    with pytest.raises(MetaCalibrationError, match="below meta.MIN_STUDIES"):
        null_i_squared([10, 10], n_trials=1_000)


def test_calibrated_p_never_returns_exactly_zero():
    """(r+1)/(B+1). A finite simulation cannot license p = 0."""
    assert calibrated_p(1.0, ICBI_N, n_trials=2_000) > 0.0


# --------------------------------------------------------------------------
# The floor curve, against the committed numbers
# --------------------------------------------------------------------------

def _committed_per_study() -> pd.DataFrame:
    from src.common.paths import RESULTS_DIR
    path = RESULTS_DIR / "2026-09-05_61ba221/icbi_coexpression_meta_per_study.parquet"
    if not path.exists():
        pytest.skip(f"committed meta table absent: {path}")
    return pd.read_parquet(path)


def test_the_status_quo_floor_reproduces_the_committed_meta_exactly():
    """A re-read that cannot reproduce the number it re-reads has not been checked."""
    per_study = _committed_per_study()
    curve = floor_curve(per_study, n_trials=2_000)
    krt8 = curve[(curve.gene == "KRT8") & (curve.patient_floor == 3)].iloc[0]
    assert krt8["k"] == 11
    assert krt8["i_squared"] == pytest.approx(0.876418, abs=1e-5)
    assert krt8["cochran_q"] == pytest.approx(80.918, abs=1e-2)
    assert krt8["verdict"] == "UNRESOLVED"
    actb = curve[(curve.gene == "ACTB") & (curve.patient_floor == 3)].iloc[0]
    assert actb["i_squared"] == pytest.approx(0.627662, abs=1e-5)
    assert actb["verdict"] == "HOLDS"


def test_every_row_carries_its_own_null():
    curve = floor_curve(_committed_per_study(), n_trials=2_000)
    stated = curve[curve.estimability == "estimated"]
    assert len(stated) == len(FLOORS_REPORTED) * 2
    assert stated["null_i_squared_median"].notna().all()
    assert stated["null_p_of_observed"].notna().all()


def test_the_influence_table_is_labelled_exploratory():
    """It names the most influential study. Nothing may read it as a licence."""
    infl = influence(_committed_per_study())
    assert infl["exploratory"].all()


# --------------------------------------------------------------------------
# The defect this job's own first draft shipped
# --------------------------------------------------------------------------

def _row(**kw) -> pd.Series:
    base = {"estimability": "estimated", "homogeneous": True, "verdict": "HOLDS"}
    return pd.Series({**base, **kw})


def test_verdict_cause_separates_the_two_routes_into_unresolved():
    over_null = _row(homogeneous=False, verdict="UNRESOLVED")
    straddles = _row(homogeneous=True, verdict="UNRESOLVED")
    assert verdict_cause(over_null) == "heterogeneity_over_calibrated_null"
    assert verdict_cause(straddles) == "homogeneous_straddles_tolerance"
    assert verdict_cause(over_null) != verdict_cause(straddles)


def test_read_verdict_does_not_call_a_changed_cause_stable():
    """The defect, as a test.

    KRT8 is UNRESOLVED at every floor -- at n>=3 because its calibrated null
    rejects homogeneity, at n>=6 because the studies AGREE the control moved and
    the interval straddles the tolerance. Those are opposite findings under one
    word. The first version of ``read_verdict`` compared labels and reported
    STABLE, which is a check that could not fail on the input it existed to catch.
    """
    frame = pd.DataFrame([
        {"gene": "KRT8", "patient_floor": 3, "estimability": "estimated",
         "homogeneous": False, "verdict": "UNRESOLVED"},
        {"gene": "KRT8", "patient_floor": 4, "estimability": "estimated",
         "homogeneous": True, "verdict": "UNRESOLVED"},
        {"gene": "KRT8", "patient_floor": 6, "estimability": "estimated",
         "homogeneous": True, "verdict": "UNRESOLVED"},
    ])
    frame["verdict_cause"] = frame.apply(verdict_cause, axis=1)
    outcome = read_verdict(frame)
    assert outcome["verdict"] != "STABLE"
    assert "heterogeneity_over_calibrated_null" in outcome["detail"]
    assert "homogeneous_straddles_tolerance" in outcome["detail"]


def test_read_verdict_reports_stable_when_the_cause_really_is_unchanged():
    """The other half: the pass path, so the test above is not vacuous."""
    frame = pd.DataFrame([
        {"gene": "ACTB", "patient_floor": f, "estimability": "estimated",
         "homogeneous": True, "verdict": "HOLDS"}
        for f in FLOORS_REPORTED
    ])
    frame["verdict_cause"] = frame.apply(verdict_cause, axis=1)
    assert read_verdict(frame)["verdict"] == "STABLE"


def test_an_intermediate_floor_cause_cannot_be_hidden_by_matching_endpoints():
    """n>=3 and n>=6 may agree while n>=4 says something else."""
    frame = pd.DataFrame([
        {"gene": "KRT8", "patient_floor": 3, "estimability": "estimated",
         "homogeneous": False, "verdict": "UNRESOLVED"},
        {"gene": "KRT8", "patient_floor": 4, "estimability": "estimated",
         "homogeneous": True, "verdict": "UNRESOLVED"},
        {"gene": "KRT8", "patient_floor": 6, "estimability": "estimated",
         "homogeneous": False, "verdict": "UNRESOLVED"},
    ])
    frame["verdict_cause"] = frame.apply(verdict_cause, axis=1)
    outcome = read_verdict(frame)
    assert outcome["verdict"] != "STABLE"
    assert "n>=4 homogeneous_straddles_tolerance" in outcome["detail"]


def test_the_calibrated_null_can_refuse_an_i_squared_below_the_old_ceiling():
    """At k=6 the calibrated null can reject an I^2 below 75%."""
    curve = floor_curve(_committed_per_study(), n_trials=40_000)
    top = curve[(curve.gene == "KRT8") & (curve.patient_floor == PATIENT_FLOOR)].iloc[0]
    assert top["i_squared"] < 0.75
    assert top["null_p_of_observed"] < NULL_ALPHA
    assert bool(top["heterogeneous_by_null"])
    assert not bool(top["homogeneous"])
    assert top["heterogeneity_alpha"] == HETEROGENEITY_ALPHA
