"""W3.7 — the baseline models, and the guards around them.

Two tests here exist because of bugs that actually happened:

- :func:`test_ph_check_returns_real_p_values` — the first implementation parsed
  ``CoxPHFitter.check_assumptions``, which returns *matplotlib axes*, not
  results. It reported "no term violates" while lifelines was printing two real
  violations to stdout. A PH check that cannot fail is worse than none.
- :func:`test_reference_levels_come_from_the_config` — letting pandas one-hot
  categoricals picks a reference alphabetically, which would silently make
  stage IV's hazard ratio relative to the wrong baseline.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.bulk.covariates import CovariateError, load_covariate_set
from src.bulk.survival import (
    PH_ALPHA,
    SurvivalError,
    _design_matrix,
    fit_cox,
    kaplan_meier_by_stage,
    power_note,
    proportional_hazards_check,
    stage_sanity_check,
)

RNG = np.random.default_rng(20260818)


@pytest.fixture
def spec():
    return load_covariate_set()


def _cohort(n=320, effect=2.2, seed=0):
    """Synthetic cohort where stage genuinely drives the hazard."""
    rng = np.random.default_rng(seed)
    stages = np.array(["I", "II", "III", "IV"])
    stage = rng.choice(stages, n, p=[0.25, 0.3, 0.28, 0.17])
    rank = pd.Series(stage).map({"I": 0, "II": 1, "III": 2, "IV": 3}).to_numpy()
    hazard = np.exp(effect * rank / 3.0)
    time = rng.exponential(1200.0 / hazard)
    censor = rng.exponential(1800.0, n)
    observed = np.minimum(time, censor)
    event = (time <= censor).astype(int)
    return pd.DataFrame(
        {
            "patient_id": [f"P{i}" for i in range(n)],
            "project": rng.choice(["TCGA-COAD", "TCGA-READ"], n),
            "stage": stage,
            "age": rng.normal(65, 11, n).round(),
            "sex": rng.choice(["Male", "Female"], n),
            "msi_status": rng.choice(["MSS", "MSI"], n, p=[0.86, 0.14]),
            "site": rng.choice(
                ["right_colon", "left_colon", "rectum", "colon_unspecified"], n
            ),
            "PFI": event,
            "PFI.time": np.maximum(observed, 1.0),
            "DSS": event,
            "DSS.time": np.maximum(observed, 1.0),
        }
    )


# ---------------------------------------------------------------------------
# The lock binds at the point of fitting
# ---------------------------------------------------------------------------


def test_fitting_against_an_unlocked_set_is_refused(spec):
    """The covariate lock has to bite where it matters — the fit, not just the
    loader."""
    proposed = dict(spec, status="proposed")
    with pytest.raises(CovariateError, match="not 'locked'"):
        fit_cox(_cohort(), proposed, endpoint="PFI")


# ---------------------------------------------------------------------------
# Design matrix
# ---------------------------------------------------------------------------


def test_reference_levels_come_from_the_config(spec):
    """stage's reference is I, so there must be no stage[I] column and there
    must be one for every other level."""
    cohort = _cohort()
    frame = _design_matrix(cohort, ["stage"], spec)
    assert "stage[I]" not in frame.columns
    assert set(frame.columns) == {"stage[II]", "stage[III]", "stage[IV]"}


def test_sex_reference_is_female_not_alphabetical_chance(spec):
    frame = _design_matrix(_cohort(), ["sex"], spec)
    assert list(frame.columns) == ["sex[Male]"]


def test_an_unexpected_level_is_refused_loudly(spec):
    """Better to stop than to silently drop a level the config never saw —
    that would quietly change the cohort."""
    cohort = _cohort()
    cohort.loc[0, "stage"] = "Stage 0"
    with pytest.raises(SurvivalError, match="not in the locked set"):
        _design_matrix(cohort, ["stage"], spec)


def test_continuous_covariates_pass_through_numerically(spec):
    frame = _design_matrix(_cohort(), ["age"], spec)
    assert frame["age"].dtype.kind == "f"


# ---------------------------------------------------------------------------
# The sanity check itself
# ---------------------------------------------------------------------------


def test_stage_comes_out_prognostic_on_a_cohort_where_it_is(spec):
    """The positive control for the whole task."""
    _, tidy = fit_cox(_cohort(), spec, endpoint="PFI")
    verdict = stage_sanity_check(tidy, endpoint="PFI")
    assert verdict["passed"]
    assert verdict["stage_IV_hazard_ratio"] > 1.5
    assert verdict["verdict"] == "stage is prognostic"


def test_the_sanity_check_fails_when_stage_carries_no_signal(spec):
    """The negative control, and the one that makes a PASS mean something. With
    stage unrelated to the hazard, the check must not report success."""
    cohort = _cohort(n=300, effect=0.0, seed=7)
    _, tidy = fit_cox(cohort, spec, endpoint="PFI")
    verdict = stage_sanity_check(tidy, endpoint="PFI")
    assert not verdict["passed"]
    assert "investigate" in verdict["verdict"].lower()


def test_the_sanity_check_refuses_a_model_without_stage():
    tidy = pd.DataFrame({"term": ["age"], "hazard_ratio": [1.0], "ci_low": [0.9],
                         "ci_high": [1.1], "p": [0.5]})
    with pytest.raises(SurvivalError, match="stage is the sanity check"):
        stage_sanity_check(tidy, endpoint="PFI")


def test_monotonicity_is_reported_but_not_required(spec):
    _, tidy = fit_cox(_cohort(), spec, endpoint="PFI")
    verdict = stage_sanity_check(tidy, endpoint="PFI")
    assert isinstance(verdict["hazard_increases_with_stage"], bool)


# ---------------------------------------------------------------------------
# Proportional hazards — the regression test
# ---------------------------------------------------------------------------


def test_ph_check_returns_real_p_values(spec):
    """REGRESSION. The first version parsed check_assumptions(), which returns
    matplotlib axes — so every term silently read as passing. Assert we get
    finite statistics for real terms, under both time transforms."""
    cohort = _cohort()
    fitter, _ = fit_cox(cohort, spec, endpoint="PFI")
    ph = proportional_hazards_check(
        fitter, cohort, spec, endpoint="PFI", context="clinical_baseline"
    )
    assert len(ph) > 0
    assert ph["p"].notna().all()
    assert np.isfinite(ph["test_statistic"]).all()
    assert set(ph["time_transform"]) == {"km", "rank"}
    assert "stage[IV]" in set(ph["term"])
    assert ph["violates_ph"].dtype == bool


def test_ph_check_flags_a_genuine_violation(spec):
    """Crossing hazards must be caught. A check that never fires is not a check.

    Females get a decreasing hazard (Weibull shape < 1), males an increasing one
    (shape > 1), so the hazard ratio reverses over follow-up — the textbook
    proportional-hazards violation.
    """
    rng = np.random.default_rng(11)
    n = 600
    cohort = _cohort(n=n, seed=11)
    sex = rng.choice(["Female", "Male"], n)
    cohort["sex"] = sex
    shape = np.where(sex == "Female", 0.6, 3.0)
    cohort["PFI.time"] = np.maximum(rng.weibull(shape) * 900.0, 1.0)
    cohort["PFI"] = 1  # uncensored, so the crossing is fully observed
    fitter, _ = fit_cox(cohort, spec, endpoint="PFI")
    ph = proportional_hazards_check(
        fitter, cohort, spec, endpoint="PFI", context="clinical_baseline"
    )
    assert ph.loc[ph["term"] == "sex[Male]", "p"].min() < PH_ALPHA
    assert ph.loc[ph["term"] == "sex[Male]", "violates_ph"].any()


# ---------------------------------------------------------------------------
# Kaplan-Meier and power
# ---------------------------------------------------------------------------


def test_km_reports_one_row_per_stage_with_a_logrank(spec):
    cohort = _cohort()
    km, logrank = kaplan_meier_by_stage(cohort, endpoint="PFI")
    assert set(km["stage"]) == {"I", "II", "III", "IV"}
    assert (km["n"] > 0).all()
    assert 0.0 <= logrank <= 1.0
    assert (km["survival_at_1y"] >= km["survival_at_3y"]).all()


def test_km_survival_falls_with_stage(spec):
    """On a cohort where stage drives the hazard, 3-year survival must fall
    from stage I to stage IV."""
    km, _ = kaplan_meier_by_stage(_cohort(), endpoint="PFI")
    by_stage = km.set_index("stage")["survival_at_3y"]
    assert by_stage["I"] > by_stage["IV"]


def test_power_note_names_the_shortfall(spec):
    cohort = _cohort(n=60)
    note = power_note(cohort, spec, "PFI", "clinical_baseline")
    assert "per df" in note
    assert ("adequate" in note) or ("UNDERPOWERED" in note)


def test_power_note_uses_the_context_df(spec):
    """DSS drops site in the clinical baseline, so it costs fewer df than PFI."""
    cohort = _cohort()
    assert "9 df" in power_note(cohort, spec, "PFI", "clinical_baseline")
    assert "6 df" in power_note(cohort, spec, "DSS", "clinical_baseline")
