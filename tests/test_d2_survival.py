"""D2's pre-specified survival-design guards, without fitting an outcome model."""

from __future__ import annotations

from copy import deepcopy

import numpy as np
import pandas as pd
import pytest

from src.bulk.covariates import load_covariate_set
from src.bulk.d2_feasibility_gate import GENE
from src.bulk.d2_survival import (
    D2_DESCRIPTIVE_CONTEXT,
    D2SurvivalError,
    d2_events_per_df,
    d2_model_spec,
    d2_spec,
    prepare_design,
    require_d2_locked,
    retained_vs_dropped,
)
from src.bulk.run_d2_survival import MODEL_PRIMARY, d2_verdict


def _design(n: int = 120) -> pd.DataFrame:
    return pd.DataFrame({
        "patient_id": [f"P{i}" for i in range(n)],
        "PFI": [1] * n,
        "GUCA2A": np.linspace(0.0, 4.0, n),
        "project": ["TCGA-COAD"] * n,
        "plate": ["PL0"] * n,
    })


def _clinical(n: int = 6) -> pd.DataFrame:
    return pd.DataFrame({
        "patient_id": [f"P{i}" for i in range(n)],
        "project": ["TCGA-COAD"] * n,
        "stage": ["I", "II", "III", "IV", "I", "II"][:n],
        "age": np.arange(50, 50 + n),
        "sex": ["Female", "Male"] * (n // 2),
        "msi_status": ["MSS"] * n,
        "site": ["left_colon"] * n,
        "PFI": [1, 0, 1, 0, 1, 0][:n],
        "PFI.time": np.arange(100, 100 + n),
        "usable_PFI": [True] * n,
    })


def _manifest(n: int = 6, *, plate_n: int = 6) -> pd.DataFrame:
    return pd.DataFrame({
        "patient_id": [f"P{i}" for i in range(n)],
        "barcode": [f"B{i}" for i in range(n)],
        "sample_type": ["01"] * n,
        "plate": [f"PL{i}" if i < plate_n else pd.NA for i in range(n)],
    })


def _purity(n: int = 6) -> pd.DataFrame:
    return pd.DataFrame({
        "patient_id": [f"P{i}" for i in range(n)],
        "method": ["absolute"] * n,
        "purity": np.linspace(0.4, 0.8, n),
        "sample_type": ["01"] * n,
    })


def _expression(n: int = 6) -> pd.DataFrame:
    return pd.DataFrame({"ENSG_GUCA2A": np.linspace(0.0, 4.0, n)},
                        index=[f"B{i}" for i in range(n)])


def _purity_config(spec):
    return next(c for c in spec["covariates"] if c["name"] == "purity")


def _primary_as_estimate(spec):
    _purity_config(spec)["source"] = "estimate_affy_extrapolated"


def _sensitivity_as_absolute(spec):
    _purity_config(spec)["sensitivity_source"] = "absolute"


def test_d2_spec_binds_endpoint_roles_purity_sources_covariates_and_strata():
    locked = d2_spec(load_covariate_set())
    assert locked["endpoints"]["PFI"]["role"] == "primary"
    assert locked["endpoints"]["DSS"]["role"] == "primary"
    assert locked["endpoints"]["OS"]["role"] == "secondary"
    purity = next(c for c in locked["covariates"] if c["name"] == "purity")
    assert purity["source"] == "absolute"
    assert purity["sensitivity_source"] == "estimate_affy_extrapolated"
    assert locked["model"]["strata"] == ["project", "plate"]


def test_d2_spec_refuses_a_changed_endpoint_role():
    spec = deepcopy(load_covariate_set())
    spec["endpoints"] = {**spec["endpoints"], "OS": {"role": "primary", "lead": False}}
    with pytest.raises(D2SurvivalError, match="OS must remain"):
        d2_spec(spec)


@pytest.mark.parametrize(("mutate", "message"), [
    (lambda s: s["endpoints"]["PFI"].update(lead=False), "PFI must remain"),
    (_primary_as_estimate, "ABSOLUTE"),
    (_sensitivity_as_absolute, "ESTIMATE"),
    (lambda s: s["contexts"]["expression_models"].update(additional_required=[]), "require plate"),
])
def test_d2_spec_forcing_inputs_trigger_every_non_role_guard(mutate, message):
    spec = deepcopy(load_covariate_set())
    mutate(spec)
    with pytest.raises(D2SurvivalError, match=message):
        d2_spec(spec)


def test_d2_prereg_must_be_locked(tmp_path):
    path = tmp_path / "d2.md"
    path.write_text("# D2\n\n**Status:** proposed\n", encoding="utf-8")
    with pytest.raises(D2SurvivalError, match="not locked"):
        require_d2_locked(path)


def test_events_per_df_counts_the_continuous_guca2a_predictor():
    stats = d2_events_per_df(_design(), d2_spec(load_covariate_set()), endpoint="PFI")
    assert stats["non_stratum_df"] == 11
    assert stats["meets_lead_floor"] is True
    assert stats["n_nonempty_strata"] == 1
    assert stats["n_event_contributing_strata"] == 1


def test_prepare_design_records_plate_exclusion_before_later_filters():
    design, attrition = prepare_design(
        _clinical(), _purity(), _manifest(plate_n=4), _expression(), "ENSG_GUCA2A",
        load_covariate_set(), endpoint="PFI",
    )
    assert attrition.iloc[0]["step"] == "clinical table"
    assert attrition.iloc[0]["n"] == 6
    plate = attrition.loc[attrition["step"] == "complete plate"].iloc[0]
    assert plate["n"] == 4
    assert plate["dropped"] == 2
    assert set(design["patient_id"]) == {"P0", "P1", "P2", "P3"}
    assert GENE in design


def test_prepare_design_separates_missing_primary_sample_from_missing_plate():
    clinical = _clinical()
    manifest = _manifest(n=5, plate_n=4)  # P5 has no primary specimen at all.
    design, attrition = prepare_design(
        clinical, _purity(), manifest, _expression(n=5), "ENSG_GUCA2A",
        load_covariate_set(), endpoint="PFI",
    )
    primary = attrition.loc[attrition["step"] == "has primary-tumour RNA sample"].iloc[0]
    plate = attrition.loc[attrition["step"] == "complete plate"].iloc[0]
    assert primary["n"] == 5
    assert primary["dropped"] == 1
    assert plate["n"] == 4
    assert plate["dropped"] == 1
    assert set(design["patient_id"]) == {"P0", "P1", "P2", "P3"}


def test_descriptive_context_omits_purity_but_keeps_the_plate_stratum():
    locked, context = d2_model_spec(load_covariate_set(), descriptive=True)
    assert context == D2_DESCRIPTIVE_CONTEXT
    assert "purity" in locked["contexts"][context]["exclude"]
    # The expression-model stratum is not a generic clinical-baseline model.
    assert locked["model"]["strata"] == ["project", "plate"]
    design, _ = prepare_design(
        _clinical(), _purity(), _manifest(plate_n=4), _expression(), "ENSG_GUCA2A",
        load_covariate_set(), endpoint="PFI", descriptive=True,
    )
    assert "purity" not in design.columns
    assert len(design) == 4


def test_retained_vs_dropped_reports_both_guca2a_distributions():
    clinical = pd.DataFrame({
        "patient_id": ["P1", "P2", "P3", "P4"],
        "usable_PFI": [True, True, True, False],
    })
    values = pd.Series([0.0, 1.0, 3.0, 8.0], index=["P1", "P2", "P3", "P4"], name="GUCA2A")
    design = pd.DataFrame({"patient_id": ["P1", "P3"]})
    got = retained_vs_dropped(clinical, values, design, endpoint="PFI", purity_method="absolute")
    by = got.set_index("membership")
    assert by.loc["retained", "n_participants"] == 2
    assert by.loc["dropped", "n_participants"] == 1
    assert by.loc["dropped", "median"] == 1.0


def _d2_stats(*, estimable=True):
    return pd.DataFrame([{
        "endpoint": "PFI", "model": MODEL_PRIMARY, "purity_method": "absolute",
        "meets_lead_floor": estimable,
    }])


def _d2_coefficients(*, pfi_hr=0.7, pfi_low=0.55, pfi_high=0.9, dss_hr=0.8,
                     estimate_hr=0.75, estimate_low=0.6, estimate_high=0.95):
    return pd.DataFrame([
        {"endpoint": "PFI", "model": MODEL_PRIMARY, "purity_method": "absolute",
         "term": "GUCA2A", "hazard_ratio": pfi_hr,
         "ci_low": pfi_low, "ci_high": pfi_high},
        {"endpoint": "DSS", "model": MODEL_PRIMARY, "purity_method": "absolute",
         "term": "GUCA2A", "hazard_ratio": dss_hr, "ci_low": 0.6, "ci_high": 1.1},
        {"endpoint": "PFI", "model": MODEL_PRIMARY,
         "purity_method": "estimate_affy_extrapolated", "term": "GUCA2A",
         "hazard_ratio": estimate_hr, "ci_low": estimate_low,
         "ci_high": estimate_high},
    ])


def _d2_ph(*, violates=False, endpoint="PFI"):
    return pd.DataFrame([{
        "endpoint": endpoint, "model": MODEL_PRIMARY, "purity_method": "absolute",
        "term": "GUCA2A", "violates_ph": violates,
    }])


def test_d2_verdict_stops_before_coefficients_when_lead_is_not_estimable():
    verdict = d2_verdict(_d2_stats(estimable=False), pd.DataFrame(), pd.DataFrame()).iloc[0]
    assert verdict["verdict"] == "NOT ESTIMABLE"


def test_d2_verdict_requires_a_purity_robust_ph_clean_lead_association():
    verdict = d2_verdict(_d2_stats(), _d2_coefficients(), _d2_ph()).iloc[0]
    assert verdict["verdict"] == "SUPPORTED ASSOCIATION"

    source_sensitive = d2_verdict(
        _d2_stats(), _d2_coefficients(estimate_hr=1.1, estimate_high=1.4), _d2_ph()
    ).iloc[0]
    assert source_sensitive["verdict"] == "PURITY-SOURCE SENSITIVE"

    ph_failed = d2_verdict(_d2_stats(), _d2_coefficients(), _d2_ph(violates=True)).iloc[0]
    assert ph_failed["verdict"] == "PH VIOLATED"


def test_a_null_interval_above_one_is_not_a_contrary_direction():
    """§6's rules are ordered: interval-includes-1 is read BEFORE the sign.

    HR 1.04 with a 95% interval of [0.88, 1.23] is a null result. Reporting it
    as DIRECTION CONTRARY turns an absence of evidence into a finding, in the
    one direction this pre-registration cannot afford.
    """
    verdict = d2_verdict(
        _d2_stats(),
        _d2_coefficients(pfi_hr=1.04, pfi_low=0.88, pfi_high=1.23, dss_hr=1.06,
                         estimate_hr=1.02, estimate_low=0.87, estimate_high=1.20),
        _d2_ph(),
    ).iloc[0]
    assert verdict["verdict"] == "NO SUPPORTED ASSOCIATION"

    contrary = d2_verdict(
        _d2_stats(),
        _d2_coefficients(pfi_hr=1.31, pfi_low=1.09, pfi_high=1.58, dss_hr=1.22,
                         estimate_hr=1.28, estimate_low=1.06, estimate_high=1.55),
        _d2_ph(),
    ).iloc[0]
    assert contrary["verdict"] == "DIRECTION CONTRARY"


def test_a_ph_violation_in_a_secondary_endpoint_does_not_silence_the_lead():
    """§6: OS is descriptive and cannot reverse the primary conclusion.

    Without an endpoint filter on the Schoenfeld table, an OS violation
    returned PH VIOLATED for the whole of D2, under a detail line naming PFI.
    """
    verdict = d2_verdict(
        _d2_stats(), _d2_coefficients(),
        _d2_ph(violates=True, endpoint="OS"),
    ).iloc[0]
    assert verdict["verdict"] == "SUPPORTED ASSOCIATION"
    assert not verdict["ph_violated"]
