"""D2's pre-specified survival-design guards, without fitting an outcome model."""

from __future__ import annotations

from copy import deepcopy

import numpy as np
import pandas as pd
import pytest

from src.bulk.covariates import load_covariate_set
from src.bulk.d2_feasibility_gate import GENE
from src.bulk.d2_survival import (
    D2SurvivalError,
    d2_events_per_df,
    d2_spec,
    prepare_design,
    require_d2_locked,
    retained_vs_dropped,
)


def _design(n: int = 120) -> pd.DataFrame:
    return pd.DataFrame({
        "patient_id": [f"P{i}" for i in range(n)],
        "PFI": [1] * n,
        "GUCA2A": np.linspace(0.0, 4.0, n),
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
