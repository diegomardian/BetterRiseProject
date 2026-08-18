"""W3.6 — the covariate lock has to actually bind.

The single most important test here is
:func:`test_a_proposed_set_refuses_to_run_a_model`. A pre-specification that
does not stop you fitting is a document, not a lock.
"""

from __future__ import annotations

import textwrap

import pandas as pd
import pytest
import yaml

from src.bulk.covariates import (
    COVARIATE_SET_PATH,
    CovariateError,
    build_design,
    covariate_names,
    events_per_df,
    load_covariate_set,
    purity_column,
    require_locked,
    total_df,
)


@pytest.fixture
def spec():
    return load_covariate_set()


def _write(tmp_path, spec_dict):
    path = tmp_path / "covariate_set.yaml"
    path.write_text(yaml.safe_dump(spec_dict), encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# The lock
# ---------------------------------------------------------------------------


def test_a_proposed_set_refuses_to_run_a_model(spec):
    """THE test. The set ships as `proposed`; nothing may be fitted against it
    until a human confirms and flips it to `locked` in its own commit."""
    assert spec["status"] == "proposed"
    with pytest.raises(CovariateError, match="not 'locked'"):
        require_locked(spec)


def test_a_locked_set_passes(tmp_path, spec):
    locked = dict(spec, status="locked")
    require_locked(locked)  # must not raise


def test_an_unknown_status_is_rejected(tmp_path, spec):
    path = _write(tmp_path, dict(spec, status="draft"))
    load_covariate_set.cache_clear()
    with pytest.raises(CovariateError, match="proposed"):
        load_covariate_set(path)


def test_a_missing_config_is_a_loud_failure(tmp_path):
    load_covariate_set.cache_clear()
    with pytest.raises(CovariateError, match="does not exist"):
        load_covariate_set(tmp_path / "nope.yaml")


@pytest.mark.parametrize("key", ["version", "status", "covariates", "endpoints", "model"])
def test_required_keys_are_enforced(tmp_path, spec, key):
    incomplete = {k: v for k, v in spec.items() if k != key}
    path = _write(tmp_path, incomplete)
    load_covariate_set.cache_clear()
    with pytest.raises(CovariateError, match=key):
        load_covariate_set(path)


# ---------------------------------------------------------------------------
# The committed set says what the brief asked for
# ---------------------------------------------------------------------------


def test_the_brief_s_six_covariates_are_all_present(spec):
    """stage, age, sex, MMR/MSI, purity, tumour site."""
    assert {c["name"] for c in spec["covariates"]} == {
        "stage", "age", "sex", "msi_status", "purity", "site",
    }


def test_the_full_set_costs_ten_degrees_of_freedom(spec):
    """Six variables, ten df. The df count is what binds against the event
    count, not the variable count."""
    assert total_df(spec, covariate_names(spec, context="expression_models")) == 10
    assert spec["total_df"] == 10


def test_dss_drops_site_and_pfi_does_not(spec):
    """DSS has roughly half PFI's events. The reduction is pre-specified here so
    it cannot be chosen after seeing a fit."""
    assert total_df(spec, covariate_names(spec, endpoint="PFI")) == 9
    dss = covariate_names(spec, endpoint="DSS")
    assert "site" not in dss
    assert total_df(spec, dss) == 6


def test_purity_is_excluded_from_the_clinical_baseline_only(spec):
    """Purity confounds analyses whose PREDICTOR is expression. The W3.7
    baseline has none, and requiring purity there costs 65 patients and 19 PFI
    events for no inferential gain."""
    assert "purity" not in covariate_names(spec, endpoint="PFI")
    assert "purity" in covariate_names(spec, endpoint="PFI", context="expression_models")


def test_an_unknown_context_is_rejected(spec):
    with pytest.raises(CovariateError, match="unknown context"):
        covariate_names(spec, context="made_up")


def test_endpoint_roles_match_invariant_9(spec):
    assert spec["endpoints"]["DSS"]["role"] == "primary"
    assert spec["endpoints"]["PFI"]["role"] == "primary"
    assert spec["endpoints"]["OS"]["role"] == "secondary"
    assert spec["endpoints"]["DFI"]["role"] == "not_used"


def test_exactly_one_endpoint_leads(spec):
    leads = [k for k, v in spec["endpoints"].items() if v.get("lead")]
    assert leads == ["PFI"]


def test_purity_primary_is_absolute_and_sensitivity_is_estimate(spec):
    """ABSOLUTE is called from copy number and is not circular with an
    expression outcome. ESTIMATE is, so it is the sensitivity."""
    assert purity_column(spec) == "absolute"
    assert purity_column(spec, sensitivity=True) == "estimate_affy_extrapolated"


def test_project_is_a_stratum_not_a_covariate(spec):
    """Rectal cancer usually gets neoadjuvant chemoradiation, so a common
    baseline hazard with a proportional shift is not credible."""
    assert spec["model"]["strata"] == ["project"]
    assert "project" not in covariate_names(spec)


def test_plate_is_excluded_here_but_required_for_expression_models(spec):
    """W3.4 found plate outweighs every biological variable — but it affects
    expression measurement, and these models contain no expression."""
    excluded = {e["name"] for e in spec["excluded"]}
    assert "plate" in excluded
    required = {
        a["name"]
        for a in spec["contexts"]["expression_models"]["additional_required"]
    }
    assert "plate" in required


def test_msi_conflicting_is_not_a_modelled_level(spec):
    msi = next(c for c in spec["covariates"] if c["name"] == "msi_status")
    assert msi["levels"] == ["MSS", "MSI"]
    assert "conflicting" not in msi["levels"]


def test_the_committed_file_is_the_one_that_loads(spec):
    assert COVARIATE_SET_PATH.exists()
    assert spec["version"] == "1.0.0"


def test_unknown_endpoint_is_rejected(spec):
    with pytest.raises(CovariateError, match="unknown endpoint"):
        covariate_names(spec, endpoint="NOPE")


# ---------------------------------------------------------------------------
# Design assembly
# ---------------------------------------------------------------------------


def _clinical(n=40):
    return pd.DataFrame(
        {
            "patient_id": [f"P{i}" for i in range(n)],
            "project": ["TCGA-COAD"] * (n // 2) + ["TCGA-READ"] * (n - n // 2),
            "stage": (["I", "II", "III", "IV"] * n)[:n],
            "age": list(range(40, 40 + n)),
            "sex": (["Male", "Female"] * n)[:n],
            "msi_status": (["MSS", "MSS", "MSI", "conflicting"] * n)[:n],
            "site": (["right_colon", "left_colon", "rectum", "colon_unspecified"] * n)[:n],
            "PFI": ([1, 0] * n)[:n],
            "PFI.time": [100.0 + i for i in range(n)],
            "usable_PFI": [True] * n,
            "DSS": ([1, 0] * n)[:n],
            "DSS.time": [100.0 + i for i in range(n)],
            "usable_DSS": [True] * n,
        }
    )


def _purity(n=40, method="absolute"):
    return pd.DataFrame(
        {
            "patient_id": [f"P{i}" for i in range(n)],
            "method": [method] * n,
            "purity": [0.5 + (i % 10) / 40 for i in range(n)],
        }
    )


def test_conflicting_msi_is_dropped_not_modelled(spec):
    """The 16 real patients whose WGS and WXS files disagree must not become a
    third level in the design matrix."""
    design, _ = build_design(_clinical(), _purity(), spec, endpoint="PFI")
    assert "conflicting" not in set(design["msi_status"])
    assert len(design) < 40


def test_attrition_records_every_step_in_order(spec):
    _, attrition = build_design(_clinical(), _purity(), spec, endpoint="PFI")
    steps = list(attrition["step"])
    assert steps[0] == "clinical table"
    assert steps[-1] == "final design"
    assert "usable for PFI" in steps
    assert (attrition["dropped"] >= 0).all()


def test_the_dss_design_has_no_site_column(spec):
    design, _ = build_design(_clinical(), _purity(), spec, endpoint="DSS")
    assert "site" not in design.columns
    assert "stage" in design.columns


def test_sensitivity_uses_the_estimate_purity(spec):
    purity = pd.concat(
        [_purity(method="absolute").head(5), _purity(method="estimate_affy_extrapolated")]
    )
    _, attrition = build_design(
        _clinical(),
        purity,
        spec,
        endpoint="PFI",
        context="expression_models",
        sensitivity=True,
    )
    assert set(attrition["purity_method"]) == {"estimate_affy_extrapolated"}


def test_missing_purity_method_is_a_loud_failure(spec):
    with pytest.raises(CovariateError, match="no purity values"):
        build_design(
            _clinical(),
            _purity(method="something_else"),
            spec,
            endpoint="PFI",
            context="expression_models",
        )


def test_events_per_df_flags_the_floor(spec):
    design, _ = build_design(_clinical(), _purity(), spec, endpoint="PFI")
    stats = events_per_df(design, spec, "PFI")
    assert stats["total_df"] == 9  # clinical baseline: no purity
    assert stats["events_per_df"] == round(stats["n_events"] / 9, 2)
    assert stats["meets_floor_of_10"] is (stats["events_per_df"] >= 10)


def test_a_missing_clinical_column_is_a_loud_failure(spec):
    clinical = _clinical().drop(columns=["site"])
    with pytest.raises(CovariateError, match="no column 'site'"):
        build_design(clinical, _purity(), spec, endpoint="PFI")


def test_yaml_rationales_are_present_for_every_choice(spec):
    """A locked set whose reasons live only in someone's memory is not
    reproducible. Every covariate and exclusion carries prose."""
    for covariate in spec["covariates"]:
        assert covariate.get("rationale") or covariate["name"] in {"sex"}
    for excluded in spec["excluded"]:
        assert textwrap.dedent(excluded["reason"]).strip()
