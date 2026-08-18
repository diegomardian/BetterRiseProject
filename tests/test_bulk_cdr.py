"""W3.5 — the curated clinical table.

Two things here are easy to get wrong in ways that never raise:

1. **TCGA's bracketed sentinels.** ``[Not Available]`` left in place becomes a
   stage level and a Cox model happily fits a coefficient for it.
2. **Global versus per-endpoint exclusion.** A patient with no ``tumor_status``
   has no DSS but a perfectly good PFI. Dropping them everywhere shrinks every
   analysis to satisfy the strictest one, and nothing complains.

Both get explicit tests, plus the arithmetic check that exclusion reasons are
mutually exclusive and account for every patient.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.bulk.cdr import (
    ENDPOINTS,
    SENTINELS,
    add_usability_flags,
    build_curated_table,
    clean_sentinels,
    cohort_reconciliation,
    endpoint_exclusion_reason,
    event_summary,
    reconciliation,
    stage_disagreement,
)


def _cdr_frame() -> pd.DataFrame:
    """Six patients, one per exclusion path plus two clean ones."""
    return pd.DataFrame(
        {
            "patient_id": [f"TCGA-AA-000{i}" for i in range(6)],
            "type": ["COAD", "COAD", "READ", "COAD", "READ", "COAD"],
            "age_at_initial_pathologic_diagnosis": [55, 61, 70, 48, 80, 66],
            "gender": ["MALE", "FEMALE", "MALE", "FEMALE", "MALE", "FEMALE"],
            "ajcc_pathologic_tumor_stage": [
                "Stage IIA", "Stage IIIB", "Stage IV", "Stage I", None, "Stage IIC",
            ],
            "treatment_outcome_first_course": [
                "Complete Remission/Response", None, "Progressive Disease", None, None, None,
            ],
            "vital_status": ["Alive", "Dead", "Dead", "Alive", "Alive", "Dead"],
            "tumor_status": [
                "TUMOR FREE", "WITH TUMOR", None, "TUMOR FREE", "TUMOR FREE", "WITH TUMOR",
            ],
            "Redaction": [None, None, None, "Redacted", None, None],
            # clean, clean, DSS undefined, redacted, zero time, missing time
            "OS": [0, 1, 1, 0, 0, 1],
            "OS.time": [400.0, 900.0, 300.0, 500.0, 0.0, np.nan],
            "DSS": [0, 1, np.nan, 0, 0, 1],
            "DSS.time": [400.0, 900.0, 300.0, 500.0, 0.0, np.nan],
            "PFI": [0, 1, 1, 0, 0, 1],
            "PFI.time": [400.0, 900.0, 300.0, 500.0, 0.0, np.nan],
            "DFI": [0, np.nan, np.nan, 0, 0, np.nan],
            "DFI.time": [400.0, np.nan, np.nan, 500.0, 0.0, np.nan],
        }
    )


# ---------------------------------------------------------------------------
# Sentinels
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("sentinel", SENTINELS)
def test_every_sentinel_becomes_nan(sentinel):
    """Left alone these become categorical levels — a "[Not Available]" stage
    group that a Cox model will happily fit a coefficient for."""
    cleaned = clean_sentinels(pd.Series(["Stage I", sentinel, "Stage IV"]))
    assert cleaned.isna().sum() == 1
    assert sentinel not in set(cleaned.dropna())


def test_clean_sentinels_leaves_numeric_columns_alone():
    numeric = pd.Series([1.0, 2.0, np.nan])
    pd.testing.assert_series_equal(clean_sentinels(numeric), numeric)


# ---------------------------------------------------------------------------
# Exclusion rules, one test per rule
# ---------------------------------------------------------------------------


def test_redacted_is_excluded_from_every_endpoint():
    row = pd.Series({"Redaction": "Redacted", "OS": 1, "OS.time": 500.0})
    assert endpoint_exclusion_reason(row, "OS") == "redacted"


def test_missing_event_indicator_is_excluded_not_treated_as_censored():
    """A missing DSS means tumour status was unknown, so "died of disease" is
    undefined. Reading it as 0 would silently add a censored patient."""
    row = pd.Series({"Redaction": None, "DSS": np.nan, "DSS.time": 500.0})
    assert endpoint_exclusion_reason(row, "DSS") == "no_event_indicator"


def test_zero_followup_is_excluded_visibly():
    """22 real COAD/READ patients have exactly 0 days. Most Cox
    implementations drop them silently; here the count is reported."""
    row = pd.Series({"Redaction": None, "OS": 0, "OS.time": 0.0})
    assert endpoint_exclusion_reason(row, "OS") == "nonpositive_followup_time"


def test_missing_followup_time_is_excluded():
    row = pd.Series({"Redaction": None, "OS": 1, "OS.time": np.nan})
    assert endpoint_exclusion_reason(row, "OS") == "no_followup_time"


def test_a_clean_patient_has_no_exclusion_reason():
    row = pd.Series({"Redaction": None, "OS": 1, "OS.time": 500.0})
    assert endpoint_exclusion_reason(row, "OS") is None


# ---------------------------------------------------------------------------
# Per-endpoint, not per-patient
# ---------------------------------------------------------------------------


def test_a_patient_missing_dss_still_contributes_to_pfi():
    """THE structural test. Patient 2 has no tumour_status so no DSS, but their
    PFI is fine. A global drop would lose them from both."""
    flagged = add_usability_flags(_cdr_frame()).set_index("patient_id")
    assert not flagged.loc["TCGA-AA-0002", "usable_DSS"]
    assert flagged.loc["TCGA-AA-0002", "exclusion_DSS"] == "no_event_indicator"
    assert flagged.loc["TCGA-AA-0002", "usable_PFI"]
    assert flagged.loc["TCGA-AA-0002", "exclusion_PFI"] is None


def test_redacted_patient_is_unusable_everywhere():
    flagged = add_usability_flags(_cdr_frame()).set_index("patient_id")
    for endpoint in ENDPOINTS:
        assert not flagged.loc["TCGA-AA-0003", f"usable_{endpoint}"]
        assert flagged.loc["TCGA-AA-0003", f"exclusion_{endpoint}"] == "redacted"


def test_exclusion_reasons_are_mutually_exclusive_and_total():
    """Every patient appears exactly once per endpoint, so the reconciliation
    adds up. If reasons could double-count, the table would overstate drops."""
    flagged = add_usability_flags(_cdr_frame())
    recon = reconciliation(flagged)
    for endpoint in ENDPOINTS:
        sub = recon.loc[recon["endpoint"] == endpoint]
        assert int(sub["n"].sum()) == len(flagged), endpoint


def test_event_summary_counts_only_usable_patients():
    flagged = add_usability_flags(_cdr_frame())
    summary = event_summary(flagged).set_index("endpoint")
    usable_pfi = int(flagged["usable_PFI"].sum())
    assert summary.loc["PFI", "n_usable"] == usable_pfi
    assert summary.loc["PFI", "n_events"] + summary.loc["PFI", "n_censored"] == usable_pfi


def test_endpoint_roles_follow_invariant_9():
    """DSS and PFI primary, OS secondary. Changing this needs a PR and two
    approvals, so it is pinned."""
    assert ENDPOINTS["DSS"][0] == "primary"
    assert ENDPOINTS["PFI"][0] == "primary"
    assert ENDPOINTS["OS"][0] == "secondary"
    assert ENDPOINTS["DFI"][0] == "not_used"


# ---------------------------------------------------------------------------
# The curated table
# ---------------------------------------------------------------------------


def test_curated_table_harmonises_stage_and_keeps_the_raw_value():
    curated = build_curated_table(add_usability_flags(_cdr_frame())).set_index("patient_id")
    assert curated.loc["TCGA-AA-0000", "stage"] == "II"
    assert curated.loc["TCGA-AA-0000", "stage_raw"] == "Stage IIA"
    assert curated.loc["TCGA-AA-0005", "stage"] == "II"  # Stage IIC
    assert pd.isna(curated.loc["TCGA-AA-0004", "stage"])


def test_curated_table_takes_site_and_msi_from_the_gdc_pull():
    """The CDR carries neither. Re-deriving them here would make a second
    source of truth for variables that already have one."""
    gdc = pd.DataFrame(
        {
            "patient_id": ["TCGA-AA-0000", "TCGA-AA-0001"],
            "site": ["right_colon", "rectum"],
            "msi_status": ["MSI", "MSS"],
            "n_distinct_calls": [1, 1],
            "stage": ["II", "III"],
        }
    )
    curated = build_curated_table(add_usability_flags(_cdr_frame()), gdc).set_index("patient_id")
    assert curated.loc["TCGA-AA-0000", "site"] == "right_colon"
    assert curated.loc["TCGA-AA-0000", "msi_status"] == "MSI"
    assert pd.isna(curated.loc["TCGA-AA-0002", "site"])


def test_stage_disagreement_is_surfaced_not_merged_over():
    gdc = pd.DataFrame(
        {
            "patient_id": ["TCGA-AA-0000", "TCGA-AA-0001"],
            "site": ["right_colon", "rectum"],
            "msi_status": ["MSI", "MSS"],
            "n_distinct_calls": [1, 1],
            "stage": ["IV", "III"],  # first one disagrees with the CDR's II
        }
    )
    curated = build_curated_table(add_usability_flags(_cdr_frame()), gdc)
    disagree = stage_disagreement(curated)
    assert list(disagree["patient_id"]) == ["TCGA-AA-0000"]
    assert disagree.iloc[0]["stage"] == "II"
    assert disagree.iloc[0]["stage_gdc"] == "IV"


def test_sex_is_normalised_from_the_cdr_uppercase():
    curated = build_curated_table(add_usability_flags(_cdr_frame()))
    assert set(curated["sex"]) == {"Male", "Female"}


def test_cohort_reconciliation_counts_both_directions():
    curated = build_curated_table(add_usability_flags(_cdr_frame()))
    recon = cohort_reconciliation(
        curated, ["TCGA-AA-0000", "TCGA-AA-0001", "TCGA-NOT-INCDR"]
    ).set_index("set")
    assert recon.loc["in both", "n"] == 2
    assert recon.loc["expression but no CDR record", "n"] == 1
    assert recon.loc["CDR but no expression", "n"] == 4
