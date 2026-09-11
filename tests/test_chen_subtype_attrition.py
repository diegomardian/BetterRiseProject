"""Gate 6 compares coverage, never arms, and never zeroes a missing value."""

from __future__ import annotations

import numpy as np
import pandas as pd

from src.reference.jobs.chen_subtype_attrition import (
    RUNG,
    compare,
    coverage_by_patient,
    summarise,
)


def _joined(rows: list[tuple[str, str, bool]]) -> pd.DataFrame:
    """(patient, polyp_type, matched) -> the shape exact_join returns."""
    frame = pd.DataFrame(
        [
            {
                "patient_id": f"Chen.{p}",
                "scRNA_biospecimen_id": f"{p}_2000001011",
                "POLYP_TYPE": t if matched else np.nan,
                "label_matched": matched,
            }
            for p, t, matched in rows
        ]
    )
    return frame


def _decomposition(values: dict[str, float | None]) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "patient_id": f"Chen.{p}", "gene": "GUCA2A", "granularity_rung": RUNG,
                "weighting": "doubly_robust",
                "estimability": "ok" if v is not None else "not_estimable",
                "compositional": 0.0, "intrinsic": v, "interaction": 0.0,
            }
            for p, v in values.items()
        ]
    )


def test_a_patient_with_any_arm_label_counts_as_covered():
    cov = coverage_by_patient(_joined([("A", "AD", True), ("B", "SER", True)]))
    assert cov["covered"].all()


def test_the_conflicted_patient_is_covered_even_though_no_arm_takes_them():
    """They are reached by the label and removed by a rule; that is not attrition."""
    cov = coverage_by_patient(
        pd.concat([_joined([("C", "AD", True)]), _joined([("C", "SER", True)])])
    )
    assert cov.loc[cov.patient_id == "Chen.C", "covered"].all()


def test_unknown_and_absent_are_distinguished_as_reasons():
    cov = coverage_by_patient(
        _joined([("D", "Unknown", True), ("E", "AD", False)])
    ).set_index("patient_id")
    assert cov.loc["Chen.D", "uncovered_reason"] == "label is Unknown"
    assert cov.loc["Chen.E", "uncovered_reason"] == "no cBioPortal record"


def test_a_not_estimable_value_is_excluded_from_the_mean_and_counted():
    """Invariant 1: averaging it in as 0.0 would move the mean toward zero."""
    cov = coverage_by_patient(_joined([("A", "AD", True), ("B", "Unknown", True)]))
    out = compare(_decomposition({"A": 2.0, "B": None}), cov).set_index("term")
    row = out.loc["intrinsic"]
    assert row["n_uncovered_not_estimable"] == 1
    assert row["mean_uncovered"] is None or pd.isna(row["mean_uncovered"])
    assert row["mean_covered"] == 2.0


def test_no_arm_contrast_appears_anywhere_in_the_output():
    cov = coverage_by_patient(
        _joined([("A", "AD", True), ("B", "SER", True), ("C", "Unknown", True)])
    )
    out = compare(_decomposition({"A": 1.0, "B": 2.0, "C": 3.0}), cov)
    assert summarise(cov, out)["arm_contrast_computed"] is False
    joined_columns = " ".join(out.columns)
    assert "AD" not in joined_columns and "SER" not in joined_columns


def test_no_threshold_is_attached_to_the_standardised_difference():
    """None was pre-committed; inventing one after seeing the table is the choice
    this project refuses elsewhere."""
    cov = coverage_by_patient(_joined([("A", "AD", True), ("B", "Unknown", True)]))
    out = compare(_decomposition({"A": 1.0, "B": 2.0}), cov)
    summary = summarise(cov, out)
    assert summary["systematic_difference_rule"] == "none pre-committed; inspected by hand"
    assert not any("pass" in str(k).lower() or "fail" in str(k).lower() for k in summary)
