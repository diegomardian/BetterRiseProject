"""Malignant vs. normal epithelium. W1.

The test that matters is `TestNormalEpitheliumValidation`. execution_plan.md §4
makes "normal epithelium not misread as tumour" the done-when for this stage, and
it is only meaningful because the CNV baseline is non-epithelial — a population
used as the reference is non-malignant by construction, so validating on it would
prove nothing.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.reference.malignancy import (
    MIN_REFERENCE_CELLS,
    REFERENCE_COMPARTMENTS,
    MalignancyError,
    call_malignancy,
    run_infercnv,
    select_cnv_reference,
    validate_normal_epithelium,
)

RNG = np.random.default_rng(20260818)


def cohort(
    n_reference: int = 200,
    n_normal_epi: int = 300,
    n_tumour_epi: int = 300,
    malignant_shift: float = 3.0,
    patient: str = "P1",
):
    """Diploid reference cells, normal epithelium that matches them, and tumour
    epithelium shifted upward in CNV score."""
    reference = RNG.normal(1.0, 0.15, n_reference)
    normal_epi = RNG.normal(1.0, 0.15, n_normal_epi)
    tumour_epi = RNG.normal(1.0 + malignant_shift * 0.15, 0.15, n_tumour_epi)

    return pd.DataFrame(
        {
            "cnv_score": np.concatenate([reference, normal_epi, tumour_epi]),
            "compartment": (
                ["immune"] * n_reference
                + ["epithelial"] * (n_normal_epi + n_tumour_epi)
            ),
            "tissue": (
                ["normal"] * (n_reference + n_normal_epi) + ["tumour"] * n_tumour_epi
            ),
            "patient_id": [patient] * (n_reference + n_normal_epi + n_tumour_epi),
        }
    )


class TestReferenceSelection:
    def test_non_epithelial_compartments_are_the_reference(self):
        """Not the matched normal epithelium — that would make the validation
        circular, since a baseline population is non-malignant by construction."""
        assert "epithelial" not in REFERENCE_COMPARTMENTS
        assert set(REFERENCE_COMPARTMENTS) == {"immune", "stromal", "endothelial"}

    def test_counts_reference_cells_per_patient(self):
        df = cohort()
        out = select_cnv_reference(df["compartment"], patient_id=df["patient_id"])
        assert int(out.loc[0, "n_reference"]) == 200
        assert int(out.loc[0, "n_epithelial"]) == 600
        assert bool(out.loc[0, "usable"])

    def test_a_patient_without_enough_reference_is_marked_unusable(self):
        df = cohort(n_reference=10)
        out = select_cnv_reference(df["compartment"], patient_id=df["patient_id"])
        assert not bool(out.loc[0, "usable"])

    def test_all_reference_compartments_count(self):
        df = pd.DataFrame(
            {
                "compartment": ["immune", "stromal", "endothelial", "epithelial"],
                "patient_id": ["P1"] * 4,
            }
        )
        out = select_cnv_reference(df["compartment"], patient_id=df["patient_id"],
                                   min_cells=3)
        assert int(out.loc[0, "n_reference"]) == 3


class TestCalls:
    def test_tumour_epithelium_is_called_malignant(self):
        df = cohort()
        calls = call_malignancy(
            df["cnv_score"], compartment=df["compartment"], patient_id=df["patient_id"]
        )
        tumour = calls[df["tissue"].to_numpy() == "tumour"]
        assert (tumour["call"].astype(str) == "malignant").mean() > 0.65

    def test_the_strict_threshold_trades_sensitivity_for_specificity(self):
        """MALIGNANT_QUANTILE = 0.99 is deliberate. At a 3-sigma separation it
        recovers about three quarters of malignant cells, and that is the right
        trade: a false MALIGNANT call moves a normal cell into the tumour arm,
        which is the direction that inflates the compositional term and the
        direction of the prior hypothesis.
        """
        df = cohort()
        strict = call_malignancy(
            df["cnv_score"], compartment=df["compartment"],
            patient_id=df["patient_id"], quantile=0.99,
        )
        loose = call_malignancy(
            df["cnv_score"], compartment=df["compartment"],
            patient_id=df["patient_id"], quantile=0.80,
        )
        is_tumour = df["tissue"].to_numpy() == "tumour"
        is_normal_epi = (df["tissue"].to_numpy() == "normal") & (
            df["compartment"].to_numpy() == "epithelial"
        )
        sensitivity = [
            float((c[is_tumour]["call"].astype(str) == "malignant").mean())
            for c in (strict, loose)
        ]
        false_positive = [
            float((c[is_normal_epi]["call"].astype(str) == "malignant").mean())
            for c in (strict, loose)
        ]
        assert sensitivity[1] > sensitivity[0]        # loose catches more
        assert false_positive[1] > false_positive[0]  # and costs specificity
        assert false_positive[0] < 0.05               # strict keeps it low

    def test_reference_cells_are_labelled_reference_not_non_malignant(self):
        """They defined the baseline; calling them non-malignant is circular."""
        df = cohort()
        calls = call_malignancy(
            df["cnv_score"], compartment=df["compartment"], patient_id=df["patient_id"]
        )
        reference = calls[calls["compartment"] == "immune"]
        assert (reference["call"].astype(str) == "reference").all()

    def test_confidence_scales_with_distance_past_the_threshold(self):
        weak = cohort(malignant_shift=1.5)
        strong = cohort(malignant_shift=6.0)
        out = []
        for df in (weak, strong):
            calls = call_malignancy(
                df["cnv_score"], compartment=df["compartment"],
                patient_id=df["patient_id"],
            )
            tumour = calls[df["tissue"].to_numpy() == "tumour"]
            out.append(float(tumour["confidence"].median()))
        assert out[1] > out[0]

    def test_threshold_is_per_patient(self):
        """A baseline pooled across patients would fold germline copy-number
        variation and per-patient capture differences into the call."""
        noisy = cohort(patient="P2")
        noisy["cnv_score"] = noisy["cnv_score"] * 3.0
        df = pd.concat([cohort(patient="P1"), noisy], ignore_index=True)
        calls = call_malignancy(
            df["cnv_score"], compartment=df["compartment"], patient_id=df["patient_id"]
        )
        thresholds = calls.groupby("patient_id", observed=True)["threshold"].first()
        assert thresholds["P2"] > thresholds["P1"]

    def test_a_patient_without_a_reference_is_not_called(self):
        """Never guessed. A cell with no baseline has no malignancy status."""
        df = cohort(n_reference=5)
        calls = call_malignancy(
            df["cnv_score"], compartment=df["compartment"], patient_id=df["patient_id"]
        )
        assert (calls["call"].astype(str) == "not_called").all()
        assert calls["confidence"].isna().all()

    def test_length_mismatch_raises(self):
        df = cohort()
        with pytest.raises(MalignancyError, match="lengths differ"):
            call_malignancy(df["cnv_score"][:10], compartment=df["compartment"],
                            patient_id=df["patient_id"])

    def test_empty_input_raises(self):
        with pytest.raises(MalignancyError, match="no cells"):
            call_malignancy([], compartment=[], patient_id=[])


class TestNormalEpitheliumValidation:
    """§4's done-when: normal epithelium must not be misread as tumour."""

    def test_clean_calls_pass(self):
        df = cohort()
        calls = call_malignancy(
            df["cnv_score"], compartment=df["compartment"], patient_id=df["patient_id"]
        )
        report = validate_normal_epithelium(calls, tissue=df["tissue"])
        assert bool(report.loc[0, "passed"])
        assert report.loc[0, "specificity"] > 0.95

    def test_normal_epithelium_miscalled_as_tumour_fails(self):
        """If this passed silently, every downstream number would be computed
        over a tumour arm contaminated with normal cells."""
        df = cohort()
        normal_epi = (df["tissue"] == "normal") & (df["compartment"] == "epithelial")
        df.loc[normal_epi, "cnv_score"] += 3.0        # make them look aneuploid
        calls = call_malignancy(
            df["cnv_score"], compartment=df["compartment"], patient_id=df["patient_id"]
        )
        report = validate_normal_epithelium(calls, tissue=df["tissue"])
        assert not bool(report.loc[0, "passed"])

    def test_reference_cells_are_excluded_from_the_check(self):
        """They are in normal tissue but are the baseline, so counting them
        would inflate specificity toward 1 automatically."""
        df = cohort()
        calls = call_malignancy(
            df["cnv_score"], compartment=df["compartment"], patient_id=df["patient_id"]
        )
        report = validate_normal_epithelium(calls, tissue=df["tissue"])
        assert int(report.loc[0, "n_normal_epithelial"]) == 300   # not 500

    def test_per_patient_rows(self):
        df = pd.concat([cohort(patient="P1"), cohort(patient="P2")], ignore_index=True)
        calls = call_malignancy(
            df["cnv_score"], compartment=df["compartment"], patient_id=df["patient_id"]
        )
        report = validate_normal_epithelium(calls, tissue=df["tissue"])
        assert set(report["patient_id"]) == {"P1", "P2"}

    def test_no_normal_epithelium_raises_rather_than_passing_vacuously(self):
        df = cohort(n_normal_epi=0)
        calls = call_malignancy(
            df["cnv_score"], compartment=df["compartment"], patient_id=df["patient_id"]
        )
        with pytest.raises(MalignancyError, match="nothing to validate"):
            validate_normal_epithelium(calls, tissue=df["tissue"])

    def test_tissue_length_mismatch_raises(self):
        df = cohort()
        calls = call_malignancy(
            df["cnv_score"], compartment=df["compartment"], patient_id=df["patient_id"]
        )
        with pytest.raises(MalignancyError, match="entries for"):
            validate_normal_epithelium(calls, tissue=["normal"])


def test_infercnv_is_an_explicit_todo():
    with pytest.raises(NotImplementedError, match="real per-patient matrices"):
        run_infercnv()


def test_min_reference_cells_is_documented_not_arbitrary():
    assert MIN_REFERENCE_CELLS >= 20
