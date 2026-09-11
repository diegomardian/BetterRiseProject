"""Gate 4 refuses the suffix-family merge that would inflate an arm."""

from __future__ import annotations

import pandas as pd
import pytest

from src.reference.jobs.chen_subtype_id_provenance import (
    IdProvenanceError,
    arm_summary,
    classify_patients,
    exact_join,
    read_labels,
)


def _crosswalk(ids: list[str]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "scRNA_biospecimen_id": ids,
            "patient_id": ["Chen_2021_Cell." + i.rsplit("_", 1)[0] for i in ids],
        }
    )


def _labels(rows: list[tuple[str, str]]) -> pd.DataFrame:
    return pd.DataFrame(
        [{"sampleId": s, "patientId": s.rsplit("_", 1)[0], "POLYP_TYPE": t} for s, t in rows]
    )


def test_two_specimens_of_one_participant_stay_distinct():
    """HTA11_866 carries two lesions; merging them would invent an arm member."""
    ids = ["HTA11_866_2000001011", "HTA11_866_3004761011"]
    joined = exact_join(_labels([(ids[0], "AD"), (ids[1], "AD")]), _crosswalk(ids))
    assert joined["label_matched"].all()
    assert joined["scRNA_biospecimen_id"].nunique() == 2
    assert joined.attrs["n_prefix_families_with_multiple_specimens"] == 1


def test_a_specimen_absent_from_the_labels_is_unmatched_not_dropped():
    ids = ["HTA11_1_2000001011", "HTA11_2_2000001011"]
    joined = exact_join(_labels([(ids[0], "AD")]), _crosswalk(ids))
    assert len(joined) == 2
    assert int(joined["label_matched"].sum()) == 1


def test_a_conflicted_patient_is_excluded_and_never_assigned():
    """Rule 1: assigning by any tiebreak is a modelling choice, not a measurement."""
    ids = ["HTA11_9_2000001011", "HTA11_9_3000001011"]
    joined = exact_join(_labels([(ids[0], "AD"), (ids[1], "SER")]), _crosswalk(ids))
    patients = classify_patients(joined)
    assert list(patients["status"]) == ["excluded_conflicting"]
    assert patients["arm"].isna().all()


def test_a_concordant_multi_lesion_patient_contributes_one_patient():
    """Rule 2: lesions are not entered twice."""
    ids = ["HTA11_7_2000001011", "HTA11_7_3000001011"]
    joined = exact_join(_labels([(ids[0], "AD"), (ids[1], "AD")]), _crosswalk(ids))
    patients = classify_patients(joined)
    assert len(patients) == 1
    assert patients.iloc[0]["n_labelled_lesions"] == 2
    assert patients.iloc[0]["arm"] == "AD"


def test_labelled_and_analysable_lesions_are_reported_separately():
    """The conflicted lesions are labelled but not analysable; both are counted."""
    ids = ["HTA11_9_2000001011", "HTA11_9_3000001011", "HTA11_5_2000001011"]
    joined = exact_join(
        _labels([(ids[0], "AD"), (ids[1], "SER"), (ids[2], "AD")]), _crosswalk(ids)
    )
    summary = arm_summary(classify_patients(joined), joined).set_index("arm")
    assert summary.loc["AD", "labelled_lesions"] == 2
    assert summary.loc["AD", "analysable_lesions"] == 1
    assert summary.loc["SER", "labelled_lesions"] == 1
    assert summary.loc["SER", "analysable_lesions"] == 0


def test_a_missing_snapshot_is_refused_rather_than_fetched_live():
    with pytest.raises(IdProvenanceError, match="run gate 3 first"):
        read_labels(pd.io.common.Path("does/not/exist.json"))
