"""MHIST's first gate must remain metadata-only and image-level."""

from __future__ import annotations

import zipfile

import pandas as pd
import pytest

from src.reference.jobs.mhist_feasibility import (
    MODELING_LICENSE_STATUS,
    MhistFeasibilityError,
    archive_image_names,
    inventory,
    read_annotations,
)


def _annotations() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "Image Name": ["MHIST_a.png", "MHIST_b.png", "MHIST_c.png", "MHIST_d.png"],
            "Majority Vote Label": ["HP", "SSA", "HP", "SSA"],
            "Number of Annotators who Selected SSA (Out of 7)": [0, 4, 3, 7],
            "Partition": ["train", "train", "test", "test"],
        }
    )


def test_inventory_is_image_level_and_refuses_a_patient_claim():
    annotations = _annotations()
    result = inventory(
        annotations.assign(ssa_votes=[0, 4, 3, 7]),
        set(annotations["Image Name"]),
        expected_images=4,
    )

    assert result["n_images"].sum() == 4
    assert set(result["evaluation_unit"]) == {"image"}
    assert result["patient_or_slide_group_key"].isna().all()
    assert not result["patient_held_out_licensed"].any()


def test_license_status_is_not_inferred_from_the_data_structure():
    """A sound split cannot override a license that restricts derivatives."""
    assert MODELING_LICENSE_STATUS == "pending_written_permission"


def test_read_annotations_rejects_a_vote_label_contradiction(tmp_path):
    annotations = _annotations()
    annotations.loc[0, "Majority Vote Label"] = "SSA"
    path = tmp_path / "annotations.csv"
    annotations.to_csv(path, index=False)

    with pytest.raises(MhistFeasibilityError, match="contradict"):
        read_annotations(path)


def test_inventory_rejects_an_archive_annotation_mismatch():
    annotations = _annotations().assign(ssa_votes=[0, 4, 3, 7])

    with pytest.raises(MhistFeasibilityError, match="disagree"):
        inventory(annotations, {"MHIST_a.png"}, expected_images=4)


def test_archive_names_reads_directory_without_decoding_pixels(tmp_path):
    path = tmp_path / "images.zip"
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("images/MHIST_a.png", b"not decoded")
        archive.writestr("images/MHIST_b.png", b"not decoded")

    assert archive_image_names(path) == {"MHIST_a.png", "MHIST_b.png"}
