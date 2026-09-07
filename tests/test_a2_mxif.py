"""A2 Step 1 keeps product and patient cardinalities honest before outcomes."""

from __future__ import annotations

import pandas as pd
import pytest

from src.reference.a2_mxif import (
    A2MxIFInventoryError,
    inventory_counts,
    participant_inventory,
    read_region_manifest,
    validate_region_manifest,
)


def _region(patient, region, tissue, kind, status="supplementary_table_1"):
    batch = f"{patient}_0000_01_01"
    slide = f"{batch}_region_{region:03d}"
    return {
        "participant_id": patient,
        "batch_name": batch,
        "slide_region": slide,
        "tissue_category": tissue,
        "broad_precancer_type": kind,
        "metadata_status": status,
        "archive_entry": (
            f"Upload_data/colon_polyp_images/{batch}_ADJ_"
            f"region_{region:03d}_downsampled.npz"
        ),
    }


def _obs():
    return pd.DataFrame(
        {
            "dataset": ["VUMC_HTAN_validation"] * 3,
            "patient_id": [
                "Chen_2021_Cell.HTA11_1",
                "Chen_2021_Cell.HTA11_1",
                "Chen_2021_Cell.HTA11_2",
            ],
            "sample_type": ["polyp", "healthy normal", "polyp"],
        }
    )


def test_archive_only_regions_cannot_silently_inherit_paper_labels():
    frame = pd.DataFrame([_region("HTA11_1", 1, "tumor", "AD", "archive_only")])
    with pytest.raises(A2MxIFInventoryError, match="may not inherit"):
        validate_region_manifest(frame)


def test_archive_and_metadata_identifiers_must_agree():
    frame = pd.DataFrame([_region("HTA11_1", 1, "tumor", "AD")])
    frame.loc[0, "slide_region"] = "HTA11_1_0000_01_01_region_002"
    with pytest.raises(A2MxIFInventoryError, match="identifier mismatch"):
        validate_region_manifest(frame)


def test_patient_not_region_is_the_crosswalk_unit():
    regions = pd.DataFrame(
        [
            _region("HTA11_1", 1, "normal", "AD"),
            _region("HTA11_1", 2, "tumor", "AD"),
            _region("HTA11_2", 1, "tumor", "SSL"),
        ]
    )
    got = participant_inventory(regions, _obs()).set_index("participant_id")
    assert len(got) == 2
    assert bool(got.loc["HTA11_1", "mxif_within_batch_pair"])
    assert bool(got.loc["HTA11_1", "avenue_a_has_both_arms"])
    assert bool(got.loc["HTA11_1", "conventional_ad_mxif_pair"])
    assert not bool(got.loc["HTA11_2", "mxif_within_batch_pair"])
    assert not bool(got.loc["HTA11_2", "avenue_a_has_both_arms"])
    assert set(got["crosswalk_level"]) == {"participant_only"}
    assert not got["specimen_exact"].any()


def test_an_unlabelled_archive_region_cannot_create_a_pair():
    unlabelled = _region("HTA11_1", 2, "", "", "archive_only")
    regions = pd.DataFrame([_region("HTA11_1", 1, "tumor", "AD"), unlabelled])
    got = participant_inventory(regions, _obs()).iloc[0]
    assert got["n_unlabeled_archive_regions"] == 1
    assert not bool(got["mxif_within_batch_pair"])


def test_normal_and_tumor_in_different_batches_do_not_create_a_pair():
    regions = pd.DataFrame(
        [
            _region("HTA11_1", 1, "normal", "AD"),
            _region("HTA11_1", 2, "tumor", "AD"),
        ]
    )
    other_batch = "HTA11_1_0000_02_02"
    regions.loc[1, "batch_name"] = other_batch
    regions.loc[1, "slide_region"] = f"{other_batch}_region_002"
    regions.loc[1, "archive_entry"] = (
        f"Upload_data/colon_polyp_images/{other_batch}_ADJ_"
        "region_002_downsampled.npz"
    )
    got = participant_inventory(regions, _obs()).iloc[0]
    assert got["n_normal_regions"] == 1
    assert got["n_tumor_regions"] == 1
    assert got["n_mxif_paired_batches"] == 0
    assert not bool(got["mxif_within_batch_pair"])


def test_committed_public_inventory_cardinalities_are_frozen():
    regions = read_region_manifest()
    both_arms = {
        "HTA11_10167", "HTA11_4255", "HTA11_6298", "HTA11_6801",
        "HTA11_7862", "HTA11_7956", "HTA11_8099", "HTA11_8622",
        "HTA11_866", "HTA11_8920", "HTA11_9408",
    }
    obs_rows = []
    for patient in sorted(set(regions["participant_id"])):
        obs_rows.append({
            "dataset": "VUMC_HTAN_validation",
            "patient_id": f"Chen_2021_Cell.{patient}",
            "sample_type": "polyp",
        })
        if patient in both_arms:
            obs_rows.append({
                "dataset": "VUMC_HTAN_validation",
                "patient_id": f"Chen_2021_Cell.{patient}",
                "sample_type": "healthy normal",
            })
    participants = participant_inventory(regions, pd.DataFrame(obs_rows))
    assert inventory_counts(regions, participants) == {
        "archive_regions": 42,
        "labeled_regions": 38,
        "unlabeled_archive_regions": 4,
        "participants": 15,
        "mxif_within_batch_pairs": 8,
        "conventional_ad_mxif_pairs": 3,
        "avenue_a_patient_matches": 15,
        "avenue_a_both_arm_patients": 11,
    }
