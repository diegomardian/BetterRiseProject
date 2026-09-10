"""UniToPatho's public-manifest gate must prevent cross-slide leakage."""

from __future__ import annotations

import pandas as pd
import pytest
import yaml

from src.reference.jobs.unitopatho_feasibility import (
    EXPECTED_CLASSES,
    EXPECTED_PATCHES,
    EXPECTED_SLIDE_GROUPS,
    MODELING_LICENSE_STATUS,
    UniToPathoFeasibilityError,
    _slide_group,
    inventory,
    read_manifest,
)


def _manifest_scale_frame() -> pd.DataFrame:
    """A compact fixture promoted to the audited size only inside the test."""
    rows = []
    for index in range(EXPECTED_PATCHES):
        split = "training" if index < 5000 else "test"
        group = f"TA.LG/slide-{index % EXPECTED_SLIDE_GROUPS}.ndpi_ROI__"
        rows.append(
            {
                "image_index": index,
                "label": "TA.LG",
                "location": f"{group}mpp0.44_reg000_crop.png",
                "slide_group": group.removesuffix(".ndpi_ROI__"),
                "split": split,
            }
        )
    # Keep every group in one split, matching the actual gate condition.
    frame = pd.DataFrame(rows)
    frame["split"] = frame["slide_group"].map(
        {
            f"TA.LG/slide-{group}": "training" if group < 150 else "test"
            for group in range(EXPECTED_SLIDE_GROUPS)
        }
    )
    return frame


def test_inventory_records_engineering_only_patient_containment():
    detail, summary = inventory(_manifest_scale_frame())

    assert detail["n_patches"].sum() == EXPECTED_PATCHES
    assert summary.loc[0, "n_slide_groups"] == EXPECTED_SLIDE_GROUPS
    assert summary.loc[0, "n_groups_crossing_supplied_splits"] == 0
    assert summary.loc[0, "patient_containment_supported_by_source"]
    assert not summary.loc[0, "molecular_endpoint_available"]
    assert not summary.loc[0, "biological_project_result_licensed"]


def test_inventory_refuses_a_slide_group_crossing_splits():
    images = _manifest_scale_frame()
    images.loc[0, "split"] = "test"

    with pytest.raises(UniToPathoFeasibilityError, match="crosses supplied splits"):
        inventory(images)


def test_source_location_must_be_a_png_with_the_roi_marker():
    assert _slide_group("HP/slide.ndpi_ROI__mpp0.44_reg000.png") == "HP/slide"
    with pytest.raises(UniToPathoFeasibilityError, match="ROI"):
        _slide_group("HP/slide.png")


def test_terms_are_not_inferred_from_an_eligible_split():
    assert MODELING_LICENSE_STATUS == "cc_by_attribution_required"
    assert EXPECTED_CLASSES == {"HP", "NORM", "TA.HG", "TA.LG", "TVA.HG", "TVA.LG"}


def test_read_manifest_parses_the_pinned_schema_without_an_image_archive(tmp_path):
    labels = sorted(EXPECTED_CLASSES)
    payload = {
        "name": "deephealth-uc2-800",
        "classes": labels,
        "images": [
            {
                "label": label,
                "location": f"{label}/{label}.ndpi_ROI__mpp0.44_reg000.png",
            }
            for label in labels
        ],
        "split": {"training": [0, 1], "validation": [2, 3], "test": [4, 5]},
    }
    path = tmp_path / "unitopath-public-800.yml"
    path.write_text(yaml.safe_dump(payload), encoding="utf-8")

    manifest, images = read_manifest(path)

    assert manifest["name"] == "deephealth-uc2-800"
    assert len(images) == len(EXPECTED_CLASSES)
    assert images["location"].is_unique
