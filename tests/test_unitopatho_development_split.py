"""The UniToPatho development split cannot leak source test slides."""

from __future__ import annotations

import pandas as pd
import pytest

from src.reference.jobs.unitopatho_development_split import (
    EXPECTED_GROUP_COUNTS,
    EXPECTED_VALIDATION_GROUPS_BY_LABEL,
    TASK_LABELS,
    UniToPathoFeasibilityError,
    build_development_split,
)


def _actual_scale_images() -> pd.DataFrame:
    source_counts = {
        "HP": {"training": 29, "validation": 2, "test": 10},
        "NORM": {"training": 12, "validation": 1, "test": 8},
        "TA.HG": {"training": 16, "validation": 1, "test": 9},
        "TA.LG": {"training": 97, "validation": 3, "test": 46},
        "TVA.HG": {"training": 11, "validation": 2, "test": 7},
        "TVA.LG": {"training": 28, "validation": 2, "test": 8},
    }
    rows = []
    index = 0
    for label, splits in source_counts.items():
        for split, count in splits.items():
            for group_number in range(count):
                rows.append(
                    {
                        "image_index": index,
                        "label": label,
                        "location": (
                            f"{label}/{split}-{group_number}.ndpi_ROI__mpp0.44_reg000.png"
                        ),
                        "slide_group": f"{label}/{split}-{group_number}",
                        "split": split,
                    }
                )
                index += 1
    return pd.DataFrame(rows)


def test_development_split_preserves_source_test_and_locked_allocations():
    assignments, summary = build_development_split(_actual_scale_images())

    assert assignments["model_split"].value_counts().to_dict() == EXPECTED_GROUP_COUNTS
    validation_counts = (
        assignments.loc[assignments["model_split"].eq("validation")]
        .groupby("label")["slide_group"]
        .nunique()
        .to_dict()
    )
    assert validation_counts == EXPECTED_VALIDATION_GROUPS_BY_LABEL
    assert assignments.loc[assignments["source_split"].eq("test"), "model_split"].eq(
        "test"
    ).all()
    assert set(summary["label"]) == set(TASK_LABELS)


def test_development_split_refuses_a_slide_with_multiple_labels():
    images = _actual_scale_images()
    duplicate = images.iloc[[0]].copy()
    duplicate["label"] = "NORM"
    images = pd.concat([images, duplicate], ignore_index=True)

    with pytest.raises(UniToPathoFeasibilityError, match="more than one label"):
        build_development_split(images)
