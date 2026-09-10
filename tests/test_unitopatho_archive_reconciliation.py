"""The UniToPatho archive must remain exactly the manifest-selected 800 set."""

from __future__ import annotations

import pandas as pd
import pytest

from src.reference.jobs.unitopatho_archive_reconciliation import (
    REQUIRED_AUXILIARY_FILES,
    reconcile,
)
from src.reference.jobs.unitopatho_feasibility import UniToPathoFeasibilityError


def _images() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "location": [
                "HP/first.ndpi_ROI__mpp0.44_reg000.png",
                "TA.LG/second.ndpi_ROI__mpp0.44_reg001.png",
            ]
        }
    )


def test_reconciliation_accepts_exact_manifest_paths_and_csvs():
    paths = [*_images()["location"], *sorted(REQUIRED_AUXILIARY_FILES)]
    detail, summary = reconcile(_images(), paths)

    assert detail.loc[0, "n_missing_paths"] == 0
    assert detail.loc[0, "n_unexpected_paths"] == 0
    assert summary.loc[0, "archive_paths_exactly_match_manifest"]
    assert not summary.loc[0, "dataset_terms_reviewed"]


def test_reconciliation_refuses_a_manifest_with_duplicate_image_locations():
    images = _images()
    images.loc[1, "location"] = images.loc[0, "location"]
    paths = [images.loc[0, "location"], *sorted(REQUIRED_AUXILIARY_FILES)]

    with pytest.raises(UniToPathoFeasibilityError, match="duplicate_rows=2"):
        reconcile(images, paths)


@pytest.mark.parametrize(
    ("paths", "expected"),
    [
        (["HP/first.ndpi_ROI__mpp0.44_reg000.png", "train.csv", "test.csv"], "missing=1"),
        (
            [
                *_images()["location"],
                *sorted(REQUIRED_AUXILIARY_FILES),
                "7000/foreign.ndpi_ROI__mpp0.44_reg000.png",
            ],
            "unexpected=1",
        ),
        ([*_images()["location"], *sorted(REQUIRED_AUXILIARY_FILES), "train.csv"], "duplicates=1"),
    ],
)
def test_reconciliation_refuses_partial_foreign_or_duplicate_archives(paths, expected):
    with pytest.raises(UniToPathoFeasibilityError, match=expected):
        reconcile(_images(), paths)
