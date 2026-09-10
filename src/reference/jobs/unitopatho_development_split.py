#!/usr/bin/env python
"""Lock UniToPatho's patient-contained task and development split before pixels.

This job reads the public 800 YAML manifest only.  It preserves the source test
set, retains every source validation slide in validation, and adds a
deterministic class-stratified subset of source-training slides until validation
contains 20 percent of each non-test class.  No image pixels are opened and no
encoder, model, score, or abstention threshold is selected.
"""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

import pandas as pd

from src.common.io import write_versioned_table
from src.common.label_provenance import (
    check_no_circular_claim,
    provenance_meta,
)
from src.common.provenance import DEFAULT_SEED
from src.reference.jobs.unitopatho_feasibility import (
    CLAIM_PROVENANCE,
    DATASET_TERMS_REVIEWED,
    EXPECTED_MANIFEST_SHA256,
    LABEL_PROVENANCE,
    UniToPathoFeasibilityError,
    read_manifest,
    sha256,
)

DEVELOPMENT_SPLIT_SALT = "unitopatho-development-v1"
VALIDATION_FRACTION = 0.20
EXPECTED_GROUP_COUNTS = {"test": 88, "training": 161, "validation": 43}
EXPECTED_VALIDATION_GROUPS_BY_LABEL = {
    "HP": 7,
    "NORM": 3,
    "TA.HG": 4,
    "TA.LG": 20,
    "TVA.HG": 3,
    "TVA.LG": 6,
}
TASK_LABELS = ("HP", "NORM", "TA.HG", "TA.LG", "TVA.HG", "TVA.LG")


def _rank(slide_group: str) -> str:
    return hashlib.sha256(f"{DEVELOPMENT_SPLIT_SALT}:{slide_group}".encode()).hexdigest()


def group_table(images: pd.DataFrame) -> pd.DataFrame:
    """Collapse manifest patches to their one patient/slide grouping unit."""
    per_group = images.groupby("slide_group", as_index=False, observed=True).agg(
        label=("label", "first"),
        source_split=("split", "first"),
        n_patches=("image_index", "size"),
        n_labels=("label", "nunique"),
        n_source_splits=("split", "nunique"),
    )
    if not per_group["n_labels"].eq(1).all():
        raise UniToPathoFeasibilityError("a source slide group has more than one label")
    if not per_group["n_source_splits"].eq(1).all():
        raise UniToPathoFeasibilityError("a source slide group crosses supplied splits")
    return per_group.drop(columns=["n_labels", "n_source_splits"])


def build_development_split(images: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Create the fixed source-test-preserving development assignment."""
    groups = group_table(images)
    if tuple(sorted(groups["label"].unique())) != tuple(sorted(TASK_LABELS)):
        raise UniToPathoFeasibilityError("source groups differ from the locked six-label task")

    groups["model_split"] = groups["source_split"]
    for label in TASK_LABELS:
        label_rows = groups[groups["label"].eq(label)]
        development_rows = label_rows[~label_rows["source_split"].eq("test")]
        target = int(-(-len(development_rows) // 5))
        source_validation = development_rows[
            development_rows["source_split"].eq("validation")
        ]
        candidates = development_rows[development_rows["source_split"].eq("training")].copy()
        needed = target - len(source_validation)
        if needed < 0:
            raise UniToPathoFeasibilityError(
                f"source validation exceeds the locked 20% target for {label}"
            )
        if needed > len(candidates):
            raise UniToPathoFeasibilityError(
                f"not enough source-training slide groups to fill validation for {label}"
            )
        selected = (
            candidates.assign(rank=candidates["slide_group"].map(_rank))
            .sort_values("rank")
            .head(needed)
        )
        groups.loc[selected.index, "model_split"] = "validation"

    counts = groups["model_split"].value_counts().to_dict()
    if counts != EXPECTED_GROUP_COUNTS:
        raise UniToPathoFeasibilityError(
            f"locked group counts changed: expected {EXPECTED_GROUP_COUNTS}, got {counts}"
        )
    validation_by_label = (
        groups.loc[groups["model_split"].eq("validation")]
        .groupby("label", observed=True)["slide_group"]
        .nunique()
        .to_dict()
    )
    if validation_by_label != EXPECTED_VALIDATION_GROUPS_BY_LABEL:
        raise UniToPathoFeasibilityError(
            "locked validation allocation changed: "
            f"expected {EXPECTED_VALIDATION_GROUPS_BY_LABEL}, got {validation_by_label}"
        )
    if not groups.loc[groups["source_split"].eq("test"), "model_split"].eq("test").all():
        raise UniToPathoFeasibilityError("the source test set must remain untouched")

    groups = groups.sort_values(["model_split", "label", "slide_group"]).reset_index(drop=True)
    summary = (
        groups.groupby(["model_split", "label"], as_index=False, observed=True)
        .agg(n_slide_groups=("slide_group", "size"), n_patches=("n_patches", "sum"))
        .sort_values(["model_split", "label"])
        .reset_index(drop=True)
    )
    return groups, summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, help="audited unitopath-public-800.yml")
    parser.add_argument("--no-write", action="store_true", help="inspect without writing artifacts")
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args()

    check_no_circular_claim(labels=LABEL_PROVENANCE, claim=CLAIM_PROVENANCE)
    manifest_path = Path(args.manifest)
    if not manifest_path.is_file():
        raise FileNotFoundError(f"UniToPatho manifest not found: {manifest_path}")
    manifest_sha256 = sha256(manifest_path)
    if manifest_sha256 != EXPECTED_MANIFEST_SHA256:
        raise UniToPathoFeasibilityError(
            "manifest SHA-256 differs from the audited public configuration: "
            f"expected {EXPECTED_MANIFEST_SHA256}, got {manifest_sha256}"
        )
    _, images = read_manifest(manifest_path)
    assignments, summary = build_development_split(images)
    print(summary.to_string(index=False))
    print("METADATA ONLY — task and patient-contained split locked; no pixels or model read.")

    if not args.no_write:
        common_meta = {
            **provenance_meta(LABEL_PROVENANCE, CLAIM_PROVENANCE),
            "roadmap": "docs/ml_abstention_roadmap.md",
            "manifest_sha256": manifest_sha256,
            "development_split_salt": DEVELOPMENT_SPLIT_SALT,
            "validation_fraction": VALIDATION_FRACTION,
            "task_labels": list(TASK_LABELS),
            "dataset_terms_reviewed": DATASET_TERMS_REVIEWED,
            "feature_extraction_authorized": False,
            "engineering_only": True,
        }
        assignment_path = write_versioned_table(
            assignments,
            "unitopatho_development_split",
            seed=DEFAULT_SEED,
            notes="Locked metadata-only six-label task and patient-contained split",
            extra_meta=common_meta,
            allow_dirty=args.allow_dirty,
        )
        summary_path = write_versioned_table(
            summary,
            "unitopatho_development_split_summary",
            seed=DEFAULT_SEED,
            notes="Locked group and patch counts; source test is untouched",
            extra_meta=common_meta,
            allow_dirty=args.allow_dirty,
        )
        print(f"wrote {assignment_path}")
        print(f"wrote {summary_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
