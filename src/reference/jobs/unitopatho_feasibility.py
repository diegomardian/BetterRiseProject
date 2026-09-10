#!/usr/bin/env python
"""Inventory UniToPatho's public manifest before downloading image pixels.

This is the WP-B data gate in ``docs/ml_abstention_roadmap.md``.  It reads a
public YAML manifest only: no H&E image is decoded, no encoder is selected,
and no model is trained.  The chosen ``800`` configuration is deliberately
kept separate from the 7,000-pixel configuration, because their supplied
splits do not agree for every overlapping slide.

    python -m src.reference.jobs.unitopatho_feasibility \
        --manifest "$HOME/Downloads/unitopath-public-800.yml"

The manifest can establish a patient-contained *engineering* split only.  It
does not provide a molecular endpoint or make diagnostic-label imitation into
a project biological result.  Dataset terms must be reviewed separately before
image download, feature extraction, or modeling.
"""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from src.common.io import write_versioned_table
from src.common.label_provenance import (
    Measurement,
    check_no_circular_claim,
    provenance_meta,
)
from src.common.provenance import DEFAULT_SEED

LABEL_PROVENANCE = Measurement(
    modality="morphology",
    assay="UniToPatho expert H&E diagnostic patch annotation",
    genes=(),
)
CLAIM_PROVENANCE = Measurement(
    modality="morphology",
    assay="UniToPatho public manifest structure and declared split",
    genes=(),
)

CONFIGURATION_NAME = "deephealth-uc2-800"
EXPECTED_MANIFEST_SHA256 = "692e85cbf6ffe5ba3bef70c19bb8d9fc544ce58c254fce4f1fd20a59cd3dfb26"
EXPECTED_PATCHES = 8669
EXPECTED_SLIDE_GROUPS = 292
EXPECTED_CLASSES = {"HP", "NORM", "TA.HG", "TA.LG", "TVA.HG", "TVA.LG"}
EXPECTED_SPLITS = {"training", "validation", "test"}
MODELING_LICENSE_STATUS = "terms_not_yet_verified"


class UniToPathoFeasibilityError(ValueError):
    """The delivered public manifest cannot support the pre-specified gate."""


def sha256(path: str | Path, block_size: int = 1024 * 1024) -> str:
    """Return the supplier-manifest SHA-256 without interpreting image data."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def _slide_group(location: str) -> str:
    """Derive the source WSI key from a documented UniToPatho patch location."""
    marker = ".ndpi_ROI__"
    if location.count(marker) != 1 or not location.endswith(".png"):
        raise UniToPathoFeasibilityError(
            "each location must be a PNG patch with exactly one '.ndpi_ROI__' marker"
        )
    return location.split(marker, maxsplit=1)[0]


def read_manifest(path: str | Path) -> tuple[dict[str, Any], pd.DataFrame]:
    """Read and validate one selected configuration's manifest fields only."""
    with Path(path).open(encoding="utf-8") as handle:
        manifest = yaml.safe_load(handle)
    if not isinstance(manifest, dict):
        raise UniToPathoFeasibilityError("manifest must contain a YAML mapping")
    required = {"name", "classes", "images", "split"}
    if not required.issubset(manifest):
        missing = sorted(required - set(manifest))
        raise UniToPathoFeasibilityError(f"manifest missing required fields: {missing}")
    if manifest["name"] != CONFIGURATION_NAME:
        raise UniToPathoFeasibilityError(
            f"expected {CONFIGURATION_NAME!r}, got {manifest['name']!r}; "
            "do not merge configurations"
        )
    declared_classes = manifest["classes"]
    if set(declared_classes) != EXPECTED_CLASSES or len(declared_classes) != len(EXPECTED_CLASSES):
        raise UniToPathoFeasibilityError(
            "manifest classes differ from the six pre-specified UniToPatho labels"
        )
    if not isinstance(manifest["images"], list) or not manifest["images"]:
        raise UniToPathoFeasibilityError("manifest images must be a nonempty list")
    if set(manifest["split"]) != EXPECTED_SPLITS:
        raise UniToPathoFeasibilityError(
            "manifest must declare exactly training, validation, and test splits"
        )

    split_for_index: dict[int, str] = {}
    for split, indices in manifest["split"].items():
        if not isinstance(indices, list) or any(not isinstance(index, int) for index in indices):
            raise UniToPathoFeasibilityError(f"{split} indices must be integers")
        for index in indices:
            if index in split_for_index:
                raise UniToPathoFeasibilityError(
                    f"image index {index} occurs in more than one supplied split"
                )
            split_for_index[index] = split

    n_images = len(manifest["images"])
    if set(split_for_index) != set(range(n_images)):
        raise UniToPathoFeasibilityError(
            "supplied splits must assign every image index exactly once"
        )

    rows: list[dict[str, object]] = []
    for index, image in enumerate(manifest["images"]):
        if not isinstance(image, dict) or set(image) != {"label", "location"}:
            raise UniToPathoFeasibilityError(
                "each image entry must contain exactly label and location"
            )
        label, location = image["label"], image["location"]
        if label not in EXPECTED_CLASSES or not isinstance(location, str):
            raise UniToPathoFeasibilityError(
                "each image must have a known label and string location"
            )
        group = _slide_group(location)
        if not group.startswith(f"{label}/"):
            raise UniToPathoFeasibilityError(
                "patch location class directory must agree with its declared label"
            )
        rows.append(
            {
                "image_index": index,
                "label": label,
                "location": location,
                "slide_group": group,
                "split": split_for_index[index],
            }
        )
    return manifest, pd.DataFrame(rows)


def inventory(images: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return split/class counts and a summary, refusing derived-slide leakage."""
    if len(images) != EXPECTED_PATCHES:
        raise UniToPathoFeasibilityError(
            f"expected {EXPECTED_PATCHES} patches in {CONFIGURATION_NAME}, got {len(images)}"
        )
    n_groups = images["slide_group"].nunique()
    if n_groups != EXPECTED_SLIDE_GROUPS:
        raise UniToPathoFeasibilityError(
            f"expected {EXPECTED_SLIDE_GROUPS} slide groups, got {n_groups}"
        )
    groups_per_split = images.groupby("slide_group", observed=True)["split"].nunique()
    if not groups_per_split.eq(1).all():
        raise UniToPathoFeasibilityError("a derived slide group crosses supplied splits")

    detail = (
        images.groupby(["split", "label"], as_index=False, observed=True)
        .agg(n_patches=("image_index", "size"), n_slide_groups=("slide_group", "nunique"))
        .sort_values(["split", "label"])
        .reset_index(drop=True)
    )
    summary = pd.DataFrame(
        [
            {
                "configuration": CONFIGURATION_NAME,
                "n_patches": len(images),
                "n_slide_groups": n_groups,
                "n_groups_crossing_supplied_splits": 0,
                "supplied_split_names": ",".join(sorted(images["split"].unique())),
                "patient_containment_supported_by_source": True,
                "molecular_endpoint_available": False,
                "biological_project_result_licensed": False,
                "modeling_license_status": MODELING_LICENSE_STATUS,
                "verdict": (
                    "ENGINEERING SPLIT STRUCTURE PASSES — "
                    "TERMS REVIEW REQUIRED BEFORE MODELING"
                ),
            }
        ]
    )
    return detail, summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, help="downloaded unitopath-public-800.yml")
    parser.add_argument(
        "--no-write", action="store_true", help="inspect without writing result artifacts"
    )
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args()

    check_no_circular_claim(labels=LABEL_PROVENANCE, claim=CLAIM_PROVENANCE)
    path = Path(args.manifest)
    if not path.is_file():
        raise FileNotFoundError(f"UniToPatho manifest not found: {path}")
    actual_sha256 = sha256(path)
    if actual_sha256 != EXPECTED_MANIFEST_SHA256:
        raise UniToPathoFeasibilityError(
            "manifest SHA-256 differs from the audited public configuration: "
            f"expected {EXPECTED_MANIFEST_SHA256}, got {actual_sha256}"
        )
    _, images = read_manifest(path)
    detail, summary = inventory(images)
    print(detail.to_string(index=False))
    print("\n" + summary.to_string(index=False))
    print("MANIFEST ONLY — no image pixels were opened and no model was selected.")

    if not args.no_write:
        common_meta = {
            **provenance_meta(LABEL_PROVENANCE, CLAIM_PROVENANCE),
            "roadmap": "docs/ml_abstention_roadmap.md",
            "manifest_sha256": actual_sha256,
            "dataset_terms_reviewed": False,
            "modeling_license_status": MODELING_LICENSE_STATUS,
            "molecular_endpoint_available": False,
            "engineering_only": True,
        }
        detail_path = write_versioned_table(
            detail,
            "unitopatho_manifest_inventory",
            seed=DEFAULT_SEED,
            notes="WP-B manifest-only inventory; no image pixels decoded and no model selected",
            extra_meta=common_meta,
            allow_dirty=args.allow_dirty,
        )
        summary_path = write_versioned_table(
            summary,
            "unitopatho_manifest_summary",
            seed=DEFAULT_SEED,
            notes="WP-B eligibility summary; engineering-only and terms remain unverified",
            extra_meta=common_meta,
            allow_dirty=args.allow_dirty,
        )
        print(f"wrote {detail_path}")
        print(f"wrote {summary_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
