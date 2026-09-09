#!/usr/bin/env python
"""Inventory MHIST before opening image pixels or selecting a model.

This is WP-A's data gate in ``docs/ml_abstention_roadmap.md``.  It verifies
the delivered archive against the supplied MD5, reads only ``annotations.csv``
and the ZIP central directory, and records whether the data can support an
image-level agreement-aware abstention benchmark.  It never decodes an image.

    python -m src.reference.jobs.mhist_feasibility \
        --annotations "$HOME/Downloads/MHIST/annotations.csv" \
        --images-zip "$HOME/Downloads/MHIST/images.zip"

Commit the job before writing a durable result.  The versioned writer refuses
a dirty tree because its git SHA must reproduce the gate.
"""

from __future__ import annotations

import argparse
import hashlib
import sys
import zipfile
from pathlib import Path

import pandas as pd

from src.common.io import write_versioned_table
from src.common.label_provenance import (
    Measurement,
    check_no_circular_claim,
    provenance_meta,
)
from src.common.provenance import DEFAULT_SEED

# Invariant 11: declared before anything is read. This job defines its
# population from morphology and claims nothing about a
# transcript programme, but "unstated" is not "none".
LABEL_PROVENANCE = Measurement(
    modality="morphology",
    assay="MHIST seven-pathologist H&E majority-vote diagnosis",
    genes=(),
)
CLAIM_PROVENANCE = Measurement(
    modality="morphology",
    assay="MHIST delivered annotation and image-inventory fields",
    genes=(),
)

EXPECTED_ARCHIVE_MD5 = "486973abb3cf7985df2034701d81e9c1"
EXPECTED_COLUMNS = (
    "Image Name",
    "Majority Vote Label",
    "Number of Annotators who Selected SSA (Out of 7)",
    "Partition",
)
EXPECTED_LABELS = {"HP", "SSA"}
EXPECTED_PARTITIONS = {"train", "test"}
EXPECTED_IMAGES = 3152
MODELING_LICENSE_STATUS = "pending_written_permission"
MODELING_LICENSE_REASON = (
    "delivered RUA permits non-commercial research use but prohibits modification and "
    "derivative works"
)


class MhistFeasibilityError(ValueError):
    """The delivered MHIST release cannot support the pre-specified gate."""


def md5(path: str | Path, block_size: int = 1024 * 1024) -> str:
    """Return a file's MD5 without interpreting its content."""
    digest = hashlib.md5()  # noqa: S324 -- required to check the supplier's MD5.
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def read_annotations(path: str | Path) -> pd.DataFrame:
    """Read and validate the four delivered annotation fields."""
    table = pd.read_csv(path)
    if tuple(table.columns) != EXPECTED_COLUMNS:
        raise MhistFeasibilityError(
            f"unexpected annotation schema: expected {list(EXPECTED_COLUMNS)}, "
            f"got {list(table.columns)}"
        )
    if table["Image Name"].isna().any() or table["Image Name"].duplicated().any():
        raise MhistFeasibilityError("image names must be present and unique")
    if set(table["Majority Vote Label"].dropna()) != EXPECTED_LABELS:
        raise MhistFeasibilityError("majority labels must be exactly HP and SSA")
    if set(table["Partition"].dropna()) != EXPECTED_PARTITIONS:
        raise MhistFeasibilityError("partitions must be exactly train and test")

    votes = pd.to_numeric(
        table["Number of Annotators who Selected SSA (Out of 7)"], errors="coerce"
    )
    if votes.isna().any() or not votes.between(0, 7).all() or not (votes % 1 == 0).all():
        raise MhistFeasibilityError("SSA votes must be integer counts from 0 through 7")
    table = table.assign(ssa_votes=votes.astype(int))
    expected_labels = table["ssa_votes"].ge(4).map({True: "SSA", False: "HP"})
    if not table["Majority Vote Label"].eq(expected_labels).all():
        raise MhistFeasibilityError("majority labels contradict the seven-reader vote count")
    return table


def archive_image_names(path: str | Path) -> set[str]:
    """Read ZIP names only; image bytes are deliberately never decoded."""
    with zipfile.ZipFile(path) as archive:
        names = {
            item.filename.removeprefix("images/")
            for item in archive.infolist()
            if not item.is_dir()
        }
    if not names or any(not name.endswith(".png") for name in names):
        raise MhistFeasibilityError("archive must contain PNG images under images/")
    return names


def inventory(
    annotations: pd.DataFrame,
    image_names: set[str],
    *,
    expected_images: int = EXPECTED_IMAGES,
) -> pd.DataFrame:
    """Return the declared-split/label inventory and refuse any mismatch."""
    annotation_names = set(annotations["Image Name"])
    if annotation_names != image_names:
        raise MhistFeasibilityError(
            "archive and annotations disagree: "
            f"archive_only={len(image_names - annotation_names)}, "
            f"annotation_only={len(annotation_names - image_names)}"
        )
    if len(annotations) != expected_images:
        raise MhistFeasibilityError(
            f"expected {expected_images} annotated images, got {len(annotations)}"
        )

    out = (
        annotations.groupby(["Partition", "Majority Vote Label"], as_index=False)
        .agg(
            n_images=("Image Name", "size"),
            min_ssa_votes=("ssa_votes", "min"),
            max_ssa_votes=("ssa_votes", "max"),
        )
        .sort_values(["Partition", "Majority Vote Label"])
        .reset_index(drop=True)
    )
    out["evaluation_unit"] = "image"
    out["patient_or_slide_group_key"] = pd.NA
    out["patient_held_out_licensed"] = False
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotations", required=True, help="delivered MHIST annotations.csv")
    parser.add_argument("--images-zip", required=True, help="delivered MHIST images.zip")
    parser.add_argument("--no-write", action="store_true", help="inspect without a result artifact")
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args()

    # Invariant 11: refuse before the first read, not at the writer.
    check_no_circular_claim(labels=LABEL_PROVENANCE, claim=CLAIM_PROVENANCE)

    annotations_path = Path(args.annotations)
    images_path = Path(args.images_zip)
    for path in (annotations_path, images_path):
        if not path.is_file():
            raise FileNotFoundError(f"MHIST input not found: {path}")
    actual_md5 = md5(images_path)
    if actual_md5 != EXPECTED_ARCHIVE_MD5:
        raise MhistFeasibilityError(
            f"images.zip MD5 mismatch: expected {EXPECTED_ARCHIVE_MD5}, got {actual_md5}"
        )

    annotations = read_annotations(annotations_path)
    result = inventory(annotations, archive_image_names(images_path))
    print(result.to_string(index=False))
    print("\nSTRUCTURAL DATA GATE PASSES FOR AN IMAGE-LEVEL BENCHMARK")
    print("No patient or slide identifier was delivered; do not claim patient-held-out evaluation.")
    print("LICENSE CLARIFICATION PENDING: the delivered RUA prohibits derivative works;")
    print("Do not embed images, train a model, or publish derived outputs without")
    print("written permission.")

    if not args.no_write:
        path = write_versioned_table(
            result,
            "mhist_feasibility",
            seed=DEFAULT_SEED,
            notes="WP-A metadata-only data gate; no pixels decoded and no model selected",
            extra_meta={
                **provenance_meta(LABEL_PROVENANCE, CLAIM_PROVENANCE),
                "roadmap": "docs/ml_abstention_roadmap.md",
                "archive_md5": actual_md5,
                "n_images": int(len(annotations)),
                "declared_partitions": sorted(annotations["Partition"].unique()),
                "has_patient_or_slide_group_key": False,
                "evaluation_unit": "image",
                "patient_held_out_licensed": False,
                "modeling_license_status": MODELING_LICENSE_STATUS,
                "modeling_license_reason": MODELING_LICENSE_REASON,
            },
            allow_dirty=args.allow_dirty,
        )
        print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
