#!/usr/bin/env python
"""Reconcile a metadata-only UniToPatho 800 archive listing to its manifest.

This is the second WP-B data gate. The input listing is produced with::

    rclone lsf -R "Meowsers:800" > unitopatho_800_listing.txt

The job reads paths only. It does not open image pixels, select an encoder, or
fit a model. It refuses any missing, extra, or duplicated file path rather than
treating a partial archive as the audited 800 configuration.
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
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
    EXPECTED_MANIFEST_SHA256,
    LABEL_PROVENANCE,
    MODELING_LICENSE_STATUS,
    UniToPathoFeasibilityError,
    read_manifest,
    sha256,
)

REQUIRED_AUXILIARY_FILES = {"train.csv", "test.csv"}


def read_listing(path: str | Path) -> list[str]:
    """Read an ``rclone lsf -R`` listing, retaining file paths only."""
    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(f"UniToPatho archive listing not found: {source}")
    paths = []
    for line in source.read_text(encoding="utf-8").splitlines():
        relative = line.strip()
        if not relative or relative.endswith("/"):
            continue
        if relative.startswith("/") or "\\" in relative or "/../" in f"/{relative}":
            raise UniToPathoFeasibilityError(
                "archive listing must use safe, root-relative POSIX paths"
            )
        paths.append(relative)
    if not paths:
        raise UniToPathoFeasibilityError("archive listing contains no file paths")
    return paths


def reconcile(images: pd.DataFrame, listing: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Require an exact one-to-one match between selected paths and listing."""
    if not images["location"].is_unique:
        duplicate_count = int(images["location"].duplicated(keep=False).sum())
        raise UniToPathoFeasibilityError(
            "manifest references a file path more than once: "
            f"duplicate_rows={duplicate_count}"
        )
    expected_images = set(images["location"])
    expected = expected_images | REQUIRED_AUXILIARY_FILES
    observed = set(listing)
    duplicates = sorted(path for path, count in Counter(listing).items() if count > 1)
    missing = sorted(expected - observed)
    unexpected = sorted(observed - expected)
    listed_images = {path for path in observed if path.endswith(".png")}
    listed_auxiliary = observed & REQUIRED_AUXILIARY_FILES

    detail = pd.DataFrame(
        [
            {
                "configuration": "deephealth-uc2-800",
                "n_expected_image_files": len(expected_images),
                "n_listed_image_files": len(listed_images),
                "n_required_auxiliary_files": len(REQUIRED_AUXILIARY_FILES),
                "n_listed_auxiliary_files": len(listed_auxiliary),
                "n_missing_paths": len(missing),
                "n_unexpected_paths": len(unexpected),
                "n_duplicate_listing_paths": len(duplicates),
            }
        ]
    )
    if missing or unexpected or duplicates:
        raise UniToPathoFeasibilityError(
            "UniToPatho archive does not exactly match the audited 800 manifest: "
            f"missing={len(missing)}, unexpected={len(unexpected)}, duplicates={len(duplicates)}"
        )

    summary = detail.assign(
        archive_paths_exactly_match_manifest=True,
        dataset_terms_reviewed=False,
        modeling_license_status=MODELING_LICENSE_STATUS,
        molecular_endpoint_available=False,
        engineering_only=True,
        verdict=(
            "ARCHIVE STRUCTURE RECONCILED — TERMS REVIEW AND LOCKED "
            "PROTOCOL REQUIRED BEFORE MODELING"
        ),
    )
    return detail, summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, help="audited unitopath-public-800.yml")
    parser.add_argument("--listing", required=True, help="metadata-only rclone lsf -R output")
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
    listing_path = Path(args.listing)
    listing = read_listing(listing_path)
    detail, summary = reconcile(images, listing)
    print(detail.to_string(index=False))
    print("\n" + summary.to_string(index=False))
    print("METADATA LISTING ONLY — no image pixels were opened and no model was selected.")

    if not args.no_write:
        common_meta = {
            **provenance_meta(LABEL_PROVENANCE, CLAIM_PROVENANCE),
            "roadmap": "docs/ml_abstention_roadmap.md",
            "manifest_sha256": manifest_sha256,
            "listing_sha256": sha256(listing_path),
            "dataset_terms_reviewed": False,
            "modeling_license_status": MODELING_LICENSE_STATUS,
            "molecular_endpoint_available": False,
            "engineering_only": True,
        }
        detail_path = write_versioned_table(
            detail,
            "unitopatho_archive_reconciliation",
            seed=DEFAULT_SEED,
            notes="WP-B metadata-only archive reconciliation; no pixels decoded or model selected",
            extra_meta=common_meta,
            allow_dirty=args.allow_dirty,
        )
        summary_path = write_versioned_table(
            summary,
            "unitopatho_archive_reconciliation_summary",
            seed=DEFAULT_SEED,
            notes="WP-B exact archive-path verdict; terms and protocol remain outstanding",
            extra_meta=common_meta,
            allow_dirty=args.allow_dirty,
        )
        print(f"wrote {detail_path}")
        print(f"wrote {summary_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
