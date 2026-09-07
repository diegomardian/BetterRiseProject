#!/usr/bin/env python
"""Build A2 Step 1: public MxIF product and Avenue-A patient crosswalk.

This job performs no image download and computes no marker outcome.  It reads
the fixed region inventory transcribed from MILWRM Supplementary Table 1 plus
the Zenodo ZIP directory, then joins patient identifiers to the cached ICBI
metadata used by Avenue A.

    python -m src.reference.jobs.a2_mxif_inventory --no-write

Remove ``--no-write`` only after committing the inventory: the versioned writer
will otherwise correctly refuse a result whose git SHA cannot reproduce it.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from src.common.io import write_versioned_table
from src.common.paths import INTERIM_DIR
from src.common.provenance import DEFAULT_SEED
from src.reference.a2_mxif import (
    MXIF_CHANNELS,
    REGION_MANIFEST,
    inventory_counts,
    participant_inventory,
    read_region_manifest,
)

ZENODO_RECORD = "https://zenodo.org/records/10557593"
ZENODO_ARCHIVE_BYTES = 4_915_599_049
ZENODO_ARCHIVE_MD5 = "9f89acceda8232b5044460e23c5671cc"
SUPPLEMENT_URL = (
    "https://media.springernature.com/original/springer-static/esm/"
    "art%3A10.1038%2Fs42003-024-06281-8/MediaObjects/42003_2024_6281_MOESM1_ESM.pdf"
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--regions", default=str(REGION_MANIFEST))
    parser.add_argument("--obs", default=str(INTERIM_DIR / "icbi_obs.parquet"))
    parser.add_argument("--no-write", action="store_true")
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args()

    regions = read_region_manifest(args.regions)
    obs_path = Path(args.obs)
    if not obs_path.exists():
        raise FileNotFoundError(
            f"ICBI metadata not found at {obs_path}; run pull_icbi_metadata first"
        )
    obs = pd.read_parquet(obs_path, columns=["dataset", "patient_id", "sample_type"])
    participants = participant_inventory(regions, obs)
    counts = inventory_counts(regions, participants)

    print("=" * 72)
    print("A2 STEP 1 — PRODUCT AND CROSSWALK ONLY; NO MARKER OUTCOMES")
    print("=" * 72)
    for key, value in counts.items():
        print(f"  {key}: {value}")
    print("\n  public product: 5.6 um/pixel, 26-channel downsampled pixel NPZ")
    print("  arrays: img [height, width, 26], ch [26], mask [height, width]")
    print("  segmentation: none; mask is a tissue mask, not a cell mask")
    print("  crosswalk: participant exact, specimen equivalence NOT established")
    print("\n" + participants.to_string(index=False))
    print("\nDECISION")
    print("  Pixel-level reuse is feasible now. The planned crypt-position/per-cell")
    print("  equivalence test is not: the public product has no cell segmentation,")
    print("  and the public conventional-adenoma within-batch paired n is 3.")
    print("  Do not choose a replacement statistic until the Synapse MANDO-product")
    print("  inventory has been obtained or its absence confirmed.")

    if not args.no_write:
        meta = {
            "a2_stage": "step_1_product_inventory",
            "zenodo_record": ZENODO_RECORD,
            "zenodo_archive_bytes": ZENODO_ARCHIVE_BYTES,
            "zenodo_archive_md5": ZENODO_ARCHIVE_MD5,
            "supplement_url": SUPPLEMENT_URL,
            "processed_resolution_um_per_pixel": 5.6,
            "processed_channels": list(MXIF_CHANNELS),
            "has_cell_segmentation": False,
            "crosswalk_level": "participant_only",
            **counts,
        }
        for name, table in (
            ("a2_mxif_regions", regions),
            ("a2_mxif_participants", participants),
        ):
            path = write_versioned_table(
                table,
                name,
                seed=DEFAULT_SEED,
                notes="A2 Step 1 metadata inventory; no marker outcomes",
                extra_meta=meta,
                allow_dirty=args.allow_dirty,
            )
            print(f"\nwrote {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
