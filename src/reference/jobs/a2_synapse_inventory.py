#!/usr/bin/env python
"""Inventory Chen Synapse entities after access is granted; download no images.

    pip install -e '.[a2]'
    synapse login --rememberMe
    python -m src.reference.jobs.a2_synapse_inventory --no-write

The first run is deliberately ``--no-write``: review the complete metadata
inventory and candidate list before making it a provenance-stamped artifact.
This job does not call ``syn.get`` without ``downloadFile=False`` and never
calls a sync/download method.
"""

from __future__ import annotations

import argparse
import json
import sys

from src.common.io import write_versioned_table
from src.common.provenance import DEFAULT_SEED
from src.reference.a2_synapse import (
    CHEN_SYNAPSE_ROOTS,
    A2SynapseError,
    login_synapse,
    walk_metadata,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-write", action="store_true")
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args()

    try:
        table = walk_metadata(login_synapse())
    except A2SynapseError as exc:
        parser.error(str(exc))
    candidates = table[table["candidate_cell_product"]].copy()
    print("=" * 72)
    print("A2 SYNAPSE INVENTORY — METADATA ONLY; NO FILE CONTENT DOWNLOADED")
    print("=" * 72)
    print(f"entities: {len(table)}")
    print(f"candidate cell-derived products for manual review: {len(candidates)}")
    print(candidates.drop(columns="annotations_json").to_string(index=False))
    print("\nDECISION GATE")
    print("  A candidate name is not evidence of a usable MANDO table. Confirm")
    print("  cell IDs, coordinates, per-cell intensities, and segmentation provenance")
    print("  before choosing a spatial statistic or downloading any file.")

    if not args.no_write:
        serialised = table.copy()
        serialised["annotations_json"] = serialised["annotations_json"].map(json.dumps)
        common = {
            "a2_stage": "synapse_metadata_inventory",
            "source_entities": list(CHEN_SYNAPSE_ROOTS),
            "file_content_downloaded": False,
            "candidate_rule": "filename or annotation keyword screen; manual verification required",
        }
        for name, frame in (
            ("a2_synapse_metadata", serialised),
            ("a2_synapse_candidate_products", serialised[serialised["candidate_cell_product"]]),
        ):
            path = write_versioned_table(
                frame,
                name,
                seed=DEFAULT_SEED,
                notes="A2 access-gated Synapse metadata inventory; no file contents downloaded",
                extra_meta=common,
                allow_dirty=args.allow_dirty,
            )
            print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
