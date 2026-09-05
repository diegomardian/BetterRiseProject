#!/usr/bin/env python
"""Pull the ICBI atlas metadata — not the atlas. W1, week 1.

    python -m src.reference.jobs.pull_icbi_metadata

execution_plan.md §8.2 says to take the metadata table before committing any
compute. No such table is published, so this reads /obs out of the 32.7 GB h5ad
over HTTP range requests and touches roughly 0.1% of the object.

Answers four live questions in one pull:

  platform mix        is there a plate-based subset, and how big (decisions
                      #8 and #14)
  per-cell depth      does it actually resolve the sparse axis-1 markers, or is
                      that just the protocol's reputation
  paired samples      the real sample size, against 36 of 62 here (#9)
  enrichment          unsorted vs CD45-sorted at atlas scale (#11)

Writes results/icbi_metadata_summary/ and caches the obs frame so reruns are free.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from src.common.io import write_versioned_table
from src.common.paths import INTERIM_DIR
from src.common.provenance import DEFAULT_SEED
from src.reference.icbi import (
    ATLAS_URL,
    depth_by_platform,
    enrichment_summary,
    epithelial_fraction,
    paired_sample_summary,
    platform_summary,
    read_atlas_obs,
)

#: Set from --allow-dirty. Default False -- this job used to pass
#: allow_dirty=True unconditionally, so it could not write a clean stamp.
ALLOW_DIRTY = False


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--url", default=ATLAS_URL,
        help="the atlas URL, or a LOCAL PATH to an already-fetched copy. "
             "Range requests exist because the object is 32.7 GB and /obs is "
             "0.1%% of it; once it is on disk that argument is gone.",
    )
    # INTERIM_DIR, not a path relative to the repo. Those are the SAME
    # directory when BRP_DATA_DIR is unset -- which is why this worked on a
    # laptop for weeks -- and different the moment it points outside the repo,
    # as it does on the cluster. The puller then wrote one place and every
    # consumer looked in another.
    parser.add_argument("--cache", default=str(INTERIM_DIR / "icbi_obs.parquet"))
    parser.add_argument("--refresh", action="store_true")
    parser.add_argument(
        "--allow-dirty", action="store_true",
        help="write from a dirty tree; the recorded sha will not reproduce it",
    )
    args = parser.parse_args()
    global ALLOW_DIRTY
    ALLOW_DIRTY = args.allow_dirty

    cache = Path(args.cache)
    if cache.exists() and not args.refresh:
        print(f"using cached obs: {cache}")
        obs = pd.read_parquet(cache)
        report = {"n_cells": len(obs), "bytes_fetched": 0, "cached": True}
    else:
        print(f"reading /obs from {args.url}")
        print("(range requests — the 32.7 GB object is not downloaded)")
        obs, report = read_atlas_obs(args.url)
        cache.parent.mkdir(parents=True, exist_ok=True)
        obs.to_parquet(cache)
        print(
            f"\nfetched {report['bytes_fetched']:,} bytes of "
            f"{report['object_size']:,} ({report['fraction_fetched']:.3%})"
        )
        if report["columns_missing"]:
            print(f"note: columns absent from this build: {report['columns_missing']}")

    print(f"\n{report['n_cells']:,} cells")

    sections = {
        "icbi_platform": ("PLATFORM MIX — is there a plate-based subset?", platform_summary),
        "icbi_depth": (
            "PER-CELL DEPTH — does plate-based actually resolve sparse markers?",
            depth_by_platform,
        ),
        "icbi_paired": (
            "PAIRED TUMOUR/NORMAL PER STUDY — the real sample size",
            paired_sample_summary,
        ),
        "icbi_enrichment": ("ENRICHMENT — unsorted vs sorted (decision #11)", enrichment_summary),
        "icbi_epithelial": ("EPITHELIAL FRACTION BY STUDY", epithelial_fraction),
    }

    outputs: dict[str, pd.DataFrame] = {}
    for name, (title, function) in sections.items():
        print("\n" + "=" * 70)
        print(title)
        print("=" * 70)
        try:
            table = function(obs)
        except Exception as exc:
            print(f"  unavailable: {exc}")
            continue
        outputs[name] = table
        print(table.head(25).to_string(index=False))
        if len(table) > 25:
            print(f"  ... {len(table) - 25} more rows in the written table")

    print("\n" + "=" * 70)
    print("READ THIS")
    print("=" * 70)
    if "icbi_depth" in outputs:
        depth = outputs["icbi_depth"]
        plate = depth[depth["plate_based"]]
        droplet = depth[~depth["plate_based"]]
        if len(plate) and len(droplet):
            ratio = float(plate["median_genes"].max() / droplet["median_genes"].max())
            print(f"  plate-based median genes/cell is {ratio:.1f}x the best droplet"
                  f" platform.")
            print("  GSE178341's epithelium sits around 2,000-3,300 genes and leaves")
            print("  49% of cells tied on axis 1. If the plate subset clears ~5,000,")
            print("  axis 1 becomes measurable there and decision #14 has an answer")
            print("  that does not require pulling axis 3 forward from week 13+.")
    if "icbi_platform" in outputs:
        plate_cells = int(
            outputs["icbi_platform"].query("plate_based")["n_cells"].sum()
        )
        print(f"\n  plate-based cells available: {plate_cells:,}")
        print("  Plate protocols have essentially no ambient soup, so this is also")
        print("  G1's fallback after decision #8 — an intrinsic signal surviving")
        print("  there is strong evidence it is not contamination (§8.2).")

    for name, table in outputs.items():
        path = write_versioned_table(
            table, name, seed=DEFAULT_SEED, allow_dirty=ALLOW_DIRTY,
            notes="ICBI CRC atlas metadata, read via HTTP range requests",
        )
        print(f"\n  {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
