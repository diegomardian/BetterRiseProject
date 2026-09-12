"""Emit the workshop-overlap and mathematical-claims matrix.

    python -m src.reference.jobs.workshop_overlap --no-write

Week 1, deliverable 3. Reads ``config/workshop_overlap.yaml`` and nothing else.
"""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from pathlib import Path

from src.common.io import write_versioned_table
from src.common.provenance import DEFAULT_SEED
from src.reference.workshop_overlap import (
    OVERLAP_PATH,
    build_overlap,
    load_overlap,
    overlap_summary,
)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--overlap", type=Path, default=OVERLAP_PATH)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--no-write", action="store_true")
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)

    frame = build_overlap(load_overlap(args.overlap))
    summary = overlap_summary(frame)

    print(f"{summary['n_rows']} rows: {summary['n_overlap_rows']} overlap, "
          f"{summary['n_math_claims']} mathematical claims")
    print(f"  already published: {summary['n_already_published']}")
    print(f"  extended:          {summary['n_extended']}")
    print(f"  new:               {summary['n_new']}")
    print(f"  contradicted:      {summary['n_contradicted']} (each needs a paper edit)")
    print(f"  novelty share of the overlap rows: {summary['novelty_share_of_overlap_rows']}")

    if args.no_write:
        print("--no-write: nothing written")
        return 0

    path = write_versioned_table(
        frame,
        "workshop_overlap",
        seed=args.seed,
        results_dir=args.results_dir,
        allow_dirty=args.allow_dirty,
        notes=(
            "What the proposed BMC Bioinformatics submission would claim against "
            "what the WMHS workshop paper already claims, plus an audit of the "
            "mathematical claims. Week 1 deliverable 3."
        ),
        extra_meta={
            **summary,
            "what_this_licenses": (
                "a precise statement of overlap for a cover letter and a related-work "
                "section, and a list of paper edits the mathematical audit requires"
            ),
            "what_this_does_not_license": (
                "a claim that the new paper is sufficiently novel; that is an "
                "editor's judgement, not a table's"
            ),
        },
    )
    print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
