"""Locate a committed result table by name.

The figure code reads every number from a versioned table. Hard-coding the
``{date}_{sha}`` directory means the figures break the moment the sweep is
re-run under a different commit — which is exactly what happens when a table is
re-derived to fix its provenance. So resolve by name, take the newest, and
*print which one was used*: the point is that the reader can trace the number to
a table, not that the path is a constant.
"""

from __future__ import annotations

import sys
from pathlib import Path

RESULTS = Path("results")


def newest(name: str) -> Path:
    """Newest ``results/*/{name}.parquet``. Exits with a message if absent."""
    matches = sorted(RESULTS.glob(f"*/{name}.parquet"))
    if not matches:
        print(
            f"no results/*/{name}.parquet — run\n"
            f"    python -m src.harness.calibration_gap\n"
            f"first (add --replicates 500 for the r500 tables).",
            file=sys.stderr,
        )
        raise SystemExit(1)
    chosen = matches[-1]
    print(f"reading {chosen}")
    return chosen
