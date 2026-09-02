"""Locate a committed result table by name.

The figure code reads every number from a versioned table. Hard-coding the
``{date}_{sha}`` directory means the figures break the moment the sweep is
re-run under a different commit — which is exactly what happens when a table is
re-derived to fix its provenance. So resolve by name, take the newest, and
*print which one was used*: the point is that the reader can trace the number to
a table, not that the path is a constant.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

RESULTS = Path("results")


def _written_at(table: Path) -> tuple[str, float]:
    """When the table was written, from its own sidecar.

    Sorting the directory names lexicographically orders correctly BY DATE and
    then arbitrarily *within* a date, because what follows the date is a commit
    sha and shas do not sort by time. Four runs on one day resolve to whichever
    sha sorts last, which is not the newest and is not reproducible from
    anything a reader can see. Every table carries a ``utc_timestamp`` in its
    sidecar; that is the field that means "when", so use it and fall back to the
    file's mtime only when the sidecar is missing.
    """
    sidecar = table.with_suffix(".meta.json")
    try:
        stamp = json.loads(sidecar.read_text()).get("utc_timestamp")
    except (OSError, ValueError):
        stamp = None
    return (stamp or "", table.stat().st_mtime)


def newest(name: str) -> Path:
    """Newest ``results/*/{name}.parquet``. Exits with a message if absent."""
    matches = list(RESULTS.glob(f"*/{name}.parquet"))
    if not matches:
        print(
            f"no results/*/{name}.parquet — run\n"
            f"    python -m src.harness.calibration_gap\n"
            f"first (add --replicates 500 for the r500 tables).",
            file=sys.stderr,
        )
        raise SystemExit(1)
    chosen = max(matches, key=_written_at)
    print(f"reading {chosen}")
    return chosen
