"""Locate a committed result table by name, for this manuscript's numbers.

Deliberately NOT a copy of ``paper/wmhs/_tables.py``. That module resolves with
``sorted(glob(...))[-1]``, which orders ``{date}_{sha7}`` directories by commit
sha rather than by time; within a single date that is arbitrary. The audit found
23 table names with ambiguous same-date runs, 9 of them holding different
content, and the name-order rule produced a wrong Crowell comparison before it
was caught.

So this resolves by the sidecar's ``utc_timestamp`` via
``src.reference.table_resolution.newest_by_time`` and prints the path it used.
Neither rule is universally correct — that is the audit's actual finding — so
when a section needs the name-order run it must say so at the call site rather
than swapping the default.
"""

from __future__ import annotations

import sys
from pathlib import Path

from src.reference.table_resolution import newest_by_time

RESULTS = Path("results")


def table(name: str) -> Path:
    """Most recently written ``results/*/{name}.parquet``.

    Exits rather than falling back to a stale run: a manuscript number that
    silently reads an older table is the failure this paper is about.
    """
    chosen = newest_by_time(RESULTS, name)
    if chosen is None:
        print(
            f"no results/*/{name}.parquet\n"
            f"  the producing command is in the `command` column of\n"
            f"  results/*/reproduction_manifest.parquet",
            file=sys.stderr,
        )
        raise SystemExit(1)
    print(f"reading {chosen}")
    return chosen
