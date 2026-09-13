"""Where a result table resolves to the wrong run, and whether it matters.

    python -m src.reference.jobs.table_resolution_audit

THE FINDING. Every job's ``newest(name)`` is ``sorted(results/*/name.parquet)[-1]``,
and a result directory is ``{date}_{sha7}``. Sorting names orders runs by commit
sha, not by time, so **within one date the choice is arbitrary**. The repo
committed six Crowell multisection runs on 2026-09-07; the lexicographically last
is the second-earliest, and the canonical table is a different one written 2h28m
later. This reproduced a wrong comparison in the Crowell reproduction job before
it was fixed; the check that would have caught it did not exist.

This job enumerates the ambiguities and, more usefully, **measures whether each
one is consequential**: ``candidate_content`` reads both candidate runs and says
whether their frames are identical. An ambiguity where the two candidates hold
the same bytes is harmless; one where they differ is a number that could change
depending on which resolver a job happened to use.

It is a finding of the same family as ledger entry L24: a resolution that is
silent about which of two things it chose. Unlike L24 it changes no committed
number -- the affected jobs that matter were checked and the Crowell one is
fixed -- but nothing prevented that, which is the point.
"""

from __future__ import annotations

import argparse
import logging
import sys
from collections.abc import Sequence
from pathlib import Path

import pandas as pd

from src.common.io import write_versioned_table
from src.common.paths import RESULTS_DIR
from src.common.provenance import DEFAULT_SEED
from src.reference.table_resolution import ambiguity_audit

log = logging.getLogger(__name__)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    rows = ambiguity_audit(RESULTS_DIR)
    frame = pd.DataFrame(rows)
    if frame.empty:
        log.info("no ambiguous same-date resolution found")
        return 0

    consequential = frame[frame["disposition"] == "differs"]
    log.info("\n%s\nAMBIGUOUS SAME-DATE TABLE RESOLUTION\n%s", "=" * 72, "=" * 72)
    log.info("%d table/date pairs where name order and time order disagree; "
             "%d hold different content.", len(frame), len(consequential))
    log.info("%s", frame[[
        "table", "date", "n_same_date_runs", "name_order_picks",
        "time_order_picks", "disposition",
    ]].to_string(index=False))

    meta = {
        "purpose": (
            "Audit the resolver every job uses to pick a committed table. "
            "newest() sorts {date}_{sha7} directory names, so within a date it "
            "orders by sha, not time. This enumerates the disagreements and "
            "measures whether the two candidates differ in content."
        ),
        "finding": (
            "a name-based resolver is silent about which of two same-date runs "
            "it read. Same family as ledger L24; it reproduced a wrong Crowell "
            "comparison before newest_by_time was introduced."
        ),
        "consequential_definition": (
            "disposition='differs' iff the two candidate frames are not exactly "
            "equal including row order; 'identical' means the ambiguity is "
            "harmless for this table because both runs hold the same bytes."
        ),
        "n_ambiguous": int(len(frame)),
        "n_consequential": int(len(consequential)),
        "fix": (
            "src/reference/table_resolution.py: newest_by_time() resolves by the "
            "sidecar utc_timestamp with an mtime fallback. Jobs that compare to "
            "a committed table use it; newest() remains for callers that rely on "
            "name order."
        ),
        "exploratory": False,
    }
    path = write_versioned_table(
        frame, "table_resolution_ambiguity", seed=DEFAULT_SEED,
        results_dir=args.results_dir, extra_meta=meta, allow_dirty=args.allow_dirty,
    )
    log.info("wrote %s (%d rows)", path, len(frame))
    return 0


if __name__ == "__main__":
    sys.exit(main())
