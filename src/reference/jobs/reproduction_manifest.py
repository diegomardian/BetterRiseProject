"""Write the audit's reproduction index.

    python -m src.reference.jobs.reproduction_manifest

One row per deliverable in ``DECISION_2026-09-11_scope_and_pivot.md`` §4, with
the module that produces it, the newest versioned result, and that result's sha,
seed and row count. A deliverable that is blocked or absent is present=False
with the reason in ``note`` rather than dropped. This is the machine-readable
form of the "reproduction instructions" the plan asks for; the prose version
belongs in the paper and can be checked against this table.
"""

from __future__ import annotations

import argparse
import logging
import sys
from collections.abc import Sequence
from pathlib import Path

from src.common.io import write_versioned_table
from src.common.paths import RESULTS_DIR
from src.common.provenance import DEFAULT_SEED
from src.reference.reproduction import build_manifest, missing_deliverables

log = logging.getLogger(__name__)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    manifest = build_manifest(RESULTS_DIR)
    log.info("\n%s\nAUDIT REPRODUCTION INDEX\n%s", "=" * 72, "=" * 72)
    log.info("%s", manifest[[
        "week", "area", "table", "present", "git_sha", "seed", "n_rows",
    ]].to_string(index=False))

    missing = missing_deliverables(manifest)
    if missing.empty:
        log.info("\n  every deliverable resolves to a committed result")
    else:
        log.info("\n  %d deliverable(s) have no committed result yet:", len(missing))
        for row in missing.itertuples():
            log.info("    week %s %-28s %s", row.week, row.table, row.note)

    meta = {
        "purpose": (
            "DECISION_2026-09-11_scope_and_pivot.md §4 week 4: package the audit "
            "for another researcher, with executable examples and reproduction "
            "instructions."
        ),
        "what_this_is": (
            "a machine-readable index. Each row names a deliverable, the module "
            "that reproduces it, and the newest versioned result's sha, seed and "
            "row count, resolved by sidecar time rather than directory name."
        ),
        "what_this_is_not": (
            "a claim that a blocked deliverable is done. guard_challenge_results "
            "and the blinded challenge carry present=False and a note naming the "
            "unassigned role."
        ),
        "n_deliverables": int(len(manifest)),
        "n_missing": int(len(missing)),
        "exploratory": False,
    }
    path = write_versioned_table(
        manifest, "reproduction_manifest", seed=DEFAULT_SEED,
        results_dir=args.results_dir, extra_meta=meta, allow_dirty=args.allow_dirty,
    )
    log.info("wrote %s (%d rows)", path, len(manifest))
    return 0


if __name__ == "__main__":
    sys.exit(main())
