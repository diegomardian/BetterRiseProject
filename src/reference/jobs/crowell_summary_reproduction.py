"""The Crowell summary, recomputed from committed derived inputs and checked.

    python -m src.reference.jobs.crowell_summary_reproduction

WHY THIS IS A REPRODUCTION RATHER THAN A RE-RUN. The Crowell h5ads were deleted
after the pre-registered run (disk protocol, md5s recorded in
``data/manifest.csv``), so what is available are the committed per-section
``crowell_feasibility`` tables. This job feeds those through the same
``per_block_did``/``aggregate`` the primary reading uses, then compares the
result to the committed ``crowell_multisection_summary`` gene by gene. The
decision record requires this be labelled **reproduction from derived inputs,
not raw-data replication**, and the sidecar says so.

The output is the summary with both columns and a ``reproduces`` flag. If the
flag is False for any gene, the reproduction failed and the committed reading is
not confirmed from derived inputs -- which is a finding, not a crash.
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
from src.reference.crowell_reproduction import compare_summaries, newest_by_time
from src.reference.jobs.crowell_multisection import aggregate, per_block_did

log = logging.getLogger(__name__)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--committed", type=Path, default=None,
                        help="committed crowell_multisection_summary to compare to")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    committed_path = args.committed or newest_by_time(
        RESULTS_DIR, "crowell_multisection_summary"
    )
    if committed_path is None:
        raise SystemExit("no results/*/crowell_multisection_summary.parquet to "
                         "compare against")
    committed = pd.read_parquet(committed_path)
    log.info("committed summary <- %s (%d genes)",
             committed_path.parent.name + "/" + committed_path.name, len(committed))

    per_block = per_block_did(RESULTS_DIR, lesion_class="adenoma")
    if per_block.empty:
        raise SystemExit(
            "no usable per-block rows. The committed crowell_feasibility "
            "tables carry the sidecars this needs; check they are present."
        )
    reproduced = aggregate(per_block)
    log.info("reproduced %d blocks x %d rows from committed feasibility tables",
             per_block["block"].nunique(), len(per_block))

    comparison = compare_summaries(reproduced, committed)
    log.info("\n%s\nREPRODUCTION OF THE CROWELL SUMMARY FROM DERIVED INPUTS\n%s",
             "=" * 72, "=" * 72)
    log.info("%s", comparison[[
        "gene", "n_blocks", "mean_did", "mean_did_committed", "mean_did_diff",
        "reproduces",
    ]].to_string(index=False))

    failed = comparison[~comparison["reproduces"]]
    if failed.empty:
        log.info("\n  All %d genes reproduce exactly. This is a reproduction from "
                 "derived inputs,\n  not a replication: the raw h5ads are not on "
                 "disk and were not read.", len(comparison))
    else:
        log.warning(
            "\n  %d gene(s) DID NOT reproduce: %s. The committed summary is not "
            "confirmed\n  from derived inputs for these.",
            len(failed), ", ".join(failed["gene"]))

    meta = {
        "purpose": (
            "DECISION_2026-09-11_scope_and_pivot.md §4 week 3: reproduce Crowell "
            "without opening new endpoints."
        ),
        "provenance_kind": (
            "reproduction from derived inputs, NOT raw-data replication. The "
            "Crowell h5ads were deleted per the prereg's disk protocol; the "
            "inputs here are the committed per-section crowell_feasibility "
            "tables and their sidecars."
        ),
        "committed_summary": str(
            committed_path.relative_to(RESULTS_DIR.parent)
            if committed_path.is_relative_to(RESULTS_DIR.parent)
            else committed_path
        ),
        "n_input_sections": int(per_block["section"].nunique()),
        "n_blocks": int(per_block["block"].nunique()),
        "all_genes_reproduce": bool(comparison["reproduces"].all()),
        "failed_genes": sorted(failed["gene"].tolist()),
        "comparison_tolerance": 1e-9,
        "why_this_is_not_a_replication": (
            "the same committed per-section tables are read back through the "
            "same DiD and Student-t aggregation. Agreement shows the committed "
            "summary is a deterministic function of committed inputs; it does "
            "not re-measure the raw CosMx data."
        ),
        "exploratory": False,
    }
    path = write_versioned_table(
        comparison, "crowell_summary_reproduction", seed=args.seed,
        results_dir=args.results_dir, extra_meta=meta, allow_dirty=args.allow_dirty,
    )
    log.info("wrote %s (%d rows)", path, len(comparison))
    # A reproduction mismatch is recorded in the artifact and logged, not a
    # process failure: the point is the table, and a job that refused to write
    # one on mismatch would leave no evidence of the mismatch.
    return 0


if __name__ == "__main__":
    sys.exit(main())
