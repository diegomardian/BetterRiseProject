"""Reproduce the DIS/VAL split from committed inputs and check it.

    python -m src.reference.jobs.disval_reproduction

The ``disval_stability`` job reads ``data/interim/icbi_obs.parquet`` for the
patient-to-collection assignment, and that interim cache is not in the
checkout. Its own committed ``disval_stability_assignment`` table carries the
same assignment, so the split can be rebuilt from derived inputs alone. This job
does that, re-runs the same ``scale_free``/``contrasts_within`` the primary
reading uses, and compares the result to the committed ``disval_stability``
table cell by cell.

Labelled a reproduction from derived inputs, not a raw-data replication: the
atlas's per-cell obs are not read, only the assignment it produced.
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
from src.reference.disval_reproduction import compare_disval
from src.reference.jobs.adenoma_decomposition_scales import STATISTICS, scale_free
from src.reference.jobs.disval_stability import (
    PRIMARY_HALVES,
    RUNG,
    WEIGHTING,
    contrasts_within,
)
from src.reference.table_resolution import newest_by_time

log = logging.getLogger(__name__)


def _load(path: Path | None, name: str) -> pd.DataFrame:
    if path is None:
        raise SystemExit(f"no results/*/{name}.parquet. Run the producing job first.")
    frame = pd.read_parquet(path)
    log.info("%s | %d rows", path.parent.name + "/" + path.name, len(frame))
    return frame


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assignment", type=Path, default=None)
    parser.add_argument("--decomposition", type=Path, default=None)
    parser.add_argument("--committed", type=Path, default=None)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    assignment_path = args.assignment or newest_by_time(
        RESULTS_DIR, "disval_stability_assignment")
    decomposition_path = args.decomposition or newest_by_time(
        RESULTS_DIR, "adenoma_decomposition")
    committed_path = args.committed or newest_by_time(
        RESULTS_DIR, "disval_stability")

    assignment = _load(assignment_path, "disval_stability_assignment")
    decomposition = _load(decomposition_path, "adenoma_decomposition")
    committed = _load(committed_path, "disval_stability")

    block = decomposition[
        (decomposition["granularity_rung"] == RUNG)
        & (decomposition["weighting"] == WEIGHTING)
    ]
    if block.empty:
        raise SystemExit(f"no {RUNG}/{WEIGHTING} rows in the decomposition")

    usable = assignment[~assignment["shared"]]
    values = scale_free(block).merge(
        usable[["patient_id", "collection"]], on="patient_id", how="inner"
    )
    if values.empty:
        raise SystemExit(
            "the assignment and the decomposition share no patients. They write "
            "the id as 'Chen_2021_Cell.HTA11_...'; if one side changed shape "
            "this is an identifier-space mismatch, not an empty cohort."
        )

    frames = []
    for label in (*PRIMARY_HALVES, "VUMC_HTAN_cohort3", "POOLED"):
        subset = values if label == "POOLED" else values[values["collection"] == label]
        if subset.empty:
            continue
        for statistic in STATISTICS:
            frame = contrasts_within(
                subset, statistic=statistic, seed=args.seed, label=label
            )
            if not frame.empty:
                frames.append(frame)
    reproduced = pd.concat(frames, ignore_index=True)

    comparison = compare_disval(reproduced, committed)
    load = comparison[comparison["statistic"] == "log_ratio"]
    log.info("\n%s\nDIS/VAL REPRODUCTION FROM DERIVED INPUTS — load-bearing\n%s",
             "=" * 72, "=" * 72)
    log.info("%s", load[[
        "half", "contrast", "n_patients", "mean", "mean_committed",
        "mean_diff", "reproduces",
    ]].to_string(index=False))

    failed = comparison[~comparison["reproduces"]]
    if failed.empty:
        log.info("\n  All %d cells reproduce exactly. Reproduction from derived "
                 "inputs,\n  not a replication: the atlas obs were not read.",
                 len(comparison))
    else:
        log.warning("\n  %d cell(s) DID NOT reproduce; first: %s",
                    len(failed),
                    failed[["half", "statistic", "contrast"]].head(3)
                    .to_dict("records"))

    meta = {
        "purpose": (
            "DECISION_2026-09-11_scope_and_pivot.md §4 week 3: reproduce DIS/VAL "
            "without opening new endpoints."
        ),
        "provenance_kind": (
            "reproduction from derived inputs, NOT raw-data replication. The job "
            "reads the committed disval_stability_assignment and the committed "
            "adenoma decomposition; the interim atlas obs are not read."
        ),
        "sources": {
            "assignment": str(assignment_path.relative_to(RESULTS_DIR.parent))
            if assignment_path.is_relative_to(RESULTS_DIR.parent) else str(assignment_path),
            "decomposition": str(decomposition_path.relative_to(RESULTS_DIR.parent))
            if decomposition_path.is_relative_to(RESULTS_DIR.parent) else str(decomposition_path),
            "committed_disval": str(committed_path.relative_to(RESULTS_DIR.parent))
            if committed_path.is_relative_to(RESULTS_DIR.parent) else str(committed_path),
        },
        "rung": RUNG,
        "weighting": WEIGHTING,
        "halves": [*PRIMARY_HALVES, "VUMC_HTAN_cohort3", "POOLED"],
        "all_cells_reproduce": bool(comparison["reproduces"].all()),
        "failed_cells": failed[["half", "statistic", "contrast"]]
        .to_dict("records") if not failed.empty else [],
        "comparison_tolerance": 1e-9,
        "what_this_is_not": (
            "an independent replication. Same lab, platform and population; the "
            "single-cohort qualifier on the adenoma result stays."
        ),
        "exploratory": False,
    }
    path = write_versioned_table(
        comparison, "disval_reproduction", seed=args.seed,
        results_dir=args.results_dir, extra_meta=meta, allow_dirty=args.allow_dirty,
    )
    log.info("wrote %s (%d rows)", path, len(comparison))
    return 0


if __name__ == "__main__":
    sys.exit(main())
