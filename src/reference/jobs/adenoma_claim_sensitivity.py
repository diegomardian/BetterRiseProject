"""The adenoma reading under every fixed choice, from committed derived inputs.

    python -m src.reference.jobs.adenoma_claim_sensitivity

A THIN LOCAL READ. No cluster, no raw data, no new endpoint: both decompositions,
the compositional estimability table and the scale-free statistics are already
committed, and this job assembles them into the artifact the second reviewer
objection asks for (``DECISION_2026-09-11_scope_and_pivot.md`` §3).

It writes three tables:

``adenoma_claim_sensitivity``
    Every contrast on every statistic, under **both** denominators and all three
    weightings, with patient counts and missingness on each row.

``estimability_attrition``
    Where patients leave: the intrinsic gate on ``n_cells_mature`` and the
    compositional gate on ``n_cells_resolved``, counted separately per gene.

``adenoma_claim_sensitivity_summary``
    The cross-block count that carries the two-block claim, per fixed-choice
    cell, so a reader can see which choices it survives under.

REPRODUCTION, NOT REPLICATION. These summaries are recomputed from derived
per-patient tables already in git. They are a **reproduction from derived
inputs**, not a replication from raw data, and the sidecar says so. Nothing here
opens a new endpoint.
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
from src.reference.claim_sensitivity import (
    contrasts_by_denominator,
    cross_block_survival,
    estimability_attrition,
)

log = logging.getLogger(__name__)


def newest(name: str) -> Path | None:
    """Newest committed table by directory name (``{date}_{sha}``), not mtime."""
    matches = sorted(RESULTS_DIR.glob(f"*/{name}.parquet"))
    return matches[-1] if matches else None


def _load(path: Path | None, name: str) -> pd.DataFrame:
    if path is None:
        raise SystemExit(f"no results/*/{name}.parquet. Run the producing job first.")
    frame = pd.read_parquet(path)
    log.info("%s | %d rows", path.parent.name + "/" + path.name, len(frame))
    return frame


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, default=None)
    parser.add_argument("--all-epithelial", type=Path, default=None)
    parser.add_argument("--compositional", type=Path, default=None)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    primary_path = args.primary or newest("adenoma_decomposition")
    epithelial_path = args.all_epithelial or newest(
        "adenoma_decomposition_all_epithelial"
    )
    compositional_path = args.compositional or newest(
        "adenoma_decomposition_compositional_estimability"
    )
    primary = _load(primary_path, "adenoma_decomposition")
    all_epithelial = _load(epithelial_path, "adenoma_decomposition_all_epithelial")
    compositional = _load(
        compositional_path, "adenoma_decomposition_compositional_estimability"
    )

    sensitivity = contrasts_by_denominator(
        primary, all_epithelial, seed=args.seed
    )
    attrition = estimability_attrition(primary, all_epithelial, compositional)
    summary = cross_block_survival(sensitivity)

    log.info("\n%s\nCROSS-BLOCK CONTRASTS EXCLUDING ZERO, by fixed-choice cell\n%s",
             "=" * 72, "=" * 72)
    shown = summary[summary["granularity_rung"].isin(["lineage", "best4"])]
    log.info("%s", shown.to_string(index=False))
    log.info("\n  `load_bearing` marks the pre-registered statistic (log_ratio). "
             "A cell where the\n  cross-block count collapses is a claim that "
             "rests on that choice.")

    log.info("\n%s\nWHERE PATIENTS LEAVE (distinct patients per rung)\n%s",
             "=" * 72, "=" * 72)
    per_rung = (
        attrition.groupby(["denominator", "granularity_rung"], observed=True)
        .agg(
            n_patients=("n_patients", "max"),
            intrinsic_missing=("n_intrinsic_missing", "sum"),
            not_estimable=("n_estimability_not_estimable", "max"),
            log_ratio_undefined=("n_log_ratio_undefined", "sum"),
        )
        .reset_index()
    )
    log.info("%s", per_rung.to_string(index=False))

    meta = {
        "purpose": (
            "DECISION_2026-09-11_scope_and_pivot.md §4 week 3: recompute the "
            "Chen adenoma summaries across all weightings, rungs, both "
            "denominators and all four scale-free statistics, with patient "
            "counts and missingness, to answer the second reviewer objection."
        ),
        "provenance_kind": (
            "reproduction from derived inputs, NOT a replication from raw data. "
            "Every row is recomputed from committed per-patient tables already "
            "in git; no new endpoint is opened."
        ),
        "sources": {
            "primary": str(primary_path.relative_to(RESULTS_DIR.parent))
            if primary_path.is_relative_to(RESULTS_DIR.parent) else str(primary_path),
            "all_epithelial": str(epithelial_path.relative_to(RESULTS_DIR.parent))
            if epithelial_path.is_relative_to(RESULTS_DIR.parent) else str(epithelial_path),
            "compositional": str(compositional_path.relative_to(RESULTS_DIR.parent))
            if compositional_path.is_relative_to(RESULTS_DIR.parent) else str(compositional_path),
        },
        "denominators": {
            "resolved": "mature cells over cells resolved to a maturity bin "
                        "(decision #14 primary)",
            "all_epithelial": "mature cells over all epithelial cells "
                              "(the denominator #14 rejected)",
        },
        "weightings": sorted(primary["weighting"].unique()),
        "statistics": sorted(sensitivity["statistic"].unique()),
        "rungs": sorted(primary["granularity_rung"].unique()),
        "not_included": (
            "detection-scale statistics (cloglog, log2_cp10k, detection) are a "
            "different estimand on a different scale and are reported in "
            "adenoma_specificity_disagreements.parquet. Including them here "
            "would compare two estimands and call it a sensitivity analysis."
        ),
        "none_is_not_zero": (
            "a None intrinsic term is dropped before any contrast and counted in "
            "estimability_attrition; it never enters a numerator as zero"
        ),
        "exploratory": False,
    }
    for frame, name in (
        (sensitivity, "adenoma_claim_sensitivity"),
        (attrition, "estimability_attrition"),
        (summary, "adenoma_claim_sensitivity_summary"),
    ):
        path = write_versioned_table(
            frame, name, seed=args.seed, results_dir=args.results_dir,
            extra_meta=meta, allow_dirty=args.allow_dirty,
        )
        log.info("wrote %s (%d rows)", path, len(frame))
    return 0


if __name__ == "__main__":
    sys.exit(main())
