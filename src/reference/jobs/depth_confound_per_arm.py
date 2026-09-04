"""The depth-confound diagnostic, persisted PER ARM. Lee/SMC.

WHY THIS EXISTS. ``depth_confound_report`` already computes, for each arm, the
correlation between depth and the maturity call, that arm's prevalence, and the
bound $\\sqrt{3p(1-p)}$ that prevalence permits. The job that writes the shipped
table throws all of it away and keeps ``worst_rho``, the maximum over arms.

That single column is why the paper's own ceiling figure could not be drawn
honestly: pairing a max-over-arms correlation with one arm's prevalence is
exactly the mispairing the paper reports making three times, so the figure could
show the bound or the observed correlations but not both. It is also why the
claim that the diagnostic "persists both arms' correlation, bound, prevalence
and depth" was not true of any committed table.

This writes one row per (patient, axis, rung, arm), so ``rho`` and
``max_attainable_rho`` in a row always come from the same cells. Nothing
downstream has to pair them, which is the only fix that actually removes the
error rather than making it less likely.

WHAT IT DOES NOT DO. It does not change ``depth_confound_report``'s
max-over-arms summary fields. Those are still computed as independent maxima and
the paper says so. Removing them is a separate change with its own review; this
one adds the rows that make them unnecessary.

    python -m src.reference.jobs.depth_confound_per_arm
"""

from __future__ import annotations

import argparse
import logging
import sys
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd

from src.common.io import write_versioned_table
from src.common.provenance import DEFAULT_SEED
from src.harness.depth_confound import depth_confound_report

log = logging.getLogger(__name__)

AXES: tuple[str, ...] = ("stem_pole", "opposite_lineage")
UNRESOLVED = "unresolved_depth"


def mature_bin(rung: str) -> str:
    """The bin this rung calls mature.

    ``RUNG_SPECS`` orders each rung's bins least- to most-mature, so the mature
    one is the last. Hard-coding "differentiated" works only for ``lineage`` and
    silently gives every other rung a prevalence of zero, which then reads as a
    bound of zero and "the check cannot fire" — a wrong answer that looks exactly
    like this paper's finding. Read it off the spec.
    """
    from src.reference.labels import RUNG_SPECS

    return RUNG_SPECS[rung].bins[-1]


def per_arm_rows(
    depth: np.ndarray,
    call: np.ndarray,
    tissue: np.ndarray,
    *,
    patient_id: str,
    axis: str,
    rung: str,
) -> list[dict]:
    """One row per arm, with each arm's rho beside its own bound."""
    report = depth_confound_report(depth, call == mature_bin(rung), tissue)
    out = []
    for arm, stats in report["per_arm"].items():
        rho = float(stats["rho_depth_vs_mature"])
        ceiling = float(stats["max_attainable_rho"])
        out.append({
            "patient_id": patient_id,
            "labeling_axis": axis,
            "granularity_rung": rung,
            "mature_bin": mature_bin(rung),
            "arm": arm,
            "n_cells": int(stats["n_cells"]),
            "median_depth": float(stats["median_depth"]),
            "prevalence": float(stats["mature_share"]),
            "rho": rho,
            "abs_rho": abs(rho),
            "max_attainable_rho": ceiling,
            # Paired within the arm, by construction. This is the column the
            # shipped table could not carry.
            "rho_vs_ceiling": (
                abs(rho) / ceiling if np.isfinite(ceiling) and ceiling > 0
                else float("nan")
            ),
            "tolerance_is_reachable": bool(
                np.isfinite(ceiling) and ceiling >= 0.20
            ),
        })
    return out


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cohort", choices=("smc", "kul3"), default="smc",
        help="which Lee cohort. Both are local. The attainable bound is "
             "analytic and cohort-independent, but which rungs are degenerate "
             "at which prevalence is not, so a second cohort is a real test.",
    )
    parser.add_argument("--raw-dir", type=Path, default=None)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    from src.common.panel import granularity_rungs
    from src.estimator.lee_io import load_lee_cohort
    from src.reference.labels import label_column

    rungs = tuple(granularity_rungs())
    log.info("loading Lee/%s …", args.cohort.upper())
    cohort = load_lee_cohort(
        args.cohort, target_genes=["GUCA2A"], axes=AXES, rungs=rungs,
        raw_dir=args.raw_dir,
    )
    depth_all = cohort.n_counts
    tissue_all = cohort.cells["tissue"]
    patient_all = cohort.cells["patient_id"]

    rows: list[dict] = []
    for patient in sorted(patient_all.unique()):
        sel = patient_all == patient
        for axis in AXES:
            for rung in rungs:
                col = label_column(axis, rung)
                if col not in cohort.labels.columns:
                    continue
                call = cohort.labels[col].reindex(cohort.cells.index)
                scored = sel & call.notna() & (call.astype(str) != UNRESOLVED)
                idx = cohort.cells.index[scored]
                t = tissue_all.loc[idx].to_numpy()
                if len(idx) < 50 or len(set(t.tolist())) < 2:
                    continue
                rows.extend(per_arm_rows(
                    depth_all.reindex(idx).to_numpy(dtype=float),
                    call.loc[idx].astype(str).to_numpy(),
                    t,
                    patient_id=str(patient), axis=axis, rung=rung,
                ))

    table = pd.DataFrame(rows)
    if not table.empty:
        table.insert(0, "study_id", cohort.study_id)
    if table.empty:
        log.error("no rows — check the cohort loaded and the labels exist")
        return 1

    unreachable = (~table["tolerance_is_reachable"]).sum()
    log.info("%d rows | %d/%d cannot reach the 0.20 tolerance at their own "
             "prevalence", len(table), unreachable, len(table))
    for rung, g in table.groupby("granularity_rung"):
        log.info("  %-16s median p %.4f  median ceiling %.3f  "
                 "reachable %d/%d", rung, g["prevalence"].median(),
                 g["max_attainable_rho"].median(),
                 int(g["tolerance_is_reachable"].sum()), len(g))

    # smc keeps the committed name; any other cohort is namespaced so it
    # cannot overwrite it, and so newest("depth_confound_per_arm") still
    # resolves to the SMC table it has always resolved to.
    name = "depth_confound_per_arm"
    if args.cohort != "smc":
        name += f"_{args.cohort}"
    path = write_versioned_table(
        table, name, seed=DEFAULT_SEED,
        results_dir=args.results_dir, allow_dirty=args.allow_dirty,
        extra_meta={
            "cohort": (
                f"{cohort.study_id} (Lee/{args.cohort.upper()}), "
                f"{patient_all.nunique()} patients, sha256-verified against "
                f"data/manifest.csv"
            ),
            "rho_tolerance": 0.20,
            "grain": "one row per (patient, labeling_axis, granularity_rung, arm)",
            "what_this_answers": (
                "Each arm's depth-maturity correlation beside THAT ARM's "
                "attainable bound sqrt(3p(1-p)). The shipped table keeps only "
                "worst_rho, the max over arms, so pairing it with a prevalence "
                "requires guessing which arm the prevalence came from."
            ),
        },
    )
    log.info("wrote %s", path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
