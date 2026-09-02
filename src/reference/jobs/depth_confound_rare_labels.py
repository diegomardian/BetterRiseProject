"""The prevalence bound, on labels rare enough for it to bite. Both Lee cohorts.

WHY THIS EXISTS. ``depth_confound_per_arm`` establishes that a rank correlation
against a binary label of prevalence $p$ cannot exceed $\\sqrt{3p(1-p)}$, so a
fixed tolerance of 0.20 is unreachable below $p = 1.3516\\%$ and the diagnostic
reports clean whatever the data does. On this project's own maturity labels that
is a theorem with nowhere to land: **no row falls in $0 < p < 1.3516\\%$**. The
58 rows that cannot reach the tolerance are all degenerate, at prevalence
exactly 0 or exactly 1, where "a threshold cannot fire" is arithmetic rather
than a finding.

The regime is empty in those labels, not in this data. Both cohorts ship a
published per-cell subtype annotation — 36 subtypes on SMC, 40 on KUL3 — and a
dozen of them sit under 1.35% prevalence, including the mature enterocyte and
goblet populations this decomposition is about. Scoring those the same way puts
several hundred rows squarely inside the band, with a large share of them
strongly confounded and none of them flagged.

WHAT MAKES THE COMPARISON FAIR. Every row here comes from
``depth_confound_report`` — the same function, the same Spearman, the same
bound, on the same QC-filtered cells and the same whole-transcriptome library
size that ``depth_confound_per_arm`` uses. The only thing that changes is which
binary label is being scored. So the two tables are directly comparable and the
difference between them is prevalence and nothing else.

WHAT IT DOES NOT CLAIM. These are the authors' published annotations, not this
project's maturity call. That is deliberate: the claim under test is about the
*statistic*, not about our labelling, and demonstrating it on labels we did not
build is the stronger version. It says nothing about whether those annotations
are correct.

    python -m src.reference.jobs.depth_confound_rare_labels
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
from src.harness.depth_confound import (
    MATURITY_DEPTH_RHO_TOLERANCE,
    depth_confound_report,
)

log = logging.getLogger(__name__)

#: Prevalence at which sqrt(3p(1-p)) crosses the 0.20 tolerance. Below this the
#: check cannot fire on any data. Solved in closed form, asserted in tests.
CROSSING_PREVALENCE = 0.0135166

#: Cells an arm needs before its correlation is worth reporting. A rank
#: correlation over a handful of cells is noise, and admitting those rows would
#: pad the blind band with rows that are uninformative for a second reason.
MIN_CELLS_PER_ARM = 100

COHORTS: tuple[tuple[str, str], ...] = (
    ("smc", "GSE132465"),
    ("kul3", "GSE144735"),
)


def rows_for_label(
    depth: np.ndarray,
    is_member: np.ndarray,
    arm: np.ndarray,
    *,
    study_id: str,
    patient_id: str,
    subtype: str,
) -> list[dict]:
    """One row per arm for a single subtype, rho beside that arm's own bound."""
    report = depth_confound_report(depth, is_member, arm)
    out = []
    for arm_name, stats in report["per_arm"].items():
        if int(stats["n_cells"]) < MIN_CELLS_PER_ARM:
            continue
        rho = float(stats["rho_depth_vs_mature"])
        prevalence = float(stats["mature_share"])
        ceiling = float(stats["max_attainable_rho"])
        out.append({
            "study_id": study_id,
            "patient_id": patient_id,
            "label_source": "author_cell_subtype",
            "subtype": subtype,
            "arm": arm_name,
            "n_cells": int(stats["n_cells"]),
            "median_depth": float(stats["median_depth"]),
            "prevalence": prevalence,
            "rho": rho,
            "abs_rho": abs(rho),
            "max_attainable_rho": ceiling,
            "rho_vs_ceiling": float(stats["rho_vs_ceiling"]),
            "tolerance_is_reachable": bool(stats["tolerance_is_reachable"]),
            # Degenerate (p exactly 0 or 1) and rare (0 < p < crossing) are
            # different findings and the paper reports them separately: a
            # threshold that cannot fire on a label with no variance is
            # arithmetic, one that cannot fire on a rare label is the blind
            # spot. Keeping them in one column would merge the two again.
            "is_degenerate": bool(prevalence in (0.0, 1.0)),
            "in_blind_band": bool(0.0 < prevalence < CROSSING_PREVALENCE),
            "flagged_by_tolerance": bool(
                np.isfinite(rho) and abs(rho) >= MATURITY_DEPTH_RHO_TOLERANCE
            ),
        })
    return out


def collect(which: str, study_id: str, raw_dir: Path | None) -> pd.DataFrame:
    from src.estimator.lee_io import load_lee_cohort

    log.info("loading %s (%s) …", which, study_id)
    cohort = load_lee_cohort(
        which, target_genes=["GUCA2A"], raw_dir=raw_dir, label_compartment=None,
    )
    cells = cohort.cells
    depth = cohort.n_counts.reindex(cells.index).to_numpy(dtype=float)
    arm = cells["tissue"].to_numpy()
    patient = cells["patient_id"]
    subtype = cells["author_cell_subtype"].astype(str)

    rows: list[dict] = []
    for pid in sorted(patient.unique()):
        sel = (patient == pid).to_numpy()
        if len(set(arm[sel].tolist())) < 2:
            continue
        for st in sorted(subtype[sel].unique()):
            rows.extend(rows_for_label(
                depth[sel], (subtype.to_numpy()[sel] == st), arm[sel],
                study_id=study_id, patient_id=str(pid), subtype=st,
            ))
    return pd.DataFrame(rows)


def summarise(table: pd.DataFrame) -> None:
    """What the table says, printed so the run itself states the finding."""
    defined = table[np.isfinite(table["rho"])]
    band = defined[defined["in_blind_band"]]
    visible = defined[defined["tolerance_is_reachable"]]

    log.info("%d rows, %d carrying a correlation", len(table), len(defined))
    log.info(
        "  blind band 0 < p < %.4f%%: %d rows | median |rho|/bound %.3f | "
        "flagged %d | above half their bound %d | max ratio %.3f",
        CROSSING_PREVALENCE * 100, len(band), band["rho_vs_ceiling"].median(),
        int(band["flagged_by_tolerance"].sum()),
        int((band["rho_vs_ceiling"] >= 0.5).sum()),
        band["rho_vs_ceiling"].max(),
    )
    log.info(
        "  can fire (p >= %.4f%%):    %d rows | median |rho|/bound %.3f | "
        "flagged %d (%.0f%%)",
        CROSSING_PREVALENCE * 100, len(visible),
        visible["rho_vs_ceiling"].median(),
        int(visible["flagged_by_tolerance"].sum()),
        100 * visible["flagged_by_tolerance"].mean() if len(visible) else 0.0,
    )
    if int(band["flagged_by_tolerance"].sum()):
        log.error("  a row in the blind band was flagged — that is impossible; "
                  "check the bound")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", type=Path, default=None)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    frames = []
    for which, study_id in COHORTS:
        try:
            frames.append(collect(which, study_id, args.raw_dir))
        except (FileNotFoundError, OSError) as exc:
            log.warning("skipping %s (%s): %s", which, study_id, exc)

    if not frames or all(f.empty for f in frames):
        log.error("no rows — check the Lee matrices are present under data/raw")
        return 1

    table = pd.concat(frames, ignore_index=True)
    for study_id, g in table.groupby("study_id"):
        log.info("%s:", study_id)
        summarise(g)
    log.info("both cohorts:")
    summarise(table)

    path = write_versioned_table(
        table, "depth_confound_rare_labels", seed=DEFAULT_SEED,
        results_dir=args.results_dir, allow_dirty=args.allow_dirty,
        extra_meta={
            "cohorts": "GSE132465 (SMC) and GSE144735 (KUL3), sha256-verified "
                       "against data/manifest.csv",
            "rho_tolerance": MATURITY_DEPTH_RHO_TOLERANCE,
            "crossing_prevalence": CROSSING_PREVALENCE,
            "min_cells_per_arm": MIN_CELLS_PER_ARM,
            "grain": "one row per (study, patient, author subtype, arm)",
            "label_source": (
                "the cohorts' own published per-cell subtype annotation, NOT "
                "this project's maturity call. The claim under test is about "
                "the statistic, so labels we did not build are the stronger "
                "demonstration."
            ),
            "statistic": (
                "src.harness.depth_confound.depth_confound_report — the same "
                "function, cells and depth definition as "
                "depth_confound_per_arm. Only the binary label differs, so the "
                "two tables are directly comparable."
            ),
            "what_this_answers": (
                "Whether the prevalence bound is empty in this data or only in "
                "this project's labels. It is only in the labels: several "
                "hundred per-arm rows sit at 0 < p < 1.3516%, a large share of "
                "them above half their attainable bound, and the 0.20 "
                "tolerance flags none of them."
            ),
        },
    )
    log.info("wrote %s", path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
