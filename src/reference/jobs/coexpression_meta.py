"""Piece 3: combine per-study coexpression readings. Never pool.

    python -m src.reference.jobs.coexpression_meta                  # newest ICBI run
    python -m src.reference.jobs.coexpression_meta --deltas <path>  # any committed table

A THIN LOCAL READ over committed per-study tables. The cluster job ends at
per-study deltas; this combines them, and it runs on a laptop from artifacts in
git. Same shape as ``run_replication``'s side-by-side.

WHAT IS POOLED, AND WHAT IS NOT. Invariant 4 forbids pooling cells across
studies; it requires pooling ESTIMATES. So each study's control shift is
estimated from its own patients, and only those per-study estimates meet.

EACH CONTROL SEPARATELY, AND BOTH MUST HOLD. ACTB and KRT8 are meta-analysed on
their own and the premise holds only if both pooled intervals sit inside the
tolerance -- mirroring the per-study rule that every control must hold. Pooling
"the worst control per study" would hide which control drives the
heterogeneity, and the heterogeneity is half the answer.

THE ORDER IS FIXED. The premise is a gate: the per-gene detection deltas are
read at the meta level ONLY if the pooled controls resolve inside tolerance. A
marker falling inside a population that has itself changed is not silencing, and
that is as true of fourteen studies as of one.

THE BASELINE THIS MUST REPRODUCE. Run against the committed three-cohort table
it should return "undecided at k = 3" -- because each of GSE178341, SMC and KUL3
came back UNRESOLVED individually. A combiner that cannot reproduce the known
answer has not been checked, so ``--deltas`` on that table is the dry run, and
it is what pre-commits this bar before the fourteen numbers exist.
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
from src.common.paths import RESULTS_DIR
from src.common.provenance import DEFAULT_SEED
from src.harness.meta import (
    MAX_I_SQUARED,
    MIN_STUDIES,
    MetaError,
    meta_analyse,
    premise_verdict,
)
from src.reference.jobs.coexpression_silencing import (
    CONTROL_LOG2_TOLERANCE,
    GENE_ROLES,
    MIN_PREMISE_PATIENTS,
)

log = logging.getLogger(__name__)

CONTROLS = tuple(g for g, role in GENE_ROLES.items() if role == "control")
TARGETS = tuple(g for g, role in GENE_ROLES.items() if role in ("target", "identity"))


def per_study_stats(
    deltas: pd.DataFrame, gene: str, column: str = "log2_cp10k_ratio"
) -> pd.DataFrame:
    """One row per study: the estimate and its standard error.

    The SE is the standard error of the mean over PATIENTS, which is the unit
    of inference (invariant 5). It is what inverse-variance weighting wants,
    and it agrees with the bootstrap interval `premise_holds` reports to within
    Monte Carlo error -- pinned by a test rather than asserted, since a meta
    layer weighting by a different quantity than the per-study check reports
    would be two answers to one question.
    """
    rows = deltas.loc[deltas["gene"] == gene, ["study_id", "patient_id", column]]
    rows = rows.dropna(subset=[column])
    if rows.empty:
        return pd.DataFrame(
            columns=["study_id", "gene", "n_patients", "estimate", "se",
                     "ci_low", "ci_high"]
        )

    out = []
    for study, block in rows.groupby("study_id", observed=True):
        # One value per patient. A study contributing many cells for one
        # patient is one patient's worth of evidence.
        per_patient = block.groupby("patient_id")[column].mean()
        n = int(per_patient.size)
        mean = float(per_patient.mean())
        se = float(per_patient.std(ddof=1) / np.sqrt(n)) if n > 1 else float("nan")
        out.append({
            "study_id": str(study), "gene": gene, "n_patients": n,
            "estimate": mean, "se": se,
            "ci_low": mean - 1.96 * se, "ci_high": mean + 1.96 * se,
            "below_patient_floor": n < MIN_PREMISE_PATIENTS,
        })
    return pd.DataFrame(out).sort_values("study_id", ignore_index=True)


def meta_control(deltas: pd.DataFrame, gene: str) -> tuple[dict, pd.DataFrame]:
    """Meta-analyse one control's log2 shift. Returns (row, per-study table)."""
    per_study = per_study_stats(deltas, gene)
    usable = per_study[~per_study["below_patient_floor"]] if not per_study.empty else per_study
    if len(usable) < MIN_STUDIES:
        return {
            "gene": gene, "k_studies": len(usable), "verdict": "UNRESOLVED",
            "detail": (
                f"{len(usable)} study/studies clear the {MIN_PREMISE_PATIENTS}-"
                f"patient floor, below the {MIN_STUDIES} a meta-analysis needs. "
                f"Report them side by side."
            ),
        }, per_study

    try:
        result = meta_analyse(usable["estimate"], usable["se"])
    except MetaError as exc:
        return {"gene": gene, "k_studies": len(usable), "verdict": "UNRESOLVED",
                "detail": str(exc)}, per_study

    verdict, detail = premise_verdict(result, CONTROL_LOG2_TOLERANCE)
    return {"gene": gene, "verdict": verdict, "detail": detail,
            "tolerance": CONTROL_LOG2_TOLERANCE, **result.as_row()}, per_study


def meta_premise(deltas: pd.DataFrame) -> tuple[str, str, pd.DataFrame, pd.DataFrame]:
    """The premise across studies. Every control must hold.

    Returns ``(verdict, detail, control_table, per_study_table)``.
    """
    rows, per_study_frames = [], []
    for gene in CONTROLS:
        row, per_study = meta_control(deltas, gene)
        rows.append(row)
        if not per_study.empty:
            per_study_frames.append(per_study)
    controls = pd.DataFrame(rows)
    per_study = (
        pd.concat(per_study_frames, ignore_index=True)
        if per_study_frames else pd.DataFrame()
    )

    verdicts = set(controls["verdict"])
    if "REFUSED" in verdicts:
        breached = controls[controls["verdict"] == "REFUSED"]["gene"].tolist()
        return "REFUSED", (
            f"control(s) {breached} moved beyond the {CONTROL_LOG2_TOLERANCE} "
            f"tolerance across studies. The arms are not comparable and the "
            f"detection reading is not licensed."
        ), controls, per_study
    if verdicts == {"HOLDS"}:
        return "HOLDS", (
            f"every control ({', '.join(CONTROLS)}) pools inside the "
            f"{CONTROL_LOG2_TOLERANCE} tolerance with tolerable heterogeneity."
        ), controls, per_study

    undecided = controls[controls["verdict"] == "UNRESOLVED"]["gene"].tolist()
    return "UNRESOLVED", (
        f"control(s) {undecided} do not resolve. Not refused and not "
        f"satisfied -- undecided, and reported as such."
    ), controls, per_study


def meta_detection(deltas: pd.DataFrame) -> pd.DataFrame:
    """Pooled detection deltas per target gene. Read ONLY if the premise holds."""
    rows = []
    for gene in TARGETS:
        per_study = per_study_stats(deltas, gene, column="delta_detect")
        if len(per_study) < MIN_STUDIES:
            rows.append({"gene": gene, "k_studies": len(per_study),
                         "note": f"below the {MIN_STUDIES}-study minimum"})
            continue
        try:
            result = meta_analyse(per_study["estimate"], per_study["se"])
        except MetaError as exc:
            rows.append({"gene": gene, "k_studies": len(per_study), "note": str(exc)})
            continue
        rows.append({"gene": gene, **result.as_row()})
    return pd.DataFrame(rows)


def newest(name: str) -> Path | None:
    matches = sorted(RESULTS_DIR.glob(f"*/{name}.parquet"))
    return matches[-1] if matches else None


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--deltas", type=Path, default=None,
        help="per-study deltas. Defaults to the newest icbi_coexpression; "
             "point it at results/2026-09-04_975cf5c/coexpression_silencing."
             "parquet for the three-cohort dry run.",
    )
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--name", default="coexpression_meta")
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    path = args.deltas or newest("icbi_coexpression")
    if path is None:
        raise SystemExit(
            "no results/*/icbi_coexpression.parquet. Run the per-study job "
            "first, or pass --deltas."
        )
    deltas = pd.read_parquet(path)
    studies = sorted(deltas["study_id"].unique())
    log.info("%s | %d studies, %d patient-gene rows",
             path.name, len(studies), len(deltas))
    for study in studies:
        block = deltas[deltas["study_id"] == study]
        log.info("  %-42s %3d patients", study, block["patient_id"].nunique())

    verdict, detail, controls, per_study = meta_premise(deltas)

    log.info("\n%s", "=" * 72)
    log.info("PER-STUDY CONTROL SHIFTS (log2, the premise's own statistic)")
    log.info("%s", "=" * 72)
    if not per_study.empty:
        show = per_study[["gene", "study_id", "n_patients", "estimate", "se"]]
        log.info("%s", show.to_string(index=False))

    log.info("\n%s", "=" * 72)
    log.info("POOLED, PER CONTROL. Each on its own; both must hold.")
    log.info("%s", "=" * 72)
    for _, row in controls.iterrows():
        log.info("  %-6s %s", row["gene"], row["verdict"])
        log.info("         %s", row["detail"])
        if "i_squared" in row and pd.notna(row.get("i_squared")):
            log.info("         I2 = %.1f%%  tau2 = %.4f  prediction "
                     "[%+.3f, %+.3f]  (I2 ceiling %.0f%%)",
                     100 * row["i_squared"], row["tau_squared"],
                     row["prediction_low"], row["prediction_high"],
                     100 * MAX_I_SQUARED)

    log.info("\n%s", "=" * 72)
    log.info("THE META PREMISE: %s", verdict)
    log.info("  %s", detail)
    log.info("%s", "=" * 72)

    detection = pd.DataFrame()
    if verdict == "HOLDS":
        detection = meta_detection(deltas)
        log.info("\nPooled detection deltas (licensed by the premise above):")
        log.info("%s", detection.to_string(index=False))
    else:
        log.info(
            "\nThe per-gene detection deltas are NOT read at the meta level.\n"
            "The premise gates them: a marker falling inside a population that\n"
            "has itself changed is not silencing, and that is as true of %d\n"
            "studies as of one. The per-study tables stand on their own.",
            len(studies),
        )

    frames = [(controls, f"{args.name}_controls"), (per_study, f"{args.name}_per_study")]
    if not detection.empty:
        frames.append((detection, f"{args.name}_detection"))
    for frame, name in frames:
        if frame.empty:
            continue
        written = write_versioned_table(
            frame, name, seed=args.seed,
            results_dir=args.results_dir, allow_dirty=args.allow_dirty,
            extra_meta={
                "source": str(path.name),
                "studies": studies,
                "meta_premise_verdict": verdict,
                "meta_premise_detail": detail,
                "control_genes": list(CONTROLS),
                "control_rule": (
                    "each control meta-analysed separately; the premise holds "
                    "only if EVERY pooled interval sits inside the tolerance"
                ),
                "tolerance_log2": CONTROL_LOG2_TOLERANCE,
                "i_squared_ceiling": MAX_I_SQUARED,
                "min_studies": MIN_STUDIES,
                "pooling": (
                    "estimates, never cells (invariant 4). Each study's shift "
                    "is estimated from its own patients."
                ),
                "detection_gated_by_premise": True,
                "exploratory": True,
                "pre_registered": False,
            },
        )
        log.info("wrote %s", written)
    return 0


if __name__ == "__main__":
    sys.exit(main())
