"""Is the 13-study heterogeneity verdict a property of the data or of the floor?

    python -m src.reference.jobs.meta_floor_sensitivity
    python -m src.reference.jobs.meta_floor_sensitivity --per-study <path>

A THIN LOCAL READ over a committed table. No cells, no cluster, no new data --
it re-reads ``results/2026-09-05_61ba221/icbi_coexpression_meta_per_study.parquet``
at three patient floors and reports what each one does to the verdict.

PRE-REGISTERED: ``docs/prereg_meta_weight_calibration.md``, committed in
``a9dc874`` before this job existed and before the floor curve was run. §6 there
states in advance that the rule will probably move KRT8 off UNRESOLVED and that
this must be read as the verdict being floor-unstable, NOT as the controls
holding.

WHAT THE FLOORS ARE, AND WHY THEY ARE NOT A SWEEP. Inverse-variance weights are
inverse chi-squares on ``n-1`` df. ``E[w]`` is finite only for ``n >= 4`` and
``Var[w]`` only for ``n >= 6``. Those two numbers, plus the status quo of 3, are
the whole grid. It is not a search over thresholds -- each floor is the smallest
``n`` at which a named moment of the weight exists, fixed by the chi-square df
condition in the pre-registration. A fourth floor chosen because it gave a
nicer answer would be the defect this job is about.

EVERY ROW CARRIES ITS OWN NULL. ``I^2`` against a fixed ceiling is meaningless
without the distribution of ``I^2`` under homogeneity at that floor's patient
counts, which at the ICBI ``n``s has median 0.269 rather than 0. So the null and
the calibrated ``p`` are emitted in the same row as the verdict, and
``check_heterogeneity_carries_its_own_null`` refuses the frame otherwise.

WHAT THIS JOB REFUSES TO DO. It does not drop a study by name. The influence
table it emits alongside is diagnostic and labelled ``exploratory`` in its own
column; it exists so that the reader can see WHERE the heterogeneity sits, not
so that anyone can remove it.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from src.common.io import write_versioned_table
from src.common.paths import RESULTS_DIR
from src.harness.meta import MAX_I_SQUARED, MetaError, meta_analyse, premise_verdict
from src.reference.jobs.coexpression_silencing import CONTROL_LOG2_TOLERANCE
from src.reference.meta_calibration import (
    DEFAULT_N_TRIALS,
    DEFAULT_SEED,
    FLOORS_REPORTED,
    MetaCalibrationError,
    calibrated_p,
    check_heterogeneity_carries_its_own_null,
    floor_label,
    null_i_squared,
    weight_inflation,
)

log = logging.getLogger(__name__)

#: The committed meta this re-reads. Named so the provenance is in the source.
DEFAULT_PER_STUDY = "2026-09-05_61ba221/icbi_coexpression_meta_per_study.parquet"

#: Significance level for the calibrated heterogeneity test. The level the fixed
#: 0.75 ceiling was implicitly standing in for and never achieved.
NULL_ALPHA = 0.05


def find_per_study(results_dir: Path) -> Path | None:
    direct = results_dir / DEFAULT_PER_STUDY
    if direct.exists():
        return direct
    matches = sorted(results_dir.glob("*/icbi_coexpression_meta_per_study.parquet"))
    return matches[-1] if matches else None


def influence(per_study: pd.DataFrame) -> pd.DataFrame:
    """Where the heterogeneity sits, per study. Diagnostic, and labelled so.

    Reports each study's share of Cochran's Q and of the fixed-effect weight
    next to the closed-form inflation ``(n-1)/(n-3)`` its ``n`` implies. A study
    whose weight share far exceeds its patient share is being over-weighted, and
    the closed form says by how much in expectation.
    """
    rows = []
    for gene, group in per_study.groupby("gene", sort=True):
        g = group.dropna(subset=["estimate", "se"])
        g = g[g["se"] > 0]
        if len(g) < 3:
            continue
        y, s = g["estimate"].to_numpy(float), g["se"].to_numpy(float)
        w = 1.0 / s ** 2
        fixed = (w * y).sum() / w.sum()
        q_i = w * (y - fixed) ** 2
        total_q = q_i.sum()
        n = g["n_patients"].to_numpy(float)
        for i, (_, study) in enumerate(g.iterrows()):
            rows.append({
                "gene": gene,
                "study_id": study["study_id"],
                "n_patients": int(n[i]),
                "estimate": float(y[i]),
                "se": float(s[i]),
                "sd_over_patients": float(s[i] * np.sqrt(n[i])),
                "weight_share": float(w[i] / w.sum()),
                "patient_share": float(n[i] / n.sum()),
                "q_contribution": float(q_i[i]),
                "q_share": float(q_i[i] / total_q) if total_q > 0 else float("nan"),
                "expected_weight_inflation": weight_inflation(int(n[i])),
                "exploratory": True,
            })
    return pd.DataFrame(rows)


def floor_curve(
    per_study: pd.DataFrame,
    tolerance: float = CONTROL_LOG2_TOLERANCE,
    floors: tuple[int, ...] = FLOORS_REPORTED,
    *,
    n_trials: int = DEFAULT_N_TRIALS,
    seed: int = DEFAULT_SEED,
) -> pd.DataFrame:
    """One row per (gene, floor), verdict and null together or not at all."""
    rows = []
    for gene, group in per_study.groupby("gene", sort=True):
        usable = group.dropna(subset=["estimate", "se"])
        usable = usable[usable["se"] > 0]
        for floor in floors:
            kept = usable[usable["n_patients"] >= floor]
            dropped = sorted(
                set(usable["study_id"]) - set(kept["study_id"])
            )
            row: dict[str, object] = {
                "gene": gene,
                "patient_floor": int(floor),
                "floor_rationale": floor_label(floor),
                "k": int(len(kept)),
                "studies_dropped": "; ".join(dropped) if dropped else "",
                "n_dropped": len(dropped),
                "min_n_patients": int(kept["n_patients"].min()) if len(kept) else 0,
                "tolerance": float(tolerance),
            }
            try:
                result = meta_analyse(kept["estimate"], kept["se"])
            except MetaError as exc:
                # Invariant 1: a meta that could not run is not a meta of zero.
                row.update({
                    "verdict": "NOT_ESTIMABLE",
                    "detail": str(exc),
                    "estimability": "not_estimable",
                })
                rows.append(row)
                continue

            counts = kept["n_patients"].to_numpy(int)
            null = null_i_squared(counts, n_trials=n_trials, seed=seed)
            verdict, detail = premise_verdict(result, tolerance)
            row.update({
                "estimability": "estimated",
                "pooled": result.pooled,
                "se": result.se,
                "ci_low": result.ci_low,
                "ci_high": result.ci_high,
                "tau_squared": result.tau_squared,
                "i_squared": result.i_squared,
                "cochran_q": result.q,
                "df": result.df,
                "ceiling": float(MAX_I_SQUARED),
                "homogeneous": bool(result.homogeneous),
                "verdict": verdict,
                "detail": detail,
                "null_p_of_observed": calibrated_p(
                    result.i_squared, counts, n_trials=n_trials, seed=seed
                ),
            })
            row.update(null.as_row())
            # The ceiling's answer and the null's answer, side by side. They are
            # not the same question: 0.75 is a fixed rule of thumb, the null is
            # what I^2 does at THESE patient counts under exact homogeneity.
            row["heterogeneous_by_ceiling"] = not bool(result.homogeneous)
            row["heterogeneous_by_null"] = bool(
                row["null_p_of_observed"] < NULL_ALPHA
            )
            row["ceiling_and_null_agree"] = bool(
                row["heterogeneous_by_ceiling"] == row["heterogeneous_by_null"]
            )
            rows.append(row)
    frame = pd.DataFrame(rows)
    check_heterogeneity_carries_its_own_null(frame)
    frame["verdict_cause"] = frame.apply(verdict_cause, axis=1)
    return frame


def verdict_cause(row: "pd.Series") -> str:
    """WHY a verdict is what it is. ``premise_verdict`` returns UNRESOLVED down
    two different routes and the label alone cannot tell them apart.

    This function exists because the first version of ``read_verdict`` compared
    labels and reported STABLE for KRT8 — which is UNRESOLVED at every floor,
    but at n>=3 because ``I^2`` exceeds the ceiling and at n>=4 and n>=6 because
    the pooled interval straddles the tolerance while the studies are
    HOMOGENEOUS. Those are opposite findings wearing one word: the first says
    the studies do not agree, the second says they agree that the control moved.
    A comparison that cannot distinguish them is a check that cannot fail.
    """
    if row["estimability"] != "estimated":
        return "not_estimable"
    if not bool(row["homogeneous"]):
        return "heterogeneity_over_ceiling"
    if row["verdict"] == "HOLDS":
        return "homogeneous_within_tolerance"
    if row["verdict"] == "REFUSED":
        return "homogeneous_beyond_tolerance"
    return "homogeneous_straddles_tolerance"


def read_verdict(frame: pd.DataFrame) -> dict[str, str]:
    """What the curve says, in the four branches §6 of the prereg fixed.

    Compares the CAUSE, never the label — see :func:`verdict_cause`.
    """
    primary = frame[frame["patient_floor"] == max(FLOORS_REPORTED)]
    weak = frame[frame["patient_floor"] == 4]
    status_quo = frame[frame["patient_floor"] == 3]

    def verdicts(f: pd.DataFrame) -> dict[str, str]:
        return dict(zip(f["gene"], f["verdict_cause"]))

    v_primary, v_weak, v_status = verdicts(primary), verdicts(weak), verdicts(status_quo)
    moved = sorted(g for g in v_status if v_status[g] != v_primary.get(g))
    disagree = sorted(g for g in v_weak if v_weak[g] != v_primary.get(g))

    def moves(gene: str) -> str:
        return f"{gene}: {v_status[gene]} -> {v_primary.get(gene)}"

    disputed = sorted(
        frame.loc[~frame["ceiling_and_null_agree"].fillna(True), "gene"].unique()
    )
    dispute_note = (
        f" And the fixed {MAX_I_SQUARED:.0%} ceiling disagrees with its own "
        f"calibrated null for {', '.join(disputed)} at one or more floors — "
        f"the ceiling calls them homogeneous where the null does not."
        if disputed else ""
    )

    if not moved:
        return {
            "verdict": "STABLE",
            "detail": (
                "no control's verdict changes CAUSE between the status-quo "
                "floor and the floor at which the weight has two moments. The "
                "heterogeneity is a property of the data, not of the weighting. "
                "The committed conclusion stands and now carries a calibration "
                "it did not have." + dispute_note
            ),
        }
    if disagree:
        return {
            "verdict": "FLOOR-UNSTABLE, AND THE FLOORS DISAGREE",
            "detail": (
                f"{'; '.join(moves(g) for g in moved)} against the status quo, "
                f"and {', '.join(disagree)} differ between n>=4 and n>=6. Report "
                f"all three floors; claim none. The instrument does not resolve "
                f"it." + dispute_note
            ),
        }
    return {
        "verdict": "FLOOR-UNSTABLE",
        "detail": (
            f"{'; '.join(moves(g) for g in moved)} once studies whose weights "
            f"are not integrable are excluded, and n>=4 and n>=6 agree. The "
            f"verdict was carried by those studies. This says the instrument "
            f"was unstable to a criterion nobody had applied; it does NOT say "
            f"the controls hold." + dispute_note
        ),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--per-study", type=Path, default=None,
                        help="committed icbi_coexpression_meta_per_study table")
    parser.add_argument("--results-dir", type=Path, default=RESULTS_DIR)
    parser.add_argument("--n-trials", type=int, default=DEFAULT_N_TRIALS)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    path = args.per_study or find_per_study(args.results_dir)
    if path is None or not Path(path).exists():
        raise SystemExit(
            f"no per-study meta table found under {args.results_dir}. "
            f"Expected {DEFAULT_PER_STUDY}; pass --per-study."
        )
    per_study = pd.read_parquet(path)
    log.info("read %s: %d rows, %d genes",
             path, len(per_study), per_study["gene"].nunique())

    log.info("\n%s\nWEIGHT INFLATION, CLOSED FORM — E[w_hat]/w_true = (n-1)/(n-3)\n%s",
             "=" * 72, "=" * 72)
    for n in sorted(per_study["n_patients"].unique()):
        inflation = weight_inflation(int(n))
        mark = "   <- no finite mean" if not np.isfinite(inflation) else ""
        log.info("  n = %2d   %s%s", int(n),
                 "inf" if not np.isfinite(inflation) else f"{inflation:.2f}x", mark)

    curve = floor_curve(per_study, n_trials=args.n_trials, seed=args.seed)
    log.info("\n%s\nTHE FLOOR CURVE — verdict and its null, per gene\n%s",
             "=" * 72, "=" * 72)
    show = ["gene", "patient_floor", "k", "i_squared", "null_i_squared_median",
            "null_p_of_observed", "heterogeneous_by_ceiling",
            "heterogeneous_by_null", "verdict", "verdict_cause"]
    log.info("%s", curve[show].to_string(index=False))

    infl = influence(per_study[~per_study["below_patient_floor"]]
                     if "below_patient_floor" in per_study.columns else per_study)
    log.info("\n%s\nWHERE THE HETEROGENEITY SITS — DIAGNOSTIC, POST-HOC\n%s",
             "=" * 72, "=" * 72)
    for gene, g in infl.groupby("gene"):
        top = g.nlargest(3, "q_share")
        log.info("  %s: %s", gene, ", ".join(
            f"{r.study_id.split('_')[0]} {r.q_share:.1%} of Q "
            f"(n={r.n_patients}, {r.weight_share:.1%} of weight)"
            for r in top.itertuples()))

    outcome = read_verdict(curve)
    log.info("\n%s\nVERDICT\n%s", "=" * 72, "=" * 72)
    log.info("  %s", outcome["verdict"])
    log.info("  %s", outcome["detail"])
    log.info(
        "\n  READ §6 OF THE PREREG BEFORE QUOTING THIS. A verdict that moves "
        "when the\n  floor moves is a statement about the INSTRUMENT. It is not "
        "evidence that\n  the controls hold, at any floor.")

    meta = {
        "prereg": "docs/prereg_meta_weight_calibration.md",
        "prereg_committed_in": "a9dc874, before this job existed",
        "source": f"{Path(path).parent.name}/{Path(path).name}",
        "what_this_is": (
            "a sensitivity of the 13-study heterogeneity verdict to the patient "
            "floor, where each floor is the smallest n at which a named moment "
            "of the inverse-variance weight exists"
        ),
        "what_this_is_not": (
            "a threshold search, and not a licence to drop a study. The three "
            "floors are fixed by the chi-square df condition (E[w] finite iff "
            "n>=4, Var[w] finite iff n>=6) and were written down before the "
            "curve was run."
        ),
        "closed_form": "E[w_hat]/w_true = (n-1)/(n-3); no finite mean at n<=3",
        "primary_floor": max(FLOORS_REPORTED),
        "floors_reported": list(FLOORS_REPORTED),
        "ceiling": float(MAX_I_SQUARED),
        "null_is_per_row": True,
        "does_not_touch": (
            "the adenoma results (avenue A, §6h, §6j) — they do not go through "
            "meta.py at all"
        ),
        "verdict": outcome,
        "exploratory": False,
        "pre_registered": True,
    }
    for frame, name in ((curve, "meta_floor_sensitivity"),
                        (infl, "meta_floor_influence")):
        log.info("wrote %s", write_versioned_table(
            frame, name, seed=args.seed, results_dir=args.results_dir,
            allow_dirty=args.allow_dirty, extra_meta=meta,
        ))
    return 0


if __name__ == "__main__":
    sys.exit(main())
