"""The interval calibration under skewed, heavy-tailed and boundary generators.

    python -m src.reference.jobs.interval_stress_calibration
    python -m src.reference.jobs.interval_stress_calibration --trials 400 \
        --seeds 2 --cohorts mlh1_methylated,adenoma_lineage   # a timed probe

A THIN LOCAL RUN. No cluster and no atlas: the cohorts are read off the same
committed patient summaries ``interval_calibration`` uses, and everything else is
simulation. It writes two tables:

``interval_stress_calibration``
    One row per (cohort, regime, method, seed): the false-positive rate of the
    interval under a true null, its binomial standard error and Wilson interval,
    and what each of the three check families detected in that cell.

``interval_stress_calibration_summary``
    The same grid collapsed over seeds, with seed-to-seed spread and each
    family's detection rate.

WHAT IT ANSWERS, AND WHAT IT CANNOT. ``interval_calibration`` measured the
percentile bootstrap's over-rejection under a **normal** per-patient effect.
The audit decision record (§3) says that result does not establish that BCa is
generally worse or that Student-t is universally calibrated; this job varies the
generator to see how far it holds. It is a statement about **these tested
generators**, not a universal claim about bootstraps.

READ THE POOLED COHORTS SEPARATELY (invariant 4). No row combines two cohorts;
each cohort carries its own n and its own cells-per-arm vector, and the summary
keeps them apart.
"""

from __future__ import annotations

import argparse
import logging
import sys
from collections.abc import Sequence
from pathlib import Path

import pandas as pd

from src.common.io import write_versioned_table
from src.common.provenance import DEFAULT_SEED
from src.reference.interval_stress import (
    CHECK_FAMILIES,
    REGIMES,
    STRESS_TAU,
    stress_calibration_table,
    summarise_stress,
)
from src.reference.jobs.interval_calibration import (
    cohort_vectors,
    newest,
    pelka_with_strata,
)
from src.reference.mlh1_arms import ARMS, arm_of

log = logging.getLogger(__name__)

#: Fixed seeds, so a cell's seed-to-seed spread is reproducible and the table
#: can be re-derived. Five is enough to see a verdict that moves with the seed;
#: the binomial standard error is the other half of the uncertainty.
SEEDS: tuple[int, ...] = (20260101, 20260102, 20260103, 20260104, 20260105)

#: Replicates per cell. 800 keeps the binomial SE near 1pp at a 10% rate, which
#: is inside the 2pp tolerance, so a cell is not called miscalibrated by noise.
DEFAULT_TRIALS: int = 800


def build_cohorts() -> dict:
    """The committed-summary cohorts, defined in one place.

    Imports the loader from ``interval_calibration`` rather than re-deriving it:
    a second definition of an arm has already cost this repo a defect, and the
    stress grid must be calibrated on exactly the cohorts the paper reports at.
    """
    pelka = pelka_with_strata()
    pelka["arm"] = pelka["mlh1_stratum"].map(arm_of)
    cohorts = {
        arm: cohort_vectors(pelka[pelka["arm"] == arm])
        for arm in ARMS if not pelka[pelka["arm"] == arm].empty
    }
    adenoma_path = newest("icbi_adenoma")
    if adenoma_path is not None:
        adenoma = pd.read_parquet(adenoma_path)
        adenoma["short_id"] = adenoma["patient_id"].astype(str)
        for rung, block in adenoma.groupby("granularity_rung", observed=True):
            cohorts[f"adenoma_{rung}"] = cohort_vectors(block)
    if not cohorts:
        raise SystemExit("no cohorts could be built from the committed summaries")
    return cohorts


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trials", type=int, default=DEFAULT_TRIALS)
    parser.add_argument(
        "--seeds", type=int, default=len(SEEDS),
        help="how many of the fixed seed list to use",
    )
    parser.add_argument(
        "--cohorts", default=None,
        help="comma-separated subset of cohort names; default all",
    )
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED,
                        help="retained for provenance; the cell seeds are fixed")
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    cohorts = build_cohorts()
    if args.cohorts:
        wanted = [c.strip() for c in args.cohorts.split(",") if c.strip()]
        missing = sorted(set(wanted) - set(cohorts))
        if missing:
            raise SystemExit(
                f"unknown cohort(s) {missing}; available: {sorted(cohorts)}"
            )
        cohorts = {k: cohorts[k] for k in wanted}
    seeds = SEEDS[: args.seeds]

    log.info("cohorts (%d), kept separate:", len(cohorts))
    for name, (n_cells, depth) in cohorts.items():
        log.info("  %-22s n=%2d patients, cells/arm median %5.0f, depth median %6.0f",
                 name, len(n_cells), float(n_cells.mean()), float(depth.mean()))
    log.info("regimes: %s", ", ".join(r.name for r in REGIMES))
    log.info("seeds: %s | trials/cell: %d | tau: %.2f", seeds, args.trials, STRESS_TAU)

    frame = stress_calibration_table(
        cohorts=cohorts, regimes=REGIMES, seeds=seeds,
        n_trials=args.trials, tau=STRESS_TAU,
    )
    summary = summarise_stress(frame)

    log.info("\n%s\nPERCENTILE FALSE-POSITIVE RATE BY REGIME (tau=%.2f)\n%s",
             "=" * 72, STRESS_TAU, "=" * 72)
    pct = summary[summary["method"] == "percentile"]
    log.info("%s", (100 * pct.pivot_table(
        index=["cohort", "n_patients"], columns="regime", values="fpr_median"
    )).round(1).to_string())
    log.info("\n%s\nSTUDENT-T FALSE-POSITIVE RATE BY REGIME (tau=%.2f)\n%s",
             "=" * 72, STRESS_TAU, "=" * 72)
    stu = summary[summary["method"] == "student_t"]
    log.info("%s", (100 * stu.pivot_table(
        index=["cohort", "n_patients"], columns="regime", values="fpr_median"
    )).round(1).to_string())
    log.info("\n%s\nWHAT EACH CHECK FAMILY DETECTED (share of seeds that fired)\n%s",
             "=" * 72, "=" * 72)
    for family in CHECK_FAMILIES:
        fired = summary[f"detection_rate_{family}"]
        log.info("  %-17s fired on %d of %d cells (max %.2f)",
                 family, int((fired > 0).sum()), len(fired),
                 float(fired.max()) if len(fired) else float("nan"))

    meta = {
        "purpose": (
            "DECISION_2026-09-11_scope_and_pivot.md §4 week 2 -- add "
            "skewed/heavy-tailed and boundary regimes to the interval "
            "calibration, quantify Monte Carlo uncertainty. Delimits the "
            "generators over which the §3a result was measured."
        ),
        "generator": (
            "the same Poisson thinning as interval_calibration "
            "(p = 1 - exp(-mu), mu = cp10k/1e4 * depth, scored through "
            "cloglog_rate and its boundary rule); only the distribution of the "
            "per-patient log fold change varies, standardised to unit variance "
            "and scaled by tau, with mean zero so every regime is a true null"
        ),
        "regimes": [
            {"name": r.name, "kind": r.kind, "cp10k": r.cp10k,
             "description": r.description}
            for r in REGIMES
        ],
        "tau": STRESS_TAU,
        "tau_source": (
            "Pelka's control genes measure tau ~ 0.2; see interval_heterogeneity "
            "from the interval_calibration run"
        ),
        "check_families": {
            "input_validation": (
                "finite values, at least two patients, non-zero spread, finite "
                "interval; ordinary hygiene, excluded from the stop rule by name"
            ),
            "original_check": (
                "the pre-audit diagnostic: measured percentile rate above the "
                "closed form by more than the tolerance. Defined for percentile "
                "only; null elsewhere"
            ),
            "audit_check": (
                "measured rate above nominal + tolerance AND above nominal on the "
                "Wilson interval"
            ),
        },
        "monte_carlo_uncertainty": (
            "binomial standard error and Wilson interval per cell; five fixed "
            "seeds per cell so seed-to-seed spread is visible"
        ),
        "what_this_does_not_establish": (
            "a universal claim about bootstrap intervals. These are the tested "
            "generators; the result is about them, per DECISION_2026-09-11 §3."
        ),
        "ownership_note": (
            "cohorts are read through interval_calibration's loader, so the "
            "stress grid runs at exactly the n and cells-per-arm the paper "
            "reports at; a second definition of an arm is what caused the "
            "15-vs-19 defect."
        ),
        "exploratory": False,
    }
    for out, name in (
        (frame, "interval_stress_calibration"),
        (summary, "interval_stress_calibration_summary"),
    ):
        path = write_versioned_table(
            out, name, seed=args.seed, results_dir=args.results_dir,
            extra_meta=meta, allow_dirty=args.allow_dirty,
        )
        log.info("wrote %s (%d rows)", path, len(out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
