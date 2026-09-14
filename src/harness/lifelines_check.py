"""Independent lifelines check of the bespoke Cox implementations.

This is a deliberately small reproduction command, not a prevalence study.
For four fixed administrative-censoring targets it fits the same generated
cohort with this repository's Breslow implementation and with lifelines'
``CoxPHFitter``. Continuous event times make Breslow and Efron tie handling
coincide here; tied-time behavior has a separate brute-force regression test.
"""

from __future__ import annotations

import argparse
import logging
import sys
import warnings
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter

from src.common.io import write_versioned_table
from src.harness.trial_survival import (
    CENSORING_TARGETS,
    censor,
    censoring_time_for,
    simulate_latent_cohort,
    stratified_cox_breslow,
    unadjusted_cox,
)

log = logging.getLogger(__name__)

FIXED_TOLERANCE = 1e-5


def compare(
    *,
    seed: int = 5,
    n_patients: int = 1500,
    targets: Sequence[float] = CENSORING_TARGETS,
    tolerance: float = FIXED_TOLERANCE,
) -> pd.DataFrame:
    """Return paired estimates and explicit validity/tolerance outcomes."""
    rows: list[dict] = []
    for target in targets:
        cohort = simulate_latent_cohort(n_patients, rng=np.random.default_rng(seed))
        trial = censor(cohort, censoring_time_for(float(target)))
        frame = trial.records.rename(columns={"time": "T", "event": "E"}).astype(
            {"treated": float}
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            independent = {
                "stratified-cox": float(
                    CoxPHFitter()
                    .fit(
                        frame[["T", "E", "treated", "stratum"]],
                        duration_col="T",
                        event_col="E",
                        strata=["stratum"],
                    )
                    .params_["treated"]
                ),
                "unadjusted-cox": float(
                    CoxPHFitter()
                    .fit(
                        frame[["T", "E", "treated"]],
                        duration_col="T",
                        event_col="E",
                    )
                    .params_["treated"]
                ),
            }
        bespoke = {
            "stratified-cox": stratified_cox_breslow(trial.records),
            "unadjusted-cox": unadjusted_cox(trial.records),
        }
        for estimator in bespoke:
            valid = bool(np.isfinite([bespoke[estimator], independent[estimator]]).all())
            difference = (
                abs(bespoke[estimator] - independent[estimator]) if valid else np.nan
            )
            rows.append(
                {
                    "censoring_target": float(target),
                    "censored_fraction": trial.censored_fraction,
                    "estimator": estimator,
                    "bespoke_estimate": bespoke[estimator],
                    "lifelines_estimate": independent[estimator],
                    "absolute_difference": difference,
                    "valid_pair": valid,
                    "within_fixed_tolerance": valid and difference <= tolerance,
                    "fixed_tolerance": tolerance,
                    "n_patients": n_patients,
                    "seed": seed,
                }
            )
    return pd.DataFrame(rows)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=5)
    parser.add_argument("--patients", type=int, default=1500)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    table = compare(seed=args.seed, n_patients=args.patients)
    log.info("\n%s", table.to_string(index=False))
    if args.results_dir is not None:
        output = write_versioned_table(
            table,
            "trial_survival_lifelines_check",
            seed=args.seed,
            results_dir=args.results_dir,
            allow_dirty=args.allow_dirty,
            extra_meta={
                "fixed_tolerance": FIXED_TOLERANCE,
                "independent_implementation": "lifelines.CoxPHFitter",
                "scope": (
                    "One fixed synthetic generator example; not evidence of prevalence "
                    "across clinical pipelines."
                ),
            },
        )
        log.info("wrote %s", output)
    if not table["within_fixed_tolerance"].all():
        log.error("at least one comparison was invalid or exceeded the fixed tolerance")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
