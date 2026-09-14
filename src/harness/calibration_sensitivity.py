"""Separate inner-bootstrap noise from patient-omission sensitivity.

Outer samples enumerate every two-patient holdout pair (45 for SMC, 15 for
KUL3). The same pair-keyed samples are reused at each inner-bootstrap budget.
Patient influence is assessed by removing every baseline pair containing that
patient; unaffected pairs are neither regenerated nor re-randomised.
"""

from __future__ import annotations

import argparse
import itertools
import logging
import sys
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd

from src.common.io import write_versioned_table
from src.common.provenance import DEFAULT_SEED
from src.harness.attenuation import SweepConfig, SweepGrid, run_sweep
from src.harness.calibration import PREREGISTERED, coverage_and_discrimination_by_count
from src.harness.calibration_gap import N_CELLS, load_cohort_arrays
from src.harness.controlled_grid import crossing_brackets

log = logging.getLogger(__name__)

INNER_BUDGETS: tuple[int, ...] = (200, 1000, 5000)
COUNT_WINDOWS: dict[str, tuple[int, ...]] = {
    "smc": (40, 50, 60, 70, 80, 90, 100, 110),
    "kul3": (200, 300, 400, 600, 800),
}
EXPECTED_PAIRS: dict[str, int] = {"smc": 45, "kul3": 15}


def holdout_pairs(patient_ids: Sequence[str]) -> tuple[tuple[str, str], ...]:
    patients = sorted(set(patient_ids))
    return tuple(itertools.combinations(patients, 2))


def _reference_config(arrays: dict) -> SweepConfig:
    tissue = np.asarray(arrays["tissue"])
    patient = np.asarray(arrays["patient_id"])
    mask = tissue == "normal"
    return SweepConfig(
        counts=np.asarray(arrays["counts"])[mask],
        cell_type=np.asarray(arrays["cell_type"])[mask].tolist(),
        patient_id=patient[mask].tolist(),
        genes=list(arrays["genes"]),
        target_gene="GUCA2A",
        mature_label="differentiated",
    )


def _one_sensitivity_run(
    config: SweepConfig,
    *,
    counts: Sequence[int],
    n_boot: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    pairs = holdout_pairs(config.patient_id)
    grid = SweepGrid(
        mature_fractions=tuple(count / N_CELLS for count in counts),
        shifts=(1.0, PREREGISTERED.detectable_shift),
        n_replicates=len(pairs),
        n_cells=N_CELLS,
        n_held_out=2,
        n_boot=n_boot,
    )
    sweep = run_sweep(
        config,
        grid,
        seed=seed,
        arms=("oracle",),
        seed_strategy="configuration",
        holdout_schedule=pairs,
    )
    rates = coverage_and_discrimination_by_count(sweep)
    return sweep, rates


def retained_pair_rows(sweep: pd.DataFrame, omitted_patient: str) -> pd.DataFrame:
    """Drop pairs containing one patient without changing any retained draw."""
    if "held_out_pair" not in sweep:
        raise ValueError("sweep must identify each scheduled pair in held_out_pair")
    contains = sweep["held_out_pair"].str.split("|").map(
        lambda pair: omitted_patient in pair
    )
    return sweep.loc[~contains].copy()


def run_inner_budget_comparison(
    arrays_by_cohort: dict[str, dict],
    *,
    budgets: Sequence[int] = INNER_BUDGETS,
    seed: int = DEFAULT_SEED,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rate_frames: list[pd.DataFrame] = []
    reference_outer: dict[str, pd.DataFrame] = {}
    for cohort, arrays in arrays_by_cohort.items():
        config = _reference_config(arrays)
        pairs = holdout_pairs(config.patient_id)
        if len(pairs) != EXPECTED_PAIRS[cohort]:
            raise ValueError(
                f"{cohort} has {len(pairs)} holdout pairs; expected {EXPECTED_PAIRS[cohort]}"
            )
        for budget in budgets:
            log.info("cohort=%s inner_bootstrap=%d pairs=%d", cohort, budget, len(pairs))
            sweep, rates = _one_sensitivity_run(
                config,
                counts=COUNT_WINDOWS[cohort],
                n_boot=int(budget),
                seed=seed,
            )
            outer = sweep[
                [
                    "frac_mature_tumour",
                    "shift",
                    "replicate",
                    "held_out_pair",
                    "intrinsic_true_parametric",
                    "intrinsic_true_realised",
                    "intrinsic_hat",
                    "n_cells_mature",
                    "seed",
                ]
            ].copy()
            # The first budget is the frozen outer-sample fingerprint. Interval
            # columns are intentionally absent: only inner resampling may differ.
            if cohort in reference_outer:
                pd.testing.assert_frame_equal(
                    reference_outer[cohort].reset_index(drop=True),
                    outer.reset_index(drop=True),
                    check_exact=True,
                )
            else:
                reference_outer[cohort] = outer
            rate_frames.append(
                rates.assign(
                    cohort=cohort,
                    pool="reference",
                    seed=seed,
                    grid="count_window",
                    binning="per_count",
                    inner_bootstrap=budget,
                    n_holdout_pairs=len(pairs),
                )
            )
    rates = pd.concat(rate_frames, ignore_index=True)
    crossings = crossing_brackets(rates)
    return rates, crossings


def run_leave_one_patient_out(
    arrays_by_cohort: dict[str, dict],
    *,
    n_boot: int = 200,
    seed: int = DEFAULT_SEED,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rate_frames: list[pd.DataFrame] = []
    influence_rows: list[dict] = []
    for cohort, arrays in arrays_by_cohort.items():
        config = _reference_config(arrays)
        all_patients = sorted(set(config.patient_id))
        baseline_sweep, baseline_rates = _one_sensitivity_run(
            config,
            counts=COUNT_WINDOWS[cohort],
            n_boot=n_boot,
            seed=seed,
        )
        baseline_labelled = baseline_rates.assign(
            cohort=cohort,
            pool="reference",
            seed=seed,
            grid="count_window",
            binning="per_count",
            inner_bootstrap=n_boot,
            n_holdout_pairs=len(holdout_pairs(config.patient_id)),
        )
        base = crossing_brackets(baseline_labelled).set_index("criterion")
        for omitted in all_patients:
            retained = retained_pair_rows(baseline_sweep, omitted)
            rates = coverage_and_discrimination_by_count(retained)
            n_retained_pairs = retained["held_out_pair"].nunique()
            labelled = rates.assign(
                cohort=cohort,
                omitted_patient=omitted,
                pool="reference",
                seed=seed,
                grid="count_window",
                binning="per_count",
                inner_bootstrap=n_boot,
                n_holdout_pairs=n_retained_pairs,
            )
            rate_frames.append(labelled)
            crossed = crossing_brackets(labelled).set_index("criterion")
            for criterion in ("coverage_and_discrimination", "coverage_only"):
                baseline_candidate = float(base.loc[criterion, "candidate"])
                omitted_candidate = float(crossed.loc[criterion, "candidate"])
                changed = not (
                    (np.isnan(baseline_candidate) and np.isnan(omitted_candidate))
                    or baseline_candidate == omitted_candidate
                )
                influence_rows.append(
                    {
                        "cohort": cohort,
                        "omitted_patient": omitted,
                        "criterion": criterion,
                        "baseline_candidate": baseline_candidate,
                        "omitted_candidate": omitted_candidate,
                        "omitted_status": crossed.loc[criterion, "status"],
                        "conclusion_changed": changed,
                    }
                )
    influence = pd.DataFrame(influence_rows)
    return pd.concat(rate_frames, ignore_index=True), influence


def influence_summary(influence: pd.DataFrame) -> pd.DataFrame:
    return (
        influence.groupby(["cohort", "criterion"], sort=True)
        .agg(
            n_patients=("omitted_patient", "count"),
            n_changing_conclusion=("conclusion_changed", "sum"),
            changing_patients=(
                "omitted_patient",
                lambda s: ",".join(
                    s[influence.loc[s.index, "conclusion_changed"]].astype(str)
                ),
            ),
        )
        .reset_index()
    )


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort", choices=("smc", "kul3", "both"), default="both")
    parser.add_argument("--budgets", type=int, nargs="+", default=list(INNER_BUDGETS))
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--skip-influence", action="store_true")
    parser.add_argument("--raw-dir", type=Path, default=None)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--suffix", default="")
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    cohorts = ("smc", "kul3") if args.cohort == "both" else (args.cohort,)
    arrays = {cohort: load_cohort_arrays(cohort, args.raw_dir) for cohort in cohorts}
    rates, crossings = run_inner_budget_comparison(
        arrays,
        budgets=args.budgets,
        seed=args.seed,
    )
    frames: list[tuple[pd.DataFrame, str]] = [
        (rates, "calibration_inner_budget_rates"),
        (crossings, "calibration_inner_budget_crossings"),
    ]
    if not args.skip_influence:
        if 200 not in args.budgets:
            raise ValueError("leave-one-patient-out comparison requires budget 200")
        omission_rates, influence = run_leave_one_patient_out(
            arrays,
            n_boot=200,
            seed=args.seed,
        )
        frames.extend(
            [
                (omission_rates, "calibration_lopo_rates"),
                (influence, "calibration_lopo_influence"),
                (influence_summary(influence), "calibration_lopo_summary"),
            ]
        )
    extra = {
        "cohorts": list(cohorts),
        "count_windows": {key: list(COUNT_WINDOWS[key]) for key in cohorts},
        "inner_bootstrap_budgets": list(args.budgets),
        "outer_samples": "Every two-patient holdout pair, enumerated once per count.",
        "uncertainty_separation": (
            "Inner-budget rows repeat fixed outer samples. Patient-influence rows "
            "drop every baseline holdout pair containing that patient without "
            "regenerating retained pairs. This is not a population confidence interval."
        ),
    }
    suffix = args.suffix or "_b" + "-".join(map(str, args.budgets))
    for frame, name in frames:
        output = write_versioned_table(
            frame,
            name + suffix,
            seed=args.seed,
            results_dir=args.results_dir,
            allow_dirty=args.allow_dirty,
            extra_meta=extra,
        )
        log.info("wrote %s", output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
