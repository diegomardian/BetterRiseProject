"""Prospectively specified control showing residual class is not performance.

The empirical full-sample mean is defined before any estimator is evaluated.
Four estimators then separate point equality, interval calibration, bias and a
legitimate departure from that reference. Both a null and non-null requested
effect are included.
"""

from __future__ import annotations

import argparse
import logging
import math
import sys
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from src.common.io import write_versioned_table

log = logging.getLogger(__name__)

REFERENCE_TOLERANCE = 1e-12
REQUESTED_EFFECTS: tuple[float, ...] = (0.0, 0.5)
KNOWN_BIAS = 0.5
Z_975 = 1.959963984540054


@dataclass(frozen=True)
class Estimate:
    point: float
    low: float
    high: float


Estimator = Callable[[np.ndarray, float], Estimate]


def _mean_interval(values: np.ndarray, centre: float, *, width_scale: float = 1.0) -> Estimate:
    if len(values) < 2 or not np.isfinite(values).all():
        return Estimate(math.nan, math.nan, math.nan)
    half_width = width_scale * Z_975 * float(values.std(ddof=1)) / math.sqrt(len(values))
    return Estimate(float(centre), float(centre - half_width), float(centre + half_width))


def _reference_mean(values: np.ndarray, reference: float) -> Estimate:
    return _mean_interval(values, reference)


def _reference_mean_narrow(values: np.ndarray, reference: float) -> Estimate:
    return _mean_interval(values, reference, width_scale=0.25)


def _known_biased(values: np.ndarray, reference: float) -> Estimate:
    return _mean_interval(values, reference + KNOWN_BIAS)


def _independent_half_mean(values: np.ndarray, _reference: float) -> Estimate:
    # Uses a prospectively selected disjoint half of the observations. It is a
    # valid estimator of the requested population mean but need not equal the
    # full-sample empirical reference.
    half = values[::2]
    return _mean_interval(half, float(half.mean()))


ESTIMATORS: dict[str, Estimator] = {
    "empirical-mean-calibrated": _reference_mean,
    "same-point-narrow-interval": _reference_mean_narrow,
    "known-bias-plus-0.5": _known_biased,
    "independent-half-mean": _independent_half_mean,
}


def run(
    *,
    seed: int,
    n_replicates: int = 2000,
    n_patients: int = 200,
    effects: Sequence[float] = REQUESTED_EFFECTS,
) -> pd.DataFrame:
    rows: list[dict] = []
    for effect in effects:
        for replicate in range(n_replicates):
            rng = np.random.default_rng([seed, int(np.float64(effect).view(np.uint64)), replicate])
            values = rng.normal(loc=effect, scale=1.0, size=n_patients)
            # Defined independently and before estimator evaluation.
            reference = float(values.mean())
            for name, estimator in ESTIMATORS.items():
                result = estimator(values, reference)
                valid = bool(np.isfinite([result.point, result.low, result.high]).all())
                rows.append(
                    {
                        "requested_effect": float(effect),
                        "replicate": replicate,
                        "n_patients": n_patients,
                        "estimator": name,
                        "empirical_reference": reference,
                        "estimate": result.point,
                        "ci_low": result.low,
                        "ci_high": result.high,
                        "valid_output": valid,
                        "residual_vs_reference": (
                            abs(result.point - reference) if valid else math.nan
                        ),
                        "difference_vs_requested": (
                            result.point - effect if valid else math.nan
                        ),
                        "covered_requested": (
                            result.low <= effect <= result.high if valid else pd.NA
                        ),
                    }
                )
    return pd.DataFrame(rows)


def summarise(runs: pd.DataFrame, *, tolerance: float = REFERENCE_TOLERANCE) -> pd.DataFrame:
    grouped = runs.groupby(["requested_effect", "estimator"], sort=True)
    out = grouped.agg(
        n_attempted=("replicate", "count"),
        n_valid=("valid_output", "sum"),
        max_residual_vs_reference=("residual_vs_reference", "max"),
        median_residual_vs_reference=("residual_vs_reference", "median"),
        bias=("difference_vs_requested", "mean"),
        rmse=(
            "difference_vs_requested",
            lambda s: float(np.sqrt(np.mean(np.square(s.dropna())))),
        ),
        interval_coverage=("covered_requested", "mean"),
    ).reset_index()
    out["n_valid"] = out["n_valid"].astype(int)
    out["n_invalid"] = out["n_attempted"] - out["n_valid"]
    out["valid_output_rate"] = out["n_valid"] / out["n_attempted"]
    out["coverage_mc_se"] = np.sqrt(
        out["interval_coverage"] * (1.0 - out["interval_coverage"]) / out["n_valid"]
    )
    out["reference_audit_outcome"] = np.select(
        [
            out["max_residual_vs_reference"] > tolerance,
            out["n_valid"] < out["n_attempted"],
        ],
        ["numerical_departure", "insufficient_valid_pairs"],
        default="numerical_equality",
    )
    return out


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=20260913)
    parser.add_argument("--replicates", type=int, default=2000)
    parser.add_argument("--patients", type=int, default=200)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    table = summarise(
        run(seed=args.seed, n_replicates=args.replicates, n_patients=args.patients)
    )
    output = write_versioned_table(
        table,
        "residual_performance_clean_control",
        seed=args.seed,
        results_dir=args.results_dir,
        allow_dirty=args.allow_dirty,
        extra_meta={
            "reference": (
                "Full-sample empirical mean, computed before estimator evaluation."
            ),
            "reference_tolerance": REFERENCE_TOLERANCE,
            "requested_effects": list(REQUESTED_EFFECTS),
            "known_bias": KNOWN_BIAS,
            "n_replicates": args.replicates,
            "n_patients": args.patients,
            "interval": "Normal 95% interval; narrow control uses 0.25 times its half-width.",
            "interpretation": (
                "Reference equality and departure are arithmetic audit outcomes, not "
                "certificates of estimator quality. Bias, RMSE and coverage are evaluated "
                "against the prospectively requested effect."
            ),
            "SYNTHETIC": "No clinical or biological data are used.",
        },
    )
    log.info("\n%s", table.to_string(index=False))
    log.info("wrote %s", output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
