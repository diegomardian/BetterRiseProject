"""Controlled comparison of the committed, extended and dense count grids.

The union is simulated once. Grid-specific tables are then subsets of that one
run, so a setting shared by two grids is the same outer draw and the same inner
bootstrap stream. This separates evaluated counts and binning from random-draw
changes that affected the historical, independently run sweeps.
"""

from __future__ import annotations

import argparse
import logging
import math
import sys
from collections.abc import Mapping, Sequence
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

from src.common.io import write_versioned_table
from src.common.provenance import DEFAULT_SEED
from src.harness.attenuation import SweepConfig, SweepGrid, run_sweep
from src.harness.calibration import (
    PREREGISTERED,
    CalibrationCriteria,
    _bin_edges,
    coverage_and_discrimination,
    coverage_and_discrimination_by_count,
)
from src.harness.calibration_gap import (
    COMMITTED_FRACTIONS,
    EXTENDED_FRACTIONS,
    N_BINS,
    N_CELLS,
    POOLS,
    SEEDS,
    _pool_mask,
    load_cohort_arrays,
)
from src.reference.jobs.cutpoint_dense_grid import DENSE_FRACTIONS

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

log = logging.getLogger(__name__)

GRID_SUBSETS: dict[str, tuple[float, ...]] = {
    "committed": COMMITTED_FRACTIONS,
    "extended": EXTENDED_FRACTIONS,
    "dense": DENSE_FRACTIONS,
}
UNION_FRACTIONS: tuple[float, ...] = tuple(
    sorted({f for values in GRID_SUBSETS.values() for f in values}, reverse=True)
)


def fixed_edges_for_union(
    fractions: Sequence[float] = UNION_FRACTIONS,
    *,
    n_cells: int = N_CELLS,
    n_bins: int = N_BINS,
) -> np.ndarray:
    """Common bin edges derived once from the union of requested counts."""
    counts = pd.Series([int(round(f * n_cells)) for f in fractions])
    return _bin_edges(counts, n_bins)


def _with_mc_uncertainty(table: pd.DataFrame) -> pd.DataFrame:
    out = table.copy()
    n = out["n_estimated"].to_numpy(dtype=float)
    for metric in ("coverage", "discrimination"):
        rate = out[metric].to_numpy(dtype=float)
        out[f"{metric}_mc_se"] = np.sqrt(rate * (1.0 - rate) / n)
        out.loc[n == 0, f"{metric}_mc_se"] = np.nan
    return out


def _select_fractions(sweep: pd.DataFrame, fractions: Sequence[float]) -> pd.DataFrame:
    requested = np.asarray(tuple(fractions), dtype=float)
    keep = np.isclose(
        sweep["frac_mature_tumour"].to_numpy(dtype=float)[:, None],
        requested[None, :],
        rtol=0.0,
        atol=1e-12,
    ).any(axis=1)
    return sweep.loc[keep].copy()


def compare_grid_subsets(
    union_sweep: pd.DataFrame,
    *,
    subsets: Mapping[str, Sequence[float]] = GRID_SUBSETS,
    fixed_bin_edges: Sequence[float] | None = None,
    criteria: CalibrationCriteria = PREREGISTERED,
) -> pd.DataFrame:
    """Three binning readings for each grid subset of one union sweep."""
    edges = (
        fixed_edges_for_union()
        if fixed_bin_edges is None
        else np.asarray(fixed_bin_edges, dtype=float)
    )
    frames: list[pd.DataFrame] = []
    for grid_name, fractions in subsets.items():
        rows = _select_fractions(union_sweep, fractions)
        if rows.empty:
            raise ValueError(f"grid subset {grid_name!r} selected no rows")
        readings = {
            "adaptive_bins": coverage_and_discrimination(rows, criteria, n_bins=N_BINS),
            "fixed_bins": coverage_and_discrimination(
                rows,
                criteria,
                n_bins=N_BINS,
                bin_edges=edges,
            ),
            "per_count": coverage_and_discrimination_by_count(rows, criteria),
        }
        for binning, table in readings.items():
            frames.append(
                _with_mc_uncertainty(table).assign(
                    grid=grid_name,
                    binning=binning,
                )
            )
    return pd.concat(frames, ignore_index=True)


def crossing_brackets(
    rates: pd.DataFrame,
    criteria: CalibrationCriteria = PREREGISTERED,
) -> pd.DataFrame:
    """Candidate crossings and adjacent evaluated counts, without endorsement."""
    keys = [
        column
        for column in (
            "cohort",
            "pool",
            "seed",
            "grid",
            "binning",
            "inner_bootstrap",
            "omitted_patient",
        )
        if column in rates
    ]
    rows: list[dict] = []
    for key, group in rates.groupby(keys, dropna=False, sort=True):
        key_values = (key,) if len(keys) == 1 else key
        base = dict(zip(keys, key_values, strict=True))
        ordered = group.sort_values("n_cells_mature").reset_index(drop=True)
        definitions = {
            "coverage_and_discrimination": (
                (ordered["coverage"] >= criteria.coverage_target)
                & (ordered["discrimination"] >= criteria.discrimination_target)
            ),
            "coverage_only": ordered["coverage"] >= criteria.coverage_target,
        }
        candidates: dict[str, float] = {}
        for criterion, qualifies in definitions.items():
            positions = np.flatnonzero(qualifies.fillna(False).to_numpy())
            if not len(positions):
                rows.append(
                    base
                    | {
                        "criterion": criterion,
                        "lower_evaluated": np.nan,
                        "upper_evaluated": np.nan,
                        "candidate": np.nan,
                        "status": "no_crossing",
                    }
                )
                candidates[criterion] = math.nan
                continue
            position = int(positions[0])
            candidate = float(ordered.loc[position, "n_cells_mature"])
            lower = (
                float(ordered.loc[position - 1, "n_cells_mature"])
                if position > 0
                else math.nan
            )
            rows.append(
                base
                | {
                    "criterion": criterion,
                    "lower_evaluated": lower,
                    "upper_evaluated": candidate,
                    "candidate": candidate,
                    "status": "bracketed" if position > 0 else "lower_bound_unobserved",
                }
            )
            candidates[criterion] = candidate
        if (
            np.isfinite(candidates["coverage_and_discrimination"])
            and candidates["coverage_and_discrimination"] == candidates["coverage_only"]
        ):
            for row in rows[-2:]:
                row["status"] = "coincident_cutpoints"
    return pd.DataFrame(rows)


def plot_crossings(rates: pd.DataFrame, output: Path) -> None:
    """One compact per-count crossing plot, faceted by cohort and pool."""
    per_count = rates[rates["binning"] == "per_count"]
    facets = list(per_count.groupby(["cohort", "pool"], sort=True))
    fig, axes = plt.subplots(
        len(facets),
        1,
        figsize=(6.5, max(2.4, 2.1 * len(facets))),
        squeeze=False,
        sharex=True,
    )
    for ax, ((cohort, pool), panel) in zip(axes[:, 0], facets, strict=True):
        # The dense grid contains the union; plotting it once avoids drawing
        # identical shared counts three times.
        panel = panel[panel["grid"] == "dense"].sort_values("n_cells_mature")
        ax.errorbar(
            panel["n_cells_mature"],
            panel["coverage"],
            yerr=1.96 * panel["coverage_mc_se"],
            marker="o",
            ms=3,
            label="coverage",
        )
        ax.errorbar(
            panel["n_cells_mature"],
            panel["discrimination"],
            yerr=1.96 * panel["discrimination_mc_se"],
            marker="s",
            ms=3,
            label="discrimination",
        )
        ax.axhline(PREREGISTERED.coverage_target, color="C0", ls="--", lw=0.8)
        ax.axhline(PREREGISTERED.discrimination_target, color="C1", ls="--", lw=0.8)
        ax.set_xscale("symlog", linthresh=5)
        ax.set_ylim(-0.03, 1.08)
        ax.set_ylabel("rate")
        ax.set_title(f"{cohort.upper()} — {pool}", fontsize=9)
    axes[-1, 0].set_xlabel("mature cells per simulated tumour arm")
    axes[0, 0].legend(frameon=False, ncol=2)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def run_controlled(
    *,
    cohorts: Sequence[str],
    n_replicates: int,
    n_boot: int,
    seeds: Sequence[int],
    raw_dir: Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run each cohort/pool/seed union once and return rates and crossings."""
    rate_frames: list[pd.DataFrame] = []
    fixed_edges = fixed_edges_for_union()
    for cohort_name in cohorts:
        arrays = load_cohort_arrays(cohort_name, raw_dir)
        tissue = pd.Series(np.asarray(arrays["tissue"]))
        for pool in POOLS:
            mask = _pool_mask(tissue, pool)
            config = SweepConfig(
                counts=np.asarray(arrays["counts"])[mask],
                cell_type=np.asarray(arrays["cell_type"])[mask].tolist(),
                patient_id=np.asarray(arrays["patient_id"])[mask].tolist(),
                genes=list(arrays["genes"]),
                target_gene="GUCA2A",
                mature_label="differentiated",
            )
            grid = SweepGrid(
                mature_fractions=UNION_FRACTIONS,
                n_replicates=n_replicates,
                n_cells=N_CELLS,
                n_boot=n_boot,
            )
            for seed in seeds:
                log.info("cohort=%s pool=%s seed=%s", cohort_name, pool, seed)
                sweep = run_sweep(
                    config,
                    grid,
                    seed=seed,
                    arms=("oracle",),
                    seed_strategy="configuration",
                )
                rate_frames.append(
                    compare_grid_subsets(
                        sweep,
                        fixed_bin_edges=fixed_edges,
                    ).assign(cohort=cohort_name, pool=pool, seed=seed)
                )
    rates = pd.concat(rate_frames, ignore_index=True)
    return rates, crossing_brackets(rates)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort", choices=("smc", "kul3", "both"), default="both")
    parser.add_argument("--replicates", type=int, default=200)
    parser.add_argument("--n-boot", type=int, default=200)
    parser.add_argument("--seeds", type=int, default=8)
    parser.add_argument("--raw-dir", type=Path, default=None)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--suffix", default="")
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    cohorts = ("smc", "kul3") if args.cohort == "both" else (args.cohort,)
    seeds = SEEDS[: args.seeds]
    rates, crossings = run_controlled(
        cohorts=cohorts,
        n_replicates=args.replicates,
        n_boot=args.n_boot,
        seeds=seeds,
        raw_dir=args.raw_dir,
    )
    extra = {
        "cohorts": list(cohorts),
        "pools": list(POOLS),
        "grid_counts": {
            name: [int(round(f * N_CELLS)) for f in values]
            for name, values in GRID_SUBSETS.items()
        },
        "fixed_bin_edges": fixed_edges_for_union().tolist(),
        "n_replicates_per_seed": args.replicates,
        "inner_bootstrap_draws": args.n_boot,
        "seeds": list(seeds),
        "seed_strategy": "configuration",
        "uncertainty": (
            "Binomial Monte Carlo standard errors conditional on each fixed cohort. "
            "These are not patient-sampling confidence intervals."
        ),
    }
    suffix = args.suffix or f"_r{args.replicates}_b{args.n_boot}"
    rate_path = write_versioned_table(
        rates,
        f"controlled_grid_rates{suffix}",
        seed=DEFAULT_SEED,
        results_dir=args.results_dir,
        allow_dirty=args.allow_dirty,
        extra_meta=extra,
    )
    crossing_path = write_versioned_table(
        crossings,
        f"controlled_grid_crossings{suffix}",
        seed=DEFAULT_SEED,
        results_dir=args.results_dir,
        allow_dirty=args.allow_dirty,
        extra_meta=extra,
    )
    plot_path = rate_path.with_name(f"controlled_grid_crossing_plot{suffix}.pdf")
    plot_crossings(rates, plot_path)
    log.info("wrote %s", rate_path)
    log.info("wrote %s", crossing_path)
    log.info("wrote %s", plot_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
