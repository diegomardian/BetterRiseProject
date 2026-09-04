"""The cutpoint calibration, re-derived from committed code.

This module exists because it did not. The calibration behind the project's
central figure was produced by a script that was never committed, stamped with a
sha pointing at a commit that did not contain it — the provenance failure
recorded in ``src/common/provenance.py``. Re-deriving the number meant writing
the sweep down. This is that file.

It answers two questions the earlier run could not separate.

**Which draw pool?** ``pseudobulk.patient_holdout`` filters eligibility on
*patient* alone, so a sweep handed the whole cohort draws both synthetic arms
from the held-out patients' cells of both tissues, pooled. The design spec lists
what stays fixed and never mentions tissue, so nobody chose this — an
implementation detail did. Both readings run here (``POOLS``) and both are
reported, rather than the one that yields a cutpoint.

**Which grid?** The committed grid varies a mature *fraction* against a fixed
2,000 cells, so the mature counts it can reach are
``{0, 20, 40, 100, 200, 400, 800}`` — nothing between 40 and 100. A cutpoint
returned from that grid cannot be distinguished from any value in the gap.
``EXTENDED_FRACTIONS`` adds nine fractions inside it and changes nothing else.

Run it::

    python -m src.harness.calibration_gap                # 50 replicates, 13 seeds
    python -m src.harness.calibration_gap --replicates 500 --seeds 8

All three tables land under ``results/{date}_{sha7}/`` with a sidecar
naming the cohort, the grids, the criteria and the committed cutpoints it is
being compared against.
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
from src.harness.attenuation import (
    DEFAULT_MATURE_FRACTIONS,
    SweepConfig,
    SweepGrid,
    run_sweep,
)
from src.harness.calibration import (
    PREREGISTERED,
    CalibrationCriteria,
    calibrate_cutpoints,
    coverage_and_discrimination,
)
from src.harness.positivity import CUTPOINTS

log = logging.getLogger(__name__)

#: Cells drawn per synthetic arm. Held fixed across the grid, which is why a
#: grid over *fractions* is also a grid over counts — and why its resolution in
#: counts is a property of the fractions somebody chose.
N_CELLS: int = 2000

#: The committed grid. Seven fractions, reaching counts
#: ``{800, 400, 200, 100, 40, 20, 0}`` at ``N_CELLS``.
COMMITTED_FRACTIONS: tuple[float, ...] = DEFAULT_MATURE_FRACTIONS

#: The nine fractions the committed grid could not reach, filling 40..100 and
#: extending a little past it: counts 150, 130, 120, 110, 90, 80, 70, 60, 50.
#: Chosen to bracket the gap, not to land anywhere in particular.
GAP_FILLING_FRACTIONS: tuple[float, ...] = (
    0.075, 0.065, 0.060, 0.055, 0.045, 0.040, 0.035, 0.030, 0.025,
)  # fmt: skip

#: The committed grid plus the nine. Nothing else differs.
EXTENDED_FRACTIONS: tuple[float, ...] = tuple(
    sorted(set(COMMITTED_FRACTIONS) | set(GAP_FILLING_FRACTIONS), reverse=True)
)

GRIDS: dict[str, tuple[float, ...]] = {
    "committed": COMMITTED_FRACTIONS,
    "extended": EXTENDED_FRACTIONS,
}

#: ``pooled`` is what the sweep code implements: eligibility filtered on patient
#: alone, so both arms draw from both tissues. ``reference`` is the second
#: reading of the same unspecified choice — reference-tissue cells only.
POOLS: tuple[str, ...] = ("pooled", "reference")

#: Seeds. The first is the project seed; the rest exist so seed-to-seed spread
#: is visible rather than assumed away, because a cutpoint stable across seeds
#: and a cutpoint stable because the grid has nowhere else to land look
#: identical in a single run.
SEEDS: tuple[int, ...] = (20260831, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 99, 12345)

TARGET_GENE: str = "GUCA2A"
AXIS: str = "stem_pole"
RUNG: str = "lineage"
MATURE_BIN: str = "differentiated"

#: Bins for the coverage/discrimination table. Log-spaced by
#: ``calibration._bin_edges``; 12 is that function's own default.
N_BINS: int = 12


def counts_reachable(fractions: Sequence[float], n_cells: int = N_CELLS) -> list[int]:
    """Mature counts a fraction grid can reach at ``n_cells``.

    The grid's resolution in the quantity the rule is *stated* in. Reporting a
    cutpoint without this is reporting the first grid point past the criteria
    rather than where the criteria first hold.
    """
    return sorted({int(round(f * n_cells)) for f in fractions})


def _pool_mask(tissue: pd.Series, pool: str) -> np.ndarray:
    """Which cells a replicate may draw from, under one reading of the pool."""
    if pool == "pooled":
        return tissue.isin(["normal", "tumour"]).to_numpy()
    if pool == "reference":
        return (tissue == "normal").to_numpy()
    raise ValueError(f"unknown pool {pool!r}; known: {list(POOLS)}")


def recovery_summary(sweep: pd.DataFrame) -> pd.DataFrame:
    """The recovery curve, per grid point, with the check that kills it.

    A recovery curve plots recovered over true. The oracle arm applies
    ``decompose`` to the sweep's *own* realised summary statistics, so it
    reproduces them exactly and the ratio reduces to realised-over-parametric
    truth — a property of the generator, carrying no information about the
    estimator. ``max_abs_residual_vs_realised`` is that one-line check made into
    a column: it is identically zero exactly when the curve cannot see the
    estimator.

    The ratio is undefined at ``shift = 1.0``, where the parametric intrinsic
    term is exactly zero by construction. That is the design's own null control,
    and it is the one grid point at which this validation statistic cannot be
    evaluated at all. Those rows are counted, not silently dropped.
    """
    rows = sweep[sweep["arm"] == "oracle"].copy()
    rows["_residual"] = (
        rows["intrinsic_hat"] - rows["intrinsic_true_realised"]
    ).abs()
    defined = np.isfinite(rows["attenuation_ratio"].to_numpy(dtype=float))
    rows["_ratio"] = np.where(defined, rows["attenuation_ratio"], np.nan)

    out = (
        rows.groupby(["frac_mature_tumour", "shift"], observed=True)
        .agg(
            n_replicates=("replicate", "count"),
            n_ratio_undefined=("_ratio", lambda s: int(s.isna().sum())),
            median_n_cells_mature=("n_cells_mature", "median"),
            ratio_median=("_ratio", "median"),
            ratio_q25=("_ratio", lambda s: s.quantile(0.25)),
            ratio_q75=("_ratio", lambda s: s.quantile(0.75)),
            max_abs_residual_vs_realised=("_residual", "max"),
        )
        .reset_index()
    )
    return out.sort_values(["shift", "frac_mature_tumour"]).reset_index(drop=True)


def run_one(
    config: SweepConfig,
    fractions: Sequence[float],
    *,
    seed: int,
    n_replicates: int,
    criteria: CalibrationCriteria = PREREGISTERED,
    n_bins: int = N_BINS,
) -> tuple[pd.DataFrame, dict, pd.DataFrame]:
    """One (pool, grid, seed) cell: sweep, then calibrate.

    Returns the per-bin table, the cutpoint row and the recovery summary.
    ``calibrate_cutpoints`` raises when no bin meets both criteria — that is
    gate criterion G4 firing, a result rather than an error, so it is caught and
    recorded as ``returned_a_cutpoint = False`` instead of aborting the run.
    """
    grid = SweepGrid(
        mature_fractions=tuple(fractions),
        n_replicates=n_replicates,
        n_cells=N_CELLS,
    )
    sweep = run_sweep(config, grid, seed=seed, arms=("oracle",))
    bins = coverage_and_discrimination(sweep, criteria, n_bins=n_bins)

    row: dict = {
        "seed": seed,
        "n_grid_points": len(fractions),
        "max_discrimination": float(bins["discrimination"].max()),
    }
    try:
        report = calibrate_cutpoints(sweep, criteria, n_bins=n_bins)
    except ValueError:
        # No n meets both targets. The routine returns no cutpoint, and that is
        # the finding — see calibrate_cutpoints' own docstring.
        row |= {"returned_a_cutpoint": False, "ok": None, "wide": None}
    else:
        row |= {
            "returned_a_cutpoint": True,
            "ok": float(report.cutpoints.ok),
            "wide": float(report.cutpoints.wide),
        }
    return bins, row, recovery_summary(sweep)


def run_calibration_gap(
    counts: np.ndarray,
    cell_type: Sequence[str],
    patient_id: Sequence[str],
    tissue: Sequence[str],
    genes: Sequence[str],
    *,
    seeds: Sequence[int] = SEEDS,
    n_replicates: int = 50,
    pools: Sequence[str] = POOLS,
    grids: dict[str, tuple[float, ...]] | None = None,
    target_gene: str = TARGET_GENE,
    mature_label: str = MATURE_BIN,
    criteria: CalibrationCriteria = PREREGISTERED,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Every (pool, grid, seed) cell. Returns ``(bins, cutpoints, recovery)``.

    Everything except pool, grid and seed is held fixed, so a difference between
    two rows is attributable to one of those three.
    """
    grids = GRIDS if grids is None else grids
    tissue = pd.Series(np.asarray(tissue))
    counts = np.asarray(counts)
    cell_type = np.asarray(cell_type)
    patient_id = np.asarray(patient_id)

    bin_frames: list[pd.DataFrame] = []
    cut_rows: list[dict] = []
    rec_frames: list[pd.DataFrame] = []

    for pool in pools:
        mask = _pool_mask(tissue, pool)
        if not mask.any():
            raise ValueError(f"pool {pool!r} selected no cells")
        config = SweepConfig(
            counts=counts[mask],
            cell_type=cell_type[mask].tolist(),
            patient_id=patient_id[mask].tolist(),
            genes=list(genes),
            target_gene=target_gene,
            mature_label=mature_label,
        )
        n_patients = len(set(config.patient_id))
        log.info("pool=%s: %d cells, %d patients", pool, int(mask.sum()), n_patients)

        for grid_name, fractions in grids.items():
            for seed in seeds:
                bins, row, rec = run_one(
                    config, fractions, seed=seed,
                    n_replicates=n_replicates, criteria=criteria,
                )
                bin_frames.append(bins.assign(pool=pool, grid=grid_name, seed=seed))
                rec_frames.append(rec.assign(pool=pool, grid=grid_name, seed=seed))
                cut_rows.append({"pool": pool, "grid": grid_name} | row)
                log.info(
                    "  pool=%s grid=%s seed=%s -> ok=%s wide=%s (max disc %.3f)",
                    pool, grid_name, seed, row["ok"], row["wide"],
                    row["max_discrimination"],
                )

    cutpoints = pd.DataFrame(cut_rows)[
        ["pool", "grid", "seed", "n_grid_points",
         "returned_a_cutpoint", "ok", "wide", "max_discrimination"]
    ]  # fmt: skip
    return (
        pd.concat(bin_frames, ignore_index=True),
        cutpoints,
        pd.concat(rec_frames, ignore_index=True),
    )


def load_cohort_arrays(which: str = "smc", raw_dir: Path | None = None):
    """One Lee cohort, as the arrays ``run_sweep`` wants.

    Raw integer counts, not CP10K: the generator realises a multiplicative shift
    by binomial thinning and Poisson augmentation, both defined on counts, and
    refuses a non-integer matrix outright.

    ``which`` selects the cohort. Both are local and sha256-verified against
    ``data/manifest.csv``. Running the same sweep on ``kul3`` is what makes the
    difference between "this rule fails its own criteria on this cohort" and
    "this rule fails its own criteria" — and the blindness check is algebraic,
    so a second cohort is a real test of whether the algebra is the whole story.
    KUL3 carries a third tissue class, ``border``, which ``_pool_mask`` excludes
    from both pools because it is neither the reference nor the diseased arm.
    """
    from src.estimator.lee_io import load_lee_cohort
    from src.reference.labels import label_column

    cohort = load_lee_cohort(
        which,
        target_genes=[TARGET_GENE],
        axes=(AXIS,),
        rungs=(RUNG,),
        raw_dir=raw_dir,
        keep_raw_counts=True,
    )
    column = label_column(AXIS, RUNG)
    labels = cohort.labels[column]

    # Cells outside the epithelial compartment carry pd.NA rather than False.
    # They are not a maturity bin and must not be pooled into one.
    labelled = labels.notna()
    index = cohort.raw_counts.index[labelled.reindex(cohort.raw_counts.index, fill_value=False)]

    # ``raw_counts`` comes back object-dtyped from the streaming parser. The
    # generator asserts a finite integer matrix, and ``np.isfinite`` raises
    # rather than returns False on object dtype — so cast here, having checked
    # the values really are integral rather than assuming it.
    raw = cohort.raw_counts.loc[index].to_numpy()
    as_float = raw.astype(np.float64)
    if not np.all(np.equal(np.mod(as_float, 1), 0)):
        raise ValueError("raw_counts carries non-integer values; expected UMI counts")

    return {
        "counts": as_float.astype(np.int64),
        "genes": list(cohort.raw_counts.columns),
        "cell_type": labels.reindex(index).astype(str).tolist(),
        "patient_id": cohort.cells.loc[index, "patient_id"].tolist(),
        "tissue": cohort.cells.loc[index, "tissue"].tolist(),
        "n_excluded_patients": len(cohort.excluded_patients),
        "study_id": cohort.study_id,
        "n_patients": int(cohort.cells.loc[index, "patient_id"].nunique()),
    }


def load_smc_arrays(raw_dir: Path | None = None):
    """Back-compatible alias. ``load_cohort_arrays("smc", ...)``."""
    return load_cohort_arrays("smc", raw_dir)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=50)
    parser.add_argument(
        "--seeds", type=int, default=len(SEEDS),
        help="how many of the committed seed list to use",
    )
    parser.add_argument(
        "--cohort", choices=("smc", "kul3"), default="smc",
        help="which Lee cohort to calibrate on. Both are local. Default smc, "
             "which is the committed sweep.",
    )
    parser.add_argument("--raw-dir", type=Path, default=None)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument(
        "--suffix", default=None,
        help="table-name suffix; defaults to _r{replicates} above 50",
    )
    parser.add_argument(
        "--allow-dirty", action="store_true",
        help="write from a dirty tree. The recorded sha will not reproduce the "
             "table, and the sidecar will say so.",
    )
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(message)s")
    seeds = SEEDS[: args.seeds]
    suffix = args.suffix
    if suffix is None:
        # smc keeps the bare names the committed tables already use; any other
        # cohort is namespaced so it cannot overwrite them.
        suffix = "" if args.cohort == "smc" else f"_{args.cohort}"
        if args.replicates != 50:
            suffix += f"_r{args.replicates}"

    log.info("loading Lee/%s …", args.cohort.upper())
    arrays = load_cohort_arrays(args.cohort, args.raw_dir)
    log.info("%d labelled epithelial cells", len(arrays["patient_id"]))

    bins, cutpoints, recovery = run_calibration_gap(
        arrays["counts"], arrays["cell_type"], arrays["patient_id"],
        arrays["tissue"], arrays["genes"],
        seeds=seeds, n_replicates=args.replicates,
    )

    meta = {
        "n_replicates": args.replicates,
        "n_boot": SweepGrid().n_boot,
        "n_bins": N_BINS,
        "n_cells_per_arm": N_CELLS,
        "n_seeds": len(seeds),
        "seeds": list(seeds),
        "cohort": (
            f"{arrays['study_id']} ({args.cohort.upper()}), "
            f"{arrays['n_patients']} matched patients, "
            f"sha256-verified against data/manifest.csv"
        ),
        "target_gene": TARGET_GENE,
        "labelling": f"label_{AXIS}_{RUNG}, mature = {MATURE_BIN!r}",
        "pools": {
            "pooled": "eligibility filtered on patient alone — both arms draw "
                      "from both tissues. What the sweep code implements.",
            "reference": "reference-tissue cells only. A second reading of the "
                         "same unspecified design choice.",
        },
        "committed_grid_targets": counts_reachable(COMMITTED_FRACTIONS),
        "extended_grid_targets": counts_reachable(EXTENDED_FRACTIONS),
        "committed_cutpoints": {"ok": CUTPOINTS.ok, "wide": CUTPOINTS.wide},
        "criteria": {
            "detectable_shift": PREREGISTERED.detectable_shift,
            "coverage_target": PREREGISTERED.coverage_target,
            "discrimination_target": PREREGISTERED.discrimination_target,
        },
        "what_this_answers": (
            "The committed grid has no point between 40 and 100 mature cells, "
            "so a cutpoint returned from it cannot be distinguished from any "
            "value in that interval. The extended grid adds nine fractions "
            "inside it and changes nothing else."
        ),
    }

    n_sweep_rows = int(recovery["n_replicates"].sum())
    max_residual = float(recovery["max_abs_residual_vs_realised"].max())
    meta |= {
        "n_oracle_sweep_rows": n_sweep_rows,
        "max_abs_residual_vs_realised": max_residual,
        "recovery_curve_note": (
            "The oracle estimate reproduces the sweep's own realised summary "
            "statistics, so max|intrinsic_hat - intrinsic_true_realised| is "
            "identically zero and the recovery ratio reduces to "
            "realised/parametric truth. The curve carries no information about "
            "the estimator. This field is that one-line check."
        ),
    }
    log.info("oracle sweep rows: %d | max |i_hat - i_realised| = %.3g",
             n_sweep_rows, max_residual)

    for df, name in ((bins, f"calibration_gap_bins{suffix}"),
                     (cutpoints, f"calibration_gap_cutpoints{suffix}"),
                     (recovery, f"calibration_gap_recovery{suffix}")):
        path = write_versioned_table(
            df, name, seed=SEEDS[0], results_dir=args.results_dir,
            extra_meta=meta, allow_dirty=args.allow_dirty,
        )
        log.info("wrote %s (%d rows)", path, len(df))

    returned = cutpoints.groupby(["pool", "grid"])["returned_a_cutpoint"].sum()
    log.info("\ncutpoints returned, by pool and grid:\n%s", returned.to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
