"""Run the diagnosis and the repair attempt. Pre-registered in
``docs/prereg_repaired_interval.md``, committed before this file produced a row.

    python -m src.harness.run_interval_repair --task both

Two studies share one driver because they share one generator and must share one
set of holdout draws:

**TASK 1, DIAGNOSIS.** What each of the reviewer's three candidate causes -- zero
inflation, skew, fixed-fraction resampling -- is worth in percentage points of
null rejection, at 50 and at 800 mature cells. Measured by substituting the
resampling population and re-measuring, in ``harness.interval_diagnosis``.

**TASK 2, REPAIR.** Six intervals (the status-quo percentile baseline and five
candidates) on the *same* generated samples, scored for null rejection at
``shift = 1.0``, and for coverage and discrimination at the pre-registered
detectable effect ``shift = 0.5``. Every rate carries its binomial Monte-Carlo
standard error, because the difference this exists to resolve -- 5% against
9.5% -- is under two standard errors at the 200 replicates the committed
diagnostic used.

WHY THE REAL GENERATOR
----------------------
Under the null the mature cells of both arms are i.i.d. draws from one pool, so
the repair study *could* be run from that pool alone, much faster. It is not.
Every replicate goes through ``pseudobulk.generate_pseudobulk`` exactly as
``papers/whms_bode/diagnose_low_counts.py`` does, so there is no question of
whether a fast path reproduces the generator -- the same code produced both.
The diagnosis, which must substitute the population, necessarily does not go
through the generator, and ``check_diagnosis_reproduces_the_generator`` is the
guard that the empirical family still lands where the generator does.

DATA
----
The Lee matrices are not in every worktree. Set ``BRP_DATA_DIR`` to a checkout
that has them; the cohort arrays are then cached to a single ``.npz`` so the
16 worker processes do not each spend two minutes parsing a 129 MB gzip.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
from collections.abc import Sequence
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Final

import numpy as np
import pandas as pd

from src.common.io import write_versioned_table
from src.harness.attenuation import _composition, _configuration_seed
from src.harness.calibration_gap import AXIS, MATURE_BIN, RUNG, TARGET_GENE, _pool_mask
from src.harness.interval_diagnosis import (
    DESIGNS,
    POOL_FAMILIES,
    attribute,
    check_design_effect_vanishes_where_the_arms_coincide,
    check_gaussian_floor_matches_its_arithmetic,
    closure_verdict,
    diagnosis_closes,
    gaussian_floor,
    null_rejection,
)
from src.harness.interval_repair import (
    CANDIDATES,
    DEFAULT_N_BOOT,
    all_intervals,
    covers,
    excludes_zero,
    mcse,
    width,
)
from src.harness.pseudobulk import generate_pseudobulk, patient_holdout

log = logging.getLogger(__name__)

COHORTS: Final[tuple[str, ...]] = ("smc", "kul3")
POOLS: Final[tuple[str, ...]] = ("pooled", "reference")

#: The counts the committed diagnostic reports, unchanged so the new rates can
#: be laid beside the old ones without a grid change confounding the comparison.
COUNTS: Final[tuple[int, ...]] = (5, 50, 100, 800)

#: The two focal counts the pre-registration attributes causes at.
FOCAL_COUNTS: Final[tuple[int, ...]] = (50, 800)

#: 1.0 is the null; 0.5 is ``calibration.PREREGISTERED.detectable_shift``.
SHIFTS: Final[tuple[float, ...]] = (1.0, 0.5)

#: Patients held out per replicate. 2 is what the committed sweep uses; 5 exists
#: so "the patient-level candidates fail" can be separated from "two clusters is
#: not a sample". KUL3 has six patients, so 5 held out leaves one to train on --
#: legal for ``patient_holdout``, and the bulk arm is not run here.
HOLDOUTS: Final[tuple[int, ...]] = (2, 5)

#: All six at the committed holdout. At holdout=5 only the candidates whose
#: behaviour can depend on the patient count, plus the percentile baseline as an
#: anchor -- running the cell-level candidates twice would double the cost to
#: re-measure a number that cannot move.
CANDIDATES_AT: Final[dict[int, tuple[str, ...]]] = {
    2: CANDIDATES,
    5: ("percentile", "patient_cluster", "pseudobulk_patient_t"),
}

#: Five independent streams. The first is the project seed and the one the
#: committed diagnostic used, so its stream is directly comparable.
SEED_STREAMS: Final[tuple[int, ...]] = (20260831, 1, 2, 3, 4)
REPLICATES_PER_STREAM: Final[int] = 400

N_CELLS: Final[int] = 2000
FRAC_MATURE_NORMAL: Final[float] = 0.40

#: Pre-registered targets. ``calibration.PREREGISTERED`` for the last two.
NULL_TARGET: Final[float] = 0.05
COVERAGE_TARGET: Final[float] = 0.90
DISCRIMINATION_TARGET: Final[float] = 0.80
MAX_ABSTENTION: Final[float] = 0.10
MAX_WIDTH_RATIO: Final[float] = 3.0


# --------------------------------------------------------------------------
# cohort arrays, loaded once and cached
# --------------------------------------------------------------------------


def cache_cohort(cohort: str, cache_dir: Path) -> Path:
    """Parse one Lee cohort once and cache the arrays the sweep needs."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / f"lee_{cohort}.npz"
    if path.exists():
        return path
    from src.harness.calibration_gap import load_cohort_arrays

    log.info("parsing Lee/%s (this is the slow part, once) …", cohort.upper())
    a = load_cohort_arrays(cohort)
    np.savez_compressed(
        path,
        counts=a["counts"],
        cell_type=np.asarray(a["cell_type"]),
        patient_id=np.asarray(a["patient_id"]),
        tissue=np.asarray(a["tissue"]),
        genes=np.asarray(a["genes"]),
        study_id=np.asarray([a["study_id"]]),
        n_patients=np.asarray([a["n_patients"]]),
    )
    log.info("cached %s (%d cells)", path, a["counts"].shape[0])
    return path


def _load_cached(path: Path, pool: str) -> dict:
    """Cohort arrays restricted to one pool reading."""
    z = np.load(path, allow_pickle=False)
    tissue = pd.Series(z["tissue"].astype(str))
    mask = _pool_mask(tissue, pool)
    return {
        "counts": z["counts"][mask],
        "cell_type": z["cell_type"].astype(str)[mask],
        "patient_id": z["patient_id"].astype(str)[mask],
        "genes": [str(g) for g in z["genes"]],
        "study_id": str(z["study_id"][0]),
        "n_patients": int(z["n_patients"][0]),
    }


# --------------------------------------------------------------------------
# one worker cell of the repair study
# --------------------------------------------------------------------------


def _replicate_seed(stream: int, count: int, shift: float, rep: int) -> int:
    """The committed diagnostic's own key, so streams line up with it."""
    return _configuration_seed(
        stream,
        n_cells=N_CELLS,
        mature_fraction=count / N_CELLS,
        shift=shift,
        replicate=rep,
    )


def repair_cell(job: dict) -> list[dict]:
    """Every candidate, every replicate, for one (cohort, pool, holdout, count,
    shift, stream). Returns one row per candidate.

    Runs in a worker process, so it takes a plain dict and returns plain dicts.
    """
    arrays = _load_cached(Path(job["cache"]), job["pool"])
    counts = arrays["counts"]
    cell_type = arrays["cell_type"]
    patient_id = arrays["patient_id"]
    genes = arrays["genes"]
    types = sorted(set(cell_type.tolist()))
    comp_n = _composition(FRAC_MATURE_NORMAL, types, MATURE_BIN)
    comp_t = _composition(job["count"] / N_CELLS, types, MATURE_BIN)
    candidates = CANDIDATES_AT[job["holdout"]]

    tally = {
        c: {"scored": 0, "rejects": 0, "covered": 0, "abstained": 0, "widths": []}
        for c in candidates
    }
    n_attempted = 0
    for rep in range(job["replicates"]):
        seed = _replicate_seed(job["stream"], job["count"], job["shift"], rep)
        _, held = patient_holdout(patient_id, n_held_out=job["holdout"], seed=seed)
        sample = generate_pseudobulk(
            counts, cell_type, patient_id, genes,
            composition_normal=comp_n, composition_tumour=comp_t,
            shift={TARGET_GENE: job["shift"]}, held_out_patients=held,
            n_cells=N_CELLS, seed=seed, mature_label=MATURE_BIN,
        )
        cells = sample.mature_expression[TARGET_GENE]
        normal, tumour = cells["normal"], cells["tumour"]
        pats = {
            arm: sample.drawn_patient_id[arm][sample.drawn_is_mature[arm]]
            for arm in ("normal", "tumour")
        }
        truth = float(sample.truth.parametric[TARGET_GENE]["normal"]["intrinsic"])
        n_attempted += 1

        intervals = all_intervals(
            normal, tumour,
            frac_mature_normal=FRAC_MATURE_NORMAL,
            seed=seed, n_boot=job["n_boot"],
            normal_patients=pats["normal"], tumour_patients=pats["tumour"],
            candidates=candidates,
        )
        for name, interval in intervals.items():
            rejected = excludes_zero(interval)
            if rejected is None:
                tally[name]["abstained"] += 1
                continue
            tally[name]["scored"] += 1
            tally[name]["rejects"] += int(rejected)
            tally[name]["covered"] += int(bool(covers(interval, truth)))
            tally[name]["widths"].append(width(interval))

    rows = []
    for name, t in tally.items():
        scored = t["scored"]
        rate = t["rejects"] / scored if scored else float("nan")
        cov = t["covered"] / scored if scored else float("nan")
        rows.append({
            "cohort": job["cohort"],
            "pool": job["pool"],
            "n_held_out": job["holdout"],
            "n_cells_mature": job["count"],
            "shift": job["shift"],
            "candidate": name,
            "seed_stream": job["stream"],
            "n_attempted": n_attempted,
            "n_scored": scored,
            "n_abstained": t["abstained"],
            "n_excludes_zero": t["rejects"],
            "n_covered": t["covered"],
            "excludes_zero_rate": rate,
            "coverage": cov,
            "median_ci_width": float(np.median(t["widths"])) if t["widths"] else float("nan"),
        })
    return rows


def diagnosis_cell(job: dict) -> list[dict]:
    """Null rejection under every pool family and design, for one cohort/pool/count.

    The empirical pools are rebuilt per replicate from the same holdout draws
    the repair study uses, so the between-holdout variation that is part of the
    phenomenon survives the substitution.
    """
    arrays = _load_cached(Path(job["cache"]), job["pool"])
    cell_type = arrays["cell_type"]
    patient_id = arrays["patient_id"]
    target_col = arrays["genes"].index(TARGET_GENE)
    target_values = arrays["counts"][:, target_col].astype(float)
    is_mature = cell_type == MATURE_BIN

    pools, seeds = [], []
    for rep in range(job["replicates"]):
        seed = _replicate_seed(job["stream"], job["count"], 1.0, rep)
        _, held = patient_holdout(patient_id, n_held_out=job["holdout"], seed=seed)
        keep = is_mature & np.isin(patient_id, list(held))
        pools.append(target_values[keep])
        seeds.append(seed)

    rows = []
    for family in POOL_FAMILIES:
        for design in DESIGNS:
            row = null_rejection(
                pools, count=job["count"], family=family, design=design,
                n_boot=job["n_boot"], seeds=seeds,
            )
            rows.append(
                row | {
                    "cohort": job["cohort"],
                    "pool": job["pool"],
                    "n_held_out": job["holdout"],
                    "seed_stream": job["stream"],
                    "closed_form_floor": gaussian_floor(
                        row["n_normal"], row["n_tumour"]
                    ),
                }
            )
    return rows


# --------------------------------------------------------------------------
# aggregation and verdicts
# --------------------------------------------------------------------------


def pool_streams(per_stream: pd.DataFrame, keys: Sequence[str]) -> pd.DataFrame:
    """Combine the seed streams into one rate, keeping the spread visible.

    Rates are recomputed from the pooled counts rather than averaged across
    streams: a mean of five ratios with different denominators is not the ratio
    of the sums, and the denominators differ whenever a stream abstains more.
    The per-stream min and max travel alongside, so a rate that is stable only
    because it was averaged is still visible as one that is not.
    """
    grouped = per_stream.groupby(list(keys), observed=True)
    out = grouped.agg(
        n_streams=("seed_stream", "nunique"),
        n_attempted=("n_attempted", "sum"),
        n_scored=("n_scored", "sum"),
        n_abstained=("n_abstained", "sum"),
        n_excludes_zero=("n_excludes_zero", "sum"),
        n_covered=("n_covered", "sum"),
        median_ci_width=("median_ci_width", "median"),
        stream_min_rate=("excludes_zero_rate", "min"),
        stream_max_rate=("excludes_zero_rate", "max"),
    ).reset_index()
    out["excludes_zero_rate"] = out["n_excludes_zero"] / out["n_scored"].replace(0, np.nan)
    out["coverage"] = out["n_covered"] / out["n_scored"].replace(0, np.nan)
    out["abstention_rate"] = out["n_abstained"] / out["n_attempted"]
    out["excludes_zero_mcse"] = [
        mcse(r, n) for r, n in zip(out["excludes_zero_rate"], out["n_scored"], strict=True)
    ]
    out["coverage_mcse"] = [
        mcse(c, n) for c, n in zip(out["coverage"], out["n_scored"], strict=True)
    ]
    return out


def repair_table(pooled: pd.DataFrame) -> pd.DataFrame:
    """One row per (cohort, pool, holdout, count, candidate): null beside power.

    The null rejection lives on the ``shift = 1.0`` rows and the coverage and
    discrimination on the ``shift = 0.5`` rows. They are joined into ONE row on
    purpose. A candidate that reaches nominal by refusing to reject anything is
    then impossible to quote without its discrimination sitting next to the
    number, which is the failure mode falsifier F2 exists to catch.
    """
    keys = ["cohort", "pool", "n_held_out", "n_cells_mature", "candidate"]
    null = pooled[pooled["shift"] == 1.0].set_index(keys)
    alt = pooled[pooled["shift"] == 0.5].set_index(keys)

    joined = pd.DataFrame(index=null.index.union(alt.index))
    joined["null_rejection"] = null["excludes_zero_rate"]
    joined["null_rejection_mcse"] = null["excludes_zero_mcse"]
    joined["null_stream_min"] = null["stream_min_rate"]
    joined["null_stream_max"] = null["stream_max_rate"]
    joined["null_abstention_rate"] = null["abstention_rate"]
    joined["null_median_width"] = null["median_ci_width"]
    joined["n_scored_null"] = null["n_scored"]
    joined["discrimination"] = alt["excludes_zero_rate"]
    joined["discrimination_mcse"] = alt["excludes_zero_mcse"]
    joined["coverage"] = alt["coverage"]
    joined["coverage_mcse"] = alt["coverage_mcse"]
    joined["alt_abstention_rate"] = alt["abstention_rate"]
    joined["alt_median_width"] = alt["median_ci_width"]
    joined["n_scored_alt"] = alt["n_scored"]
    joined = joined.reset_index()

    baseline = joined[joined["candidate"] == "percentile"].set_index(
        ["cohort", "pool", "n_held_out", "n_cells_mature"]
    )
    idx = pd.MultiIndex.from_frame(
        joined[["cohort", "pool", "n_held_out", "n_cells_mature"]]
    )
    joined["width_ratio_vs_percentile"] = (
        joined["null_median_width"].to_numpy()
        / baseline["null_median_width"].reindex(idx).to_numpy()
    )
    joined["baseline_discrimination"] = baseline["discrimination"].reindex(idx).to_numpy()
    joined["verdict"] = [verdict(row) for _, row in joined.iterrows()]
    joined["interval_object"] = "within_patient_over_cells"
    return joined.sort_values(
        ["cohort", "pool", "n_held_out", "n_cells_mature", "candidate"]
    ).reset_index(drop=True)


def verdict(row) -> str:
    """The pre-registered three-part pass criterion, and the hollow passes.

    Order matters. A candidate is only called ``pass`` when it clears the null
    target AND answers often enough AND still discriminates. The two hollow
    passes are named rather than collapsed into ``fail``, because "reaches
    nominal by abstaining" and "never reaches nominal" are different findings
    about a method, and the first is the one that looks like success in a
    summary table.
    """
    rate, se = row["null_rejection"], row["null_rejection_mcse"]
    if not np.isfinite(rate):
        return "all_abstained"
    calibrated = rate <= NULL_TARGET + 2 * se
    if not calibrated:
        return "fails_null"
    if row["null_abstention_rate"] > MAX_ABSTENTION:
        return "passes_by_abstention"
    width_ratio = row["width_ratio_vs_percentile"]
    lost_power = (
        np.isfinite(row["discrimination"])
        and row["discrimination"] < DISCRIMINATION_TARGET
        and np.isfinite(row["baseline_discrimination"])
        and row["baseline_discrimination"] >= DISCRIMINATION_TARGET
    )
    if (np.isfinite(width_ratio) and width_ratio > MAX_WIDTH_RATIO) or lost_power:
        return "passes_by_width"
    if not np.isfinite(row["discrimination"]) or row["discrimination"] < DISCRIMINATION_TARGET:
        return "calibrated_but_underpowered"
    return "pass"


class ReproductionError(ValueError):
    """The diagnosis does not reproduce the generator it claims to describe."""


def check_diagnosis_reproduces_the_generator(
    diagnosis: pd.DataFrame, repair: pd.DataFrame, *, tolerance_mcse: float = 4.0
) -> None:
    """The ``empirical`` family under ``fixed_fraction`` IS the committed sweep.

    The diagnosis substitutes the resampling population, so it does not go
    through ``generate_pseudobulk``. Its baseline family therefore has to be
    checked against the thing it is standing in for -- otherwise every
    percentage point it attributes is attributed within a simulation that may
    not be the one under discussion.

    Compared against the percentile candidate's null rejection at the same
    cohort, pool and count, at ``n_held_out = 2``. More than four combined
    Monte-Carlo standard errors apart and the substitution is not faithful.
    """
    base = diagnosis[
        (diagnosis["family"] == "empirical")
        & (diagnosis["design"] == "fixed_fraction")
        & (diagnosis["n_held_out"] == 2)
        if "n_held_out" in diagnosis
        else (diagnosis["family"] == "empirical")
        & (diagnosis["design"] == "fixed_fraction")
    ][["cohort", "pool", "n_cells_mature", "null_rejection", "null_rejection_mcse"]]
    ref = repair[(repair["candidate"] == "percentile") & (repair["n_held_out"] == 2)][
        ["cohort", "pool", "n_cells_mature", "null_rejection", "null_rejection_mcse"]
    ]
    merged = base.merge(
        ref, on=["cohort", "pool", "n_cells_mature"],
        how="inner", suffixes=("_diagnosis", "_generator"),
    )
    bad = []
    for _, r in merged.iterrows():
        a, b = r["null_rejection_diagnosis"], r["null_rejection_generator"]
        if not (np.isfinite(a) and np.isfinite(b)):
            continue
        se = float(np.hypot(r["null_rejection_mcse_diagnosis"],
                            r["null_rejection_mcse_generator"]))
        if se <= 0:
            continue
        if abs(a - b) > tolerance_mcse * se:
            bad.append(
                f"{r['cohort']}/{r['pool']} n={r['n_cells_mature']}: "
                f"diagnosis {a:.4f} vs generator {b:.4f} "
                f"({abs(a - b) / se:.1f} MCSE)"
            )
    if not len(merged):
        raise ReproductionError(
            "no (cohort, pool, count) cell is shared between the diagnosis and "
            "the generator run, so the substitution was never checked against "
            "the thing it stands in for"
        )
    if bad:
        raise ReproductionError(
            "the diagnosis' empirical baseline does not reproduce the generator "
            "it stands in for, so nothing it attributes is an attribution about "
            "the committed sweep: " + "; ".join(bad)
        )


# --------------------------------------------------------------------------
# driver
# --------------------------------------------------------------------------


def _jobs(cache: dict[str, Path], *, replicates: int, n_boot: int, task: str) -> list[dict]:
    jobs: list[dict] = []
    for cohort in COHORTS:
        for pool in POOLS:
            for stream in SEED_STREAMS:
                if task in ("repair", "both"):
                    for holdout in HOLDOUTS:
                        for count in COUNTS:
                            for shift in SHIFTS:
                                jobs.append({
                                    "kind": "repair", "cohort": cohort, "pool": pool,
                                    "holdout": holdout, "count": count, "shift": shift,
                                    "stream": stream, "replicates": replicates,
                                    "n_boot": n_boot, "cache": str(cache[cohort]),
                                })
                if task in ("diagnosis", "both"):
                    for count in FOCAL_COUNTS:
                        jobs.append({
                            "kind": "diagnosis", "cohort": cohort, "pool": pool,
                            "holdout": 2, "count": count, "stream": stream,
                            "replicates": replicates, "n_boot": n_boot,
                            "cache": str(cache[cohort]),
                        })
    return jobs


def _run_job(job: dict) -> list[dict]:
    return repair_cell(job) if job["kind"] == "repair" else diagnosis_cell(job)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task", choices=("diagnosis", "repair", "both"), default="both")
    parser.add_argument("--replicates", type=int, default=REPLICATES_PER_STREAM,
                        help="per seed stream; five streams are run")
    parser.add_argument("--n-boot", type=int, default=DEFAULT_N_BOOT)
    parser.add_argument("--workers", type=int, default=min(os.cpu_count() or 4, 12))
    parser.add_argument("--cache-dir", type=Path, default=Path(".interval_cache"))
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--suffix", default="")
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    cache = {c: cache_cohort(c, args.cache_dir) for c in COHORTS}
    jobs = _jobs(cache, replicates=args.replicates, n_boot=args.n_boot, task=args.task)
    log.info("%d worker cells over %d processes", len(jobs), args.workers)

    rows: list[dict] = []
    done = 0
    with ProcessPoolExecutor(max_workers=args.workers) as pool_exec:
        futures = [pool_exec.submit(_run_job, j) for j in jobs]
        for fut in as_completed(futures):
            rows.extend(fut.result())
            done += 1
            if done % 25 == 0 or done == len(jobs):
                log.info("  %d/%d cells", done, len(jobs))

    frame = pd.DataFrame(rows)
    repair_rows = frame[frame["candidate"].notna()] if "candidate" in frame else pd.DataFrame()
    diag_rows = frame[frame["family"].notna()] if "family" in frame else pd.DataFrame()

    meta = {
        "prereg": "docs/prereg_repaired_interval.md",
        "interval_object": "within_patient_over_cells (src/harness/interval.py)",
        "not_this_object": (
            "the over-patients percentile bootstrap analysed in "
            "src/reference/interval_calibration.py -- a different estimand"
        ),
        "pool_source": "lee_raw",
        "cohorts": {"smc": "GSE132465, 10 matched patients",
                    "kul3": "GSE144735, 6 matched patients"},
        "target_gene": TARGET_GENE,
        "labelling": f"label_{AXIS}_{RUNG}, mature = {MATURE_BIN!r}",
        "candidates": list(CANDIDATES),
        "candidates_at_holdout": {str(k): list(v) for k, v in CANDIDATES_AT.items()},
        "counts": list(COUNTS),
        "shifts": list(SHIFTS),
        "holdouts": list(HOLDOUTS),
        "seed_streams": list(SEED_STREAMS),
        "replicates_per_stream": args.replicates,
        "replicates_total": args.replicates * len(SEED_STREAMS),
        "n_boot": args.n_boot,
        "n_cells_per_arm": N_CELLS,
        "frac_mature_normal": FRAC_MATURE_NORMAL,
        "targets": {"null_rejection": NULL_TARGET, "coverage": COVERAGE_TARGET,
                    "discrimination": DISCRIMINATION_TARGET,
                    "max_abstention": MAX_ABSTENTION,
                    "max_width_ratio_vs_baseline": MAX_WIDTH_RATIO},
        "post_hoc": False,
    }

    written = []
    if len(repair_rows):
        per_stream = repair_rows.drop(columns=[c for c in ("family", "design") if c in repair_rows])
        pooled = pool_streams(
            per_stream,
            ["cohort", "pool", "n_held_out", "n_cells_mature", "shift", "candidate"],
        )
        table = repair_table(pooled)
        for df, name in ((per_stream, "interval_repair_by_seed"),
                         (pooled, "interval_repair_pooled"),
                         (table, "interval_repair_rates")):
            written.append(write_versioned_table(
                df, f"{name}{args.suffix}", seed=SEED_STREAMS[0],
                results_dir=args.results_dir, extra_meta=meta,
                allow_dirty=args.allow_dirty,
            ))
        log.info("\n%s", table[[
            "cohort", "pool", "n_held_out", "n_cells_mature", "candidate",
            "null_rejection", "null_rejection_mcse", "discrimination", "verdict",
        ]].to_string(index=False))

    if len(diag_rows):
        keys = ["cohort", "pool", "n_cells_mature", "family", "design"]
        agg = diag_rows.groupby(keys, observed=True).agg(
            n_replicates=("n_replicates", "sum"),
            n_scored=("n_scored", "sum"),
            n_pool_undefined=("n_pool_undefined", "sum"),
            n_normal=("n_normal", "first"),
            n_tumour=("n_tumour", "first"),
            median_ci_width=("median_ci_width", "median"),
            pool_mean=("pool_mean", "mean"),
            pool_var=("pool_var", "mean"),
            pool_skew=("pool_skew", "mean"),
            pool_form=("pool_form", "first"),
            skew_matched=("skew_matched", "all"),
            closed_form_floor=("closed_form_floor", "first"),
            _rejects=("null_rejection", lambda s: np.nan),
        ).reset_index()
        # rates recomputed from pooled counts, not averaged across streams
        counts_by = diag_rows.assign(
            _r=diag_rows["null_rejection"] * diag_rows["n_scored"]
        ).groupby(keys, observed=True)["_r"].sum().reset_index()
        agg = agg.drop(columns=["_rejects"]).merge(counts_by, on=keys)
        agg["null_rejection"] = agg["_r"] / agg["n_scored"].replace(0, np.nan)
        agg["null_rejection_mcse"] = [
            mcse(r, n) for r, n in zip(agg["null_rejection"], agg["n_scored"], strict=True)
        ]
        agg = agg.drop(columns=["_r"])
        agg["interval_object"] = "within_patient_over_cells"

        check_gaussian_floor_matches_its_arithmetic(agg.to_dict("records"))
        check_design_effect_vanishes_where_the_arms_coincide(agg.to_dict("records"))

        attributions = []
        for (cohort, pool_name), block in agg.groupby(["cohort", "pool"], observed=True):
            for count in FOCAL_COUNTS:
                a = attribute(block.to_dict("records"), count=count)
                a |= {"cohort": cohort, "pool": pool_name,
                      "closure_verdict": closure_verdict(a),
                      "diagnosis_closes": diagnosis_closes(a)}
                attributions.append(a)
        attribution = pd.DataFrame(attributions)

        for df, name in ((agg, "interval_diagnosis_rates"),
                         (attribution, "interval_diagnosis_attribution")):
            written.append(write_versioned_table(
                df, f"{name}{args.suffix}", seed=SEED_STREAMS[0],
                results_dir=args.results_dir, extra_meta=meta,
                allow_dirty=args.allow_dirty,
            ))
        log.info("\n%s", attribution.to_string(index=False))

    if len(repair_rows) and len(diag_rows):
        check_diagnosis_reproduces_the_generator(agg, table)

    for path in written:
        log.info("wrote %s", path)
    print(json.dumps({"tables": [str(p) for p in written]}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
