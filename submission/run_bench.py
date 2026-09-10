#!/usr/bin/env python
"""Run the competitor benchmark and its three extensions.

    python -m submission.run_bench bench            # the committed six-world run
    python -m submission.run_bench sweep            # 1. detectable-effect sweep
    python -m submission.run_bench misspecification # 2. NB / zero-inflated arms
    python -m submission.run_bench interval-curve   # 3. width ratio, n = 2..200
    python -m submission.run_bench all              # all four

WHERE EACH THING IS WRITTEN, AND WHY IT IS NOT ONE PLACE
--------------------------------------------------------
``bench`` writes to ``submission/results/`` — deliberately NOT to the
project's ``results/``, which is versioned by git sha and governed by the
frozen schema. That was the original decision and it stands.

The three extensions write to ``results/{date}_{sha7}/`` through
``common.io.write_versioned_table``, because they are the tables a paper cites
and a cited table has to carry the sha that reproduces it.
``write_versioned_table`` REFUSES a dirty tree, untracked code included, so a
table can only be produced from committed code. Nothing here is a result about
colorectal cancer either way, and every sidecar says so.

Reproducibility: bit-identical across processes. World seeds are CRC32, not
``hash()`` (see ``bench.world_seed``); every entry point takes ``--seed`` and
defaults to the same pinned value.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from src.common.io import write_versioned_table
from src.common.paths import RESULTS_DIR
from src.common.provenance import provenance_record
from submission.bench import (
    BENCH_WORLDS,
    MEAN_NORMAL,
    N_CELLS,
    refusal_table,
    run_bench,
    sensitivity_where_estimable,
)
from submission.expression_models import EXPRESSION_MODELS
from submission.extensions import (
    CLOSED_FORM_AGREEMENT_TOLERANCE,
    COMMITTED_DETECTABLE_SHIFT,
    CURVE_N_MAX,
    CURVE_N_MIN,
    DETECTABLE_SHIFT_GRID,
    INVARIANT_WORLDS,
    MISSPECIFICATION_SEEDS,
    annihilated_invariance_proof,
    check_annihilated_refusal_is_model_invariant,
    check_closed_form_tracks_simulation,
    check_committed_value_is_in_grid,
    check_seeds_separate_model_from_stream,
    check_step_matches_grid,
    check_sweep_brackets_the_step,
    detectable_shift_sweep,
    interval_width_curve,
    misspecification_sweep,
    model_effect_vs_seed_noise,
)

SEED = 20260829
OUT = Path(__file__).parent / "results"

#: The committed simulation the closed-form curve is cross-checked against.
#: Pinned by sha, not globbed for the newest: a cross-check whose reference
#: moves when somebody reruns a job is not a cross-check.
SIMULATION_REFERENCE = "2026-09-06_a0483ae/interval_calibration.parquet"

SYNTHETIC_NOTE = (
    "Every number here comes from simulated cells with analytically known truth, "
    "or from a closed form with no data in it. Nothing in this table is a result "
    "about colorectal cancer."
)


def _banner(title: str) -> None:
    print("=" * 96)
    print(title)
    print("=" * 96)


# ---------------------------------------------------------------------------


def cmd_bench(args: argparse.Namespace) -> int:
    """The original six-world benchmark, unchanged."""
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"{len(BENCH_WORLDS)} worlds x {args.n_replicates} replicates, seed {args.seed}")
    print(f"{N_CELLS:,} cells per arm, mean_normal={MEAN_NORMAL}\n")

    bench, skipped = run_bench(seed=args.seed, n_replicates=args.n_replicates)
    refusals = refusal_table(bench)
    sensitivity = sensitivity_where_estimable(bench)

    for name, frame in (
        ("bench_raw", bench),
        ("refusal_table", refusals),
        ("sensitivity_where_estimable", sensitivity),
    ):
        frame.to_parquet(OUT / f"{name}.parquet", index=False)

    _banner(
        "REFUSAL — where the intrinsic estimand DOES NOT EXIST (annihilated: 0 mature\n"
        "tumour cells). false_confidence_rate counts numbers RETURNED, not numbers wrong."
    )
    print(refusals.to_string(index=False))
    print()
    _banner(
        "SENSITIVITY — where it DOES exist and a real intrinsic effect is present.\n"
        "The counterweight: refusing always scores 0 here."
    )
    print(sensitivity.to_string(index=False))

    if skipped:
        print("\nNOT RUN, with reasons (never omitted silently):")
        for name, why in sorted(skipped.items()):
            print(f"  {name}: {why}")

    meta = provenance_record(seed=args.seed, notes="Standalone competitor benchmark.")
    meta.update(
        {
            "n_replicates": args.n_replicates,
            "n_cells_per_arm": N_CELLS,
            "mean_normal": MEAN_NORMAL,
            "worlds": {w.name: w.why for w in BENCH_WORLDS},
            "skipped_methods": skipped,
            "SYNTHETIC": SYNTHETIC_NOTE,
        }
    )
    (OUT / "bench.meta.json").write_text(json.dumps(meta, indent=2))
    print(f"\nwrote 3 tables + bench.meta.json to {OUT}")
    return 0


def cmd_sweep(args: argparse.Namespace) -> int:
    """1. Sweep VarianceGateKitagawaMethod.DETECTABLE_SHIFT."""
    grid = DETECTABLE_SHIFT_GRID
    print(
        f"{len(grid)} detectable-effect values x {len(BENCH_WORLDS)} worlds x "
        f"{args.n_replicates} replicates, seed {args.seed}"
    )
    print(f"committed value: s_detect = {COMMITTED_DETECTABLE_SHIFT}\n")

    sweep, steps = detectable_shift_sweep(
        seed=args.seed, n_replicates=args.n_replicates, grid=grid
    )
    check_committed_value_is_in_grid(sweep)
    check_sweep_brackets_the_step(sweep)
    check_step_matches_grid(
        sweep, steps, seed=args.seed, n_replicates=args.n_replicates, grid=grid
    )

    wide = sweep.pivot_table(
        index="s_detect", columns="world", values="width_gate_refusals"
    )
    _banner(
        "1. WIDTH-GATE REFUSALS per world, sweeping the gate's own detectable effect.\n"
        f"Committed value {COMMITTED_DETECTABLE_SHIFT}, which is FAR BELOW the whole live\n"
        "range of the parameter — which is why it never binds. The gate is flat at 0\n"
        "until ~0.86, saturates by ~0.955, and by ~0.96 it refuses in worlds with ~100\n"
        "mature cells. See the step table for how narrow the usable band is."
    )
    print(wide.to_string())
    print()
    _banner(
        "THE STEP, exactly. s* is the per-replicate s_detect at which the gate flips;\n"
        "the width gate refuses iff s_detect > s*. step_width is max(s*) - min(s*)."
    )
    print(steps.to_string(index=False))

    written = [
        write_versioned_table(
            sweep,
            "width_gate_detectable_shift_sweep",
            seed=args.seed,
            notes=(
                "Width-gate refusals per world across a sweep of "
                "VarianceGateKitagawaMethod.DETECTABLE_SHIFT. SYNTHETIC."
            ),
            extra_meta={
                "n_replicates": int(args.n_replicates),
                "grid": [float(s) for s in grid],
                "committed_detectable_shift": float(COMMITTED_DETECTABLE_SHIFT),
                "n_cells_per_arm": N_CELLS,
                "mean_normal": MEAN_NORMAL,
                "generator": "poisson",
                "guards_passed": [
                    "check_committed_value_is_in_grid",
                    "check_sweep_brackets_the_step",
                    "check_step_matches_grid",
                ],
                "SYNTHETIC": SYNTHETIC_NOTE,
            },
            allow_dirty=args.allow_dirty,
        ),
        write_versioned_table(
            steps,
            "width_gate_step_location",
            seed=args.seed,
            notes=(
                "Exact per-world location of the width gate's step, from the "
                "closed-form critical s_detect per replicate. SYNTHETIC."
            ),
            extra_meta={
                "n_replicates": int(args.n_replicates),
                "committed_detectable_shift": float(COMMITTED_DETECTABLE_SHIFT),
                "SYNTHETIC": SYNTHETIC_NOTE,
            },
            allow_dirty=args.allow_dirty,
        ),
    ]
    for path in written:
        print(f"\nwrote {path}")
    return 0


def cmd_misspecification(args: argparse.Namespace) -> int:
    """2. The committed gates against negative-binomial and zero-inflated arms."""
    seeds = MISSPECIFICATION_SEEDS
    print(
        f"{len(EXPRESSION_MODELS)} count models x {len(BENCH_WORLDS)} worlds x "
        f"{args.n_replicates} replicates x {len(seeds)} seeds {list(seeds)}"
    )
    for model in EXPRESSION_MODELS:
        print(
            f"  {model.name:20s} kind={model.kind:8s} "
            f"Var(mu={MEAN_NORMAL:g}) = {model.theoretical_variance(MEAN_NORMAL):8.1f}"
        )
    print()

    sweep = misspecification_sweep(seeds=seeds, n_replicates=args.n_replicates)
    check_seeds_separate_model_from_stream(sweep)
    check_annihilated_refusal_is_model_invariant(sweep)
    separation = model_effect_vs_seed_noise(sweep)

    wide = sweep.pivot_table(
        index=["world", "expression_model"],
        columns=["gate", "seed"],
        values="n_refused",
    ).reindex(columns=["count_gate", "width_gate", "no_gate"], level=0)
    _banner(
        "2. REFUSALS under a MISSPECIFIED generator. Committed method classes,\n"
        "unchanged. Every model has the same marginal mean; only the variance moves.\n"
        "One column per seed, because at one seed the model and the RNG stream offset\n"
        "are perfectly confounded."
    )
    print(wide.to_string())
    print()
    _banner(
        "IS IT THE MODEL OR THE SEED? Changing the count model shifts the stream\n"
        "before the maturity draw, so each model draws a different realisation of the\n"
        "SAME Binomial. A gate whose model spread sits inside its seed spread has not\n"
        "responded to the generator at all."
    )
    print(separation.to_string(index=False))
    print()
    _banner(
        f"INVARIANCE in {list(INVARIANT_WORLDS)} — a theorem, not a measurement.\n"
        "The expression distribution is sampled ZERO times there."
    )
    print(annihilated_invariance_proof())

    path = write_versioned_table(
        sweep,
        "misspecified_generator_refusals",
        seed=args.seed,
        notes=(
            "Committed gate classes against Poisson, negative-binomial and "
            "zero-inflated mature-cell arms. SYNTHETIC."
        ),
        extra_meta={
            "n_replicates": int(args.n_replicates),
            "sweep_seeds": [int(x) for x in seeds],
            "n_cells_per_arm": N_CELLS,
            "mean_normal": MEAN_NORMAL,
            "generators": {
                m.name: {
                    "kind": m.kind,
                    "dispersion": m.dispersion,
                    "zero_inflation": m.zero_inflation,
                    "theoretical_variance_at_mu": m.theoretical_variance(MEAN_NORMAL),
                }
                for m in EXPRESSION_MODELS
            },
            "invariant_worlds": list(INVARIANT_WORLDS),
            "guards_passed": [
                "check_seeds_separate_model_from_stream",
                "check_annihilated_refusal_is_model_invariant",
            ],
            "SYNTHETIC": SYNTHETIC_NOTE,
        },
        allow_dirty=args.allow_dirty,
    )
    separation_path = write_versioned_table(
        separation,
        "misspecified_generator_model_vs_seed",
        seed=args.seed,
        notes=(
            "Per (world, gate): spread of refusals across count models versus "
            "across seeds. Says which differences in the refusal table can be "
            "attributed to the generator at all. SYNTHETIC."
        ),
        extra_meta={
            "sweep_seeds": [int(x) for x in seeds],
            "n_replicates": int(args.n_replicates),
            "SYNTHETIC": SYNTHETIC_NOTE,
        },
        allow_dirty=args.allow_dirty,
    )
    print(f"\nwrote {path}")
    print(f"wrote {separation_path}")
    return 0


def cmd_interval_curve(args: argparse.Namespace) -> int:
    """3. width_ratio and expected_false_positive_rate over the whole domain."""
    reference = RESULTS_DIR / SIMULATION_REFERENCE
    simulated = pd.read_parquet(reference) if reference.exists() else None
    if simulated is None:
        print(f"NOT CROSS-CHECKED: {reference} is absent. Closed form only.")

    curve = interval_width_curve(
        n_min=args.n_min, n_max=args.n_max, simulated=simulated
    )
    if simulated is not None:
        check_closed_form_tracks_simulation(curve)

    _banner(
        "3. THE PERCENTILE BOOTSTRAP'S MISCALIBRATION, as a function of n ALONE.\n"
        "No data enters either column. At n=2 a nominal 5% test rejects ~40% of the time."
    )
    shown = [2, 3, 4, 5, 6, 8, 10, 15, 19, 20, 30, 44, 60, 100, 150, 200]
    view = curve[curve["n_patients"].isin(shown)][
        [
            "n_patients",
            "width_ratio_vs_t",
            "closed_form_false_positive_rate",
            "inflation_factor",
            "verdict",
            "simulated_false_positive_rate_median",
            "abs_diff_vs_median",
            "n_simulation_cells",
        ]
    ]
    print(view.to_string(index=False))

    checked = curve[curve["n_simulation_cells"] > 0]
    if len(checked):
        print(
            f"\nCROSS-CHECK vs {SIMULATION_REFERENCE} at n = "
            f"{sorted(int(n) for n in checked['n_patients'])}:"
        )
        print(
            f"  max |closed - simulated median| = "
            f"{checked['abs_diff_vs_median'].max():.4f} "
            f"({100 * checked['abs_diff_vs_median'].max():.2f} pp), "
            f"tolerance {CLOSED_FORM_AGREEMENT_TOLERANCE}"
        )
        print(
            f"  max |closed - worst single cell| = "
            f"{checked['abs_diff_worst_cell'].max():.4f} "
            f"({100 * checked['abs_diff_worst_cell'].max():.2f} pp) — the closed form "
            f"is a floor; skew pushes individual cells above it."
        )

    calibrated = curve.loc[curve["verdict"] == "CALIBRATED", "n_patients"]
    if len(calibrated):
        print(
            f"\n  first n at which the percentile bootstrap is CALIBRATED by the "
            f"project's own tolerance: n = {int(calibrated.min())}"
        )

    path = write_versioned_table(
        curve,
        "interval_width_curve",
        seed=args.seed,
        notes=(
            "width_ratio(n) and expected_false_positive_rate(n) over n=2..200, "
            "cross-checked against the committed simulation. CLOSED FORM: no data "
            "enters either column."
        ),
        extra_meta={
            "n_min": int(args.n_min),
            "n_max": int(args.n_max),
            "alpha": 0.05,
            "simulation_reference": SIMULATION_REFERENCE,
            "cross_checked": simulated is not None,
            "agreement_tolerance": CLOSED_FORM_AGREEMENT_TOLERANCE,
            "guards_passed": (
                ["check_closed_form_tracks_simulation"] if simulated is not None else []
            ),
            "SYNTHETIC": SYNTHETIC_NOTE,
        },
        allow_dirty=args.allow_dirty,
    )
    print(f"\nwrote {path}")
    return 0


COMMANDS = {
    "bench": cmd_bench,
    "sweep": cmd_sweep,
    "misspecification": cmd_misspecification,
    "interval-curve": cmd_interval_curve,
}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "command", nargs="?", default="bench", choices=[*COMMANDS, "all"]
    )
    parser.add_argument("--n-replicates", type=int, default=200)
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--n-min", type=int, default=CURVE_N_MIN)
    parser.add_argument("--n-max", type=int, default=CURVE_N_MAX)
    parser.add_argument(
        "--allow-dirty",
        action="store_true",
        help=(
            "write a versioned table from a dirty tree. For scratch runs only: "
            "the recorded sha would not reproduce the table."
        ),
    )
    args = parser.parse_args(argv)

    pd.set_option("display.width", 220)
    pd.set_option("display.max_rows", 400)

    names = list(COMMANDS) if args.command == "all" else [args.command]
    for name in names:
        if len(names) > 1:
            print(f"\n\n{'#' * 96}\n### {name}\n{'#' * 96}")
        rc = COMMANDS[name](args)
        if rc:
            return rc
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
