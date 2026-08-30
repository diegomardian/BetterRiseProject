#!/usr/bin/env python
"""Run the competitor benchmark and write its three tables.

    python -m submission.run_bench

Writes to ``submission/results/`` — deliberately NOT to the project's
``results/``, which is versioned by git sha and governed by the frozen schema.
Nothing here is a result about colorectal cancer; every number comes from
simulated cells whose truth is known in closed form, and the tables say so.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from src.common.provenance import provenance_record
from submission.bench import (
    BENCH_WORLDS,
    MEAN_NORMAL,
    N_CELLS,
    refusal_table,
    run_bench,
    sensitivity_where_estimable,
)

SEED = 20260829
OUT = Path(__file__).parent / "results"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-replicates", type=int, default=200)
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args()

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

    pd.set_option("display.width", 200)
    print("=" * 96)
    print("REFUSAL — where the intrinsic estimand DOES NOT EXIST (annihilated: 0 mature")
    print("tumour cells). false_confidence_rate counts numbers RETURNED, not numbers wrong.")
    print("=" * 96)
    print(refusals.to_string(index=False))
    print()
    print("=" * 96)
    print("SENSITIVITY — where it DOES exist and a real intrinsic effect is present.")
    print("The counterweight: refusing always scores 0 here.")
    print("=" * 96)
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
            "SYNTHETIC": (
                "Every number here comes from simulated cells with analytically known "
                "truth. Nothing in this directory is a result about colorectal cancer."
            ),
        }
    )
    (OUT / "bench.meta.json").write_text(json.dumps(meta, indent=2))
    print(f"\nwrote 3 tables + bench.meta.json to {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
