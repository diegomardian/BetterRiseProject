"""The cutpoint calibration on a grid that can reach where adenoma lives.

    python -m src.reference.jobs.cutpoint_dense_grid --cohort smc
    python -m src.reference.jobs.cutpoint_dense_grid --cohort kul3 --replicates 200

WHY THIS EXISTS. ``docs/W2_HANDOFF_TO_NEXT_AGENT.md`` §5 item 6 is *"recalibrate
cutpoints on a denser 5-50 grid"*, listed as **not started**, and avenue A is
what makes it urgent: the adenoma ``best4`` rung has a median of **30 mature
cells per arm** and a range of 16-63, and the existing grids cannot resolve
anything between 0 and 20.

    committed  ->  0, 20, 40, 100, 200, 400, 800
    extended   ->  0, 20, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 150, 200, 400, 800

``ok`` is not the problem — it calibrates to 90-100 on SMC and 400-800 on KUL3,
far above anything adenoma reaches. **``wide`` is the problem**, because ``wide``
is the boundary between ``wide_interval`` (an intrinsic term is written) and
``not_estimable`` (it is ``None``, invariant 1), and at a median of 30 cells
that boundary decides most of the rung. The committed runs put ``wide`` at
40-45 on SMC and 65-100 on KUL3 against a provisional 20 — three values that
disagree, none of them measured on a grid with a point between 0 and 20.

THE GRID KEEPS THE HIGH POINTS, AND THAT IS NOT COSMETIC. A dense grid covering
only 5-50 would put ``ok``'s crossing outside the grid, ``calibrate_cutpoints``
would find no qualifying bin, and it would raise the G4 message -- reporting
non-identifiability that was really a grid too short to reach the answer. That
is defect #2 in ``docs/HANDOFF.md`` §3, *"calibration grid: grid cannot reach
the crossing"*, and re-introducing it while fixing its sibling would be
characteristic. So ``DENSE_FRACTIONS`` is the extended grid **plus** 5..50 by
fives, never a replacement.

OWNERSHIP. ``src/harness/`` is W2's under CONTRIBUTING §2 and **nothing in it is
edited here.** ``run_calibration_gap`` already takes ``grids`` as a parameter;
this passes a different one and reuses W2's sweep, criteria and cutpoint rule
unchanged. The number this produces is evidence for W2's decision, not a
substitute for it -- ``positivity.CUTPOINTS`` is theirs to move.
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
from src.harness.calibration import PREREGISTERED
from src.harness.calibration_gap import (
    EXTENDED_FRACTIONS,
    N_CELLS,
    SEEDS,
    counts_reachable,
    load_cohort_arrays,
    run_calibration_gap,
)
from src.harness.positivity import CUTPOINTS

log = logging.getLogger(__name__)

#: 5 to 50 in fives, as fractions of the fixed 2,000-cell draw. This is the
#: window `best4` actually occupies (16-63 mature cells per arm, median 30).
DENSE_COUNTS: tuple[int, ...] = (5, 10, 15, 20, 25, 30, 35, 40, 45, 50)

#: The extended grid PLUS the dense window. Never a replacement -- see the
#: module docstring on why a short grid manufactures a G4.
DENSE_FRACTIONS: tuple[float, ...] = tuple(
    sorted(set(EXTENDED_FRACTIONS) | {c / N_CELLS for c in DENSE_COUNTS},
           reverse=True)
)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort", choices=("smc", "kul3"), default="smc")
    parser.add_argument("--replicates", type=int, default=200)
    parser.add_argument("--seeds", type=int, default=8)
    parser.add_argument("--raw-dir", type=Path, default=None)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    log.info("grid resolution, in mature cells:")
    log.info("  extended : %s", counts_reachable(EXTENDED_FRACTIONS))
    log.info("  dense    : %s", counts_reachable(DENSE_FRACTIONS))
    log.info("  adenoma best4 occupies 16-63, median 30 -- the extended grid "
             "has NO point below 20")

    log.info("\nloading Lee/%s ...", args.cohort.upper())
    arrays = load_cohort_arrays(args.cohort, args.raw_dir)
    bins, cutpoints, recovery = run_calibration_gap(
        arrays["counts"], arrays["cell_type"], arrays["patient_id"],
        arrays["tissue"], arrays["genes"],
        seeds=SEEDS[: args.seeds], n_replicates=args.replicates,
        grids={"dense": DENSE_FRACTIONS}, criteria=PREREGISTERED,
    )

    log.info("\n%s\nCUTPOINTS ON THE DENSE GRID — %s\n%s",
             "=" * 72, args.cohort.upper(), "=" * 72)
    summary = cutpoints.groupby(["pool", "grid"]).agg(
        seeds=("seed", "nunique"),
        returned=("returned_a_cutpoint", "sum"),
        ok_median=("ok", "median"), ok_min=("ok", "min"), ok_max=("ok", "max"),
        wide_median=("wide", "median"), wide_min=("wide", "min"),
        wide_max=("wide", "max"),
        max_discrimination=("max_discrimination", "max"),
    )
    log.info("%s", summary.to_string())
    log.info("\n  in use (provisional): ok=%d wide=%d", CUTPOINTS.ok, CUTPOINTS.wide)
    log.info(
        "\n  READ `wide` HERE, NOT `ok`. `ok` calibrates far above anything the "
        "adenoma\n  cohort reaches, so it does not decide that reading. `wide` "
        "is the boundary\n  between an intrinsic term being written and being "
        "None, and at a median of\n  30 mature cells it decides most of the "
        "`best4` rung.")

    meta = {
        "purpose": (
            "W2_HANDOFF_TO_NEXT_AGENT.md §5 item 6 -- recalibrate cutpoints on "
            "a denser 5-50 grid, listed as not started. Avenue A is what makes "
            "it urgent: adenoma best4 runs 16-63 mature cells per arm and the "
            "existing grids have no point between 0 and 20."
        ),
        "cohort": arrays.get("study_id", args.cohort),
        "n_patients": arrays.get("n_patients"),
        "grid_counts_dense": counts_reachable(DENSE_FRACTIONS),
        "grid_counts_extended": counts_reachable(EXTENDED_FRACTIONS),
        "why_the_high_points_stay": (
            "a grid covering only 5-50 would put ok's crossing outside it, and "
            "calibrate_cutpoints would raise its G4 message for a grid too "
            "short rather than an estimator that cannot resolve the effect -- "
            "HANDOFF §3 defect 2, re-introduced while fixing its sibling"
        ),
        "criteria": {
            "detectable_shift": PREREGISTERED.detectable_shift,
            "coverage_target": PREREGISTERED.coverage_target,
            "discrimination_target": PREREGISTERED.discrimination_target,
        },
        "provisional_cutpoints": {"ok": CUTPOINTS.ok, "wide": CUTPOINTS.wide,
                                  "source": CUTPOINTS.source},
        "ownership": (
            "src/harness/ is W2's and nothing in it is edited here. "
            "run_calibration_gap already takes `grids`; this passes a different "
            "one and reuses W2's sweep, criteria and cutpoint rule unchanged. "
            "positivity.CUTPOINTS is W2's to move."
        ),
        "n_replicates": args.replicates,
        "seeds": list(SEEDS[: args.seeds]),
        "exploratory": False,
    }
    suffix = f"_{args.cohort}_r{args.replicates}"
    for frame, name in ((cutpoints, f"cutpoint_dense_grid_cutpoints{suffix}"),
                        (bins, f"cutpoint_dense_grid_bins{suffix}"),
                        (recovery, f"cutpoint_dense_grid_recovery{suffix}")):
        if isinstance(frame, pd.DataFrame) and frame.empty:
            continue
        log.info("wrote %s", write_versioned_table(
            frame, name, seed=DEFAULT_SEED, results_dir=args.results_dir,
            allow_dirty=args.allow_dirty, extra_meta=meta,
        ))
    return 0


if __name__ == "__main__":
    sys.exit(main())
