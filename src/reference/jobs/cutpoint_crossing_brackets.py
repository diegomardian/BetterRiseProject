"""The cutpoint calibration crossing, reported as a bracket. Week 2, audit.

    python -m src.reference.jobs.cutpoint_crossing_brackets
    python -m src.reference.jobs.cutpoint_crossing_brackets \
        --smc-bins results/2026-09-12_c49ea0f/cutpoint_dense_grid_bins_smc_r200.parquet \
        --kul3-bins results/2026-09-12_c49ea0f/cutpoint_dense_grid_bins_kul3_r200.parquet

THE DELIVERABLE. ``DECISION_2026-09-11_scope_and_pivot.md`` §4, week 2, asks for
``cutpoint_crossing_brackets.parquet``: report calibration crossings as
**brackets, not precise thresholds**, and keep SMC and KUL3 separate. This turns
the committed per-bin calibration into that table. It reads artifacts already in
git, re-runs no sweep, and keeps the two cohorts as separate rows rather than
pooling them.

WHAT THE BRACKET ADDS. The committed ``cutpoint_dense_grid`` table reports a
single number per pool and seed -- the smallest bin whose median mature count met
the criteria. The crossing is somewhere between that bin and the one below it,
and the point estimate cannot say where. The bracket is that interval. It is
also two brackets, not one: ``wide`` (coverage only) and ``ok`` (coverage and
discrimination) are separate questions and are reported separately. On the
``pooled`` draw pool ``ok`` never holds while ``wide`` does, so the point rule --
which requires ``ok`` -- writes ``wide = None``. The bracket table recovers it.

Read the summary table for the verdict. ``n_not_identifiable`` counts seeds that
could not bracket at all; it is the corrected answer for the pooled draw and it
is a count, not a wide threshold.
"""

from __future__ import annotations

import argparse
import logging
import sys
from collections.abc import Sequence
from pathlib import Path

import pandas as pd

from src.common.io import write_versioned_table
from src.common.paths import RESULTS_DIR
from src.common.provenance import DEFAULT_SEED
from src.harness.calibration import PREREGISTERED
from src.reference.cutpoint_brackets import brackets_by_seed, summarise_brackets

log = logging.getLogger(__name__)

#: The two Lee cohorts, kept separate end to end (invariant 4 is about pooling
#: datasets; a table can carry both as long as no row combines them).
COHORTS: tuple[str, ...] = ("smc", "kul3")


def newest_bins(cohort: str) -> Path | None:
    """The most recent committed dense-grid bins table for one cohort.

    Sorted by path, not by mtime: the directory name is ``{date}_{sha7}`` so
    lexicographic order is chronological, and mtime is not reproducible across a
    checkout. See HANDOFF §4 for the ``paper/*/_tables.py`` version that resolves
    by mtime and should not be copied.
    """
    matches = sorted(RESULTS_DIR.glob(f"*/cutpoint_dense_grid_bins_{cohort}_r*.parquet"))
    return matches[-1] if matches else None


def load_bins(paths: dict[str, Path]) -> pd.DataFrame:
    """Read one bins table per cohort and tag each row with its cohort."""
    frames = []
    for cohort, path in paths.items():
        frame = pd.read_parquet(path)
        frames.append(frame.assign(cohort=cohort))
    bins = pd.concat(frames, ignore_index=True)
    expected = {"pool", "grid", "seed", "n_cells_mature", "verdict"}
    missing = expected - set(bins.columns)
    if missing:
        raise SystemExit(
            f"bins table(s) missing {sorted(missing)}; these are not "
            f"coverage_and_discrimination tables. Pass --smc-bins/--kul3-bins."
        )
    return bins


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--smc-bins", type=Path, default=None)
    parser.add_argument("--kul3-bins", type=Path, default=None)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    inputs: dict[str, Path] = {}
    for cohort, override in (("smc", args.smc_bins), ("kul3", args.kul3_bins)):
        path = override or newest_bins(cohort)
        if path is None:
            raise SystemExit(
                f"no results/*/cutpoint_dense_grid_bins_{cohort}_r*.parquet. "
                f"Run src.reference.jobs.cutpoint_dense_grid --cohort {cohort} "
                f"first, or pass --{cohort}-bins."
            )
        inputs[cohort] = path
        log.info("cohort %s <- %s", cohort, path.relative_to(RESULTS_DIR.parent)
                 if path.is_relative_to(RESULTS_DIR.parent) else path)

    bins = load_bins(inputs)
    per_seed = brackets_by_seed(bins)
    summary = summarise_brackets(per_seed)

    log.info("\ncrossing brackets, by cohort / pool / criterion:\n%s",
             summary.to_string(index=False))

    # The finding the point table cannot state: where the point rule returned no
    # cutpoint, is there a coverage crossing anyway?
    hidden_by_ok = summary[
        (summary["criterion"] == "ok") & (summary["n_with_upper"] == 0)
    ][["cohort", "pool"]]
    if not hidden_by_ok.empty:
        log.info(
            "\n`ok` not identifiable on %s; a `wide` crossing still exists there "
            "(the point rule, which requires `ok`, reports wide=None for these).",
            ", ".join(f"{r.cohort}/{r.pool}" for r in hidden_by_ok.itertuples()),
        )

    meta = {
        "purpose": (
            "DECISION_2026-09-11_scope_and_pivot.md §4 week 2 -- report "
            "calibration crossings as brackets, not precise thresholds. Derived "
            "from the committed per-bin calibration; no sweep is re-run."
        ),
        "source_tables": {
            cohort: str(path.relative_to(RESULTS_DIR.parent))
            if path.is_relative_to(RESULTS_DIR.parent) else str(path)
            for cohort, path in inputs.items()
        },
        "criteria": {
            "wide": "coverage >= coverage_target",
            "ok": "coverage >= coverage_target AND discrimination >= discrimination_target",
            "coverage_target": PREREGISTERED.coverage_target,
            "discrimination_target": PREREGISTERED.discrimination_target,
            "detectable_shift": PREREGISTERED.detectable_shift,
        },
        "bracket_definition": (
            "lower = n_cells_mature median of the highest bin that fails the "
            "criterion below the crossing; upper = median of the first bin that "
            "meets it. Endpoints are bin medians, the calibration's own "
            "resolution; no interpolation to integer counts is attempted. "
            "status='not_identifiable' means no bin met the criterion; "
            "status='below_first_bin' means the smallest tested bin already met "
            "it, so the crossing is at or below the grid rather than absent."
        ),
        "cohorts_kept_separate": list(COHORTS),
        "n_input_bins": int(len(bins)),
    }
    for frame, name in (
        (per_seed, "cutpoint_crossing_brackets"),
        (summary, "cutpoint_crossing_brackets_summary"),
    ):
        path = write_versioned_table(
            frame, name, seed=DEFAULT_SEED, results_dir=args.results_dir,
            extra_meta=meta, allow_dirty=args.allow_dirty,
        )
        log.info("wrote %s (%d rows)", path, len(frame))
    return 0


if __name__ == "__main__":
    sys.exit(main())
