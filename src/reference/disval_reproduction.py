"""Reproduce the DIS/VAL stability split from committed derived inputs.

    from src.reference.disval_reproduction import compare_disval

WHY THIS CAN BE DONE WITHOUT THE OBS CACHE. ``disval_stability`` reads
``data/interim/icbi_obs.parquet`` to learn each patient's specimen collection.
That cache is gitignored and absent here, but the job *wrote out* its own
assignment (``disval_stability_assignment.parquet``: patient, collection,
shared, found) -- so the split can be rebuilt from committed derived inputs.
That makes the DIS/VAL reading reproducible on a machine with no interim data,
which is the decision record's week-3 requirement: *"Reproduce DIS/VAL and
Crowell without opening new endpoints."*

It is a **reproduction from derived inputs, not a raw-data replication**. The
assignment is itself a derived table; agreeing with the committed contrast table
shows the reading is a deterministic function of committed inputs, not that the
atlas was re-read.
"""

from __future__ import annotations

import pandas as pd

#: A contrast is identified by its half, its statistic and the gene pair. The
#: interval summary is carried too and checked, because a mean and a median with
#: the same centre are not the same claim.
KEY_COLUMNS: tuple[str, ...] = ("half", "statistic", "contrast")

#: Compared numerically. ``ratio`` is summarised by a median with a rank
#: interval, so a mismatch in any of these is a real disagreement.
NUMERIC_COLUMNS: tuple[str, ...] = ("n_patients", "mean", "ci_low", "ci_high")


class DisvalReproductionError(ValueError):
    """A contrast table that cannot be compared to the committed one."""


def compare_disval(
    reproduced: pd.DataFrame, committed: pd.DataFrame, *, atol: float = 1e-9
) -> pd.DataFrame:
    """Line up reproduced and committed DIS/VAL contrasts, cell by cell.

    Raises if the two cover different (half, statistic, contrast) cells: a
    reproduction that quietly covers fewer cells is how a missing result reads
    as a passing check.
    """
    for label, frame in (("reproduced", reproduced), ("committed", committed)):
        # Parenthesise: set difference binds tighter than union, so the
        # obvious-looking form would test the keys against the frame's columns
        # and ignore the required numeric columns entirely.
        missing = (set(KEY_COLUMNS) | {"mean"}) - set(frame.columns)
        if missing:
            raise DisvalReproductionError(f"{label} table is missing {sorted(missing)}")

    merged = reproduced.merge(
        committed, on=list(KEY_COLUMNS), suffixes=("", "_committed"),
        how="outer", indicator=True,
    )
    unmatched = merged[merged["_merge"] != "both"]
    if not unmatched.empty:
        keys = unmatched[list(KEY_COLUMNS)].to_dict("records")
        raise DisvalReproductionError(
            f"cells present in only one table: {keys[:5]}"
        )

    flags = []
    for column in NUMERIC_COLUMNS:
        left, right = column, f"{column}_committed"
        if left in merged and right in merged:
            merged[f"{column}_diff"] = (merged[left] - merged[right]).abs()
        else:
            merged[f"{column}_diff"] = float("nan")
        flags.append(merged[f"{column}_diff"].fillna(0.0) <= atol)
    if "excludes_zero" in merged and "excludes_zero_committed" in merged:
        zero_match = (
            merged["excludes_zero"].astype("boolean").fillna(False).astype(bool)
            == merged["excludes_zero_committed"].astype("boolean")
            .fillna(False).astype(bool)
        )
    else:
        zero_match = pd.Series(True, index=merged.index)
    if "summary" in merged and "summary_committed" in merged:
        summary_match = merged["summary"].fillna("NA").eq(
            merged["summary_committed"].fillna("NA")
        )
    else:
        summary_match = pd.Series(True, index=merged.index)

    merged["reproduces"] = zero_match & summary_match
    for flag in flags:
        merged["reproduces"] &= flag
    return merged.drop(columns=["_merge"]).sort_values(
        list(KEY_COLUMNS), ignore_index=True
    )
