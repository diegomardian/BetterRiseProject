"""Reproduce the Crowell summary from committed derived inputs, and check it.

    from src.reference.crowell_reproduction import compare_summaries

WHY A SEPARATE ARTIFACT. The Crowell h5ads were deleted after the run, per the
pre-registration's disk protocol, and re-downloadable from recorded md5s. That
means the summary can only be *reproduced from derived inputs* -- the committed
per-section ``crowell_feasibility`` tables -- not replicated from raw data. The
audit decision record names this table explicitly and requires that distinction
be labelled rather than elided (``DECISION_2026-09-11_scope_and_pivot.md`` §4,
week 3: *"Label as reproduction from derived inputs, not raw-data replication."*)

This module does not re-implement the DiD. It calls the same
``per_block_did``/``aggregate`` the primary job uses and reports whether the
result equals the committed summary, gene by gene. A reproduction that skipped
the comparison would be a re-run, and a re-run that was never checked is the
kind of green light this project is about.
"""

from __future__ import annotations

import pandas as pd

#: Absolute tolerance for a reproduction match. The aggregate is a mean and a
#: Student-t interval over at most seven block values, so an exact equality is
#: expected; the tolerance only absorbs float ordering, never a real difference.
MATCH_ATOL: float = 1e-9


class ReproductionError(ValueError):
    """A summary frame that cannot be compared to the committed one."""


def compare_summaries(
    reproduced: pd.DataFrame, committed: pd.DataFrame, *,
    atol: float = MATCH_ATOL,
) -> pd.DataFrame:
    """Line up the reproduced and committed per-gene summaries, gene by gene.

    Returns one row per gene with both values and an explicit ``reproduces``
    flag. A gene present in one and not the other is a failure, not a missing
    value: the reproduction must account for the same genes the committed
    summary does.
    """
    for label, frame in (("reproduced", reproduced), ("committed", committed)):
        missing = {"gene", "mean_did", "n_blocks"} - set(frame.columns)
        if missing:
            raise ReproductionError(f"{label} summary is missing {sorted(missing)}")

    merged = reproduced.merge(
        committed,
        on="gene", how="outer", suffixes=("", "_committed"), indicator=True,
    )
    unmatched = merged[merged["_merge"] != "both"]
    if not unmatched.empty:
        raise ReproductionError(
            f"genes present in only one summary: {sorted(unmatched['gene'])}. "
            f"A reproduction must cover the same genes as the committed table."
        )

    merged["mean_did_diff"] = (
        merged["mean_did"] - merged["mean_did_committed"]
    ).abs()
    for column in ("ci_low", "ci_high"):
        left, right = f"{column}", f"{column}_committed"
        if left in merged and right in merged:
            merged[f"{column}_diff"] = (merged[left] - merged[right]).abs()
        else:
            merged[f"{column}_diff"] = float("nan")

    n_blocks_match = merged["n_blocks"] == merged["n_blocks_committed"]
    interval_match = _interval_verdict(merged, "interval")
    zero_match = _interval_verdict(merged, "excludes_zero")
    merged["reproduces"] = (
        (merged["mean_did_diff"] <= atol)
        & (merged["ci_low_diff"].fillna(0.0) <= atol)
        & (merged["ci_high_diff"].fillna(0.0) <= atol)
        & n_blocks_match
        & interval_match
        & zero_match
    )
    return merged.drop(columns=["_merge"]).sort_values("gene", ignore_index=True)


def _interval_verdict(merged: pd.DataFrame, column: str) -> pd.Series:
    """Whether a categorical column matches, treating both-present-None as equal.

    ``interval`` is the string ``"refused"`` below MIN_STUDIES and None is never
    used; ``excludes_zero`` is a real boolean. A column absent on either side is
    ignored rather than failing the reproduction, because an older committed
    table may not carry it.
    """
    left, right = column, f"{column}_committed"
    if left not in merged or right not in merged:
        return pd.Series(True, index=merged.index)
    return merged[left].fillna("NA").astype(str).eq(
        merged[right].fillna("NA").astype(str)
    )
