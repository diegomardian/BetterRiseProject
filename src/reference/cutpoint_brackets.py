"""Calibration crossings as brackets, on the scale the calibration measured.

    from src.reference.cutpoint_brackets import brackets_by_seed, summarise_brackets

WHY A BRACKET AND NOT A THRESHOLD. ``calibrate_cutpoints`` returns the smallest
*bin median* mature-cell count whose bin meets the pre-registered criteria. That
number is a point on a coarse grid: the crossing it reports is somewhere between
the last bin that fails and the first bin that passes, and the point estimate
does not say where. The decision record for the audit
(``DECISION_2026-09-11_scope_and_pivot.md`` §4, week 2) asks for the crossing as
a **bracket**, precisely so that the resolution of the measurement travels with
the number. This module computes that bracket from the committed per-bin table;
it does not re-run the sweep and it does not interpolate.

THE SECOND CRITERION IS NOT THE FIRST ONE. ``wide`` and ``ok`` are different
questions and this module keeps them apart:

``wide`` — coverage >= ``coverage_target``. The interval contains the truth.
``ok``   — coverage >= target **and** discrimination >= ``discrimination_target``.
           The interval also excludes zero at the detectable effect.

``calibrate_cutpoints`` returns **no** report when no bin is ``ok``, so on a pool
where coverage is fine but discrimination never clears, the committed cutpoints
table records ``wide = None`` even though a coverage crossing exists. That is
not a property of the data; it is the point rule declining to emit half an
answer. On the project's ``pooled`` draw pool both Lee cohorts are exactly this
case. Reporting the two criteria separately recovers a crossing that the point
rule discards, and the recovered crossing is a real one.

``None`` IS NOT A ZERO. A pool with no bin meeting a criterion yields
``status="not_identifiable"`` and a null bracket. A criterion met by the first
bin yields ``status="below_first_bin"`` and a null *lower* endpoint, because the
crossing is at or below the smallest count tested rather than absent. The two
nulls mean different things and neither is a number.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from src.harness.calibration import PREREGISTERED, CalibrationCriteria

#: The two questions a calibration crossing can answer. Kept as a tuple so the
#: order is fixed and a table always carries both.
CRITERIA: tuple[str, ...] = ("wide", "ok")

#: Which per-bin verdicts satisfy each criterion. The mapping is read off
#: ``calibration._verdict``: ``ok`` is both targets met, ``wide_interval`` is
#: coverage met with discrimination below target, ``not_estimable`` and
#: ``all_abstained`` meet neither. A verdict not named here is an error, not a
#: silent pass — a new verdict string must be classified before it can count.
_MEETS: dict[str, frozenset[str]] = {
    "wide": frozenset({"ok", "wide_interval"}),
    "ok": frozenset({"ok"}),
}


class BracketError(ValueError):
    """A bins table that cannot be bracketed without guessing."""


def criterion_met(verdict: str, criterion: str) -> bool:
    """Whether one bin's verdict satisfies ``criterion``.

    ``all_abstained`` meets nothing, and this is the invariant-1 case at bin
    level: every replicate in that bin declined to answer, so coverage and
    discrimination are NaN. ``NaN >= target`` is False, which is the behaviour we
    want, but relying on that coercion here would make the rule depend on the
    comparison rather than on the verdict. It is explicit.
    """
    if criterion not in _MEETS:
        raise BracketError(
            f"unknown criterion {criterion!r}; known: {list(_MEETS)}"
        )
    if verdict == "all_abstained":
        return False
    known = set().union(*_MEETS.values()) | {"not_estimable", "all_abstained"}
    if verdict not in known:
        raise BracketError(
            f"unclassified verdict {verdict!r}. A new verdict string must be "
            f"assigned to a criterion (or explicitly excluded) before it can be "
            f"bracketed; treating it as a failure by default is how a real "
            f"crossing gets hidden."
        )
    return verdict in _MEETS[criterion]


def crossing_bracket(
    bins: pd.DataFrame,
    criterion: str,
    criteria: CalibrationCriteria = PREREGISTERED,
) -> dict[str, object]:
    """Bracket one (pool, grid, seed) crossing for one criterion.

    ``bins`` is the per-bin frame for a single sweep: one row per bin, carrying
    ``n_cells_mature`` and ``verdict``. Any other columns are ignored. Rows are
    sorted by ``n_cells_mature`` here rather than trusted to arrive sorted.

    Returns a flat dict with ``status``, the two endpoints and enough context to
    audit the bracket. The endpoints are the **bin medians** the calibration
    itself used; no interpolation to integer counts is attempted, because the
    calibration never measured at that resolution.
    """
    missing = {"n_cells_mature", "verdict"} - set(bins.columns)
    if missing:
        raise BracketError(
            f"bins table is missing {sorted(missing)}; it is not a "
            f"coverage_and_discrimination table."
        )
    if bins.empty:
        raise BracketError("bins table is empty; there is no crossing to bracket")

    ordered = bins.sort_values("n_cells_mature", kind="stable").reset_index(drop=True)
    meets = [criterion_met(v, criterion) for v in ordered["verdict"]]
    counts = ordered["n_cells_mature"].astype(float)

    out: dict[str, object] = {
        "criterion": criterion,
        "n_bins": int(len(ordered)),
        "n_bins_meeting": int(sum(meets)),
        "coverage_target": float(criteria.coverage_target),
        "discrimination_target": float(criteria.discrimination_target),
        "detectable_shift": float(criteria.detectable_shift),
    }

    if not any(meets):
        out |= {
            "status": "not_identifiable",
            "lower_n_cells_mature": np.nan,
            "upper_n_cells_mature": np.nan,
            "point_cutpoint": pd.NA,
            "first_meeting_bin": pd.NA,
            "max_coverage": float(ordered["coverage"].max())
            if "coverage" in ordered
            else np.nan,
            "max_discrimination": float(ordered["discrimination"].max())
            if "discrimination" in ordered
            else np.nan,
        }
        return out

    first = meets.index(True)
    upper = float(counts.iloc[first])
    out |= {
        "first_meeting_bin": int(first),
        "upper_n_cells_mature": upper,
        # The point estimate calibrate_cutpoints would report, so a reader can
        # see the threshold this bracket replaces.
        "point_cutpoint": int(upper),
    }
    if first == 0:
        out |= {
            "status": "below_first_bin",
            "lower_n_cells_mature": np.nan,
        }
    else:
        out |= {
            "status": "bracketed",
            "lower_n_cells_mature": float(counts.iloc[first - 1]),
        }
    if "coverage" in ordered:
        out["max_coverage"] = float(ordered["coverage"].max())
    if "discrimination" in ordered:
        out["max_discrimination"] = float(ordered["discrimination"].max())
    return out


def brackets_by_seed(
    bins: pd.DataFrame,
    *,
    criteria: CalibrationCriteria = PREREGISTERED,
    group_columns: tuple[str, ...] = ("cohort", "pool", "grid", "seed"),
) -> pd.DataFrame:
    """One row per (group, criterion): the crossing bracket for that sweep.

    ``bins`` carries the group columns and the per-bin ``n_cells_mature`` /
    ``verdict``. Every group is bracketed for both criteria, so the two answers
    stay beside each other and one can never be read as the other.
    """
    missing = set(group_columns) - set(bins.columns)
    if missing:
        raise BracketError(f"bins table is missing group columns {sorted(missing)}")

    rows: list[dict] = []
    for key, frame in bins.groupby(list(group_columns), observed=True, sort=True):
        key_tuple = key if isinstance(key, tuple) else (key,)
        base = dict(zip(group_columns, key_tuple, strict=True))
        for criterion in CRITERIA:
            rows.append(base | crossing_bracket(frame, criterion, criteria))
    return pd.DataFrame(rows)


def summarise_brackets(
    per_seed: pd.DataFrame,
    *,
    group_columns: tuple[str, ...] = ("cohort", "pool", "grid", "criterion"),
) -> pd.DataFrame:
    """Monte Carlo uncertainty of each bracket across seeds.

    The per-seed brackets are the measurement; this collapses them to one row
    per group and reports **how many seeds could not bracket at all** as a real
    count rather than folding it into a threshold. ``n_not_identifiable`` is the
    headline for a pooled draw: it is not a wide cutpoint, it is no cutpoint.
    """
    # Booleans up front so the aggregation below is a plain column reduction
    # rather than a lambda that would have to know about NaNs.
    frame = per_seed.assign(
        _has_upper=per_seed["upper_n_cells_mature"].notna(),
        _below_first=per_seed["status"].eq("below_first_bin"),
    )
    summary = (
        frame.groupby(list(group_columns), observed=True, sort=True)
        .agg(
            n_seeds=("seed", "nunique"),
            n_with_upper=("_has_upper", "sum"),
            n_below_first_bin=("_below_first", "sum"),
            lower_min=("lower_n_cells_mature", "min"),
            lower_max=("lower_n_cells_mature", "max"),
            upper_min=("upper_n_cells_mature", "min"),
            upper_max=("upper_n_cells_mature", "max"),
            point_min=("point_cutpoint", "min"),
            point_max=("point_cutpoint", "max"),
        )
        .reset_index()
    )
    summary["n_not_identifiable"] = summary["n_seeds"] - summary["n_with_upper"]
    # Stable means every seed that produced a bracket produced the same one; a
    # pool where no seed produced one is not "stable", it is absent, so the
    # comparison is guarded on there being at least one bracket.
    both_sides_fixed = summary["lower_min"].eq(summary["lower_max"]) | summary[
        "lower_min"
    ].isna()
    summary["bracket_stable"] = (
        summary["n_with_upper"].gt(0)
        & summary["upper_min"].eq(summary["upper_max"])
        & both_sides_fixed
    )
    ordered = [*group_columns, "n_seeds", "n_with_upper", "n_not_identifiable",
               "n_below_first_bin", "lower_min", "lower_max", "upper_min",
               "upper_max", "point_min", "point_max", "bracket_stable"]
    return summary[ordered]
