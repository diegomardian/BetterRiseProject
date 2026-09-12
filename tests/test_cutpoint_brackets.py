"""Forcing inputs for the cutpoint crossing bracket.

Every test here exists to make a specific claim in ``cutpoint_brackets`` fail if
it stops holding. The module's whole job is to keep three distinctions that are
easy to collapse: a crossing that is absent is not one that is wide, the first
bin already passing is not the same as no bin passing, and ``wide`` and ``ok``
are different questions even though the point rule emits only one number.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.reference.cutpoint_brackets import (
    BracketError,
    brackets_by_seed,
    criterion_met,
    crossing_bracket,
    summarise_brackets,
)


def _bins(rows: list[tuple[float, str]], **extra) -> pd.DataFrame:
    """Rows of (n_cells_mature, verdict), with plausible rates filled in."""
    frame = pd.DataFrame(rows, columns=["n_cells_mature", "verdict"])
    # A rate that is consistent with the verdict, so the max-rate context fields
    # are never the thing under test.
    coverage = {"ok": 0.96, "wide_interval": 0.93, "not_estimable": 0.75,
                "all_abstained": np.nan}
    discrimination = {"ok": 0.85, "wide_interval": 0.55,
                      "not_estimable": 0.55, "all_abstained": np.nan}
    frame["coverage"] = [coverage[v] for v in frame["verdict"]]
    frame["discrimination"] = [discrimination[v] for v in frame["verdict"]]
    for key, value in extra.items():
        frame[key] = value
    return frame


def test_criterion_met_by_verdict():
    assert criterion_met("ok", "wide")
    assert criterion_met("wide_interval", "wide")
    assert not criterion_met("not_estimable", "wide")
    # `ok` is the stricter question: coverage alone is not enough.
    assert not criterion_met("wide_interval", "ok")
    assert criterion_met("ok", "ok")


def test_all_abstained_meets_neither_criterion():
    """A bin nobody answered is not a bin that passed.

    This is the invariant-1 case at bin level: coverage is NaN and every
    comparison against it is False. The module must return False deliberately
    rather than inherit the coercion.
    """
    assert not criterion_met("all_abstained", "wide")
    assert not criterion_met("all_abstained", "ok")


def test_unknown_verdict_is_refused_not_treated_as_failure():
    with pytest.raises(BracketError, match="unclassified verdict"):
        criterion_met("sort_of_ok", "wide")


def test_a_clean_crossing_is_bracketed_between_two_bins():
    frame = _bins([
        (10.0, "not_estimable"),
        (20.0, "not_estimable"),
        (30.0, "wide_interval"),
        (40.0, "ok"),
    ])
    wide = crossing_bracket(frame, "wide")
    assert wide["status"] == "bracketed"
    assert wide["lower_n_cells_mature"] == 20.0
    assert wide["upper_n_cells_mature"] == 30.0
    assert wide["point_cutpoint"] == 30

    # The stricter criterion crosses later, and its bracket uses the wide bin as
    # its failing lower endpoint. Two questions, two brackets.
    ok = crossing_bracket(frame, "ok")
    assert ok["status"] == "bracketed"
    assert ok["lower_n_cells_mature"] == 30.0
    assert ok["upper_n_cells_mature"] == 40.0
    assert ok["point_cutpoint"] == 40


def test_rows_out_of_order_are_sorted_before_bracketing():
    frame = _bins([
        (40.0, "ok"),
        (10.0, "not_estimable"),
        (30.0, "wide_interval"),
        (20.0, "not_estimable"),
    ])
    ok = crossing_bracket(frame, "ok")
    assert (ok["lower_n_cells_mature"], ok["upper_n_cells_mature"]) == (30.0, 40.0)


def test_not_identifiable_is_null_not_a_wide_threshold():
    frame = _bins([(10.0, "not_estimable"), (20.0, "not_estimable")])
    out = crossing_bracket(frame, "ok")
    assert out["status"] == "not_identifiable"
    assert np.isnan(out["lower_n_cells_mature"])
    assert np.isnan(out["upper_n_cells_mature"])
    # Not 0, and not the largest count. A number here would be invariant 1's
    # failure mode: an unestimable crossing reported as a value.
    assert pd.isna(out["point_cutpoint"])
    assert out["n_bins_meeting"] == 0


def test_crossing_below_the_first_bin_has_no_lower_endpoint():
    frame = _bins([(5.0, "ok"), (10.0, "ok")])
    out = crossing_bracket(frame, "wide")
    assert out["status"] == "below_first_bin"
    assert np.isnan(out["lower_n_cells_mature"])
    assert out["upper_n_cells_mature"] == 5.0


def test_abstained_bin_is_a_failing_lower_endpoint():
    """A bin that abstained everywhere is the last failure below the crossing.

    It must sit on the lower side of the bracket, not be skipped and not be read
    as passing. Drawn from the real pattern: the smallest bins abstain, the
    crossing is above them.
    """
    frame = _bins([(2.5, "all_abstained"), (10.0, "wide_interval")])
    out = crossing_bracket(frame, "wide")
    assert out["status"] == "bracketed"
    assert out["lower_n_cells_mature"] == 2.5
    assert out["upper_n_cells_mature"] == 10.0


def test_missing_columns_are_refused():
    with pytest.raises(BracketError, match="missing"):
        crossing_bracket(pd.DataFrame({"n_cells_mature": [1.0]}), "wide")


def _per_seed(cohort: str, pool: str, criterion: str, seed: int,
              status: str, lower, upper) -> dict:
    return {
        "cohort": cohort, "pool": pool, "grid": "dense", "seed": seed,
        "criterion": criterion, "status": status,
        "lower_n_cells_mature": lower, "upper_n_cells_mature": upper,
        "point_cutpoint": int(upper) if upper == upper else pd.NA,
    }


def test_summary_separates_identifiable_from_absent():
    rows = [
        _per_seed("smc", "pooled", "wide", 1, "bracketed", 42.5, 70.0),
        _per_seed("smc", "pooled", "wide", 2, "bracketed", 42.5, 70.0),
        _per_seed("smc", "pooled", "ok", 1, "not_identifiable", np.nan, np.nan),
        _per_seed("smc", "pooled", "ok", 2, "not_identifiable", np.nan, np.nan),
    ]
    summary = summarise_brackets(pd.DataFrame(rows))
    wide = summary[(summary.criterion == "wide")].iloc[0]
    ok = summary[(summary.criterion == "ok")].iloc[0]
    assert wide["n_with_upper"] == 2
    assert wide["n_not_identifiable"] == 0
    assert wide["bracket_stable"]
    assert ok["n_with_upper"] == 0
    assert ok["n_not_identifiable"] == 2
    # No bracket is not "stable". It is absent, and the flag says so.
    assert not ok["bracket_stable"]


def test_summary_does_not_count_not_identifiable_as_a_wide_cutpoint():
    """The point rule's blind spot, as a table row.

    `calibrate_cutpoints` requires `ok`, so a pool with a coverage crossing and
    no `ok` crossing gets `wide=None` from the point table. If the bracket
    summary folded `not_identifiable` into a threshold, it would repeat that
    mistake in the other direction.
    """
    rows = [
        _per_seed("kul3", "pooled", "wide", s, "bracketed", 42.5, 70.0)
        for s in (1, 2)
    ] + [
        _per_seed("kul3", "pooled", "ok", s, "not_identifiable", np.nan, np.nan)
        for s in (1, 2)
    ]
    summary = summarise_brackets(pd.DataFrame(rows))
    assert summary[summary.criterion == "ok"]["n_with_upper"].eq(0).all()
    assert summary[summary.criterion == "wide"]["upper_max"].eq(70.0).all()


def test_brackets_by_seed_emits_both_criteria_per_group():
    bins = pd.concat([
        _bins(
            [(10.0, "not_estimable"), (20.0, "wide_interval")],
            cohort="smc", pool="pooled", grid="dense", seed=1,
        ),
        _bins(
            [(10.0, "not_estimable"), (20.0, "not_estimable")],
            cohort="smc", pool="pooled", grid="dense", seed=2,
        ),
    ], ignore_index=True)
    out = brackets_by_seed(bins)
    assert set(out["criterion"]) == {"wide", "ok"}
    assert len(out) == 4
    seed1_wide = out[(out.seed == 1) & (out.criterion == "wide")].iloc[0]
    assert seed1_wide["status"] == "bracketed"
    assert (seed1_wide["lower_n_cells_mature"], seed1_wide["upper_n_cells_mature"]) == (
        10.0, 20.0,
    )
