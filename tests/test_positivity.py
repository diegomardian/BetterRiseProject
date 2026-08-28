"""The positivity cutpoints, both arms. W2.

The intrinsic rule is W4's import and must not move by accident, so it is pinned
here rather than only exercised through the harness. The compositional rule is
decision #22 and its tests carry the one thing the decision record got wrong.
"""

from __future__ import annotations

import pandas as pd
import pytest

from src.harness.positivity import (
    COMPOSITIONAL_CUTPOINTS,
    CUTPOINTS,
    classify_compositional_estimability,
    classify_counts_frame,
    classify_estimability,
    estimability_verdicts,
)

# ---------------------------------------------------------------------------
# The rule itself
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("n", "expected"),
    [
        (0, "not_estimable"),
        (19, "not_estimable"),
        (20, "wide_interval"),
        (49, "wide_interval"),
        (50, "ok"),
        (5000, "ok"),
    ],
)
def test_compositional_bands_at_their_boundaries(n, expected):
    assert classify_compositional_estimability(n) == expected


def test_compositional_rejects_a_negative_count():
    with pytest.raises(ValueError, match="negative"):
        classify_compositional_estimability(-1)


def test_the_committed_numbers_are_the_committed_numbers():
    """Decision #22 fixed 50/20 in public on issue #36, before it was applied.

    A recalibration is allowed; changing them without one is the thing
    pre-commitment exists to prevent, so it should fail a test rather than pass
    review.
    """
    assert (COMPOSITIONAL_CUTPOINTS.ok, COMPOSITIONAL_CUTPOINTS.wide) == (50, 20)
    assert "#22" in COMPOSITIONAL_CUTPOINTS.source


def test_the_intrinsic_rule_is_unchanged_by_the_second_arm():
    """W4 imports ``classify_estimability``; adding an arm must not move it."""
    assert (CUTPOINTS.ok, CUTPOINTS.wide) == (50, 20)
    assert classify_estimability(50) == "ok"
    assert classify_estimability(20) == "wide_interval"
    assert classify_estimability(19) == "not_estimable"


# ---------------------------------------------------------------------------
# The two arms together
# ---------------------------------------------------------------------------


def test_verdicts_report_both_arms_separately():
    """Not one folded value: the two say different things about the same row."""
    verdicts = estimability_verdicts(n_cells_mature=3, n_cells_resolved=800)
    assert verdicts["estimability"] == "not_estimable"
    assert verdicts["compositional_estimability"] == "ok"


def test_mature_cells_exceeding_resolved_cells_raises():
    """The guard that can actually fire — a caller passing the wrong column.

    ``n_cells_epithelial`` as the denominator, or two rungs mixed, both land
    here. Compare with the guard family this repo keeps finding: a check that
    cannot fail is worse than no check.
    """
    with pytest.raises(ValueError, match="subset of"):
        estimability_verdicts(n_cells_mature=60, n_cells_resolved=40)


@pytest.mark.parametrize("n_mature", [50, 51, 120, 3547])
@pytest.mark.parametrize("extra_unresolved", [0, 1, 40, 900])
def test_an_ok_intrinsic_arm_forces_an_ok_compositional_arm(n_mature, extra_unresolved):
    """The structural fact, and the correction to decision #22's own caveat.

    Mature cells are a subset of resolved cells, so ``n_cells_mature >= 50``
    implies ``n_cells_resolved >= 50``. The compositional gate can only ever bind
    where the intrinsic gate has already flagged the row.

    The issue-#36 comment that recorded the decision claimed the opposite — that
    it binds "on rows where the intrinsic arm is already ok ... and nowhere
    else". This test is here so that reading cannot come back.
    """
    verdicts = estimability_verdicts(
        n_cells_mature=n_mature, n_cells_resolved=n_mature + extra_unresolved
    )
    assert verdicts["estimability"] == "ok"
    assert verdicts["compositional_estimability"] == "ok"


def test_the_compositional_gate_binds_only_below_the_intrinsic_one():
    """Counted over the whole feasible grid, not argued from the definition."""
    binds_with_intrinsic_ok = 0
    binds_with_intrinsic_flagged = 0
    for resolved in range(0, 200):
        for mature in range(0, resolved + 1):
            verdicts = estimability_verdicts(
                n_cells_mature=mature, n_cells_resolved=resolved
            )
            if verdicts["compositional_estimability"] == "ok":
                continue
            if verdicts["estimability"] == "ok":
                binds_with_intrinsic_ok += 1
            else:
                binds_with_intrinsic_flagged += 1
    assert binds_with_intrinsic_ok == 0
    assert binds_with_intrinsic_flagged > 0


# ---------------------------------------------------------------------------
# The frame helper — what the gate memo reads
# ---------------------------------------------------------------------------


def _counts_frame() -> pd.DataFrame:
    """A ``mature_cell_counts()``-shaped frame, keys and all."""
    return pd.DataFrame(
        {
            "patient_id": ["C119", "C165", "C124", "C001"],
            "tissue": ["tumour", "normal", "tumour", "tumour"],
            "labeling_axis": ["stem_pole"] * 4,
            "granularity_rung": ["lineage"] * 4,
            "n_cells_mature": [0, 21, 53, 600],
            "n_cells_resolved": [0, 26, 55, 900],
        }
    )


def test_frame_helper_adds_both_verdicts_and_keeps_the_keys():
    out = classify_counts_frame(_counts_frame())
    assert list(out["estimability"]) == [
        "not_estimable",
        "wide_interval",
        "ok",
        "ok",
    ]
    assert list(out["compositional_estimability"]) == [
        "not_estimable",
        "wide_interval",
        "ok",
        "ok",
    ]
    for key in ("patient_id", "tissue", "labeling_axis", "granularity_rung"):
        assert key in out.columns


def test_frame_helper_does_not_mutate_its_input():
    counts = _counts_frame()
    classify_counts_frame(counts)
    assert "estimability" not in counts.columns


def test_frame_helper_rejects_a_frame_that_is_not_mature_cell_counts():
    counts = _counts_frame()
    counts.loc[0, "n_cells_mature"] = counts.loc[0, "n_cells_resolved"] + 1
    with pytest.raises(ValueError, match="subset of resolved"):
        classify_counts_frame(counts)


def test_frame_helper_names_the_column_it_is_missing():
    counts = _counts_frame().drop(columns=["n_cells_resolved"])
    with pytest.raises(ValueError, match="n_cells_resolved"):
        classify_counts_frame(counts)
