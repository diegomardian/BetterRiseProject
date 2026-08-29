"""SUPERSEDED. The one labeller is ``src.reference.labels.assign_labels``.

**Nothing in the pipeline calls this module any more** (2026-08-28). It is kept
importable so existing tests and notebooks still resolve, and kept documented
rather than deleted so the record of why it was retired travels with the code.

It thresholds a cohort-wide quantile of an inverted marker score, with no depth
handling and no ``unresolved`` category. On Lee that made "mature" mean "sampled
zero stem markers" — 32% of epithelial cells, 4.7x shallower than the rest — and
since normal epithelium is 4.3x shallower than tumour, it manufactured a 46-point
compositional loss in the hypothesised direction out of dropout. Issue #44, gate
memo section 11.

``assign_labels`` thins marker counts to a common depth, marks cells below the
floor ``unresolved_depth`` (not scored, never counted as immature), takes cut
points from each patient's own normal arm, and gates ``best4`` on its actual
markers. Under it the same cohort's gap falls to 25 points and the two finest
rungs stop being identical. Decision #13, answered.

**Do not route around this by calling ``label_cohort`` again.** Two labellers is
the condition that produced the bug; one of them being better is not a fix.

--- original docstring ---

Per-cell maturity labels for structurally different labelling axes. W4.

Builds directly on ``src.reference.assert_no_target_leakage`` (CLAUDE.md
invariant 2): before scoring any cell against an axis's marker genes, this
module checks those markers against the run's target-gene set. It cannot see
the caller's labels either — call it, don't assume it, same warning
``src/reference/README.md`` gives W1.

PROVISIONAL, in the same spirit as ``src/harness/positivity.py``'s cutpoints:
the quantile thresholds below are a reasonable, citable first pass — the
``best4`` threshold of 0.95 comes directly from README design decision 3
("BEST4+ program... under 5% of epithelium even in healthy colon"), the rest
are a coarser-to-finer monotonic schedule matching how the granularity rungs
are described (epithelial = broadest, best4 = narrowest). They are not
calibrated against real Lee or Pelka data. Swap them for whatever W1/W2 lands
on once real cells are in hand; this file's job is the plumbing (leakage
guard, one column per axis/rung, never overwriting another), not the final
number.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Final

import pandas as pd

from src.common.panel import axis_genes, granularity_rungs
from src.reference.signature import assert_no_target_leakage

#: rung -> quantile of the axis score above which a cell counts as "mature".
#: best4 = 0.95 is grounded in README design decision 3; the rest are a
#: monotonic coarse-to-fine schedule, not independently calibrated.
PROVISIONAL_RUNG_QUANTILES: Final[dict[str, float]] = {
    "epithelial": 0.50,
    "lineage": 0.65,
    "crypt_position": 0.85,
    "best4": 0.95,
}

#: Axes whose marker genes point AWAY from maturity: axis 1 (stem_pole) scores
#: high near the stem, axis 2 (opposite_lineage) scores high for the goblet
#: program. Both are inverted before thresholding so "high score" means
#: "consistent with the mature/absorptive call" for every transcript-based
#: axis. Chromatin and spatial (axis 3) carry no marker genes yet — week 13+
#: — and axis_score() raises for them rather than guessing.
_INVERTED_AXES: Final[tuple[str, ...]] = ("stem_pole", "opposite_lineage")


def axis_score(
    expression: pd.DataFrame,
    axis: str,
    *,
    target_genes: Sequence[str],
) -> pd.Series:
    """Per-cell score for one labelling axis: mean expression of its markers.

    Guards invariant 2 before touching the data: the axis's marker genes must
    not overlap ``target_genes``, the panel genes under test for this run —
    see the MUC2/TFF3 collision, docs/open_decisions.md #1. A run testing
    either of those targets may not call this with ``axis="opposite_lineage"``.
    """
    genes = axis_genes(axis)
    assert_no_target_leakage(genes, target_genes, context=f"axis '{axis}' markers")
    usable = [g for g in genes if g in expression.columns]
    if not usable:
        raise ValueError(
            f"none of axis {axis!r}'s marker genes are in the expression matrix "
            f"— chromatin and spatial axes carry no transcript markers "
            f"(week 13+), and every other axis should have at least one"
        )
    return expression.loc[:, usable].mean(axis=1)


def classify_maturity(
    expression: pd.DataFrame,
    axis: str,
    rung: str,
    *,
    target_genes: Sequence[str],
) -> pd.Series:
    """Boolean per-cell maturity call for one axis at one granularity rung.

    Thresholds at ``PROVISIONAL_RUNG_QUANTILES[rung]`` of the (possibly
    inverted, see ``_INVERTED_AXES``) axis score. This is a real cutoff you
    can run today, not a placeholder that raises — but it is provisional, and
    swapping it out does not change any caller's code, only the constants
    above.
    """
    if rung not in PROVISIONAL_RUNG_QUANTILES:
        raise KeyError(
            f"unknown rung {rung!r}; frozen rungs are {sorted(PROVISIONAL_RUNG_QUANTILES)}"
        )
    score = axis_score(expression, axis, target_genes=target_genes)
    if axis in _INVERTED_AXES:
        score = -score
    threshold = score.quantile(PROVISIONAL_RUNG_QUANTILES[rung])
    return score >= threshold


def label_cohort(
    expression: pd.DataFrame,
    *,
    patient_id: Sequence[str],
    tissue: Sequence[str],
    target_genes: Sequence[str],
    axes: Sequence[str] = ("stem_pole", "opposite_lineage"),
    rungs: Sequence[str] | None = None,
) -> pd.DataFrame:
    """One maturity column per (axis, rung), alongside patient_id and tissue.

    Never overwrites another axis/rung's column — each gets its own
    ``mature__{axis}__{rung}`` column, per the W1/W4 task table's "labels in
    separate columns, never overwriting each other." Only axes 1 and 2
    (transcript-based) are usable here; axis 3 is week 13+.
    """
    if len(patient_id) != len(expression) or len(tissue) != len(expression):
        raise ValueError("patient_id and tissue must have one entry per cell")

    out = pd.DataFrame({"patient_id": list(patient_id), "tissue": list(tissue)})
    for axis in axes:
        for rung in rungs or granularity_rungs():
            out[f"mature__{axis}__{rung}"] = classify_maturity(
                expression, axis, rung, target_genes=target_genes
            ).to_numpy()
    return out
