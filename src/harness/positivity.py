"""Positivity cutpoints — where the third segment comes from. W2 owns these.

If a tumour has essentially no mature cells left, "how much does each mature
cell make" is not a hard question, it is an undefined one.

There are **two** cutpoints, on two different counts, because there are two arms:

- the **intrinsic** arm gates on ``n_cells_mature`` — ``CUTPOINTS``,
  ``classify_estimability``. PROVISIONAL (execution_plan.md §4, W2 week 5):
  replaced by values derived from where CI width crosses a stated threshold on
  the simulation harness — derived, not chosen. When W2 recalibrates, edit
  ``CUTPOINTS`` here and nothing else; W4 imports the function rather than
  reimplementing the rule.
- the **compositional** arm gates on ``n_cells_resolved`` —
  ``COMPOSITIONAL_CUTPOINTS``, ``classify_compositional_estimability``.
  Pre-committed as decision #22, not provisional.

They are reported separately and must not be folded into one verdict; see
``estimability_verdicts``.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING, Literal, NamedTuple

if TYPE_CHECKING:  # pandas is not a runtime dependency of the cutpoints
    import pandas as pd


class Cutpoints(NamedTuple):
    """Mature-cell counts at which the intrinsic term changes status."""

    ok: int  # n >= ok            -> contributes to the intrinsic estimate
    wide: int  # wide <= n < ok    -> wide-interval flag, sensitivity analysis
    # n < wide                     -> not estimable. Undefined, NOT zero.
    source: str


PROVISIONAL = Cutpoints(ok=50, wide=20, source="provisional (execution_plan.md §4)")

#: Swap this for the calibrated object at week 5. Keep the provisional one
#: around so the gate memo can show both.
CUTPOINTS: Cutpoints = PROVISIONAL

#: If more than this share of patients fall below the ``wide`` threshold,
#: non-identifiability stops being a caveat and becomes the headline result
#: (gate criterion G4).
NON_IDENTIFIABILITY_HEADLINE_FRACTION = 0.50

Estimability = Literal["ok", "wide_interval", "not_estimable"]


def classify_estimability(n_cells_mature: int, cutpoints: Cutpoints = CUTPOINTS) -> Estimability:
    """Map a mature-cell count to the schema's ``estimability`` value.

    >>> classify_estimability(120)
    'ok'
    >>> classify_estimability(30)
    'wide_interval'
    >>> classify_estimability(3)
    'not_estimable'

    A ``not_estimable`` verdict means the caller MUST write ``None`` into
    ``intrinsic``, not ``0.0`` (CLAUDE.md invariant 1). The compositional term
    is still estimable in that case — do not drop the row.
    """
    if n_cells_mature < 0:
        raise ValueError(f"n_cells_mature={n_cells_mature} is negative")
    if n_cells_mature >= cutpoints.ok:
        return "ok"
    if n_cells_mature >= cutpoints.wide:
        return "wide_interval"
    return "not_estimable"


# ---------------------------------------------------------------------------
# The compositional arm — decision #22, pre-committed 2026-08-27
# ---------------------------------------------------------------------------
#
# ``classify_estimability`` gates the INTRINSIC arm on ``n_cells_mature``: too
# few mature cells and "how much does each mature cell make" is undefined. There
# was no matching gate on the COMPOSITIONAL arm, so a ``mature_fraction`` of
# 0.92 computed on 9% of the epithelium read identically to one computed on 90%
# (issue #36, raised by W1 as decision #20).
#
# The quantity is ``n_cells_resolved``, not ``unresolved_fraction``. The interval
# on a proportion is driven by the count in its denominator, not by the share
# excluded: 40-of-60 and 400-of-600 share an ``unresolved_fraction`` and have
# very different precision. Same reason the intrinsic cutpoint is a count.
#
# The numbers are the intrinsic arm's, unchanged. Symmetry is the honest default
# rather than a new invention, and the threshold was fixed in public on issue #36
# before it was applied to anything.

#: Decision #22. Same shape and same numbers as the intrinsic rule, applied to
#: the count in the denominator of ``mature_fraction``. Recalibrate this
#: separately from ``CUTPOINTS`` only with a written reason — the two answer
#: different questions and the symmetry is a default, not a finding.
COMPOSITIONAL_CUTPOINTS: Cutpoints = Cutpoints(
    ok=50, wide=20, source="decision #22, pre-committed 2026-08-27 (issue #36)"
)


def classify_compositional_estimability(
    n_cells_resolved: int, cutpoints: Cutpoints = COMPOSITIONAL_CUTPOINTS
) -> Estimability:
    """Map a resolved-epithelium count to an estimability value for the fraction.

    >>> classify_compositional_estimability(120)
    'ok'
    >>> classify_compositional_estimability(30)
    'wide_interval'
    >>> classify_compositional_estimability(3)
    'not_estimable'

    ``n_cells_resolved`` is ``n_cells_epithelial - n_cells_unresolved`` from
    :func:`src.reference.labels.mature_cell_counts` — the denominator
    ``mature_fraction`` is actually computed on. A cell that could not be
    labelled is not a cell measured to be immature (open decision #14), so it is
    excluded from that denominator rather than counted as a failure to mature.
    This cutpoint is what stops the exclusion from being silent.

    ``not_estimable`` here means the **compositional** term is undefined for the
    row. That is a different claim from the intrinsic ``not_estimable`` that
    :func:`classify_estimability` returns, and CLAUDE.md invariant 1 governs
    both: undefined is ``None``, never ``0.0``.
    """
    if n_cells_resolved < 0:
        raise ValueError(f"n_cells_resolved={n_cells_resolved} is negative")
    if n_cells_resolved >= cutpoints.ok:
        return "ok"
    if n_cells_resolved >= cutpoints.wide:
        return "wide_interval"
    return "not_estimable"


def estimability_verdicts(
    n_cells_mature: int,
    n_cells_resolved: int,
    *,
    intrinsic_cutpoints: Cutpoints = CUTPOINTS,
    compositional_cutpoints: Cutpoints = COMPOSITIONAL_CUTPOINTS,
) -> dict[str, object]:
    """Both arms' verdicts for one (patient, tissue, axis, rung) row.

    WHY THIS IS NOT ONE VALUE
    -------------------------
    ``src/schema.py`` is frozen and ``estimability`` is a single enum, so it
    cannot say "intrinsic ok, compositional wide". Until the gate decides it
    needs both on disk per row — a frozen-schema PR, two approvals — the
    compositional verdict travels in harness tables and the gate memo, and this
    function is where it comes from. Do not fold the two into one column by
    taking the worse of them: "the fraction is imprecise" and "there are too few
    mature cells to ask about expression" are different findings, and the second
    is this project's contribution.

    WHAT THE SECOND GATE ACTUALLY BINDS ON — measured, and not what was
    predicted when the decision was recorded
    ------------------------------------------------------------------
    Mature cells are a **subset** of resolved cells, so
    ``n_cells_mature <= n_cells_resolved`` always holds. It follows that
    ``n_cells_mature >= 50`` implies ``n_cells_resolved >= 50``: wherever the
    intrinsic arm is ``ok``, the compositional arm is ``ok`` too, necessarily.
    The compositional gate can therefore only ever bind on rows the intrinsic
    gate has **already** flagged.

    The pre-commitment on issue #36 asserted the inverse — that it binds "on rows
    where the intrinsic arm is already ok ... and nowhere else". That is wrong,
    and it is recorded here rather than quietly fixed, because it changes what
    the decision is worth. Counted on W1's 928-row
    ``mature_cell_counts_full.parquet``, over the 696 rows outside the epithelial
    rung (whose ``mature_fraction`` is 1.0 by construction):

    - rows where the compositional gate binds and the intrinsic arm is ``ok``:
      **0**
    - rows where it binds and the intrinsic arm is not ``ok``: **108**

    So the rule adds a second, independent reason to distrust 108 rows that were
    already flagged, and rescues nothing. It does not reach the "middle of the
    range" exposure issue #36 raised — enough mature cells to clear ``ok``, too
    few resolved cells for the fraction to carry much — because on a **count**
    cutpoint that set is structurally empty. Reaching it would need a cutpoint on
    ``unresolved_fraction``, which decision #22 **declines** with the number
    attached: 4 rows in 2 patients pass both count gates with more than half the
    epithelium unresolved, none above 60%. #22 lists what would reopen it.

    Raises ``ValueError`` if ``n_cells_mature > n_cells_resolved``. That cannot
    happen for a row built by ``mature_cell_counts`` and means the caller passed
    the wrong column — ``n_cells_epithelial`` as the denominator, or counts from
    two different rungs. This project has lost four separate days to a statistic
    computed over the wrong population; a guard that can actually fire is cheap.
    """
    if n_cells_mature > n_cells_resolved:
        raise ValueError(
            f"n_cells_mature={n_cells_mature} exceeds "
            f"n_cells_resolved={n_cells_resolved}. Mature cells are a subset of "
            f"resolved cells, so this row cannot come from mature_cell_counts() "
            f"— check that the denominator is n_cells_resolved and that both "
            f"counts are from the same (axis, rung)."
        )
    return {
        "n_cells_mature": int(n_cells_mature),
        "n_cells_resolved": int(n_cells_resolved),
        "estimability": classify_estimability(n_cells_mature, intrinsic_cutpoints),
        "compositional_estimability": classify_compositional_estimability(
            n_cells_resolved, compositional_cutpoints
        ),
        "intrinsic_threshold": intrinsic_cutpoints.wide,
        "compositional_threshold": compositional_cutpoints.wide,
    }


def classify_counts_frame(
    counts: pd.DataFrame,
    *,
    intrinsic_cutpoints: Cutpoints = CUTPOINTS,
    compositional_cutpoints: Cutpoints = COMPOSITIONAL_CUTPOINTS,
) -> pd.DataFrame:
    """Both verdicts for every row of a ``mature_cell_counts()`` frame.

    Returns a copy with ``estimability`` and ``compositional_estimability``
    added. The frame keeps its own keys — ``patient_id``, ``tissue``,
    ``labeling_axis``, ``granularity_rung`` — so nothing here needs to know how
    many rungs there are.

    The **epithelial rung is not dropped**, because dropping rows is not this
    function's call. Read its output knowing that rung's ``mature_fraction`` is
    1.0 by construction — a denominator choice, still open on W1's side — so its
    compositional verdict describes the precision of a quantity that cannot vary.
    """
    required = {"n_cells_mature", "n_cells_resolved"}
    missing = sorted(required - set(counts.columns))
    if missing:
        raise ValueError(
            f"counts is missing column(s) {missing}. Pass the frame from "
            f"src.reference.labels.mature_cell_counts()."
        )
    out = counts.copy()
    mature = out["n_cells_mature"].astype("int64")
    resolved = out["n_cells_resolved"].astype("int64")
    bad = mature > resolved
    if bool(bad.any()):
        example = out.loc[bad].head(3)
        raise ValueError(
            f"{int(bad.sum())} row(s) have n_cells_mature > n_cells_resolved, "
            f"e.g.\n{example}\nMature cells are a subset of resolved cells, so "
            f"this frame did not come from mature_cell_counts()."
        )
    out["estimability"] = [classify_estimability(int(n), intrinsic_cutpoints) for n in mature]
    out["compositional_estimability"] = [
        classify_compositional_estimability(int(n), compositional_cutpoints) for n in resolved
    ]
    return out


def gate_g4_verdict(
    n_cells_mature: Sequence[int],
    cutpoints: Cutpoints = CUTPOINTS,
    *,
    n_unmatched_patients: int,
) -> dict:
    """Gate criterion G4: are fewer than 50% of patients below the threshold?

    G4'S POPULATION IS MATCHED PATIENTS ONLY
    ----------------------------------------
    Decided 2026-08-23 in answer to issue #9, before any decomposition result
    existed — which is the only time it can be decided honestly.

    ``n_cells_mature`` must contain **one entry per matched patient**: a patient
    with both a tumour and a normal arm. ``n_unmatched_patients`` is the count of
    patients excluded for having no normal arm, and it is **required**, not
    defaulted, because a default is exactly how the wrong population gets used
    without anyone choosing it.

    Why they must not be mixed. On GSE178341 only 36 of 62 patients have a
    matched normal. The compositional term is Δ(mature fraction) against the
    patient's *own* normal, so the other 26 contribute to neither arm. Feeding
    them in as ``n_cells_mature = 0`` would add 26 guaranteed-below rows to 36
    real ones and push ``fraction_below`` toward the 50% line on its own. G4's
    failure consequence is pre-committed — non-identifiability becomes the
    headline result — so a **cohort-design fact**, how many patients had a normal
    sample taken, would be reported as a **positivity finding** about mature-cell
    depletion. Those are different claims and only one of them is about biology.

    The unmatched count is still returned, so the coverage fact travels with the
    verdict instead of being lost. It simply never enters the fraction.

    This is the same failure shape W4 found twice in one afternoon — a cutoff
    computed over a mixed population, producing a plausible wrong number that
    leans toward the hypothesis. It was latent here rather than active, because
    nothing had called it with real data yet.

    Returns the numbers the week-5 gate needs, with a pre-committed consequence
    attached so it is decided against the criterion and not against the mood in
    the room.
    """
    n_cells_mature = list(n_cells_mature)
    n = len(n_cells_mature)
    if n == 0:
        raise ValueError("no matched patients supplied")
    if n_unmatched_patients < 0:
        raise ValueError(f"n_unmatched_patients={n_unmatched_patients} is negative")

    below = sum(1 for c in n_cells_mature if c < cutpoints.wide)
    frac = below / n
    passes = frac < NON_IDENTIFIABILITY_HEADLINE_FRACTION
    total = n + n_unmatched_patients
    return {
        "population": "matched_only",
        "n_patients": n,
        "n_unmatched_excluded": n_unmatched_patients,
        "n_patients_in_cohort": total,
        "matched_fraction": n / total if total else float("nan"),
        "n_below_threshold": below,
        "fraction_below": frac,
        "threshold": cutpoints.wide,
        "passes": passes,
        "consequence": (
            "Proceed."
            if passes
            else "Non-identifiability finding with diagnostics becomes the headline "
            "result, not a caveat. This is a real paper (execution_plan.md §5, G4)."
        ),
    }
