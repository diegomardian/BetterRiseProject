"""Positivity cutpoints — where the third segment comes from. W2 owns these.

If a tumour has essentially no mature cells left, "how much does each mature
cell make" is not a hard question, it is an undefined one.

The values below are PROVISIONAL (execution_plan.md §4, W2 week 5). They are
replaced by cutpoints derived from where CI width crosses a stated threshold on
the simulation harness — derived, not chosen. When W2 recalibrates, edit
``CUTPOINTS`` here and nothing else: W4 imports this function rather than
reimplementing the rule.
"""

from __future__ import annotations

from typing import Literal, NamedTuple


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


def gate_g4_verdict(n_cells_mature: list[int], cutpoints: Cutpoints = CUTPOINTS) -> dict:
    """Gate criterion G4: are fewer than 50% of patients below the threshold?

    Returns the numbers the week-5 gate needs, with a pre-committed consequence
    attached so it is decided against the criterion and not against the mood in
    the room.
    """
    n = len(n_cells_mature)
    if n == 0:
        raise ValueError("no patients supplied")
    below = sum(1 for c in n_cells_mature if c < cutpoints.wide)
    frac = below / n
    passes = frac < NON_IDENTIFIABILITY_HEADLINE_FRACTION
    return {
        "n_patients": n,
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
