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

from collections.abc import Sequence
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
