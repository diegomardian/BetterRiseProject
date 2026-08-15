"""Kitagawa decomposition. W4.

Kitagawa (1955) demographic standardisation. NOT regression Oaxaca-Blinder,
which decomposes by regression coefficients and answers a different question.

    compositional = Δ(mature fraction) × normal per-cell mean
    intrinsic     = tumour mature fraction × Δ(per-cell mean)
    interaction   = Δ(mature fraction) × Δ(per-cell mean)

The split is not unique: normal-weighted and tumour-weighted give different
answers and the difference lives in the interaction term. Report both plus
doubly-robust. CLAUDE.md invariant 7: the interaction term is reported
separately and never folded into either arm.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

Weighting = Literal["normal", "tumour", "doubly_robust"]


@dataclass(frozen=True)
class Decomposition:
    """One (patient, gene, rung, axis, weighting) result, pre-schema.

    ``intrinsic`` is None when the estimate is not identifiable. It is never
    0.0 in that case — CLAUDE.md invariant 1. ``src.schema.validate_results``
    will reject the frame if this is violated, but the honest value starts here.
    """

    compositional: float | None
    intrinsic: float | None
    interaction: float | None
    weighting: Weighting
    n_cells_mature: int


def decompose(
    frac_mature_normal: float,
    frac_mature_tumour: float,
    mean_normal: float,
    mean_tumour: float,
    *,
    n_cells_mature: int,
    weighting: Weighting = "normal",
) -> Decomposition:
    """Two-term Kitagawa split with the interaction reported separately.

    The identity holds exactly for the normal weighting::

        total = compositional + intrinsic + interaction

    Estimability is NOT decided here — call
    ``src.harness.classify_estimability(n_cells_mature)`` and null out the
    intrinsic term yourself, so the decision is visible at the call site rather
    than buried in the arithmetic.

    W4 — the scalar identity below is correct and unit-tested; the real work is
    the per-patient version over the AnnData, the doubly-robust reweighting
    (weeks 3-4) and the patient-level bootstrap (weeks 4-5).
    """
    d_frac = frac_mature_tumour - frac_mature_normal
    d_mean = mean_tumour - mean_normal

    if weighting == "normal":
        compositional = d_frac * mean_normal
        intrinsic = frac_mature_normal * d_mean
    elif weighting == "tumour":
        compositional = d_frac * mean_tumour
        intrinsic = frac_mature_tumour * d_mean
    elif weighting == "doubly_robust":
        raise NotImplementedError(
            "W4 — doubly-robust reweighting, weeks 3-4. Report agreement with "
            "the plain version quantified."
        )
    else:
        raise ValueError(f"unknown weighting {weighting!r}")

    # Sign convention: normal-weighted puts the cross term positive, tumour-
    # weighted has already absorbed it, so it is reported as its negation.
    interaction = d_frac * d_mean * (1.0 if weighting == "normal" else -1.0)

    return Decomposition(
        compositional=compositional,
        intrinsic=intrinsic,
        interaction=interaction,
        weighting=weighting,
        n_cells_mature=n_cells_mature,
    )


def bootstrap_over_patients(*args, **kwargs):
    """Patient-level bootstrap. CLAUDE.md invariant 5 — patients, not cells.

    W4 — unimplemented, weeks 4-5. Resampling cells inflates n by roughly the
    number of cells per patient and produces intervals that are wrong by
    an order of magnitude in the direction that flatters the result.
    """
    raise NotImplementedError("W4 — see src/estimator/README.md, weeks 4-5.")
