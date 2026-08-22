"""What "known ground truth" means. W2.

The harness applies a multiplicative shift ``s`` to the per-cell mean of target
genes, in mature cells only. That choice is what makes the truth analytic rather
than simulated::

    m_n = mu              normal per-cell mean in mature cells
    m_t = mu * s          tumour per-cell mean, after silencing
    dm  = mu * (s - 1)

and the Kitagawa terms follow in closed form. See :func:`analytic_terms`.

TWO TRUTHS, BOTH RECORDED
-------------------------
*Parametric* truth is what we set: the composition and shift we asked for.
*Realised* truth is what the drawn cells actually have — the empirical fractions
and means in the sample that got generated.

They differ by sampling noise, and the difference is not a nuisance. Recovery
measured against realised truth isolates estimator bias. Recovery measured
against parametric truth carries sampling noise as well. Report only one and a
sampling artefact reads as estimator bias, or vice versa. So both are recorded
on every sample, and the attenuation curve is computed against both.

A NOTE ON THE ESTIMATOR DEFINITION
----------------------------------
The boxed definition in execution_plan.md §4 reads::

    compositional = Δ(mature fraction) × normal per-cell mean
    intrinsic     = tumour mature fraction × Δ(per-cell mean)
    interaction   = Δ(mature fraction) × Δ(per-cell mean)

Those three do not sum to the total. They mix the two weightings: the
compositional term is the normal-weighted one and the intrinsic term is the
tumour-weighted one. Each weighting closes exactly on its own —

    normal:  Δf·m_n + f_n·Δm + Δf·Δm        = f_t·m_t − f_n·m_n
    tumour:  Δf·m_t + f_t·Δm − Δf·Δm        = f_t·m_t − f_n·m_n

— and ``src.estimator.kitagawa.decompose`` implements them coherently, one
weighting at a time. The prose is a loose gloss on "the split is not unique",
not a single triple. This module follows the implementation, and the identity is
asserted in tests for both weightings.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import numpy as np

Weighting = Literal["normal", "tumour"]

#: The identity must close to about this, in units of the total.
IDENTITY_TOL = 1e-9


def analytic_terms(
    frac_mature_normal: float,
    frac_mature_tumour: float,
    mean_normal: float,
    shift: float,
    *,
    weighting: Weighting = "normal",
) -> dict[str, float]:
    """Closed-form Kitagawa terms for a multiplicative shift.

    ``shift=1.0`` is the null: no silencing. It must return an intrinsic term of
    exactly ``0.0`` — not a small number — which is why the arithmetic is
    written as ``(shift - 1.0)`` and not as a difference of two computed means.

    Returns keys ``compositional``, ``intrinsic``, ``interaction``, ``total``.
    The first three sum to the fourth, exactly, for either weighting.
    """
    if shift < 0:
        raise ValueError(f"shift={shift} is negative; a multiplicative shift cannot be")
    for name, f in (("normal", frac_mature_normal), ("tumour", frac_mature_tumour)):
        if not 0.0 <= f <= 1.0:
            raise ValueError(f"frac_mature_{name}={f} outside [0, 1]")

    f_n, f_t, mu = frac_mature_normal, frac_mature_tumour, mean_normal
    d_frac = f_t - f_n
    d_mean = mu * (shift - 1.0)

    if weighting == "normal":
        compositional = d_frac * mu
        intrinsic = f_n * d_mean
        interaction = d_frac * d_mean
    elif weighting == "tumour":
        compositional = d_frac * (mu * shift)
        intrinsic = f_t * d_mean
        interaction = -d_frac * d_mean
    else:
        raise ValueError(f"unknown weighting {weighting!r}")

    total = f_t * (mu * shift) - f_n * mu
    return {
        "compositional": compositional,
        "intrinsic": intrinsic,
        "interaction": interaction,
        "total": total,
    }


def identity_residual(terms: dict[str, float]) -> float:
    """How far the three terms are from summing to the total. Should be ~0."""
    return abs(
        terms["compositional"] + terms["intrinsic"] + terms["interaction"] - terms["total"]
    )


def assert_identity_closes(terms: dict[str, float], tol: float = IDENTITY_TOL) -> None:
    """The truth must be self-consistent before it is used to judge anything."""
    residual = identity_residual(terms)
    scale = max(abs(terms["total"]), 1.0)
    if residual / scale > tol:
        raise AssertionError(
            f"Kitagawa identity does not close: residual={residual:.3e} "
            f"(relative {residual / scale:.3e} > {tol:.1e}). The ground truth is "
            f"internally inconsistent — nothing measured against it means anything."
        )


@dataclass(frozen=True)
class GroundTruth:
    """What the generator knows and the estimator must recover.

    ``parametric`` is what we asked for; ``realised`` is what the drawn cells
    actually have. Both carry the four Kitagawa keys per weighting.
    """

    #: cell type -> fraction, in the normal sample
    composition_normal: dict[str, float]
    #: cell type -> fraction, in the tumour sample
    composition_tumour: dict[str, float]
    #: gene -> multiplicative shift applied in mature cells. 1.0 means untouched.
    shift: dict[str, float]
    #: realised mature-cell count. Drives the estimability verdict.
    n_cells_mature: int
    #: patients contributing cells to this sample. Holdout is by patient.
    patient_ids: tuple[str, ...]
    seed: int
    mature_label: str = "mature_colonocyte"

    #: gene -> weighting -> terms. Filled by the generator.
    parametric: dict[str, dict[str, dict[str, float]]] = field(default_factory=dict)
    realised: dict[str, dict[str, dict[str, float]]] = field(default_factory=dict)
    #: gene -> the four scalars ``kitagawa.decompose()`` takes, measured on the
    #: cells actually drawn. This is what an estimator with perfect access to
    #: the single-cell data would see, and feeding it to ``decompose()`` is how
    #: the harness tests W4's arithmetic rather than restating its own truth.
    realised_stats: dict[str, dict[str, float]] = field(default_factory=dict)

    def __post_init__(self) -> None:
        for label, comp in (
            ("normal", self.composition_normal),
            ("tumour", self.composition_tumour),
        ):
            total = sum(comp.values())
            if not np.isclose(total, 1.0, atol=1e-6):
                raise ValueError(f"{label} composition sums to {total}, not 1.0")
            if any(v < 0 for v in comp.values()):
                raise ValueError(f"{label} composition has a negative fraction")
        if self.mature_label not in self.composition_normal:
            raise ValueError(
                f"mature_label {self.mature_label!r} is not a cell type in the "
                f"normal composition: {sorted(self.composition_normal)}"
            )
        if self.n_cells_mature < 0:
            raise ValueError(f"n_cells_mature={self.n_cells_mature} is negative")
        if any(s < 0 for s in self.shift.values()):
            raise ValueError("shift has a negative factor")

    @property
    def is_null(self) -> bool:
        """True when no gene was shifted. The intrinsic term must recover zero."""
        return all(s == 1.0 for s in self.shift.values())

    def expected_estimability(self) -> str:
        """What ``classify_estimability`` should say about this sample.

        Imported lazily so ``truth`` stays importable without the rest of the
        harness, and so the cutpoint swap at week 5 propagates here too.
        """
        from src.harness.positivity import classify_estimability

        return classify_estimability(self.n_cells_mature)


def parametric_truth(
    gt: GroundTruth,
    mean_normal: dict[str, float],
    *,
    weightings: tuple[Weighting, ...] = ("normal", "tumour"),
) -> dict[str, dict[str, dict[str, float]]]:
    """Analytic terms per (gene, weighting), from what we asked for.

    ``mean_normal`` is the per-cell mean of each gene in mature cells of the
    normal sample — the ``mu`` in the module docstring.
    """
    f_n = gt.composition_normal[gt.mature_label]
    f_t = gt.composition_tumour[gt.mature_label]
    out: dict[str, dict[str, dict[str, float]]] = {}
    for gene, mu in mean_normal.items():
        s = gt.shift.get(gene, 1.0)
        out[gene] = {}
        for w in weightings:
            terms = analytic_terms(f_n, f_t, mu, s, weighting=w)
            assert_identity_closes(terms)
            out[gene][w] = terms
    return out


def realised_truth(
    counts_normal: np.ndarray,
    counts_tumour: np.ndarray,
    is_mature_normal: np.ndarray,
    is_mature_tumour: np.ndarray,
    genes: list[str],
    *,
    weightings: tuple[Weighting, ...] = ("normal", "tumour"),
) -> tuple[dict[str, dict[str, dict[str, float]]], dict[str, dict[str, float]]]:
    """Terms computed from the cells actually drawn, not from what we asked for.

    Parameters are cells x genes count matrices and boolean mature masks. The
    fractions and means here are empirical, so these terms carry sampling noise
    — that is the point of having them alongside the parametric ones.

    Returns ``(terms, stats)``. ``stats`` carries the four scalars
    ``kitagawa.decompose()`` takes, so a caller can run the real estimator on
    them rather than reading these terms back out — which would only restate
    the harness's own arithmetic.
    """
    if counts_normal.shape[1] != len(genes) or counts_tumour.shape[1] != len(genes):
        raise ValueError("count matrices and gene list disagree on width")

    f_n = float(is_mature_normal.mean()) if is_mature_normal.size else 0.0
    f_t = float(is_mature_tumour.mean()) if is_mature_tumour.size else 0.0

    stats: dict[str, dict[str, float]] = {}
    out: dict[str, dict[str, dict[str, float]]] = {}
    for j, gene in enumerate(genes):
        m_n = float(counts_normal[is_mature_normal, j].mean()) if is_mature_normal.any() else 0.0
        m_t = float(counts_tumour[is_mature_tumour, j].mean()) if is_mature_tumour.any() else 0.0
        d_frac, d_mean = f_t - f_n, m_t - m_n
        stats[gene] = {
            "frac_mature_normal": f_n,
            "frac_mature_tumour": f_t,
            "mean_normal": m_n,
            "mean_tumour": m_t,
        }
        out[gene] = {}
        for w in weightings:
            if w == "normal":
                terms = {
                    "compositional": d_frac * m_n,
                    "intrinsic": f_n * d_mean,
                    "interaction": d_frac * d_mean,
                }
            else:
                terms = {
                    "compositional": d_frac * m_t,
                    "intrinsic": f_t * d_mean,
                    "interaction": -d_frac * d_mean,
                }
            terms["total"] = f_t * m_t - f_n * m_n
            assert_identity_closes(terms)
            out[gene][w] = terms
    return out, stats
