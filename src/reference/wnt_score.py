"""The Wnt-target score, and the partial correlation that makes it mean something.

INVARIANT 8 IS THE WHOLE DESIGN OF THIS FILE. It says, in CLAUDE.md:

    CTNNB1 / TCF7L2 transcript level is not Wnt activity. Use a target
    signature (AXIN2, NKD1, RNF43, NOTUM, TCF7); drop ASCL2/LGR5 when the stem
    axis is in play.

Both halves are load-bearing here. The first is why ``SIGNATURE`` is five
downstream targets rather than the pathway's own transcription factors — a
cell's CTNNB1 mRNA says nothing about whether beta-catenin is in the nucleus.
The second is sharper than it looks: ASCL2 and LGR5 are Wnt targets *and* they
are two of the five markers in the ``stem_pole`` labelling axis. Scoring them
would put the labeller's own markers inside the predictor and make any
correlation with maturity true by construction. The stem axis is always in play
in this project, so they are never in this signature.

WHY A PARTIAL CORRELATION AND NOT A CORRELATION. The question is whether Wnt
activity tracks differentiation output *within* cells that are already called
mature. Two things would produce that association without any Wnt involvement:

*Residual maturity.* The mature bin is a bin, not a point. Cells inside it still
vary in differentiation state, Wnt targets fall with maturity, and GUCA2A rises
with it. An unconditioned correlation recovers that gradient and reports it as
Wnt. This is the single most likely false positive in the design.

*Library depth.* Two per-cell expression scores measured in the same cell are
correlated through library size even after CP10K normalisation, because
detection is depth-dependent. Nothing about that is biological.

So both are conditioned on, and the housekeeping genes are scored through the
identical path as the empirical check on whether the conditioning worked —
whatever ACTB and KRT8 return is the floor, and a target gene's correlation is
interpretable only as the amount by which it exceeds it.

SPEARMAN RATHER THAN PEARSON, because per-cell counts are zero-inflated and
heavy-tailed, so a Pearson correlation over them is dominated by the handful of
cells that fired. That is the same argument that put this project on detection
rather than on means in the first place.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

#: The five Wnt target genes of invariant 8. Not CTNNB1, not TCF7L2, and
#: emphatically not ASCL2 or LGR5 — see the module docstring.
SIGNATURE: tuple[str, ...] = ("AXIN2", "NKD1", "RNF43", "NOTUM", "TCF7")

#: Genes invariant 8 names and this signature must NEVER contain. Kept as data
#: so the guard can say which rule was broken rather than only that one was.
FORBIDDEN: dict[str, str] = {
    "CTNNB1": "transcript level is not pathway activity (invariant 8)",
    "TCF7L2": "transcript level is not pathway activity (invariant 8)",
    "ASCL2": "a stem_pole axis marker; the stem axis is in play (invariant 8)",
    "LGR5": "a stem_pole axis marker; the stem axis is in play (invariant 8)",
}

#: A patient needs this many mature cells before a within-patient correlation
#: over cells is worth computing. Below it the correlation is a statement about
#: a handful of cells and its own standard error swamps it.
MIN_CELLS_FOR_CORRELATION = 30

#: How high the Wnt/stem correlation may run before the score is judged to be
#: measuring maturity rather than Wnt. Invariant 8 requires this correlation be
#: REPORTED when the stem axis is in play; this is the threshold at which
#: reporting turns into withholding the mechanism reading.
WNT_STEM_CORRELATION_CEILING = 0.70


class WntScoreError(ValueError):
    """A Wnt score that cannot be read as one."""


def assert_no_signature_leakage(
    panel_genes: object, axis_genes: object, signature: object = SIGNATURE
) -> None:
    """The signature must share no gene with the panel or the labelling axes.

    THE SAME FAILURE MODE AS INVARIANT 2, ONE OBJECT OVER. ``build_signature()``
    asserts that target genes never appear in the labels, because a silenced
    mature cell must not be readable as an absent mature cell. Here the risk is
    that the *predictor* contains a label marker, which would make the
    correlation with maturity true by construction and the whole test circular.

    Checked rather than assumed, and the failing input is committed.
    """
    signature = tuple(str(g) for g in signature)
    panel = {str(g) for g in panel_genes}
    axes = {str(g) for g in axis_genes}

    forbidden = sorted(set(signature) & set(FORBIDDEN))
    if forbidden:
        raise WntScoreError(
            "invariant 8 names these and the signature contains them: "
            + "; ".join(f"{g} — {FORBIDDEN[g]}" for g in forbidden)
        )
    on_panel = sorted(set(signature) & panel)
    if on_panel:
        raise WntScoreError(
            f"{on_panel} are on the frozen panel AND in the Wnt signature. A "
            f"predictor that contains a scored gene correlates with it for "
            f"arithmetic reasons."
        )
    in_axes = sorted(set(signature) & axes)
    if in_axes:
        raise WntScoreError(
            f"{in_axes} are labelling-axis markers AND in the Wnt signature. "
            f"The labeller places cells using them, so a correlation with "
            f"maturity would be true by construction — this is invariant 2's "
            f"failure mode with the predictor and the label swapped."
        )


def wnt_score(
    counts, gene_names, *, depth, signature: tuple[str, ...] = SIGNATURE
) -> tuple[np.ndarray, dict]:
    """Per-cell Wnt-target score: mean CP10K over the signature genes present.

    Returns ``(score, report)``. The report names which signature genes were
    found, because a score built from two of five is a different quantity from
    one built from five and a reader must be able to tell.

    Genes absent from the matrix are dropped and named, never imputed. Fewer
    than two present raises: a "signature" of one gene is that gene.
    """
    symbols = np.asarray([str(g) for g in gene_names])
    depth = np.asarray(depth, dtype=float)

    located: list[tuple[str, int]] = []
    missing: list[str] = []
    for gene in signature:
        hit = np.flatnonzero(symbols == gene)
        if hit.size:
            located.append((gene, int(hit[0])))
        else:
            missing.append(gene)

    if len(located) < 2:
        raise WntScoreError(
            f"only {len(located)} of {len(signature)} signature genes are in "
            f"the matrix ({missing} absent). A signature of one gene is that "
            f"gene, and this project has been bitten four times by identifier "
            f"spaces -- suspect symbols vs Ensembl before concluding the genes "
            f"are not measured."
        )

    safe_depth = np.where(depth > 0, depth, np.nan)
    columns = []
    for _, index in located:
        column = counts[:, index]
        values = (np.asarray(column.todense()).ravel()
                  if hasattr(column, "todense") else np.asarray(column))
        columns.append(values.astype(float) / safe_depth * 1e4)
    score = np.nanmean(np.vstack(columns), axis=0)
    return score, {
        "signature": list(signature),
        "genes_present": [g for g, _ in located],
        "genes_absent": missing,
        "n_genes_used": len(located),
    }


def _rank(values: np.ndarray) -> np.ndarray:
    """Average ranks, so ties do not create structure."""
    return pd.Series(np.asarray(values, dtype=float)).rank(
        method="average").to_numpy()


def partial_spearman(
    x: np.ndarray, y: np.ndarray, conditioners: np.ndarray
) -> float:
    """Spearman correlation of ``x`` and ``y`` with ``conditioners`` removed.

    Rank-transform everything, regress the two ranked variables on the ranked
    conditioners by least squares, and correlate the residuals. This is the
    standard construction and it is written out here rather than imported so
    that what is being conditioned on is visible at the call site.

    Returns ``nan`` when either residual is constant — a correlation with a
    constant is undefined, not zero.
    """
    x, y = _rank(x), _rank(y)
    z = np.asarray(conditioners, dtype=float)
    if z.ndim == 1:
        z = z[:, None]
    z = np.column_stack([_rank(z[:, j]) for j in range(z.shape[1])])
    design = np.column_stack([np.ones(len(x)), z])

    keep = np.isfinite(x) & np.isfinite(y) & np.all(np.isfinite(design), axis=1)
    if keep.sum() < MIN_CELLS_FOR_CORRELATION:
        return float("nan")
    x, y, design = x[keep], y[keep], design[keep]

    coef_x, *_ = np.linalg.lstsq(design, x, rcond=None)
    coef_y, *_ = np.linalg.lstsq(design, y, rcond=None)
    rx, ry = x - design @ coef_x, y - design @ coef_y
    if rx.std() == 0 or ry.std() == 0:
        return float("nan")
    return float(np.corrcoef(rx, ry)[0, 1])


def wnt_stem_verdict(correlation: float) -> dict:
    """Is the Wnt score just the maturity score? Invariant 8's required report.

    A high correlation does not mean the score is wrong — Wnt targets really do
    fall with differentiation. It means this test cannot separate the two, and
    the mechanism reading is withheld rather than reported with a caveat.
    """
    if not np.isfinite(correlation):
        return {"verdict": "UNDEFINED",
                "detail": "the Wnt/stem correlation could not be computed"}
    if abs(correlation) >= WNT_STEM_CORRELATION_CEILING:
        return {
            "verdict": "WNT SCORE TRACKS MATURITY",
            "detail": (
                f"|r| = {abs(correlation):.3f} against a ceiling of "
                f"{WNT_STEM_CORRELATION_CEILING}. The score and the labeller's "
                f"own maturity axis are close to the same quantity here, so a "
                f"correlation with differentiation output cannot be attributed "
                f"to Wnt. The mechanism reading is WITHHELD, not caveated."
            ),
        }
    return {
        "verdict": "SEPARABLE",
        "detail": (
            f"|r| = {abs(correlation):.3f}, below {WNT_STEM_CORRELATION_CEILING}. "
            f"The Wnt score carries variation the maturity axis does not, which "
            f"is what the partial correlation needs in order to mean anything."
        ),
    }
