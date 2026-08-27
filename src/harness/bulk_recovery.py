"""Recovering the intrinsic term from bulk. W2.

READ THIS BEFORE USING ANYTHING HERE
------------------------------------
CLAUDE.md invariant 6 forbids cell-type-specific expression imputation from
bulk. This module does exactly that. It is not a loophole — **measuring how
badly it fails is the entire point of §2.2**, and the failure curve is a
publishable object. Nothing in here may produce a number that reaches a results
table as an estimate. It exists to be shown inadequate, quantitatively, instead
of being argued about.

Every function is prefixed ``attenuated_`` so a call site cannot pretend it did
not know.

WHY IT IS ATTENUATED
--------------------
Bulk gives ``Σ_c f_c · μ_c`` for a gene: one number, many cell types. Recovering
``μ_mature`` means subtracting every other type's contribution, so the estimate
inherits the error in every fraction and every reference profile. As the mature
fraction falls, the mature contribution shrinks toward the noise in a
subtraction of much larger numbers, and the estimate degrades non-linearly. The
published band is ×0.6–0.8; the sweep measures ours.

The bias is **directional**: it shrinks the intrinsic term and leaves the
compositional term intact, which pushes toward "compositional" — also our prior
hypothesis. A method that confirms your expectation for methodological reasons
is the worst kind, and this is the specific mechanism by which it would happen.

ON INVARIANT 2
--------------
The reference matrix used to estimate *fractions* excludes target genes, or the
deconvolution would be informed by the gene under test. The target gene's own
per-cell-type profile is taken separately, from **training** patients — a prior
on where the gene lives, learned from cells outside the sample being estimated.
That is not leakage; it is what a reference is.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd


class ReferenceSeamError(ValueError):
    """``reference_profiles`` was called without saying which matrix it is.

    Not a leak — nothing has leaked yet. It means the call cannot be checked
    against invariant 2 because the caller did not state its intent, and the
    default would have been the unsafe one.
    """


def reference_profiles(
    counts: np.ndarray,
    cell_type: Sequence[str],
    genes: Sequence[str],
    *,
    exclude_genes: Sequence[str] | None = None,
    include_targets: bool = False,
) -> pd.DataFrame:
    """Mean profile per cell type: genes x cell types.

    THE CALLER MUST SAY WHICH OF THE TWO MATRICES IT WANTS
    ------------------------------------------------------
    This function builds both, and they have opposite requirements under
    invariant 2:

    - **the fractions reference** — pass ``exclude_genes=[target, ...]``. The
      deconvolution must not be informed by the gene under test.
    - **the target gene's own per-cell-type profile** — pass
      ``include_targets=True``. That is not leakage, it is what a reference is:
      a prior on where the gene lives, learned from *training* patients outside
      the sample being estimated. See the module docstring.

    Neither is a default, and passing neither raises
    :class:`ReferenceSeamError`. ``exclude_genes=()`` used to be the signature's
    default, which meant a caller who simply forgot got the target genes in
    their deconvolution reference and no complaint — invariant 2 violated by
    omission. ``LeeCohort.raw_counts`` deliberately carries the panel so that
    W2's generator can apply a multiplicative shift to it, so the unsafe call
    was one keystroke away and nothing would have said so.

    An **empty** ``exclude_genes`` is refused for the same reason
    ``build_signature`` refuses an empty target set: it silently disables the
    invariant for exactly the genes the panel is built on.
    """
    if include_targets and exclude_genes:
        raise ReferenceSeamError(
            f"include_targets=True and exclude_genes={list(exclude_genes)} "
            f"contradict each other. The target-gene profile keeps the targets; "
            f"the fractions reference drops them. Build them in two calls."
        )
    if not include_targets:
        if exclude_genes is None:
            raise ReferenceSeamError(
                "reference_profiles needs to know which matrix this is. Pass "
                "exclude_genes=[target] for the reference that estimates "
                "fractions (CLAUDE.md invariant 2 — the deconvolution must not "
                "be informed by the gene under test), or include_targets=True "
                "for the target gene's own per-cell-type profile. There is no "
                "safe default: the counts frames in this project carry the "
                "panel on purpose."
            )
        if len(list(exclude_genes)) == 0:
            raise ReferenceSeamError(
                "exclude_genes is empty, which excludes nothing and disables "
                "invariant 2 for exactly the genes the panel is built on. Pass "
                "the target genes, or include_targets=True if you meant the "
                "target profile."
            )

    exclude_genes = () if exclude_genes is None else exclude_genes
    cell_type = np.asarray(cell_type)
    genes = list(genes)
    drop = set(exclude_genes)
    keep = [j for j, g in enumerate(genes) if g not in drop]
    kept_genes = [genes[j] for j in keep]

    types = sorted(set(cell_type.tolist()))
    profile = np.zeros((len(kept_genes), len(types)))
    for k, t in enumerate(types):
        rows = cell_type == t
        if rows.any():
            profile[:, k] = counts[np.ix_(rows, keep)].mean(axis=0)
    return pd.DataFrame(profile, index=kept_genes, columns=types)


def attenuated_mature_mean(
    bulk: np.ndarray,
    genes: Sequence[str],
    *,
    gene: str,
    fractions: pd.Series,
    n_cells: int,
    target_profile: pd.Series,
    mature_label: str = "mature_colonocyte",
) -> float:
    """Per-cell mean of ``gene`` in mature cells, backed out of bulk.

    ``bulk[gene] = n_cells · Σ_c f_c · μ_c``, so::

        μ_mature = (bulk[gene]/n_cells − Σ_{c≠mature} f_c·μ_c) / f_mature

    where the ``μ_c`` for the non-mature types come from ``target_profile``,
    fitted on training patients. Every term on the right carries error, and
    ``f_mature`` in the denominator is what makes the whole thing blow up as the
    mature compartment empties.

    Returns ``nan`` when ``f_mature`` is zero — undefined, not zero. The caller
    must not turn that into 0.0 (CLAUDE.md invariant 1).
    """
    genes = list(genes)
    if gene not in genes:
        raise KeyError(f"{gene!r} is not in the bulk gene list")
    if mature_label not in fractions.index:
        raise KeyError(f"{mature_label!r} is not among the estimated fractions")

    f_mature = float(fractions[mature_label])
    if f_mature <= 0.0:
        return float("nan")

    per_cell_total = float(bulk[genes.index(gene)]) / n_cells
    other = sum(
        float(fractions[c]) * float(target_profile.get(c, 0.0))
        for c in fractions.index
        if c != mature_label
    )
    return (per_cell_total - other) / f_mature


def attenuation_ratio(recovered: float | None, truth: float | None) -> float:
    """``recovered / truth``, the quantity the §2.2 curve plots.

    1.0 is perfect recovery, 0.6–0.8 the published band, 0.0 total attenuation.

    Returns ``nan`` when the ratio is undefined: a zero truth (the null rows,
    which are read for absolute error instead), or a ``None`` estimate. ``None``
    is what the bulk arm returns once the mature compartment is too empty to
    ask — the case the third segment exists for — so it reaches this function
    routinely and must not be coerced to a number.
    """
    if recovered is None or truth is None:
        return float("nan")
    if truth == 0.0 or not np.isfinite(truth) or not np.isfinite(recovered):
        return float("nan")
    return recovered / truth
