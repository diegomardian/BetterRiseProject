"""The adapter every deconvolution method implements. W2, weeks 3-4.

The bake-off compares five methods written in two languages by five groups with
five opinions about input format. One protocol, five thin adapters, and the
ranking code never learns which is which.

Staged, per the week-1 decision: NNLS, nu-SVR and MuSiC first — all installable
from the pinned env, no external accounts. CIBERSORTx needs a Docker image and a
registered token; BayesPrism is R and slow. Those two land only if they land,
and :func:`available_methods` reports what actually ran so a skipped method is
visible in the results rather than silently absent.
"""

from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np
import pandas as pd


class Deconvolver(ABC):
    """Estimate cell-type fractions from bulk, given a reference S matrix.

    Fractions only. No cell-type-specific expression recovery — CLAUDE.md
    invariant 6. Bulk gives fractions at r ~ 0.92; intrinsic estimates from bulk
    come back attenuated x0.6-0.8 in the direction of our prior hypothesis,
    which is the worst possible direction for a bias to point.
    """

    #: Short name used in results tables and rankings.
    name: str = "abstract"
    #: False when the method needs something the pinned env cannot provide.
    requires_external: bool = False

    @abstractmethod
    def fit_predict(self, bulk: np.ndarray, signature: pd.DataFrame) -> pd.Series:
        """Return fractions indexed by cell type, non-negative and summing to 1.

        Parameters
        ----------
        bulk:
            (n_genes,) expression for one sample, on the signature's gene index.
        signature:
            genes x cell types, from ``src.reference.build_signature``.
        """

    def is_available(self) -> tuple[bool, str]:
        """(usable, reason). Overridden by adapters with external dependencies."""
        return True, "available"

    def _check_inputs(self, bulk: np.ndarray, signature: pd.DataFrame) -> None:
        if bulk.ndim != 1:
            raise ValueError(f"bulk must be 1-D, got shape {bulk.shape}")
        if bulk.shape[0] != signature.shape[0]:
            raise ValueError(
                f"bulk has {bulk.shape[0]} genes, signature has {signature.shape[0]}. "
                f"Both must be on the shared gene index — that is what makes "
                f"integration a join."
            )
        if signature.isna().to_numpy().any():
            raise ValueError("signature contains NaN")

    @staticmethod
    def _normalise(fractions: pd.Series) -> pd.Series:
        """Clip negatives and renormalise to sum to 1."""
        clipped = fractions.clip(lower=0.0)
        total = clipped.sum()
        if total <= 0:
            # Degenerate fit. Uniform is a more honest answer than all-zero,
            # and the bake-off metrics will punish it appropriately.
            return pd.Series(1.0 / len(clipped), index=clipped.index)
        return clipped / total


def available_methods(methods: list[Deconvolver]) -> tuple[list[Deconvolver], dict[str, str]]:
    """Split a method list into (usable, {skipped_name: reason}).

    The skip reasons go into the bake-off sidecar. "CIBERSORTx was not run
    because no token was configured" is a materially different statement from
    "CIBERSORTx was not run", and only one of them is reportable.
    """
    usable, skipped = [], {}
    for m in methods:
        ok, reason = m.is_available()
        (usable.append(m) if ok else skipped.update({m.name: reason}))
    return usable, skipped
