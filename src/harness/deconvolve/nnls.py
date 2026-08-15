"""Plain non-negative least squares. The baseline the others must beat.

No regularisation, no feature selection, no marker weighting. If a method with
all three cannot beat this, the extra machinery is not earning its place, and
that is a reportable result rather than an awkward one.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import nnls

from src.harness.deconvolve.base import Deconvolver


class NNLSDeconvolver(Deconvolver):
    name = "nnls"
    requires_external = False

    def __init__(self, *, normalise_columns: bool = True) -> None:
        #: Scale each signature column to unit sum, so a cell type with a larger
        #: library does not absorb signal purely for being bigger.
        self.normalise_columns = normalise_columns

    def fit_predict(self, bulk: np.ndarray, signature: pd.DataFrame) -> pd.Series:
        self._check_inputs(bulk, signature)

        s = signature.to_numpy(dtype=float)
        if self.normalise_columns:
            col_sums = s.sum(axis=0)
            col_sums[col_sums == 0] = 1.0
            s = s / col_sums

        coefs, _ = nnls(s, np.asarray(bulk, dtype=float))
        return self._normalise(pd.Series(coefs, index=signature.columns))
