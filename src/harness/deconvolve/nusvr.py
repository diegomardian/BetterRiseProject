"""nu-Support Vector Regression. The CIBERSORT-family approach. W2, week 3-4.

Robustness comes from high dimensionality — this is the method behind
execution_plan.md §2.1 error #4, which is why the signature must be 500-2000
genes and not the 11-gene panel. The bake-off quantifies that claim rather than
arguing it.

Standard recipe: fit nu-SVR at several nu, take the fit with the lowest RMSE,
clip negative coefficients and renormalise.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from src.harness.deconvolve.base import Deconvolver

#: The three values the original CIBERSORT sweeps over.
DEFAULT_NU: tuple[float, ...] = (0.25, 0.50, 0.75)


class NuSVRDeconvolver(Deconvolver):
    name = "nusvr"
    requires_external = False

    def __init__(self, *, nu_values: tuple[float, ...] = DEFAULT_NU, C: float = 1.0) -> None:
        self.nu_values = nu_values
        self.C = C

    def is_available(self) -> tuple[bool, str]:
        try:
            import sklearn  # noqa: F401
        except ImportError:
            return False, "scikit-learn not installed"
        return True, "available"

    def fit_predict(self, bulk: np.ndarray, signature: pd.DataFrame) -> pd.Series:
        from sklearn.svm import NuSVR

        self._check_inputs(bulk, signature)

        # Standardise both sides — nu-SVR is scale-sensitive and the signature
        # columns differ in library size by orders of magnitude.
        x = signature.to_numpy(dtype=float)
        y = np.asarray(bulk, dtype=float)
        x_z = (x - x.mean()) / (x.std() or 1.0)
        y_z = (y - y.mean()) / (y.std() or 1.0)

        best_coefs, best_rmse = None, np.inf
        for nu in self.nu_values:
            model = NuSVR(nu=nu, C=self.C, kernel="linear").fit(x_z, y_z)
            coefs = model.coef_.ravel()
            rmse = float(np.sqrt(np.mean((x_z @ coefs - y_z) ** 2)))
            if rmse < best_rmse:
                best_coefs, best_rmse = coefs, rmse

        return self._normalise(pd.Series(best_coefs, index=signature.columns))
