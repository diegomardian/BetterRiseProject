"""Deconvolution adapters for the bake-off. W2, weeks 3-4.

Staged by week-1 decision:

===============  ========  ===========================================
Method           Stage     Note
===============  ========  ===========================================
``nnls``         now       baseline, scipy only
``nusvr``        wk 3      scikit-learn, in the pinned env
``music``        wk 3-4    R; called via Rscript subprocess, not rpy2
``cibersortx``   if avail  Docker image + registered token
``bayesprism``   if avail  R, slow
===============  ========  ===========================================

Anything not run is reported by name and reason in the bake-off sidecar. A
missing method must be visible, not absent.
"""

from src.harness.deconvolve.base import Deconvolver, available_methods
from src.harness.deconvolve.nnls import NNLSDeconvolver

__all__ = ["Deconvolver", "NNLSDeconvolver", "available_methods"]
