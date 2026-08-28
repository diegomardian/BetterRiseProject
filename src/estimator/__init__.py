"""W4 — Estimator & replication. Owns the Kitagawa decomposition.

Develops on the Lee cohorts (GSE132465 SMC, GSE144735 KUL3) so it is not queued
behind W1, and delivers independent replication as a by-product. Blocked by its
own cohort only.

Cross-checks against cacoa and QuasiMed. NOT CoCoA-diff — it assumes cell
fractions are not a mediator, which assumes away the compositional arm.

Multiple testing: Benjamini-Hochberg within tier, reported separately for each
term.

Takes its estimability rule from ``src.harness.classify_estimability`` rather
than reimplementing thresholds — W2 recalibrates those at week 5.

Two seams worth knowing before using this package:

- ``weighting="doubly_robust"`` is a **pooled-reference** split, not AIPW, and
  its three terms do not sum to the total — ``compositional + intrinsic`` does,
  with ``interaction`` disclosing the cross term the two arms absorbed. Use
  ``ADDITIVE_WEIGHTINGS`` / ``identity_residual`` rather than summing columns.
  docs/open_decisions.md #9.
- ``qc_flags`` needs a ``compartment`` column and will not default it away.
  Pooled-across-compartment MAD bounds cut SMC's tumour epithelium 29.6 points
  harder than its normal arm. docs/open_decisions.md #12.
"""

from src.estimator.ingest import differential_retention, qc_flags
from src.estimator.kitagawa import (
    ADDITIVE_WEIGHTINGS,
    Decomposition,
    attach_intrinsic_ci,
    bootstrap_over_patients,
    decompose,
    decompose_cohort,
    identity_residual,
)

__all__ = [
    "ADDITIVE_WEIGHTINGS",
    "Decomposition",
    "attach_intrinsic_ci",
    "bootstrap_over_patients",
    "decompose",
    "decompose_cohort",
    "differential_retention",
    "identity_residual",
    "qc_flags",
]
