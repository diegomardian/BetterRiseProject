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
"""

from src.estimator.kitagawa import Decomposition, decompose

__all__ = ["Decomposition", "decompose"]
