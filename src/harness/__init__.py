"""W2 — Method. The simulation harness, the bake-off, and the gate.

The only place true ground truth exists. W2 adjudicates the week-5 gate: G3 is
entirely W2's, and G2 and G4 are read against W2's cutpoints.

Produces:
  - pseudobulk generator: held-out patients mixed at known fractions with known
    per-cell shifts, including near-zero mature-cell edge cases
  - deconvolution bake-off: nu-SVR, CIBERSORTx, MuSiC, BayesPrism, NNLS baseline
  - the attenuation curve (§2.2) — intrinsic recovery vs. mature-cell fraction,
    a publishable object on its own
  - calibrated positivity cutpoints, replacing the provisional ones
  - negative controls: housekeeping genes, within-patient label permutation

Blocked by: W1's five-patient pilot, end of week 2. Nothing else.

``classify_estimability`` is re-exported deliberately — W4 must use W2's rule
rather than reimplementing the thresholds.
"""

from src.harness.positivity import (
    CUTPOINTS,
    Cutpoints,
    classify_estimability,
    gate_g4_verdict,
)

__all__ = ["CUTPOINTS", "Cutpoints", "classify_estimability", "gate_g4_verdict"]
