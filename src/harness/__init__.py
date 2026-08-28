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
rather than reimplementing the thresholds. So is its compositional twin,
``classify_compositional_estimability`` (decision #22): the two arms are gated
on different counts and must not be folded into one verdict.
"""

from src.harness.controls import (
    HOUSEKEEPING_GENES,
    housekeeping_panel,
    permute_labels_within_patient,
)
from src.harness.positivity import (
    COMPOSITIONAL_CUTPOINTS,
    CUTPOINTS,
    Cutpoints,
    classify_compositional_estimability,
    classify_counts_frame,
    classify_estimability,
    estimability_verdicts,
    gate_g4_verdict,
)
from src.harness.pseudobulk import (
    PseudobulkSample,
    generate_pseudobulk,
    patient_holdout,
)
from src.harness.results import write_harness_table
from src.harness.truth import GroundTruth, analytic_terms, assert_identity_closes

__all__ = [
    "COMPOSITIONAL_CUTPOINTS",
    "CUTPOINTS",
    "HOUSEKEEPING_GENES",
    "Cutpoints",
    "GroundTruth",
    "PseudobulkSample",
    "analytic_terms",
    "assert_identity_closes",
    "classify_compositional_estimability",
    "classify_counts_frame",
    "classify_estimability",
    "estimability_verdicts",
    "gate_g4_verdict",
    "generate_pseudobulk",
    "housekeeping_panel",
    "patient_holdout",
    "permute_labels_within_patient",
    "write_harness_table",
]
