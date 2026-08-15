"""Pseudobulk generator with known ground truth. W2, weeks 2-3.

The harness is the only place a true (compositional, intrinsic) split exists.
Everything downstream is measured against what this module generates.

Design requirement, from execution_plan.md §4: it must be able to generate
arbitrary (composition, intrinsic) pairs INCLUDING near-zero mature-cell edge
cases. Those edge cases are not an afterthought — they are how the third
segment gets validated. A generator that only produces comfortable mixtures
validates only the easy half.

Write the design spec (week 1) before this file. All four review it.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class GroundTruth:
    """What the generator knows and the estimator must recover."""

    fractions: dict[str, float]  # cell type -> fraction of the mixture
    per_cell_shift: dict[str, float]  # gene -> log-fold shift applied within cells
    n_cells_mature: int  # drives the expected estimability verdict
    seed: int

    def __post_init__(self) -> None:
        total = sum(self.fractions.values())
        if not np.isclose(total, 1.0, atol=1e-6):
            raise ValueError(f"fractions sum to {total}, not 1.0")


def generate_pseudobulk(*args, **kwargs):
    """Mix held-out patients at known fractions with known per-cell shifts.

    W2 — unimplemented. Signature is yours to design; the spec is due week 1 and
    the implementation weeks 2-3.

    Non-negotiables when you write it:
      - hold patients out, not cells (the patient is the unit of inference)
      - the seed is an argument and lands in the GroundTruth record
      - mature-cell counts spanning 0 to well above the ``ok`` cutpoint
      - a passthrough case: zero shift must return zero intrinsic, not noise
    """
    raise NotImplementedError("W2 — see src/harness/README.md, weeks 2-3.")


def negative_control_permutation(*args, **kwargs):
    """Within-patient label permutation. Should destroy BOTH terms.

    W2 — unimplemented. If either term survives permutation, the estimator is
    reading something other than the labels.
    """
    raise NotImplementedError("W2 — negative controls, see src/harness/README.md.")
