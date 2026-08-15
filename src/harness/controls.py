"""Negative controls. W2 owns these.

Two of them, and they fail in different directions:

- **Within-patient label permutation** should destroy BOTH terms. If either
  survives, the estimator is reading something other than the labels — batch,
  depth, patient identity, anything.
- **Housekeeping genes** should show NEITHER term. They are expressed in every
  compartment at similar levels, so there is no composition effect to find and
  nothing to silence.

A harness with no negative controls cannot distinguish "the estimator works"
from "the estimator returns numbers".
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

#: Stably expressed across colonic compartments. Neither term should appear.
#: Deliberately excludes anything on the panel or on a labelling axis — check
#: with ``src.common.panel`` before adding to this list.
HOUSEKEEPING_GENES: tuple[str, ...] = (
    "ACTB",
    "GAPDH",
    "TBP",
    "RPL13A",
    "PPIA",
    "B2M",
    "HPRT1",
    "SDHA",
    "YWHAZ",
    "UBC",
)


def housekeeping_panel(exclude: Sequence[str] = ()) -> list[str]:
    """The housekeeping negative-control set, minus anything the caller excludes."""
    drop = set(exclude)
    return [g for g in HOUSEKEEPING_GENES if g not in drop]


def assert_housekeeping_are_not_panel_genes() -> None:
    """Guard against a housekeeping gene that is also a target or a label.

    Called from the tests. If this ever fires, the control is not a control.
    """
    from src.common.panel import axis_genes, load_axes, panel_genes

    label_genes: set[str] = set()
    for axis in load_axes()["axes"]:
        label_genes |= set(axis_genes(axis))
    overlap = (set(HOUSEKEEPING_GENES) & set(panel_genes())) | (
        set(HOUSEKEEPING_GENES) & label_genes
    )
    if overlap:
        raise AssertionError(
            f"housekeeping control genes overlap the panel or the labelling axes: "
            f"{sorted(overlap)}. A control that is also a target is not a control."
        )


def permute_labels_within_patient(
    cell_type: Sequence[str],
    patient_id: Sequence[str],
    *,
    seed: int,
) -> np.ndarray:
    """Shuffle cell-type labels within each patient, preserving per-patient counts.

    Within-patient, so the permutation destroys the label-expression association
    while leaving patient composition and every batch effect exactly as it was.
    A between-patient shuffle would also scramble composition and would not
    isolate what we want to isolate.
    """
    cell_type = np.asarray(cell_type)
    patient_id = np.asarray(patient_id)
    if cell_type.size != patient_id.size:
        raise ValueError(
            f"cell_type has {cell_type.size} entries, patient_id has {patient_id.size}"
        )

    rng = np.random.default_rng(seed)
    out = cell_type.copy()
    for patient in np.unique(patient_id):
        rows = np.flatnonzero(patient_id == patient)
        out[rows] = rng.permutation(cell_type[rows])
    return out


def permutation_preserved_counts(
    original: Sequence[str],
    permuted: Sequence[str],
    patient_id: Sequence[str],
) -> bool:
    """True if the permutation kept every per-patient cell-type count intact.

    A permutation that changed the counts would change composition, which would
    make a surviving compositional term unsurprising rather than alarming.
    """
    original, permuted, patient_id = (
        np.asarray(original),
        np.asarray(permuted),
        np.asarray(patient_id),
    )
    for patient in np.unique(patient_id):
        rows = patient_id == patient
        a = np.unique(original[rows], return_counts=True)
        b = np.unique(permuted[rows], return_counts=True)
        if not (np.array_equal(a[0], b[0]) and np.array_equal(a[1], b[1])):
            return False
    return True
