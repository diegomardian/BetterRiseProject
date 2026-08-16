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
import pandas as pd

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


def run_negative_controls(
    counts: np.ndarray,
    cell_type: Sequence[str],
    patient_id: Sequence[str],
    genes: Sequence[str],
    *,
    target_gene: str,
    housekeeping: Sequence[str],
    composition_normal: dict[str, float],
    composition_tumour: dict[str, float],
    n_cells: int = 2000,
    n_replicates: int = 20,
    seed: int,
    weighting: str = "normal",
    mature_label: str = "mature_colonocyte",
) -> pd.DataFrame:
    """Run both negative controls end to end. Emits ``CONTROLS_COLUMNS``.

    Three arms, and the contrast between them is the whole value:

    ``target``
        The real gene, genuinely shifted. The positive reference — if this does
        not show an intrinsic term, the other two prove nothing.
    ``housekeeping``
        Unshifted, and flat across cell types. Both terms should be ~0. Signal
        here is a bug in the harness or the estimator, never biology.
    ``permuted``
        Cell-type labels shuffled **within patient**, so composition and every
        batch effect stay intact while the label-expression association is
        destroyed. Both terms must collapse. If either survives, the estimator
        is reading something other than the labels.

    A control that has never been run is not a control, which is why this
    returns a table rather than an assertion — the numbers go in the gate memo.
    """
    from src.harness.pseudobulk import generate_pseudobulk, patient_holdout

    cell_type = np.asarray(cell_type)
    patient_id = np.asarray(patient_id)
    genes = list(genes)

    rows: list[dict] = []
    for rep in range(n_replicates):
        rep_seed = seed + rep
        _, held = patient_holdout(patient_id, n_held_out=2, seed=rep_seed)

        arms = {
            "target": (cell_type, {target_gene: 0.5}),
            "housekeeping": (cell_type, {g: 1.0 for g in housekeeping if g in genes}),
            "permuted": (
                permute_labels_within_patient(cell_type, patient_id, seed=rep_seed),
                {target_gene: 0.5},
            ),
        }
        for control, (labels, shift) in arms.items():
            if not shift:
                continue
            sample = generate_pseudobulk(
                counts, labels, patient_id, genes,
                composition_normal=composition_normal,
                composition_tumour=composition_tumour,
                shift=shift, held_out_patients=held, n_cells=n_cells,
                seed=rep_seed, mature_label=mature_label,
            )
            for gene, terms in sample.truth.realised.items():
                for term, value in terms[weighting].items():
                    if term == "total":
                        continue
                    rows.append(
                        {
                            "control": control,
                            "gene": gene,
                            "term": term,
                            "weighting": weighting,
                            "value": value,
                            "ci_low": None,
                            "ci_high": None,
                            "seed": rep_seed,
                        }
                    )
    return pd.DataFrame(rows)


def summarise_negative_controls(controls: pd.DataFrame) -> pd.DataFrame:
    """Median absolute term per (control, term), scaled by the target arm.

    ``relative_to_target`` is the number to read: the housekeeping and permuted
    arms should be a small fraction of the target arm. An arm that comes back
    near 1.0 has not been controlled at all.
    """
    summary = (
        controls.assign(abs_value=controls["value"].abs())
        .groupby(["control", "term"])["abs_value"]
        .median()
        .reset_index()
        .rename(columns={"abs_value": "median_abs"})
    )
    target = summary[summary["control"] == "target"].set_index("term")["median_abs"]
    summary["relative_to_target"] = [
        row.median_abs / target[row.term] if target.get(row.term) else float("nan")
        for row in summary.itertuples()
    ]
    return pd.DataFrame(summary)


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
