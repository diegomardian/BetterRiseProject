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
        The **same sample** as ``target`` — signal genuinely placed in the true
        mature cells — re-estimated under a shuffled mature mask.
    ``label_blind``
        What an estimator that ignores labels entirely would report: the same
        arithmetic with each arm's mean taken over *all* drawn cells. This is
        the reference the ``permuted`` arm is read against.

    HOW TO READ THE PERMUTED ARM — IT DOES NOT GO TO ZERO
    -----------------------------------------------------
    "Permutation destroys both terms" is the usual shorthand and it is wrong
    here, for a reason worth stating rather than tuning around.

    Silencing 40% of the cells moves the mean of *any* random subset, so after
    shuffling, ``Δ(per-cell mean)`` is not zero — it is the whole-sample
    difference, diluted. Analytically the permuted intrinsic term converges on
    ``f · (mean over all tumour cells − mean over all normal cells)``, which is
    exactly the ``label_blind`` arm. On this cohort that lands near 0.6 of the
    target arm, and a test asserting "≈ 0" would fail against correct code.

    So the control is **permuted ≈ label_blind**: under shuffled labels the
    estimator extracts nothing the labels were carrying, and returns precisely
    what a label-blind estimator returns. That is a sharper statement than
    "small", and it is the one that would actually catch an estimator reading
    batch, depth or patient identity instead of the labels — any of those would
    push the permuted arm *away* from label-blind.

    WHY THE PERMUTATION IS APPLIED AT ESTIMATION, NOT AT GENERATION
    ---------------------------------------------------------------
    Shuffling before generating would apply the silencing to whichever cells the
    shuffle called mature — a real effect on a random subset, not a null. The
    first version did that and came back at 23% of target, which read as a
    partial estimator failure and was a mis-specified control.

    WHAT THIS CONTROL CANNOT TEST
    -----------------------------
    The **compositional** term. The generator *imposes* the mature fraction when
    it draws cells, so a within-sample shuffle preserves the count by
    construction and Δ(mature fraction) is unchanged. Its survival is a property
    of the harness, not a finding about the estimator.

    A control that has never been run is not a control, which is why this
    returns a table rather than an assertion — the numbers go in the gate memo.
    """
    from src.estimator.kitagawa import decompose
    from src.harness.pseudobulk import generate_pseudobulk, patient_holdout

    cell_type = np.asarray(cell_type)
    patient_id = np.asarray(patient_id)
    genes = list(genes)

    def _emit(rows, control, gene, terms, rep_seed):
        for term, value in terms.items():
            if term == "total":
                continue
            rows.append(
                {
                    "control": control, "gene": gene, "term": term,
                    "weighting": weighting, "value": value,
                    "ci_low": None, "ci_high": None, "seed": rep_seed,
                }
            )

    rows: list[dict] = []
    for rep in range(n_replicates):
        rep_seed = seed + rep
        _, held = patient_holdout(patient_id, n_held_out=2, seed=rep_seed)
        shared = dict(
            composition_normal=composition_normal,
            composition_tumour=composition_tumour,
            held_out_patients=held, n_cells=n_cells,
            seed=rep_seed, mature_label=mature_label,
        )

        target_sample = generate_pseudobulk(
            counts, cell_type, patient_id, genes, shift={target_gene: 0.5}, **shared
        )
        for gene, terms in target_sample.truth.realised.items():
            _emit(rows, "target", gene, terms[weighting], rep_seed)

        hk_shift = {g: 1.0 for g in housekeeping if g in genes}
        if hk_shift:
            hk_sample = generate_pseudobulk(
                counts, cell_type, patient_id, genes, shift=hk_shift, **shared
            )
            for gene, terms in hk_sample.truth.realised.items():
                _emit(rows, "housekeeping", gene, terms[weighting], rep_seed)

        # Re-estimate the TARGET sample under a shuffled mature mask.
        rng = np.random.default_rng(rep_seed)
        for gene, arms in target_sample.drawn_expression.items():
            masks = {
                arm: rng.permutation(target_sample.drawn_is_mature[arm])
                for arm in ("normal", "tumour")
            }
            f_n = float(masks["normal"].mean())
            f_t = float(masks["tumour"].mean())
            m_n = float(arms["normal"][masks["normal"]].mean()) if masks["normal"].any() else 0.0
            m_t = float(arms["tumour"][masks["tumour"]].mean()) if masks["tumour"].any() else 0.0
            d = decompose(
                f_n, f_t, m_n, m_t,
                n_cells_mature=int(masks["tumour"].sum()), weighting=weighting,
            )
            _emit(
                rows, "permuted", gene,
                {
                    "compositional": d.compositional,
                    "intrinsic": d.intrinsic,
                    "interaction": d.interaction,
                },
                rep_seed,
            )

            # The reference the permuted arm is read against: every cell
            # counted, labels ignored entirely.
            blind = decompose(
                f_n, f_t,
                float(arms["normal"].mean()),
                float(arms["tumour"].mean()),
                n_cells_mature=int(masks["tumour"].sum()), weighting=weighting,
            )
            _emit(
                rows, "label_blind", gene,
                {
                    "compositional": blind.compositional,
                    "intrinsic": blind.intrinsic,
                    "interaction": blind.interaction,
                },
                rep_seed,
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
