"""Pseudobulk generator with known ground truth. W2, weeks 2-3.

The harness is the only place a true (compositional, intrinsic) split exists.
Everything downstream is measured against what this module generates.

HOW THE SHIFT REACHES INTEGER COUNTS
------------------------------------
The shift is defined multiplicatively **on the mean**: a factor ``s`` means the
mature-cell mean of that gene is ``s`` times what it was. Real cells give us
counts, not means, so the factor is realised on counts by a mechanism that
preserves the expectation exactly:

===========  ===========================================  ==================
``s``        mechanism                                    ``E[new]``
===========  ===========================================  ==================
``== 1``     untouched                                    ``counts``
``< 1``      ``Binomial(counts, s)``                      ``s * counts``
``> 1``      ``counts + Poisson(counts * (s - 1))``        ``s * counts``
===========  ===========================================  ==================

``s == 1`` is a genuine no-op rather than a resample, so the null is exactly
null and not merely small. That property is asserted in the tests and it is what
lets the ``s = 1.0`` row of the attenuation sweep serve as a control.

HOLDOUT IS BY PATIENT
---------------------
CLAUDE.md invariant 5 is about the bootstrap, but it applies here for the same
reason: cells within a patient are not independent draws. :func:`patient_holdout`
splits on ``patient_id`` and :func:`generate_pseudobulk` refuses to draw from a
patient outside the requested set.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np

from src.harness.truth import GroundTruth, parametric_truth, realised_truth


@dataclass(frozen=True)
class PseudobulkSample:
    """One matched (normal, tumour) pseudobulk pair with its truth attached."""

    bulk_normal: np.ndarray  # (n_genes,) summed counts
    bulk_tumour: np.ndarray  # (n_genes,)
    genes: tuple[str, ...]
    truth: GroundTruth

    @property
    def depth_normal(self) -> int:
        return int(self.bulk_normal.sum())

    @property
    def depth_tumour(self) -> int:
        return int(self.bulk_tumour.sum())


def patient_holdout(
    patient_id: Sequence[str],
    *,
    n_held_out: int,
    seed: int,
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    """Split patients into (train, held_out). Deterministic given the seed.

    Returns patient IDs, not cell indices — the whole point is that the split
    happens at the level of the unit of inference.
    """
    patients = sorted(set(patient_id))
    if n_held_out < 1 or n_held_out >= len(patients):
        raise ValueError(
            f"n_held_out={n_held_out} must be in [1, {len(patients) - 1}] "
            f"for {len(patients)} patients"
        )
    rng = np.random.default_rng(seed)
    order = rng.permutation(len(patients))
    held = tuple(sorted(patients[i] for i in order[:n_held_out]))
    train = tuple(sorted(patients[i] for i in order[n_held_out:]))
    return train, held


def _draw_cells(
    rng: np.random.Generator,
    cell_type: np.ndarray,
    eligible: np.ndarray,
    composition: Mapping[str, float],
    n_cells: int,
) -> np.ndarray:
    """Indices of ``n_cells`` cells drawn to match ``composition``.

    Cell types are filled to their integer quota with the remainder allocated to
    the largest fractions, so the realised composition is as close to the
    requested one as an integer count allows. A type with no eligible cells
    contributes nothing — the realised fraction then differs from the parametric
    one, which is exactly why realised truth is recorded separately.
    """
    quotas = {t: int(np.floor(f * n_cells)) for t, f in composition.items()}
    remainder = n_cells - sum(quotas.values())
    if remainder > 0:
        by_frac = sorted(composition, key=lambda t: composition[t], reverse=True)
        for t in by_frac[:remainder]:
            quotas[t] += 1

    picked: list[np.ndarray] = []
    for t, q in quotas.items():
        if q <= 0:
            continue
        pool = np.flatnonzero(eligible & (cell_type == t))
        if pool.size == 0:
            continue  # realised composition will differ; recorded, not patched
        picked.append(rng.choice(pool, size=q, replace=True))
    if not picked:
        return np.empty(0, dtype=int)
    return np.concatenate(picked)


def _apply_shift(
    rng: np.random.Generator,
    counts: np.ndarray,
    mature_mask: np.ndarray,
    gene_idx: Mapping[str, int],
    shift: Mapping[str, float],
) -> np.ndarray:
    """Multiplicative shift on the mean, realised on counts. See module docstring."""
    out = counts.copy()
    if not mature_mask.any():
        return out
    rows = np.flatnonzero(mature_mask)
    for gene, s in shift.items():
        if s == 1.0:
            continue  # exact null — never resample
        j = gene_idx.get(gene)
        if j is None:
            raise KeyError(f"shift names gene {gene!r}, which is not in the matrix")
        block = out[rows, j]
        if s < 1.0:
            out[rows, j] = rng.binomial(block.astype(np.int64), s)
        else:
            out[rows, j] = block + rng.poisson(block * (s - 1.0))
    return out


def generate_pseudobulk(
    counts: np.ndarray,
    cell_type: Sequence[str],
    patient_id: Sequence[str],
    genes: Sequence[str],
    *,
    composition_normal: Mapping[str, float],
    composition_tumour: Mapping[str, float],
    shift: Mapping[str, float],
    held_out_patients: Sequence[str],
    n_cells: int,
    seed: int,
    mature_label: str = "mature_colonocyte",
) -> PseudobulkSample:
    """Mix held-out patients at known fractions with a known per-cell shift.

    Parameters
    ----------
    counts:
        cells x genes, raw integer counts.
    composition_normal, composition_tumour:
        cell type -> fraction, each summing to 1. The mature fraction in the
        tumour composition is the knob the attenuation sweep turns; take it to
        zero deliberately.
    shift:
        gene -> multiplicative factor applied in mature cells of the TUMOUR
        sample only. ``1.0`` leaves the gene untouched. An all-ones shift is the
        null and must recover an intrinsic term of exactly zero.
    held_out_patients:
        Only cells from these patients are drawn. Enforced, not assumed.

    Returns
    -------
    PseudobulkSample
        with both parametric and realised truth attached.
    """
    counts = np.asarray(counts)
    cell_type = np.asarray(cell_type)
    patient_id = np.asarray(patient_id)
    genes = list(genes)

    if counts.ndim != 2:
        raise ValueError(f"counts must be 2-D (cells x genes), got shape {counts.shape}")
    if counts.shape[0] != cell_type.size or counts.shape[0] != patient_id.size:
        raise ValueError(
            f"counts has {counts.shape[0]} cells but cell_type has {cell_type.size} "
            f"and patient_id has {patient_id.size}"
        )
    if counts.shape[1] != len(genes):
        raise ValueError(f"counts has {counts.shape[1]} genes but {len(genes)} names given")
    if n_cells < 1:
        raise ValueError(f"n_cells={n_cells} must be positive")

    held = set(held_out_patients)
    unknown = held - set(patient_id.tolist())
    if unknown:
        raise ValueError(f"held_out_patients not present in the data: {sorted(unknown)}")
    eligible = np.isin(patient_id, list(held))
    if not eligible.any():
        raise ValueError("no cells belong to the held-out patients")

    rng = np.random.default_rng(seed)
    gene_idx = {g: j for j, g in enumerate(genes)}

    idx_n = _draw_cells(rng, cell_type, eligible, composition_normal, n_cells)
    idx_t = _draw_cells(rng, cell_type, eligible, composition_tumour, n_cells)

    cells_n = counts[idx_n]
    cells_t = counts[idx_t]
    mature_n = cell_type[idx_n] == mature_label
    mature_t = cell_type[idx_t] == mature_label

    cells_t = _apply_shift(rng, cells_t, mature_t, gene_idx, shift)

    shifted_genes = [g for g in shift]
    mean_normal = {
        g: float(cells_n[mature_n, gene_idx[g]].mean()) if mature_n.any() else 0.0
        for g in shifted_genes
    }

    truth = GroundTruth(
        composition_normal=dict(composition_normal),
        composition_tumour=dict(composition_tumour),
        shift=dict(shift),
        n_cells_mature=int(mature_t.sum()),
        patient_ids=tuple(sorted(held)),
        seed=seed,
        mature_label=mature_label,
    )
    # Frozen dataclass, but these are the two fields the generator fills.
    object.__setattr__(truth, "parametric", parametric_truth(truth, mean_normal))
    # realised_truth works on a matrix whose columns ARE the genes it is asked
    # about, so slice down to the shifted set rather than passing the full width.
    cols = [gene_idx[g] for g in shifted_genes]
    object.__setattr__(
        truth,
        "realised",
        realised_truth(
            cells_n[:, cols], cells_t[:, cols], mature_n, mature_t, shifted_genes
        ),
    )

    return PseudobulkSample(
        bulk_normal=cells_n.sum(axis=0),
        bulk_tumour=cells_t.sum(axis=0),
        genes=tuple(genes),
        truth=truth,
    )
