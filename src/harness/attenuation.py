"""The §2.2 attenuation sweep. W2, weeks 4-5.

A grid, not a point: how well is the intrinsic term recovered, as a function of
mature-cell fraction? Expected shape is extremes separable, middle band not —
but the shape is a result, so nothing here asserts it.

TWO ARMS, AND THE DIFFERENCE BETWEEN THEM IS THE CURVE
------------------------------------------------------
``oracle``
    ``kitagawa.decompose()`` on summary statistics computed from the cells
    themselves. This is what W4 does on single-cell data, and it is the
    reliable half. It should track truth closely; when it does not, the
    estimator is the problem (**gate criterion G3**).

``bulk``
    Deconvolve the pseudobulk for fractions, then back the mature-cell mean out
    of bulk via ``bulk_recovery.attenuated_mature_mean``, then decompose. This
    is the thing CLAUDE.md invariant 6 forbids using for results. It is run
    here **only** to measure how far it can be pushed, which is the §2.2
    calibration curve and a publishable object on its own.

Reporting the bulk arm without the oracle arm beside it would leave estimator
error and bulk attenuation confounded, and they have different consequences: one
means fix the estimator, the other means do not use bulk for this.

A NOTE ON IMPORTING W4
----------------------
This module imports ``src.estimator.kitagawa.decompose``. Workstreams do not
normally reach into each other, but the harness's entire job is to run the
estimator against known truth — there is no version of that which does not call
it. The dependency is one-way and it is on W4's public API.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from src.estimator.kitagawa import decompose
from src.harness.bulk_recovery import (
    attenuated_mature_mean,
    attenuation_ratio,
    reference_profiles,
)
from src.harness.deconvolve.base import Deconvolver
from src.harness.deconvolve.nnls import NNLSDeconvolver
from src.harness.positivity import classify_estimability
from src.harness.pseudobulk import generate_pseudobulk, patient_holdout

#: Pre-registered in docs/harness_design_spec.md §4, before any sweep was run.
#: 1.0 is the null and must recover zero.
DEFAULT_SHIFTS: tuple[float, ...] = (1.0, 0.8, 0.5, 0.25)

#: Log-spaced down to a compartment that is essentially gone.
DEFAULT_MATURE_FRACTIONS: tuple[float, ...] = (0.40, 0.20, 0.10, 0.05, 0.02, 0.01, 0.0)


@dataclass(frozen=True)
class SweepGrid:
    """What varies and what is held fixed. Everything not here is fixed."""

    mature_fractions: tuple[float, ...] = DEFAULT_MATURE_FRACTIONS
    shifts: tuple[float, ...] = DEFAULT_SHIFTS
    n_replicates: int = 50
    #: Held fixed across the grid — see the design spec's "held fixed" list.
    n_cells: int = 2000
    frac_mature_normal: float = 0.40
    n_held_out: int = 2

    def points(self) -> list[tuple[int, float, float]]:
        """``(grid_id, mature_fraction, shift)`` for every cell of the grid."""
        return [
            (i, f, s)
            for i, (f, s) in enumerate(
                (f, s) for f in self.mature_fractions for s in self.shifts
            )
        ]

    def __post_init__(self) -> None:
        if self.n_replicates < 1:
            raise ValueError(f"n_replicates={self.n_replicates} must be positive")
        if 1.0 not in self.shifts:
            raise ValueError(
                "the grid must include shift=1.0 — it is the null arm, and "
                "without it there is nothing to check the sweep against"
            )


@dataclass
class SweepConfig:
    """The cohort the sweep draws from, plus the gene under test."""

    counts: np.ndarray
    cell_type: Sequence[str]
    patient_id: Sequence[str]
    genes: Sequence[str]
    target_gene: str
    mature_label: str = "mature_colonocyte"
    deconvolver: Deconvolver = field(default_factory=NNLSDeconvolver)


def _composition(mature_frac: float, types: list[str], mature_label: str) -> dict[str, float]:
    others = [t for t in types if t != mature_label]
    rest = (1.0 - mature_frac) / len(others) if others else 0.0
    return {mature_label: mature_frac} | {t: rest for t in others}


def _oracle_arm(sample, target_gene: str, weighting: str) -> dict[str, float]:
    """Run W4's ``decompose()`` on the empirical cell-level summary statistics.

    This must go through the real estimator. Reading
    ``sample.truth.realised`` back out and calling it an estimate would compare
    the harness's arithmetic against itself and return a ratio of exactly 1.0
    for reasons having nothing to do with whether the estimator works.
    """
    stats = sample.truth.realised_stats[target_gene]
    d = decompose(
        stats["frac_mature_normal"],
        stats["frac_mature_tumour"],
        stats["mean_normal"],
        stats["mean_tumour"],
        n_cells_mature=sample.truth.n_cells_mature,
        weighting=weighting,
    )
    return {
        "compositional_hat": d.compositional,
        "intrinsic_hat": d.intrinsic,
        "interaction_hat": d.interaction,
    }


def run_sweep(
    config: SweepConfig,
    grid: SweepGrid,
    *,
    seed: int,
    weighting: str = "normal",
) -> pd.DataFrame:
    """Run the grid and return an ``attenuation``-shaped frame.

    One row per (grid point, replicate, arm). Both truths are carried on every
    row — parametric and realised — because recovery against realised isolates
    estimator bias while recovery against parametric also carries sampling
    noise, and reporting one without the other confuses the two.
    """
    counts = np.asarray(config.counts)
    cell_type = np.asarray(config.cell_type)
    patient_id = np.asarray(config.patient_id)
    genes = list(config.genes)
    types = sorted(set(cell_type.tolist()))
    target = config.target_gene

    if target not in genes:
        raise KeyError(f"target_gene {target!r} is not in the gene list")
    if config.mature_label not in types:
        raise KeyError(f"mature_label {config.mature_label!r} is not a cell type")

    rows: list[dict] = []
    for grid_id, mature_frac, shift in grid.points():
        comp_t = _composition(mature_frac, types, config.mature_label)
        comp_n = _composition(grid.frac_mature_normal, types, config.mature_label)

        for rep in range(grid.n_replicates):
            rep_seed = seed + grid_id * 10_000 + rep
            train, held = patient_holdout(
                patient_id, n_held_out=grid.n_held_out, seed=rep_seed
            )
            sample = generate_pseudobulk(
                counts, cell_type, patient_id, genes,
                composition_normal=comp_n, composition_tumour=comp_t,
                shift={target: shift}, held_out_patients=held,
                n_cells=grid.n_cells, seed=rep_seed,
                mature_label=config.mature_label,
            )
            truth_p = sample.truth.parametric[target][weighting]
            truth_r = sample.truth.realised[target][weighting]
            n_mature = sample.truth.n_cells_mature
            estimability = classify_estimability(n_mature)

            base = {
                "grid_id": grid_id,
                "replicate": rep,
                "gene": target,
                "weighting": weighting,
                "frac_mature_tumour": mature_frac,
                "shift": shift,
                "n_cells_mature": n_mature,
                "estimability": estimability,
                "compositional_true_parametric": truth_p["compositional"],
                "intrinsic_true_parametric": truth_p["intrinsic"],
                "compositional_true_realised": truth_r["compositional"],
                "intrinsic_true_realised": truth_r["intrinsic"],
                "ci_low": None,
                "ci_high": None,
                "seed": rep_seed,
            }

            # Both arms are scored against the PARAMETRIC truth — the known
            # split we asked for, which is what §2.2's question is about. The
            # realised truth is carried alongside so estimator bias and
            # sampling noise can be separated afterwards.
            for arm, terms in (
                ("oracle", _oracle_arm(sample, target, weighting)),
                ("bulk", _bulk_arm(sample, config, grid, train, weighting)),
            ):
                rows.append(
                    base
                    | {"arm": arm}
                    | terms
                    | {
                        "attenuation_ratio": attenuation_ratio(
                            terms["intrinsic_hat"], truth_p["intrinsic"]
                        )
                    }
                )

    return pd.DataFrame(rows)


def _bulk_arm(
    sample, config: SweepConfig, grid: SweepGrid, train, weighting
) -> dict[str, float | None]:
    """Deconvolve, back the mature mean out of bulk, decompose. Attenuated."""
    counts = np.asarray(config.counts)
    cell_type = np.asarray(config.cell_type)
    patient_id = np.asarray(config.patient_id)
    genes = list(config.genes)
    target = config.target_gene

    train_rows = np.isin(patient_id, list(train))
    # Fractions come from a signature WITHOUT the target gene (invariant 2).
    signature = reference_profiles(
        counts[train_rows], cell_type[train_rows], genes, exclude_genes=[target]
    )
    # The target gene's own cell-type profile, from training patients only.
    target_profile = reference_profiles(
        counts[train_rows], cell_type[train_rows], genes
    ).loc[target]

    keep = [j for j, g in enumerate(genes) if g != target]
    try:
        frac_t = config.deconvolver.fit_predict(sample.bulk_tumour[keep], signature)
        frac_n = config.deconvolver.fit_predict(sample.bulk_normal[keep], signature)
    except ValueError:
        return {
            "compositional_hat": None, "intrinsic_hat": None, "interaction_hat": None
        }

    mean_t = attenuated_mature_mean(
        sample.bulk_tumour, genes, gene=target, fractions=frac_t,
        n_cells=grid.n_cells, target_profile=target_profile,
        mature_label=config.mature_label,
    )
    mean_n = attenuated_mature_mean(
        sample.bulk_normal, genes, gene=target, fractions=frac_n,
        n_cells=grid.n_cells, target_profile=target_profile,
        mature_label=config.mature_label,
    )
    if not (np.isfinite(mean_t) and np.isfinite(mean_n)):
        # Undefined, not zero. The compositional term survives; the intrinsic
        # one does not (CLAUDE.md invariant 1).
        return {
            "compositional_hat": None, "intrinsic_hat": None, "interaction_hat": None
        }

    d = decompose(
        float(frac_n[config.mature_label]),
        float(frac_t[config.mature_label]),
        mean_n,
        mean_t,
        n_cells_mature=sample.truth.n_cells_mature,
        weighting=weighting,
    )
    return {
        "compositional_hat": d.compositional,
        "intrinsic_hat": d.intrinsic,
        "interaction_hat": d.interaction,
    }


def summarise_sweep(sweep: pd.DataFrame) -> pd.DataFrame:
    """Collapse replicates to one row per (arm, mature fraction, shift).

    This is the plottable curve: median attenuation ratio and its spread,
    against mature-cell fraction, per arm.
    """
    grouped = sweep.groupby(["arm", "frac_mature_tumour", "shift"], dropna=False)
    out = grouped.agg(
        n_replicates=("replicate", "count"),
        n_cells_mature_median=("n_cells_mature", "median"),
        intrinsic_true=("intrinsic_true_realised", "median"),
        intrinsic_hat_median=("intrinsic_hat", "median"),
        attenuation_median=("attenuation_ratio", "median"),
        attenuation_q25=("attenuation_ratio", lambda s: s.quantile(0.25)),
        attenuation_q75=("attenuation_ratio", lambda s: s.quantile(0.75)),
        n_not_estimable=("estimability", lambda s: int((s == "not_estimable").sum())),
        n_hat_undefined=("intrinsic_hat", lambda s: int(s.isna().sum())),
    ).reset_index()
    return out


def null_arm_recovers_zero(sweep: pd.DataFrame, *, tol: float = 1e-12) -> bool:
    """The shift=1.0 rows must carry a *parametric* intrinsic term of exactly zero.

    Check this before reading any other number in a sweep. If it is False the
    generator is not producing the null it claims to, and every attenuation
    ratio downstream is measured against a moving target.

    Note which truth this reads. The **parametric** null is exact by
    construction — ``analytic_terms`` computes ``(shift - 1.0)``, so 1.0 gives a
    hard zero. The **realised** null is not, and must not be asserted to be:
    the normal and tumour samples are different draws of cells, so their
    empirical means differ by sampling noise even when nothing was silenced.
    Use :func:`null_arm_noise_ratio` to check that noise is small relative to a
    real effect, which is the meaningful version of the question.
    """
    null = sweep[sweep["shift"] == 1.0]
    if null.empty:
        raise ValueError("sweep has no shift=1.0 rows to check")
    return bool(null["intrinsic_true_parametric"].abs().max() <= tol)


def null_arm_noise_ratio(sweep: pd.DataFrame, *, effect_shift: float = 0.5) -> float:
    """Sampling noise in the null, as a fraction of a real effect.

    Median ``|realised intrinsic|`` at shift=1.0 over the same at
    ``effect_shift``. Small is good; approaching 1.0 means the sweep cannot
    distinguish silencing from the noise floor at this cell count, which is a
    finding about the design rather than about the estimator.
    """
    oracle = sweep[sweep["arm"] == "oracle"]
    noise = oracle.loc[oracle["shift"] == 1.0, "intrinsic_true_realised"].abs().median()
    effect = oracle.loc[
        oracle["shift"] == effect_shift, "intrinsic_true_realised"
    ].abs().median()
    if not effect or not np.isfinite(effect):
        raise ValueError(f"no usable shift={effect_shift} rows to compare against")
    return float(noise / effect)


def bulk_overconfidence(sweep: pd.DataFrame) -> pd.DataFrame:
    """Rows where bulk returned a number and the truth was 'not estimable'.

    This is the project's thesis, measured: *every existing method returns a
    number; none flags that the intrinsic estimate is meaningless in a tumour
    with no mature cells left.* Deconvolution assigns a non-zero mature fraction
    even to a sample with zero mature cells, so the division in
    ``attenuated_mature_mean`` goes through and produces a confident, wrong
    answer where the honest output is ``None``.

    Bulk cannot count cells, so it cannot apply the positivity rule at all. That
    is not a tuning problem — it is the structural reason the third segment
    needs single-cell data to exist.
    """
    bulk = sweep[sweep["arm"] == "bulk"]
    return bulk[
        (bulk["estimability"] == "not_estimable") & bulk["intrinsic_hat"].notna()
    ]
