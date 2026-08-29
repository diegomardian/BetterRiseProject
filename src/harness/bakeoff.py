"""Deconvolution bake-off. W2, weeks 3-4.

Ranks methods on **fraction recovery** against harness ground truth — fractions
only, never cell-type-specific expression (CLAUDE.md invariant 6). The question
is which deconvolver to trust for the mature-colonocyte fraction that Stage 4
depends on.

Skipped methods are reported by name and reason rather than quietly omitted.
"CIBERSORTx was not run because no token was configured" and "CIBERSORTx was not
run" are different statements and only one of them is reportable.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

from src.harness.deconvolve.base import Deconvolver, available_methods
from src.harness.pseudobulk import PseudobulkSample
from src.reference.signature import assert_no_target_leakage


def _bulk_on_signature(
    sample: PseudobulkSample, signature: pd.DataFrame, arm: str
) -> np.ndarray:
    """Project a sample's bulk vector onto the signature's gene index."""
    bulk = sample.bulk_tumour if arm == "tumour" else sample.bulk_normal
    series = pd.Series(bulk, index=list(sample.genes))
    missing = [g for g in signature.index if g not in series.index]
    if missing:
        raise ValueError(
            f"{len(missing)} signature gene(s) are absent from the sample, e.g. "
            f"{missing[:5]}. Both must sit on the shared gene index — that is "
            f"what makes integration a join and not a negotiation."
        )
    return series.reindex(signature.index).to_numpy(dtype=float)


def marker_ranked_genes(signature: pd.DataFrame) -> list[str]:
    """Signature genes ordered by cell-type specificity, most specific first.

    Score is ``max_c(profile) / mean_c(profile)``: a gene concentrated in one
    cell type scores high, a gene expressed evenly scores ~1. Used by
    :func:`signature_width_comparison` so that a narrow signature gets the
    *best* k genes rather than k random ones — otherwise the comparison is a
    strawman and proves nothing about width.
    """
    values = signature.to_numpy(dtype=float)
    means = values.mean(axis=1)
    means[means == 0] = np.finfo(float).eps
    specificity = values.max(axis=1) / means
    order = np.argsort(-specificity)
    return [signature.index[i] for i in order]


def run_bakeoff(
    samples: Sequence[PseudobulkSample],
    signature: pd.DataFrame,
    methods: Sequence[Deconvolver],
    *,
    seed: int,
    arm: str = "tumour",
    target_genes: Iterable[str] = (),
    alias_map: Mapping[str, str] | None = None,
    mature_label: str = "mature_colonocyte",
) -> tuple[pd.DataFrame, dict[str, str]]:
    """Run every available method over every sample. Returns ``(table, skipped)``.

    ``table`` is ``BAKEOFF_COLUMNS``-shaped: one row per
    (method, sample, cell type). ``skipped`` maps a method name to why it did
    not run — put it in the results sidecar.

    ``target_genes`` is asserted absent from the signature (invariant 2). Pass
    it; the check is cheap and the failure it prevents is a classifier that
    cannot detect the phenomenon it was built to detect.

    ``alias_map`` maps target SYMBOLS onto the identifier space the signature is
    indexed on, and is required whenever those differ. This call site was one of
    the four that issue #35 found comparing symbols against Ensembl ids — an
    intersection of two disjoint namespaces, empty whatever the data, reported as
    a pass. W1 fixed the guard so it now REFUSES that comparison rather than
    passing it; this passes the map so the check can actually run. A real
    Ensembl-indexed S matrix without a map raises ``LeakageGuardError``, which is
    the correct outcome and not a regression.
    """
    if target_genes:
        assert_no_target_leakage(
            signature.index,
            target_genes,
            context="the bake-off signature",
            alias_map=alias_map,
        )
    usable, skipped = available_methods(list(methods))
    if not usable:
        raise ValueError(f"no deconvolution method is available: {skipped}")

    rows: list[dict] = []
    for sample_id, sample in enumerate(samples):
        bulk = _bulk_on_signature(sample, signature, arm)
        truth = (
            sample.truth.composition_tumour
            if arm == "tumour"
            else sample.truth.composition_normal
        )
        for method in usable:
            hat = method.fit_predict(bulk, signature)
            for cell_type in signature.columns:
                rows.append(
                    {
                        "method": method.name,
                        "sample_id": sample_id,
                        "cell_type": cell_type,
                        "fraction_true": float(truth.get(cell_type, 0.0)),
                        "fraction_hat": float(hat.get(cell_type, 0.0)),
                        "n_signature_genes": int(signature.shape[0]),
                        "seed": seed,
                    }
                )
    table = pd.DataFrame(rows)
    table.attrs["mature_label"] = mature_label
    return table, skipped


def rank_methods(
    bakeoff: pd.DataFrame, *, by: str = "rmse", mature_label: str | None = None
) -> pd.DataFrame:
    """Rank on Pearson r and RMSE against harness truth, best first.

    ``mature_label`` additionally reports the error on that one cell type. The
    overall RMSE averages across compartments, and a method can look respectable
    there while being poor on precisely the fraction this project needs.
    """
    groups = ["method", "n_signature_genes"]
    out = []
    for key, g in bakeoff.groupby(groups):
        true, hat = g["fraction_true"].to_numpy(), g["fraction_hat"].to_numpy()
        r = float(np.corrcoef(true, hat)[0, 1]) if true.std() and hat.std() else float("nan")
        row = {
            "method": key[0],
            "n_signature_genes": key[1],
            "r": r,
            "rmse": float(np.sqrt(np.mean((hat - true) ** 2))),
            "bias": float(np.mean(hat - true)),
            "n_samples": int(g["sample_id"].nunique()),
        }
        if mature_label is not None:
            m = g[g["cell_type"] == mature_label]
            row["rmse_mature"] = (
                float(np.sqrt(np.mean((m["fraction_hat"] - m["fraction_true"]) ** 2)))
                if len(m)
                else float("nan")
            )
        out.append(row)
    ascending = by in {"rmse", "rmse_mature", "bias"}
    return pd.DataFrame(out).sort_values(by, ascending=ascending).reset_index(drop=True)


def signature_width_comparison(
    samples: Sequence[PseudobulkSample],
    signature: pd.DataFrame,
    method: Deconvolver,
    widths: Sequence[int],
    *,
    seed: int,
    target_genes: Iterable[str] = (),
    alias_map: Mapping[str, str] | None = None,
    mature_label: str = "mature_colonocyte",
) -> pd.DataFrame:
    """Fraction recovery as a function of signature width. Settles §2.1 error #4.

    The Executive Brief proposed deconvolving on the 11-gene panel. That claim
    cannot be tested literally: the panel genes are targets, and invariant 2
    keeps them out of any reference matrix. So the narrow arm here is the
    **best 11 non-target markers** by :func:`marker_ranked_genes`, not the panel
    and not 11 random genes.

    That is a stronger test than the original, and worth a line in the write-up:
    if even the most discriminative 11 genes underperform 500, the claim is
    settled by dimensionality rather than by which genes were picked.
    """
    ranked = marker_ranked_genes(signature)
    frames = []
    for width in widths:
        if width > len(ranked):
            raise ValueError(
                f"width={width} exceeds the {len(ranked)} genes in the signature"
            )
        narrowed = signature.loc[ranked[:width]]
        table, _ = run_bakeoff(
            samples, narrowed, [method], seed=seed,
            arm="tumour", target_genes=target_genes, alias_map=alias_map,
            mature_label=mature_label,
        )
        frames.append(table)
    return pd.concat(frames, ignore_index=True)
