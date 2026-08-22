"""Reference (S) matrix construction. W1.

CLAUDE.md invariant 2 is enforced here: target genes never appear in labels or
the reference matrix. If GUCA2A is in the reference, a silenced mature cell is
readable as an absent mature cell, and the classifier cannot detect the
phenomenon it was built to detect. That is Executive-Brief error #1 and it is
the reason this function asserts rather than warns.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from typing import Any

import numpy as np
import pandas as pd

#: execution_plan.md §2.1 error #4 — nu-SVR robustness comes from high
#: dimensionality. The 11-gene panel is for interpretation, never deconvolution.
MIN_SIGNATURE_GENES = 500
MAX_SIGNATURE_GENES = 2000

#: Bulk CRC is 30-60% non-epithelial. Without these columns, stromal signal is
#: absorbed arbitrarily into the epithelial ones — the CMS4 failure mode.
REQUIRED_NON_EPITHELIAL = ("stromal", "immune", "endothelial")


class LeakageError(AssertionError):
    """A target gene reached the reference matrix or the labels."""


def assert_no_target_leakage(
    genes: Iterable[str],
    target_genes: Iterable[str],
    *,
    context: str,
) -> None:
    """Hard stop if any target gene appears in ``genes``. CLAUDE.md invariant 2.

    Call this from label construction too, not only from build_signature.
    """
    targets = set(target_genes)
    leaked = sorted(targets & set(genes))
    if leaked:
        raise LeakageError(
            f"{len(leaked)} target gene(s) leaked into {context}: {leaked}. "
            f"Target genes never appear in labels or the reference matrix "
            f"(CLAUDE.md invariant 2) — a silenced mature cell must not be "
            f"readable as an absent mature cell."
        )


def build_signature(
    expression: pd.DataFrame,
    cell_type: Sequence[str],
    *,
    target_genes: Iterable[str],
    gene_index: Sequence[str],
    n_genes: int = 1000,
    require_non_epithelial: bool = True,
) -> pd.DataFrame:
    """Build an S matrix: genes (rows) x cell types (columns), on the fixed index.

    Parameters
    ----------
    expression:
        cells x genes, log-normalised.
    cell_type:
        Per-cell label at one granularity rung. Built from a labelling axis that
        excludes ``target_genes`` — assert that upstream with
        :func:`assert_no_target_leakage`, this function cannot see your labels.
    target_genes:
        The panel genes under test for this run. Excluded unconditionally.
    gene_index:
        The fixed shared index (config/gene_index/). W3 emits bulk on the same
        index so integration is a join, not a negotiation.

    Notes
    -----
    Marker selection itself is W1's call and is left unimplemented here — the
    scaffold owns the guard rails, not the biology. Fill in ``_select_markers``.
    """
    targets = set(target_genes)
    if not targets:
        raise ValueError(
            "target_genes is empty. Pass the panel genes under test — an empty "
            "target set silently disables invariant 2."
        )
    if not MIN_SIGNATURE_GENES <= n_genes <= MAX_SIGNATURE_GENES:
        raise ValueError(
            f"n_genes={n_genes} outside [{MIN_SIGNATURE_GENES}, {MAX_SIGNATURE_GENES}]. "
            f"The 11-gene panel is for interpretation, not deconvolution "
            f"(execution_plan.md §2.1 error #4)."
        )
    if len(cell_type) != len(expression):
        raise ValueError(f"cell_type has {len(cell_type)} entries for {len(expression)} cells")

    assert_no_target_leakage(gene_index, targets, context="the shared gene index")

    usable = [g for g in expression.columns if g in set(gene_index)]
    assert_no_target_leakage(usable, targets, context="the reference gene pool")

    types = sorted(set(cell_type))
    if require_non_epithelial:
        missing = [c for c in REQUIRED_NON_EPITHELIAL if not any(c in t.lower() for t in types)]
        if missing:
            raise ValueError(
                f"reference is missing compartment(s) {missing}. Bulk CRC is "
                f"30-60% non-epithelial; without these columns stromal signal is "
                f"absorbed arbitrarily (execution_plan.md §2.1 error #3). Pass "
                f"require_non_epithelial=False only for single-cell-only work."
            )

    markers = _select_markers(expression.loc[:, usable], cell_type, n_genes=n_genes)
    assert_no_target_leakage(markers, targets, context="the selected marker set")

    profiles = (
        expression.loc[:, markers]
        .assign(_ct=list(cell_type))
        .groupby("_ct", observed=True)
        .mean()
        .T
    )
    profiles.index.name = "gene"
    profiles.columns.name = "cell_type"

    # Reindex onto the fixed index so every S matrix is joinable with every other.
    profiles = profiles.reindex([g for g in gene_index if g in profiles.index])
    assert_no_target_leakage(profiles.index, targets, context="the emitted S matrix")
    return profiles


#: A gene detected in fewer than this share of a cell type's cells is not a
#: marker for it, however large the fold change looks. Guards against a handful
#: of cells with a big count setting the reference profile for a whole column.
MIN_DETECTION_RATE = 0.10

#: Pseudocount for the log fold change, so a gene absent from the rest of the
#: cells does not produce an infinite score and monopolise the ranking.
LFC_PSEUDOCOUNT = 1e-3


def _select_markers(
    expression: pd.DataFrame,
    cell_type: Sequence[str],
    *,
    n_genes: int,
) -> list[str]:
    """Pick `n_genes` discriminative markers: one-vs-rest, quota per cell type.

    Why this and not something else
    -------------------------------
    Deconvolution recovers cell fractions by asking which mixture of reference
    columns explains a bulk profile, so the signature has to *separate the
    columns*. That makes one-vs-rest differential expression the natural
    criterion — a gene earns its place by being high in one column and low in the
    others — and it is what CIBERSORT-style signatures use.

    Three choices inside it are worth stating, because W2's bake-off will be
    interpreted against them (execution_plan.md §4):

    1. **A per-cell-type quota, not a global ranking.** A global top-N is
       dominated by the most abundant and most transcriptionally distinct
       compartments, so the rare columns get few genes and their fractions
       become the least identifiable. In this cohort epithelium and T cells
       would swamp endothelium and mast cells. Each type gets an equal quota
       first; only leftovers are filled by global rank.
    2. **A detection floor** (`MIN_DETECTION_RATE`). Log fold change alone
       rewards a gene detected in three cells of a type with a large count. Such
       a gene is noise in the reference profile and worse in the bulk, where it
       is diluted. A gene must be seen in at least 10% of the type's cells.
    3. **Mean of log, not log of mean.** `expression` arrives log-normalised;
       averaging on that scale weights cells equally rather than letting the
       deepest cells set the profile.

    What it does not do: no HVG pre-filter, no variance stabilisation, no
    marker-list prior. Those would each be defensible, and each would need its
    own justification against the bake-off. This is the plain version, chosen so
    that when W2 compares deconvolution methods the signature is not itself a
    confounder.

    Determinism: ties break on gene name, so the same input gives the same
    signature. Invariant 10 depends on it.
    """
    labels = pd.Series(list(cell_type), index=expression.index, name="cell_type")
    means = expression.groupby(labels, observed=True).mean()
    rates = expression.gt(0).groupby(labels, observed=True).mean()
    return select_markers_from_aggregates(means, rates, n_genes=n_genes)


def select_markers_from_aggregates(
    means: pd.DataFrame,
    rates: pd.DataFrame,
    *,
    n_genes: int,
) -> list[str]:
    """The ranking itself, given per-type means and detection rates.

    Split out so the dense and sparse paths share one implementation — the sparse
    path computes the same two aggregates without ever materialising a
    cells-by-genes array.
    """
    if n_genes < MIN_SIGNATURE_GENES:
        raise ValueError(
            f"n_genes={n_genes} is below the {MIN_SIGNATURE_GENES} minimum. "
            f"nu-SVR robustness comes from high dimensionality; the 11-gene panel "
            f"is for interpretation, not deconvolution (execution_plan.md §2.1 "
            f"error #4)."
        )
    types = sorted(means.index)
    if len(types) < 2:
        raise ValueError(
            f"need at least two cell types to select discriminative markers, got "
            f"{types}. A one-column reference cannot separate anything."
        )

    # One-vs-rest on the log scale: the type's mean against the mean of the
    # other types, each type weighted equally so a large type does not define
    # "the rest" on its own.
    scores: dict[str, pd.Series] = {}
    for cell_type_name in types:
        rest = means.drop(index=cell_type_name).mean(axis=0)
        lfc = means.loc[cell_type_name] - rest + LFC_PSEUDOCOUNT
        eligible = rates.loc[cell_type_name] >= MIN_DETECTION_RATE
        scores[cell_type_name] = lfc.where(eligible, other=-np.inf)

    quota = max(1, n_genes // len(types))
    chosen: list[str] = []
    seen: set[str] = set()
    for cell_type_name in types:
        ranked = scores[cell_type_name]
        ranked = ranked[np.isfinite(ranked)]
        # Sort by score, then by name, so ties are resolved identically every run.
        order = sorted(ranked.index, key=lambda g: (-ranked[g], g))
        for gene in order:
            if len(chosen) >= n_genes:
                break
            if gene not in seen:
                chosen.append(gene)
                seen.add(gene)
            if sum(1 for g in chosen if g in set(order[:quota])) >= quota:
                break

    if len(chosen) < n_genes:
        # Fill the remainder by best score across any type, still deterministic.
        best = pd.DataFrame(scores).max(axis=1)
        best = best[np.isfinite(best)]
        for gene in sorted(best.index, key=lambda g: (-best[g], g)):
            if len(chosen) >= n_genes:
                break
            if gene not in seen:
                chosen.append(gene)
                seen.add(gene)

    if len(chosen) < min(n_genes, MIN_SIGNATURE_GENES):
        raise ValueError(
            f"only {len(chosen)} of the requested {n_genes} genes cleared the "
            f"detection floor across {len(types)} cell types, below the "
            f"{MIN_SIGNATURE_GENES} minimum. "
            f"nu-SVR robustness comes from high dimensionality "
            f"(execution_plan.md §2.1 error #4) — lower MIN_DETECTION_RATE or "
            f"widen the gene pool rather than shipping a thin signature."
        )
    return chosen


def _group_aggregates(
    matrix: Any,
    gene_names: Sequence[str],
    cell_type: Sequence[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per-cell-type mean expression and detection rate, without densifying.

    A pilot-sized matrix is 25,959 x 43,113, which is 8.3 GB dense in float64 and
    simply will not allocate. Both aggregates are group sums, so an indicator
    matrix does the work sparsely: ``G @ X`` for sums, ``G @ (X > 0)`` for
    detection counts, where ``G`` is types-by-cells.

    Returns ``(means, rates)``, both types x genes.
    """
    import scipy.sparse as sp

    labels = pd.Series([str(c) for c in cell_type], name="cell_type")
    types = sorted(labels.unique())
    if matrix.shape[0] != len(labels):
        raise ValueError(
            f"matrix has {matrix.shape[0]} cells but cell_type has {len(labels)}"
        )

    position = {t: i for i, t in enumerate(types)}
    rows = labels.map(position).to_numpy()
    indicator = sp.csr_matrix(
        (np.ones(len(rows)), (rows, np.arange(len(rows)))),
        shape=(len(types), len(rows)),
    )
    per_type = np.asarray(indicator.sum(axis=1)).ravel().reshape(-1, 1)

    sparse_matrix = matrix if sp.issparse(matrix) else sp.csr_matrix(matrix)
    sums = np.asarray((indicator @ sparse_matrix).todense())
    detected = np.asarray((indicator @ (sparse_matrix > 0)).todense())

    means = pd.DataFrame(sums / per_type, index=types, columns=list(gene_names))
    rates = pd.DataFrame(detected / per_type, index=types, columns=list(gene_names))
    return means, rates


def normalise_sparse(matrix: Any, *, target_sum: float = 1e4) -> Any:
    """CP10K then log1p, preserving sparsity.

    Both steps map zero to zero — row scaling is diagonal, and log1p(0) is 0 — so
    the matrix never densifies.
    """
    import scipy.sparse as sp

    sparse_matrix = sp.csr_matrix(matrix, dtype=np.float32)
    totals = np.asarray(sparse_matrix.sum(axis=1)).ravel()
    scale = np.divide(
        target_sum, totals, out=np.zeros_like(totals), where=totals > 0
    )
    scaled = sp.diags(scale) @ sparse_matrix
    scaled.data = np.log1p(scaled.data)
    return scaled.tocsr()


def build_signature_sparse(
    matrix: Any,
    gene_names: Sequence[str],
    cell_type: Sequence[str],
    *,
    target_genes: Iterable[str],
    gene_index: Sequence[str],
    n_genes: int = 1000,
    require_non_epithelial: bool = True,
    already_normalised: bool = False,
) -> pd.DataFrame:
    """:func:`build_signature` for a sparse matrix. Same guards, no densification.

    The dense path allocates cells x genes, which is 8.3 GB at pilot scale and
    grows with the cohort. This aggregates sparsely and only materialises the
    selected markers — 800 genes rather than 43,113.

    Every guard from the dense path runs, in the same order: the target set may
    not be empty, `n_genes` stays inside [500, 2000], the shared gene index and
    the reference pool are checked for target leakage, the non-epithelial
    compartments must be present, and the emitted matrix is checked again.
    """
    targets = set(target_genes)
    if not targets:
        raise ValueError(
            "target_genes is empty. Pass the panel genes under test — an empty "
            "target set silently disables invariant 2."
        )
    if not MIN_SIGNATURE_GENES <= n_genes <= MAX_SIGNATURE_GENES:
        raise ValueError(
            f"n_genes={n_genes} outside [{MIN_SIGNATURE_GENES}, {MAX_SIGNATURE_GENES}]. "
            f"The 11-gene panel is for interpretation, not deconvolution "
            f"(execution_plan.md §2.1 error #4)."
        )
    if len(cell_type) != matrix.shape[0]:
        raise ValueError(
            f"cell_type has {len(cell_type)} entries for {matrix.shape[0]} cells"
        )

    assert_no_target_leakage(gene_index, targets, context="the shared gene index")

    names = [str(g) for g in gene_names]
    index_set = set(gene_index)
    keep = np.array([g in index_set for g in names], dtype=bool)
    if not keep.any():
        raise ValueError(
            "no gene in the matrix appears on the shared gene index. Check the "
            "identifier form — the index is unversioned Ensembl ids "
            "(open decision #3)."
        )
    usable = [g for g, k in zip(names, keep, strict=True) if k]
    assert_no_target_leakage(usable, targets, context="the reference gene pool")

    types = sorted({str(c) for c in cell_type})
    if require_non_epithelial:
        missing = [c for c in REQUIRED_NON_EPITHELIAL if not any(c in t.lower() for t in types)]
        if missing:
            raise ValueError(
                f"reference is missing compartment(s) {missing}. Bulk CRC is "
                f"30-60% non-epithelial; without these columns stromal signal is "
                f"absorbed arbitrarily (execution_plan.md §2.1 error #3)."
            )

    subset = matrix[:, keep]
    if not already_normalised:
        subset = normalise_sparse(subset)

    means, rates = _group_aggregates(subset, usable, cell_type)
    markers = select_markers_from_aggregates(means, rates, n_genes=n_genes)
    assert_no_target_leakage(markers, targets, context="the selected marker set")

    profiles = means.loc[:, markers].T
    profiles.index.name = "gene"
    profiles.columns.name = "cell_type"
    profiles = profiles.reindex([g for g in gene_index if g in profiles.index])
    assert_no_target_leakage(profiles.index, targets, context="the emitted S matrix")
    return profiles
