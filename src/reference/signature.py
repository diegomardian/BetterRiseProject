"""Reference (S) matrix construction. W1.

CLAUDE.md invariant 2 is enforced here: target genes never appear in labels or
the reference matrix. If GUCA2A is in the reference, a silenced mature cell is
readable as an absent mature cell, and the classifier cannot detect the
phenomenon it was built to detect. That is Executive-Brief error #1 and it is
the reason this function asserts rather than warns.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence

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
    if n_genes < MIN_SIGNATURE_GENES:
        raise ValueError(
            f"n_genes={n_genes} is below the {MIN_SIGNATURE_GENES} minimum. "
            f"nu-SVR robustness comes from high dimensionality; the 11-gene panel "
            f"is for interpretation, not deconvolution (execution_plan.md §2.1 "
            f"error #4)."
        )

    labels = pd.Series(list(cell_type), index=expression.index, name="cell_type")
    types = sorted(labels.unique())
    if len(types) < 2:
        raise ValueError(
            f"need at least two cell types to select discriminative markers, got "
            f"{types}. A one-column reference cannot separate anything."
        )

    detected = expression.gt(0)
    means = expression.groupby(labels, observed=True).mean()
    rates = detected.groupby(labels, observed=True).mean()

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
