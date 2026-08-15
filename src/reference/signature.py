"""Reference (S) matrix construction. W1.

CLAUDE.md invariant 2 is enforced here: target genes never appear in labels or
the reference matrix. If GUCA2A is in the reference, a silenced mature cell is
readable as an absent mature cell, and the classifier cannot detect the
phenomenon it was built to detect. That is Executive-Brief error #1 and it is
the reason this function asserts rather than warns.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence

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


def _select_markers(
    expression: pd.DataFrame,
    cell_type: Sequence[str],
    *,
    n_genes: int,
) -> list[str]:
    """Pick n_genes discriminative markers. W1 — unimplemented.

    Suggested starting point: per-cell-type one-vs-rest differential expression,
    take the top-k per type, union, then cap at n_genes. Whatever you choose,
    write down why in the docstring — the bake-off in W2 will be interpreted
    against this choice.
    """
    raise NotImplementedError(
        "W1 owns marker selection — see src/reference/README.md, week 4-5."
    )
