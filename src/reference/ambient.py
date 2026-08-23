"""Ambient RNA: measuring it, and correcting it. W1, week 2.

GSE178341 as deposited has had its empty droplets removed (docs/open_decisions.md
#8), so the usual route to a soup profile — average the cell-free droplets — is
unavailable. CellBender cannot run at all. SoupX runs only if handed a profile
estimated some other way, and DecontX infers contamination from cluster
structure.

That makes the *measurement* in this module more important than the correction.
:func:`contamination_fraction` estimates ambient contamination **without any
empty droplets**, using genes that cannot be expressed in the cells you are
looking at. If a colonocyte shows haemoglobin counts, no colonocyte made them —
they came from the soup, and their size relative to the soup's own composition
is the contamination rate.

Why that matters here specifically. Both surviving correction methods infer the
ambient profile from the cells themselves, which risks absorbing the very signal
this project measures: a rare mature cell's genuine low-level GUCA2A is exactly
what a cluster-based method may reassign to "contamination from elsewhere". An
estimate built from impossible genes is independent of that machinery, so it can
audit the correction rather than inherit its assumptions. It also feeds gate
criterion G1 directly.

Matrix orientation is **cells x genes**, as in ``ingest.py`` and ``qc.py``.
"""

from __future__ import annotations

from typing import Any, Final

import numpy as np
import pandas as pd

from src.common.panel import panel_genes
from src.reference.signature import assert_no_target_leakage

#: Genes that a cell of the given compartment cannot transcribe. Residual counts
#: in these are ambient by construction.
#:
#: Chosen to be lineage-exclusive and highly expressed in *some* population, so
#: they are well represented in the soup and give a stable ratio. Haemoglobin is
#: the workhorse — erythrocytes are abundant, lyse readily, and no epithelial
#: cell makes HBB.
#:
#: None of these may be panel genes; asserted at call time (CLAUDE.md invariant
#: 2). Using a target gene to measure contamination would be circular.
#: The epithelial set is deliberately broad. On the pilot, the original six
#: (haemoglobin, PTPRC, IGKC, IGHG1) held under 0.2% of the soup in 10 of 23
#: samples, because their share is dominated by how many erythrocytes and plasma
#: cells that particular dissociation happened to contain — so those samples came
#: back unestimable. The pan-leukocyte transcripts below are present in every
#: sample and abundant, which stabilises the denominator.
#:
#: Chosen for being **constitutively** haematopoietic. Deliberately excluded:
#: CD74 and HLA-DRA (interferon-inducible in epithelium), and the fibroblast
#: collagens (COL1A1/COL1A2/COL3A1), which are abundant and tempting but can be
#: switched on in partial-EMT tumour cells — exactly the state this project
#: studies, so using them would put the measurement inside the phenomenon.
IMPOSSIBLE_GENES: Final[dict[str, frozenset[str]]] = {
    "epithelial": frozenset(
        {
            # erythroid
            "HBB", "HBA1", "HBA2",
            # plasma / B
            "IGKC", "IGHG1", "MS4A1",
            # pan-leukocyte, constitutive and highly expressed
            "PTPRC", "SRGN", "LAPTM5", "CORO1A", "CD52", "ARHGDIB",
            # T
            "CD3D", "CD3E", "CD2",
            # myeloid
            "TYROBP", "FCER1G", "AIF1",
        }
    ),
    "immune": frozenset({"HBB", "HBA1", "HBA2", "EPCAM", "KRT8", "KRT18", "KRT19"}),
    "stromal": frozenset({"HBB", "HBA1", "HBA2", "EPCAM", "PTPRC", "KRT8"}),
    "endothelial": frozenset({"HBB", "HBA1", "HBA2", "EPCAM", "PTPRC", "KRT8"}),
}


class AmbientError(ValueError):
    """Contamination could not be estimated from the input given."""


def _row_sums(matrix: Any) -> np.ndarray:
    return np.asarray(matrix.sum(axis=1)).ravel()


def _gene_positions(gene_names: Any, wanted: Any) -> np.ndarray:
    names = list(gene_names)
    lookup = {str(n): i for i, n in enumerate(names)}
    return np.array([lookup[g] for g in wanted if g in lookup], dtype=int)


def soup_profile_from_cells(matrix: Any, gene_names: Any) -> pd.Series:
    """Per-gene share of total counts, pooled over all cells.

    A stand-in for the empty-droplet profile. The soup is largely mRNA from
    cells that lysed during dissociation, so the tissue's overall composition
    approximates it — imperfectly, because it also contains the real biological
    signal we are trying to measure.

    This is what you pass to SoupX's ``setSoupProfile()`` when there is no raw
    droplet table to call ``estimateSoup()`` on. Compute it **per sample**: the
    soup is a property of one dissociation, not of the study.

    Pass the **whole sample, every compartment** — never a single population.
    The impossible-gene estimator divides the observed rate by this profile's
    share of those genes, and if the profile is built only from the cells being
    tested, both terms are the same number and the ratio degenerates to 1.0.
    The erythrocytes and leukocytes that actually transcribe HBB and PTPRC are
    what put those genes in the profile at all.

    Because the pooled average includes the tested cells, the estimate is biased
    slightly **upward** in proportion to how much of the sample they represent —
    they dilute the profile at exactly the genes the ratio depends on. Negligible
    when the tested population is a minority; worth stating when it is not.
    """
    totals = np.asarray(matrix.sum(axis=0)).ravel().astype(float)
    grand = totals.sum()
    if grand <= 0:
        raise AmbientError("matrix has no counts; cannot build a soup profile")
    return pd.Series(totals / grand, index=list(gene_names), name="soup_fraction")


def contamination_fraction(
    matrix: Any,
    gene_names: Any,
    *,
    cell_mask: Any,
    impossible: Any = None,
    compartment: str = "epithelial",
    soup: pd.Series | None = None,
    min_soup_share: float = 0.002,
    max_plausible: float = 0.95,
) -> float:
    """Ambient contamination in the masked cells, from impossible genes alone.

    The estimator, which is the one SoupX uses for its marker-based route::

        observed = counts of impossible genes in these cells / their total counts
        expected = share of those genes in the soup, if the cell were pure soup
        rho      = observed / expected

    A cell of this compartment cannot make these transcripts, so everything
    observed is ambient. Dividing by the soup's own share of them converts that
    into "what fraction of this cell's counts came from soup".

    Parameters
    ----------
    cell_mask:
        Boolean, selecting cells known to belong to `compartment`. Get it from a
        malignancy or lineage call — not from the panel.
    soup:
        Soup profile over the same genes. Defaults to
        :func:`soup_profile_from_cells` on the whole matrix, which is the
        approximation available without empty droplets. The default is
        deliberately the **whole matrix**, not the masked cells — see that
        function's note on why a single-compartment profile degenerates.
    min_soup_share:
        Refuse to estimate when the impossible genes hold less than this share of
        the soup — they cannot discriminate if they are barely in it.
    max_plausible:
        Refuse when the raw ratio reaches this. A ratio at or above 1 says the
        masked cells carry as much impossible-gene signal as pure soup would
        give, which is not contamination but a **violated assumption**: doublets,
        or cells misassigned to this compartment. Observed on C107_T_1_1_0, which
        returned exactly 1.000 with 543 epithelial cells in a 2,737-cell sample.
        Clipping that to 1.0 and calling it "100% ambient" would be a fabricated
        number; refusing is the honest answer.

    Returns a fraction in [0, 1]. Report it per sample and log it next to the
    correction — it is the number that says whether the correction was
    necessary and whether it went far enough.
    """
    genes = frozenset(impossible) if impossible is not None else IMPOSSIBLE_GENES.get(compartment)
    if not genes:
        raise AmbientError(
            f"no impossible-gene set for compartment {compartment!r}; "
            f"known: {sorted(IMPOSSIBLE_GENES)}"
        )

    # Using a target gene to measure contamination would be circular.
    assert_no_target_leakage(genes, panel_genes(), context="the impossible-gene set")

    mask = np.asarray(cell_mask, dtype=bool)
    if mask.shape[0] != matrix.shape[0]:
        raise AmbientError(f"cell_mask has {mask.shape[0]} entries for {matrix.shape[0]} cells")
    if not mask.any():
        raise AmbientError("cell_mask selects no cells")

    positions = _gene_positions(gene_names, sorted(genes))
    if positions.size == 0:
        raise AmbientError(
            f"none of the impossible genes {sorted(genes)} are in the matrix. "
            f"Check gene naming — symbols vs Ensembl IDs (open decision #3)."
        )

    if soup is None:
        soup = soup_profile_from_cells(matrix, gene_names)
    expected = float(np.asarray(soup)[positions].sum())
    if expected <= 0:
        raise AmbientError(
            "the impossible genes have zero share of the soup, so contamination "
            "is not identifiable from them. Pick genes that are abundant "
            "somewhere in the tissue."
        )
    if expected < min_soup_share:
        raise AmbientError(
            f"the impossible genes hold only {expected:.4%} of the soup, below the "
            f"{min_soup_share:.4%} floor. This happens when the masked population "
            f"dominates the sample: the pooled profile is then mostly those same "
            f"cells, so observed and expected converge and the ratio runs to 1.0 "
            f"regardless of the truth. Estimate the soup where the source "
            f"populations are actually present, or report this sample as not "
            f"estimable — do not report the 1.0."
        )

    selected = matrix[mask, :]
    observed_counts = float(np.asarray(selected[:, positions].sum()).ravel().sum())
    total_counts = float(_row_sums(selected).sum())
    if total_counts <= 0:
        raise AmbientError("masked cells have no counts")

    ratio = (observed_counts / total_counts) / expected
    if ratio >= max_plausible:
        raise AmbientError(
            f"raw ratio {ratio:.3f} means the masked cells hold as much "
            f"impossible-gene signal as pure soup would give. That is not a "
            f"contamination level — it says the impossible-gene assumption is "
            f"violated for these cells (doublets, or a wrong compartment call). "
            f"Report this sample as not estimable rather than as ~100% ambient."
        )
    return float(max(ratio, 0.0))


def contamination_by_sample(
    matrix: Any,
    gene_names: Any,
    *,
    sample_id: Any,
    cell_mask: Any,
    compartment: str = "epithelial",
    min_soup_share: float = 0.002,
) -> pd.DataFrame:
    """Per-sample contamination estimate. One row per sample.

    Per sample, not pooled: the soup is a property of a single dissociation, and
    a study-wide number would average away exactly the variation that decides
    which samples to distrust.
    """
    samples = pd.Series(list(sample_id), name="sample_id")
    if len(samples) != matrix.shape[0]:
        raise AmbientError(f"sample_id has {len(samples)} entries for {matrix.shape[0]} cells")
    mask = np.asarray(cell_mask, dtype=bool)

    rows: list[dict[str, Any]] = []
    for sample, index in samples.groupby(samples).groups.items():
        positions = np.asarray(index, dtype=int)
        in_sample = np.zeros(matrix.shape[0], dtype=bool)
        in_sample[positions] = True
        selected = in_sample & mask
        if not selected.any():
            rows.append(
                {"sample_id": sample, "n_cells": 0, "mask_share": 0.0,
                 "soup_share": np.nan, "contamination": np.nan}
            )
            continue
        # Soup estimated within the sample — see the docstring above.
        soup = soup_profile_from_cells(matrix[in_sample, :], gene_names)
        genes = frozenset(IMPOSSIBLE_GENES.get(compartment, ()))
        share = float(np.asarray(soup)[_gene_positions(gene_names, sorted(genes))].sum())
        row = {
            "sample_id": sample,
            "n_cells": int(selected.sum()),
            "mask_share": float(selected.sum() / max(int(in_sample.sum()), 1)),
            "soup_share": share,
        }
        try:
            row["contamination"] = contamination_fraction(
                matrix,
                gene_names,
                cell_mask=selected,
                compartment=compartment,
                soup=soup,
                min_soup_share=min_soup_share,
            )
        except AmbientError:
            # Not estimable in this sample. None, never a number — the same
            # distinction CLAUDE.md invariant 1 makes for the intrinsic term.
            row["contamination"] = np.nan
        rows.append(row)
    return pd.DataFrame(rows).sort_values("sample_id", ignore_index=True)


# ---------------------------------------------------------------------------
# Correction
# ---------------------------------------------------------------------------


def run_soupx(*args: Any, **kwargs: Any) -> Any:
    """SoupX in degraded mode. W1, week 2.

    No raw droplet table exists for GSE178341, so the standard
    ``estimateSoup()`` route is unavailable. Construct the channel with
    ``calcSoupProfile = FALSE`` and hand it a profile from
    :func:`soup_profile_from_cells`, computed **per sample**, via
    ``setSoupProfile()``. ``autoEstCont()`` then derives the contamination
    fraction from clusters and marker genes without touching empty droplets.

    Cross-check whatever ``autoEstCont()`` returns against
    :func:`contamination_fraction` — they estimate the same quantity by
    different routes, and disagreement is informative.
    """
    raise NotImplementedError(
        "W1 — SoupX needs the real matrices and a per-sample soup profile. "
        "See docs/open_decisions.md #11."
    )


def run_decontx(*args: Any, **kwargs: Any) -> Any:
    """DecontX (``bioconductor-celda``). W1, week 2. **The second method.**

    execution_plan.md line 167 originally asked for "SoupX and CellBender, both,
    compared". That was restated on 2026-08-22: CellBender requires unfiltered
    droplets and open decision #8 established that none exist in any public
    source for this deposit — GEO ran dropletUtils upstream. DecontX takes its
    place and the deliverable is unchanged, two methods compared.

    Needs no empty droplets — it models each cell as a mixture of its own
    cluster's distribution and contamination from the others.

    Caveat to carry into the write-up: because it defines contamination as
    counts resembling other clusters, it can absorb genuine low-level expression
    of a marker in a rare population. That is precisely this project's signal,
    which is why :func:`contamination_fraction` and gate criterion G1 exist as
    independent checks rather than confirmations.
    """
    raise NotImplementedError(
        "W1 — DecontX needs the real matrices and cluster assignments. "
        "bioconductor-celda is already pinned in env/w1_reference.yml."
    )
