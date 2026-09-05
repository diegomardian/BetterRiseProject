"""CLAUDE.md invariant 2 — target genes never reach the labels or the reference.

Executive Brief error #1: GUCA2A appeared in the reference matrix and was also
the target. A silenced mature cell must not be readable as an absent mature
cell, or the classifier cannot detect the phenomenon it was built to detect.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.common.panel import panel_genes, tier_genes
from src.reference.signature import (
    MAX_SIGNATURE_GENES,
    MIN_SIGNATURE_GENES,
    LeakageError,
    assert_no_target_leakage,
    build_signature,
)

TARGETS = tier_genes("A")  # GUCA2A, GUCA2B, OTOP2, CA7


def _clean_index(n=1200):
    return [f"GENE{i:05d}" for i in range(n)]


def _expression(index, n_cells=50):
    rng = np.random.default_rng(0)
    return pd.DataFrame(rng.normal(size=(n_cells, len(index))), columns=index)


def test_leakage_check_passes_on_a_clean_index():
    assert_no_target_leakage(_clean_index(), TARGETS, context="test")


def test_leakage_check_names_the_offending_genes():
    index = _clean_index() + ["GUCA2A"]
    with pytest.raises(LeakageError, match="GUCA2A"):
        assert_no_target_leakage(index, TARGETS, context="the reference gene pool")


def _mixed_cell_types(n=50):
    """A reference needs the non-epithelial compartments (§2.1 error #3)."""
    return [["stromal", "immune", "endothelial", "epithelial"][i % 4] for i in range(n)]


def test_build_signature_filters_a_target_out_of_the_gene_index():
    """Open decision #12. The shared index carries panel genes ON PURPOSE — they
    are W3's outcome variables, and decision #2 requires 23/23 panel coverage in
    the intersection. So build_signature filters them out of the reference pool
    rather than refusing the index. Invariant 2 is about what reaches the
    reference MATRIX, and it still holds: CA7 is not in the result.
    """
    index = _clean_index() + ["CA7"]
    signature = build_signature(
        _expression(index),
        _mixed_cell_types(),
        target_genes=TARGETS,
        gene_index=index,
    )
    assert "CA7" not in set(signature.index)
    assert not set(TARGETS) & set(signature.index)


def test_a_target_in_the_expression_matrix_still_never_reaches_the_signature():
    """The case that matters. Cells genuinely express CA7 — that is the whole
    point — and it is on the index. Neither fact may put it in the reference.
    """
    index = _clean_index() + ["CA7"]
    expression = _expression(index)
    expression["CA7"] = 100.0  # loudly expressed, top of any marker ranking
    signature = build_signature(
        expression, _mixed_cell_types(), target_genes=TARGETS, gene_index=index
    )
    assert "CA7" not in set(signature.index)


def test_build_signature_refuses_an_empty_target_set():
    """An empty target set silently disables the invariant. Loud failure instead."""
    index = _clean_index()
    with pytest.raises(ValueError, match="target_genes is empty"):
        build_signature(
            _expression(index), ["stromal"] * 50, target_genes=[], gene_index=index
        )


@pytest.mark.parametrize("n_genes", [11, 100, MIN_SIGNATURE_GENES - 1, MAX_SIGNATURE_GENES + 1])
def test_signature_size_bounds_are_enforced(n_genes):
    """nu-SVR robustness comes from high dimensionality. The 11-gene panel is
    for interpretation, not deconvolution (§2.1 error #4)."""
    index = _clean_index()
    with pytest.raises(ValueError, match="interpretation, not deconvolution"):
        build_signature(
            _expression(index),
            ["stromal"] * 50,
            target_genes=TARGETS,
            gene_index=index,
            n_genes=n_genes,
        )


def test_reference_must_carry_non_epithelial_compartments():
    """Bulk CRC is 30-60% non-epithelial. Without these columns stromal signal
    is absorbed arbitrarily — the CMS4 failure mode (§2.1 error #3)."""
    index = _clean_index()
    with pytest.raises(ValueError, match="missing compartment"):
        build_signature(
            _expression(index),
            ["mature_colonocyte", "stem"] * 25,
            target_genes=TARGETS,
            gene_index=index,
        )


def test_a_signature_is_produced_once_every_guard_passes():
    """The scaffold owned the guard rails; W1 has now supplied the biology.

    Previously this asserted we landed on a NotImplementedError. Marker
    selection is implemented, so the assertion becomes: all four leakage checks
    run, the compartment requirement is met, and a usable signature comes back.
    """
    index = _clean_index()
    labels = ["mature_colonocyte", "stem", "stromal", "immune", "endothelial"] * 10
    signature = build_signature(
        _expression(index),
        labels,
        target_genes=TARGETS,
        gene_index=index,
    )
    assert MIN_SIGNATURE_GENES <= len(signature) <= MAX_SIGNATURE_GENES
    assert set(signature.columns) == set(labels)
    # The emitted matrix is still checked against the target set.
    assert_no_target_leakage(signature.index, TARGETS, context="the emitted S matrix")


def test_whole_panel_can_be_held_out_at_once():
    """The panel is largely one co-regulated program — holding out just the one
    gene under test is not enough. Passing the full panel must work."""
    index = _clean_index()
    assert_no_target_leakage(index, panel_genes(), context="test")
    with pytest.raises(LeakageError):
        assert_no_target_leakage(index + ["CDX2"], panel_genes(), context="test")


class TestIdentifierSpaceHoles:
    """Found reviewing PR #33. The guard refused symbol-vs-unversioned-Ensembl
    but not the other two Ensembl forms this project actually handles."""

    ALIAS = {"GUCA2A": "ENSG00000197273", "GUCA2B": "ENSG00000044012",
             "OTOP2": "ENSG00000183034", "CA7": "ENSG00000168748"}

    @pytest.mark.parametrize("index", [
        ["ENSG00000197273.5"],        # versioned — how TCGA STAR counts arrive
        ["ENSG00000243485.5_4"],      # CellRanger-suffixed — the deposit's raw feature_id
        ["ENSG00000197273", "GUCA2B"],  # mixed — catches one target, misses the other
    ])
    def test_it_refuses_rather_than_passing_vacuously(self, index):
        from src.reference.signature import LeakageGuardError, assert_no_target_leakage
        with pytest.raises(LeakageGuardError, match="cannot check invariant 2"):
            assert_no_target_leakage(index, TARGETS, context="test")

    def test_a_mixed_index_is_checked_in_BOTH_forms_once_translated(self):
        """The subtle half: a mixed index can hold a target under either form,
        and checking only the translated form misses the symbol-written ones."""
        from src.reference.signature import LeakageError, assert_no_target_leakage
        with pytest.raises(LeakageError, match="invariant 2"):
            assert_no_target_leakage(
                ["ENSG00000197273", "NOT_A_TARGET"], TARGETS,
                context="test", alias_map=self.ALIAS,
            )
        with pytest.raises(LeakageError, match="invariant 2"):
            assert_no_target_leakage(
                ["GUCA2B", "ENSG00000999999"], TARGETS,
                context="test", alias_map=self.ALIAS,
            )

    def test_a_clean_ensembl_index_still_passes(self):
        from src.reference.signature import assert_no_target_leakage
        assert_no_target_leakage(
            ["ENSG00000999999"], TARGETS, context="test", alias_map=self.ALIAS
        )


# ---------------------------------------------------------------------------
# A2: the linear-scale profile. Selection must not move with it.


def _tiny_cohort(n_cells: int = 400, n_genes: int = 1200, seed: int = 5):
    """Raw counts with compartment structure and real within-type dispersion.

    Dispersion matters: the two profile scales differ by Jensen's gap, which is
    zero when every cell of a type is identical. A fixture without it would let
    a broken linear path pass.
    """
    import numpy as np
    import scipy.sparse as sp

    rng = np.random.default_rng(seed)
    types = ["differentiated", "stem_like", "immune", "stromal", "endothelial"]
    rows, labels = [], []
    for j, t in enumerate(types):
        base = rng.gamma(2.0, 3.0, n_genes)
        base[j * 200:(j + 1) * 200] *= 20.0
        for _ in range(n_cells // len(types)):
            rows.append(rng.poisson(base * rng.lognormal(0.0, 0.8)))
            labels.append(t)
    counts = sp.csr_matrix(np.vstack(rows).astype(np.float32))
    genes = [f"ENSG{i:08d}" for i in range(n_genes)]
    return counts, genes, labels


def test_the_linear_profile_selects_exactly_the_same_markers():
    """The one property that makes a rebuilt matrix attributable.

    If emitting linearly also changed WHICH genes are chosen, a 2.0.0 matrix
    would differ from 1.0.0 in two ways at once and no downstream difference
    could be assigned to either.
    """
    from src.reference.signature import build_signature_sparse

    counts, genes, labels = _tiny_cohort()
    common = dict(
        gene_names=genes, cell_type=labels, target_genes=["ENSG99999999"],
        gene_index=genes, n_genes=500,
    )
    log_profile = build_signature_sparse(counts, **common, profile_scale="log1p")
    linear_profile = build_signature_sparse(counts, **common, profile_scale="linear")

    assert list(log_profile.index) == list(linear_profile.index), (
        "the two scales chose different markers; selection must not move"
    )
    assert list(log_profile.columns) == list(linear_profile.columns)


def test_the_linear_profile_is_not_merely_expm1_of_the_log_one():
    """expm1(mean(log1p(x))) is a GEOMETRIC mean, biased low by Jensen.

    That is the repair `--linearise-reference` applies to the committed
    matrices, and it is an approximation. A properly rebuilt linear profile is
    the arithmetic mean and must exceed it.
    """
    import numpy as np

    from src.reference.signature import build_signature_sparse

    counts, genes, labels = _tiny_cohort()
    common = dict(
        gene_names=genes, cell_type=labels, target_genes=["ENSG99999999"],
        gene_index=genes, n_genes=500,
    )
    log_profile = build_signature_sparse(counts, **common, profile_scale="log1p")
    linear_profile = build_signature_sparse(counts, **common, profile_scale="linear")

    approximation = np.expm1(log_profile)
    exceeds = (linear_profile.to_numpy() >= approximation.to_numpy() - 1e-6)
    assert exceeds.mean() > 0.95, (
        f"the arithmetic mean falls below the geometric one on "
        f"{100 * (1 - exceeds.mean()):.1f}% of entries; Jensen says it cannot"
    )
    relative = (
        (linear_profile - approximation).abs() / linear_profile.replace(0, np.nan)
    ).stack().median()
    assert relative > 0.05, (
        f"expm1 of the log profile is within {relative:.2%} of the real linear "
        f"one, so the rebuild would not be worth doing -- re-derive A1's +0.166 "
        f"before relying on it"
    )


def test_the_log_profile_is_unchanged_by_the_new_parameter():
    """1.0.0 must stay byte-identical, or every table built on it moves."""
    from src.reference.signature import build_signature_sparse

    counts, genes, labels = _tiny_cohort()
    common = dict(
        gene_names=genes, cell_type=labels, target_genes=["ENSG99999999"],
        gene_index=genes, n_genes=500,
    )
    default = build_signature_sparse(counts, **common)
    explicit = build_signature_sparse(counts, **common, profile_scale="log1p")
    pd.testing.assert_frame_equal(default, explicit)


def test_an_unknown_profile_scale_is_refused():
    from src.reference.signature import build_signature_sparse

    counts, genes, labels = _tiny_cohort()
    with pytest.raises(ValueError, match="profile_scale must be"):
        build_signature_sparse(
            counts, gene_names=genes, cell_type=labels,
            target_genes=["ENSG99999999"], gene_index=genes, n_genes=500,
            profile_scale="sqrt",
        )


def test_linear_is_refused_when_the_caller_already_normalised():
    """Renormalising an already-normalised matrix would rescale twice."""
    from src.reference.signature import build_signature_sparse

    counts, genes, labels = _tiny_cohort()
    with pytest.raises(ValueError, match="needs raw counts"):
        build_signature_sparse(
            counts, gene_names=genes, cell_type=labels,
            target_genes=["ENSG99999999"], gene_index=genes, n_genes=500,
            already_normalised=True, profile_scale="linear",
        )


def test_the_linear_profile_still_excludes_targets():
    """Invariant 2 does not weaken because the scale changed."""
    from src.reference.signature import build_signature_sparse

    counts, genes, labels = _tiny_cohort()
    target = genes[10]
    profile = build_signature_sparse(
        counts, gene_names=genes, cell_type=labels, target_genes=[target],
        gene_index=genes, n_genes=500, profile_scale="linear",
    )
    assert target not in set(profile.index)
