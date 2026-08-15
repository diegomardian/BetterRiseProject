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


def test_build_signature_rejects_a_target_in_the_gene_index():
    index = _clean_index() + ["CA7"]
    with pytest.raises(LeakageError, match="invariant 2"):
        build_signature(
            _expression(index),
            ["stromal"] * 50,
            target_genes=TARGETS,
            gene_index=index,
        )


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


def test_marker_selection_is_reached_only_after_every_guard_passes():
    """The scaffold owns the guard rails; W1 owns the biology. With a clean
    index and all compartments present, we should land on W1's TODO."""
    index = _clean_index()
    labels = ["mature_colonocyte", "stem", "stromal", "immune", "endothelial"] * 10
    with pytest.raises(NotImplementedError, match="W1 owns marker selection"):
        build_signature(
            _expression(index),
            labels,
            target_genes=TARGETS,
            gene_index=index,
        )


def test_whole_panel_can_be_held_out_at_once():
    """The panel is largely one co-regulated program — holding out just the one
    gene under test is not enough. Passing the full panel must work."""
    index = _clean_index()
    assert_no_target_leakage(index, panel_genes(), context="test")
    with pytest.raises(LeakageError):
        assert_no_target_leakage(index + ["CDX2"], panel_genes(), context="test")
