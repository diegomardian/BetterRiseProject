"""Marker selection and the S matrix. W1, weeks 4-5.

The design choices under test are the ones W2's bake-off will be interpreted
against: a per-cell-type quota rather than a global ranking, a detection floor,
and deterministic tie-breaking.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.reference.signature import (
    MAX_SIGNATURE_GENES,
    MIN_DETECTION_RATE,
    MIN_SIGNATURE_GENES,
    _select_markers,
    build_signature,
)

RNG = np.random.default_rng(20260820)

TYPES = ("mature_colonocyte", "stem", "stromal", "immune", "endothelial")
N_PER_TYPE = 60
N_GENES = 1500


def cohort(n_marker_genes: int = 60):
    """Each cell type has its own private marker block, plus shared background."""
    genes = [f"G{i:05d}" for i in range(N_GENES)]
    matrix = RNG.gamma(1.0, 0.5, size=(len(TYPES) * N_PER_TYPE, N_GENES))
    labels: list[str] = []
    private: dict[str, list[str]] = {}
    for index, cell_type in enumerate(TYPES):
        rows = slice(index * N_PER_TYPE, (index + 1) * N_PER_TYPE)
        block = slice(index * n_marker_genes, (index + 1) * n_marker_genes)
        matrix[rows, block] += 6.0
        private[cell_type] = genes[block]
        labels += [cell_type] * N_PER_TYPE
    return pd.DataFrame(matrix, columns=genes), labels, private


class TestMarkerSelection:
    def test_returns_the_requested_number(self):
        expression, labels, _ = cohort()
        markers = _select_markers(expression, labels, n_genes=800)
        assert len(markers) == 800

    def test_every_cell_type_gets_representation(self):
        """A global top-N is dominated by the most distinct compartments, and the
        rare columns end up with almost no genes — which is exactly where
        fraction estimates then become unidentifiable."""
        expression, labels, private = cohort()
        markers = set(_select_markers(expression, labels, n_genes=600))
        for cell_type, genes in private.items():
            hit = len(markers & set(genes))
            assert hit > 0, f"{cell_type} contributed no markers"

    def test_representation_is_roughly_balanced_across_types(self):
        expression, labels, private = cohort()
        markers = set(_select_markers(expression, labels, n_genes=500))
        counts = [len(markers & set(genes)) for genes in private.values()]
        assert min(counts) >= 0.4 * max(counts), counts

    def test_the_true_markers_are_preferred_over_background(self):
        expression, labels, private = cohort()
        markers = set(_select_markers(expression, labels, n_genes=500))
        true_markers = set().union(*private.values())
        assert len(markers & true_markers) / len(markers) > 0.5

    def test_a_gene_below_the_detection_floor_is_excluded(self):
        """Log fold change alone rewards a gene seen in three cells with a large
        count. That is noise in the reference and worse in the bulk."""
        expression, labels, _ = cohort()
        rare = "G00000"
        expression[rare] = 0.0
        # Present in only two cells of one type, but enormous there.
        expression.loc[expression.index[:2], rare] = 500.0
        markers = _select_markers(expression, labels, n_genes=600)
        assert rare not in markers

    def test_a_gene_just_above_the_floor_is_eligible(self):
        expression, labels, _ = cohort()
        gene = "G00001"
        expression[gene] = 0.0
        n_on = int(np.ceil(MIN_DETECTION_RATE * N_PER_TYPE)) + 2
        expression.loc[expression.index[:n_on], gene] = 50.0
        assert gene in _select_markers(expression, labels, n_genes=900)

    def test_selection_is_deterministic(self):
        """Invariant 10 depends on it: the same input must give the same
        signature, ties included."""
        expression, labels, _ = cohort()
        first = _select_markers(expression, labels, n_genes=700)
        second = _select_markers(expression, labels, n_genes=700)
        assert first == second

    def test_no_duplicates(self):
        expression, labels, _ = cohort()
        markers = _select_markers(expression, labels, n_genes=800)
        assert len(markers) == len(set(markers))

    def test_a_single_cell_type_is_refused(self):
        expression, labels, _ = cohort()
        with pytest.raises(ValueError, match="at least two cell types"):
            _select_markers(expression, ["only_one"] * len(labels), n_genes=600)

    def test_requesting_fewer_than_the_minimum_is_refused(self):
        """A distinct failure from 'not enough genes cleared the floor' — the
        first message conflated them."""
        expression, labels, _ = cohort()
        with pytest.raises(ValueError, match="below the .* minimum"):
            _select_markers(expression, labels, n_genes=11)

    def test_a_thin_signature_raises_rather_than_shipping(self):
        """nu-SVR robustness comes from high dimensionality (§2.1 error #4).
        Too few eligible genes must fail loudly, not quietly ship."""
        expression, labels, _ = cohort()
        sparse = expression.copy()
        sparse.iloc[:, 100:] = 0.0        # only 100 genes ever detected
        with pytest.raises(ValueError, match="below the"):
            _select_markers(sparse, labels, n_genes=600)


class TestSMatrix:
    def _index(self, expression):
        return list(expression.columns)

    def test_shape_and_columns(self):
        expression, labels, _ = cohort()
        signature = build_signature(
            expression, labels,
            target_genes=["GUCA2A", "MLH1"], gene_index=self._index(expression),
            n_genes=800,
        )
        assert set(signature.columns) == set(TYPES)
        assert MIN_SIGNATURE_GENES <= len(signature) <= MAX_SIGNATURE_GENES

    def test_the_profile_is_the_per_type_mean(self):
        expression, labels, private = cohort()
        signature = build_signature(
            expression, labels,
            target_genes=["GUCA2A"], gene_index=self._index(expression), n_genes=600,
        )
        # A type's own private markers should be highest in its own column.
        for cell_type, genes in private.items():
            present = [g for g in genes if g in signature.index]
            if not present:
                continue
            column = signature.loc[present].mean()
            assert column.idxmax() == cell_type

    def test_rows_follow_the_shared_gene_index_order(self):
        """§3.4: integration is a join, not a negotiation. Order must be the
        index's, not the selection's."""
        expression, labels, _ = cohort()
        index = self._index(expression)
        signature = build_signature(
            expression, labels,
            target_genes=["GUCA2A"], gene_index=index, n_genes=600,
        )
        positions = [index.index(g) for g in signature.index]
        assert positions == sorted(positions)

    def test_missing_a_compartment_is_refused(self):
        """Bulk CRC is 30-60% non-epithelial; omitting these columns is the CMS4
        failure mode (§2.1 error #3)."""
        expression, labels, _ = cohort()
        epithelial_only = ["mature_colonocyte" if x != "stem" else "stem" for x in labels]
        with pytest.raises(ValueError, match="missing compartment"):
            build_signature(
                expression, epithelial_only,
                target_genes=["GUCA2A"], gene_index=self._index(expression),
                n_genes=600,
            )
