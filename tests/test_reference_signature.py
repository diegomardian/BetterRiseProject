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


class TestSparsePath:
    """The dense path allocates cells x genes — 8.3 GB at pilot scale, which is
    where the first S matrix attempt died. The sparse path must give the same
    answer without ever materialising that array."""

    #: The S matrix is Ensembl-keyed while the panel is written as symbols, so
    #: the invariant-2 guard needs the translation or it compares two spaces and
    #: passes without testing anything (issue #35).
    ALIAS = {"GUCA2A": "ENSG00000197273"}

    def _sparse_cohort(self, n_cells=600, n_genes=2000):
        from scipy import sparse

        rng = np.random.default_rng(7)
        genes = [f"ENSG{i:08d}" for i in range(n_genes)]
        types = list(TYPES)
        labels = [types[i % len(types)] for i in range(n_cells)]

        dense = rng.poisson(0.4, size=(n_cells, n_genes)).astype(np.float32)
        block = n_genes // len(types)
        for index, cell_type in enumerate(types):
            rows = [i for i, t in enumerate(labels) if t == cell_type]
            dense[np.ix_(rows, range(index * block, (index + 1) * block))] += 15.0
        return sparse.csr_matrix(dense), genes, labels

    def test_it_produces_a_signature_without_densifying(self):
        from src.reference.signature import build_signature_sparse

        matrix, genes, labels = self._sparse_cohort()
        signature = build_signature_sparse(
            matrix, genes, labels,
            target_genes=["GUCA2A"], gene_index=genes, n_genes=600,
            alias_map=self.ALIAS,
        )
        assert MIN_SIGNATURE_GENES <= len(signature) <= MAX_SIGNATURE_GENES
        assert set(signature.columns) == set(TYPES)

    def test_it_matches_the_dense_path(self):
        """Same aggregates, so the same markers and the same profiles."""
        from src.reference.signature import build_signature_sparse, normalise_sparse

        matrix, genes, labels = self._sparse_cohort()
        normalised = normalise_sparse(matrix)
        dense = pd.DataFrame(np.asarray(normalised.todense()), columns=genes)

        sparse_signature = build_signature_sparse(
            matrix, genes, labels, target_genes=["GUCA2A"],
            gene_index=genes, n_genes=600, alias_map=self.ALIAS,
        )
        dense_signature = build_signature(
            dense, labels, target_genes=[self.ALIAS["GUCA2A"]],
            gene_index=genes, n_genes=600,
        )
        assert list(sparse_signature.index) == list(dense_signature.index)
        np.testing.assert_allclose(
            sparse_signature.to_numpy(), dense_signature.to_numpy(), rtol=1e-4
        )

    def test_normalisation_preserves_sparsity(self):
        """CP10K is a diagonal scaling and log1p(0) is 0, so nothing densifies."""
        from src.reference.signature import normalise_sparse

        matrix, _, _ = self._sparse_cohort()
        out = normalise_sparse(matrix)
        assert out.nnz == matrix.nnz

    def test_every_guard_still_runs(self):
        from src.reference.signature import (
            LeakageError,
            LeakageGuardError,
            build_signature_sparse,
        )

        matrix, genes, labels = self._sparse_cohort()
        # Symbol targets against an Ensembl index must REFUSE, not pass.
        with pytest.raises(LeakageGuardError, match="cannot check invariant 2"):
            build_signature_sparse(
                matrix, genes, labels, target_genes=["GUCA2A"],
                gene_index=genes, n_genes=600,
            )
        # Translated, it catches a real leak in Ensembl space.
        with pytest.raises(LeakageError, match="invariant 2"):
            build_signature_sparse(
                matrix, genes, labels, target_genes=["GUCA2A"],
                gene_index=genes, n_genes=600,
                alias_map={"GUCA2A": genes[0]},
            )
        with pytest.raises(ValueError, match="target_genes is empty"):
            build_signature_sparse(
                matrix, genes, labels, target_genes=[], gene_index=genes, n_genes=600
            )
        with pytest.raises(LeakageError, match="invariant 2"):
            build_signature_sparse(
                matrix, genes, labels, target_genes=[genes[0]],
                gene_index=genes, n_genes=600,
            )
        epithelial_only = ["mature_colonocyte"] * len(labels)
        with pytest.raises(ValueError, match="missing compartment"):
            build_signature_sparse(
                matrix, genes, epithelial_only, target_genes=["GUCA2A"],
                gene_index=genes, n_genes=600, alias_map=self.ALIAS,
            )

    def test_genes_off_the_shared_index_are_dropped(self):
        """W3 joins on this index; anything not on it cannot participate."""
        from src.reference.signature import build_signature_sparse

        matrix, genes, labels = self._sparse_cohort()
        index = genes[:1200]
        signature = build_signature_sparse(
            matrix, genes, labels, target_genes=["GUCA2A"],
            gene_index=index, n_genes=600, alias_map=self.ALIAS,
        )
        assert set(signature.index) <= set(index)

    def test_a_matrix_sharing_nothing_with_the_index_raises(self):
        from src.reference.signature import build_signature_sparse

        matrix, genes, labels = self._sparse_cohort()
        with pytest.raises(ValueError, match="unversioned Ensembl"):
            build_signature_sparse(
                matrix, genes, labels, target_genes=["GUCA2A"],
                gene_index=["NOTHING_MATCHES"], n_genes=600,
            )
