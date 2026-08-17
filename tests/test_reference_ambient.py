"""Ambient contamination measured without empty droplets. W1.

GSE178341 ships no empty droplets (docs/open_decisions.md #8), so the estimator
under test here is the method-independent replacement: contamination read off
genes the cells cannot possibly transcribe. The recovery tests inject soup at a
known rate and check it comes back.

Note the fixture shape. The sample must contain the populations that *do* make
the impossible genes — erythrocytes, leukocytes — because they are what puts
those genes into the pooled soup profile. A single-compartment fixture makes the
observed rate and the expected rate the same number and the estimator degenerates
to 1.0.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from src.common.panel import panel_genes
from src.reference.ambient import (
    IMPOSSIBLE_GENES,
    AmbientError,
    contamination_by_sample,
    contamination_fraction,
    run_decontx,
    run_soupx,
    soup_profile_from_cells,
)
from src.reference.signature import LeakageError

#: HBB/HBA1/PTPRC are impossible in epithelium; the rest is generic signal.
GENES = ["HBB", "HBA1", "PTPRC", "EPCAM", "KRT8", "GENE1", "GENE2", "GENE3"]
IMPOSSIBLE_IDX = [0, 1, 2]
NATIVE_IDX = [3, 4, 5, 6, 7]
IMPOSSIBLE = {"HBB", "HBA1", "PTPRC"}

#: The true ambient profile. Dominated by the impossible genes, as it would be
#: with abundant lysing erythrocytes and leukocytes.
SOUP = np.zeros(len(GENES))
SOUP[IMPOSSIBLE_IDX] = [0.35, 0.20, 0.15]
SOUP[NATIVE_IDX] = [0.10, 0.08, 0.06, 0.04, 0.02]

#: Native epithelial expression — exactly zero at the impossible genes.
NATIVE_EPI = np.zeros(len(GENES))
NATIVE_EPI[NATIVE_IDX] = [0.4, 0.3, 0.1, 0.1, 0.1]


def _misassigned(n_epi: int = 200, n_other: int = 4000, depth: int = 20_000):
    """A sample whose "epithelial" cells are in fact immune cells.

    The compartment call is wrong, so the masked cells genuinely transcribe the
    impossible genes and the raw ratio goes to ~1. This is what C107_T_1_1_0's
    exact 1.000 was really telling us — a violated assumption, not 100% ambient.
    """
    matrix, mask = build(0.05, n_epi=n_epi, n_other=n_other, depth=depth)
    matrix[mask] = np.rint(np.full(len(GENES), depth) * SOUP).astype(np.int64)
    return matrix, mask


def build(rho: float, n_epi: int = 200, n_other: int = 4000, depth: int = 20_000):
    """A mixed sample: epithelium contaminated at `rho`, plus the cells that
    actually make the impossible genes and so feed the soup.

    Returns ``(matrix, epithelial_mask)``.
    """
    epi_profile = (1 - rho) * NATIVE_EPI + rho * SOUP
    epi = np.rint(np.outer(np.full(n_epi, depth), epi_profile)).astype(np.int64)
    other = np.rint(np.outer(np.full(n_other, depth), SOUP)).astype(np.int64)

    matrix = np.vstack([epi, other])
    mask = np.zeros(matrix.shape[0], dtype=bool)
    mask[:n_epi] = True
    return matrix, mask


class TestSoupProfile:
    def test_profile_sums_to_one(self):
        matrix, _ = build(0.1)
        profile = soup_profile_from_cells(matrix, GENES)
        assert profile.sum() == pytest.approx(1.0)
        assert list(profile.index) == GENES

    def test_profile_recovers_the_true_soup(self):
        """Pooled over all compartments, the tissue average approximates the soup."""
        matrix, _ = build(0.05)
        profile = soup_profile_from_cells(matrix, GENES)
        assert np.allclose(profile.to_numpy(), SOUP, atol=0.02)

    def test_empty_matrix_raises(self):
        with pytest.raises(AmbientError, match="no counts"):
            soup_profile_from_cells(np.zeros((5, len(GENES))), GENES)

    def test_sparse_matches_dense(self):
        matrix, _ = build(0.2)
        dense = soup_profile_from_cells(matrix, GENES)
        spars = soup_profile_from_cells(sparse.csr_matrix(matrix), GENES)
        pd.testing.assert_series_equal(dense, spars)


class TestContaminationRecovery:
    """Inject a known contamination rate; get it back."""

    @pytest.mark.parametrize("rho", [0.02, 0.05, 0.10, 0.25, 0.40])
    def test_known_rate_is_recovered(self, rho):
        matrix, mask = build(rho)
        estimate = contamination_fraction(
            matrix, GENES, cell_mask=mask, impossible=IMPOSSIBLE
        )
        assert estimate == pytest.approx(rho, abs=0.02)

    def test_clean_data_estimates_near_zero(self):
        matrix, mask = build(0.0)
        estimate = contamination_fraction(
            matrix, GENES, cell_mask=mask, impossible=IMPOSSIBLE
        )
        assert estimate < 0.01

    def _estimate(self, rho: float, **build_kwargs) -> float:
        matrix, mask = build(rho, **build_kwargs)
        return contamination_fraction(matrix, GENES, cell_mask=mask, impossible=IMPOSSIBLE)

    def test_estimate_increases_with_contamination(self):
        assert self._estimate(0.30) > self._estimate(0.05)

    def test_bias_is_upward_when_the_tested_population_dominates(self):
        """Documented property: pooling includes the tested cells, diluting the
        denominator at exactly the genes the ratio depends on."""
        minority = self._estimate(0.1, n_epi=100, n_other=8000)
        majority = self._estimate(0.1, n_epi=4000, n_other=1000)
        assert majority > minority
        assert minority == pytest.approx(0.1, abs=0.01)

    def test_result_is_bounded(self):
        matrix, mask = build(0.4)
        estimate = contamination_fraction(matrix, GENES, cell_mask=mask, impossible={"HBB"})
        assert 0.0 <= estimate <= 1.0

    def test_sparse_input_agrees(self):
        matrix, mask = build(0.15)
        dense = contamination_fraction(matrix, GENES, cell_mask=mask, impossible=IMPOSSIBLE)
        spars = contamination_fraction(
            sparse.csr_matrix(matrix), GENES, cell_mask=mask, impossible=IMPOSSIBLE
        )
        assert dense == pytest.approx(spars, abs=1e-9)


class TestInvariantTwo:
    """Using a target gene to measure contamination would be circular."""

    def test_panel_gene_in_the_impossible_set_raises(self):
        matrix, mask = build(0.1)
        with pytest.raises(LeakageError, match="GUCA2A"):
            contamination_fraction(
                matrix, GENES, cell_mask=mask, impossible={"HBB", "GUCA2A"}
            )

    def test_the_builtin_sets_are_panel_clean(self):
        targets = set(panel_genes())
        for compartment, genes in IMPOSSIBLE_GENES.items():
            assert not (targets & genes), f"{compartment} set collides with the panel"


class TestGuards:
    def test_unknown_compartment_raises(self):
        matrix, mask = build(0.1)
        with pytest.raises(AmbientError, match="no impossible-gene set"):
            contamination_fraction(matrix, GENES, cell_mask=mask, compartment="glia")

    def test_mask_length_mismatch_raises(self):
        matrix, _ = build(0.1)
        with pytest.raises(AmbientError, match="entries for"):
            contamination_fraction(
                matrix, GENES, cell_mask=np.ones(7, dtype=bool), impossible={"HBB"}
            )

    def test_empty_mask_raises(self):
        matrix, mask = build(0.1)
        with pytest.raises(AmbientError, match="selects no cells"):
            contamination_fraction(
                matrix, GENES, cell_mask=np.zeros_like(mask), impossible={"HBB"}
            )

    def test_genes_absent_from_the_matrix_raise_with_a_naming_hint(self):
        matrix, mask = build(0.1)
        with pytest.raises(AmbientError, match="symbols vs Ensembl"):
            contamination_fraction(
                matrix, GENES, cell_mask=mask, impossible={"ENSG00000244734"}
            )

    def test_cells_expressing_an_impossible_gene_refuse_rather_than_saturate(self):
        """The C107_T_1_1_0 case: it returned exactly 1.000.

        A ratio at or above 1 says the masked cells carry as much impossible-gene
        signal as pure soup would. That is not contamination, it is a violated
        assumption — doublets, or cells misassigned to this compartment. Here the
        violation is simulated by giving the "epithelial" cells real HBB.
        """
        matrix, mask = _misassigned()
        with pytest.raises(AmbientError, match="assumption is violated"):
            contamination_fraction(matrix, GENES, cell_mask=mask, impossible=IMPOSSIBLE)

    def test_the_guard_can_be_relaxed_deliberately(self):
        matrix, mask = _misassigned()
        value = contamination_fraction(
            matrix, GENES, cell_mask=mask, impossible=IMPOSSIBLE, max_plausible=10.0
        )
        assert value > 0.9

    def test_a_healthy_mixture_passes_the_guard(self):
        matrix, mask = build(0.05)
        assert contamination_fraction(
            matrix, GENES, cell_mask=mask, impossible=IMPOSSIBLE
        ) == pytest.approx(0.05, abs=0.02)

    def test_impossible_genes_absent_from_the_soup_raise(self):
        """Not identifiable if the marker genes are nowhere in the tissue."""
        matrix, mask = build(0.1)
        matrix[:, IMPOSSIBLE_IDX] = 0
        with pytest.raises(AmbientError, match="not identifiable"):
            contamination_fraction(
                matrix, GENES, cell_mask=mask, impossible={"HBB", "HBA1"}
            )


class TestPerSample:
    def _two_samples(self, rho_a: float, rho_b: float):
        mat_a, mask_a = build(rho_a, n_epi=200, n_other=2000)
        mat_b, mask_b = build(rho_b, n_epi=200, n_other=2000)
        matrix = np.vstack([mat_a, mat_b])
        mask = np.concatenate([mask_a, mask_b])
        sample_id = ["s_clean"] * len(mask_a) + ["s_dirty"] * len(mask_b)
        return matrix, mask, sample_id

    def test_degenerate_sample_reports_nan_not_a_number(self):
        """Per-sample, an unestimable sample must come back NaN — the same
        None-is-not-zero distinction invariant 1 makes for the intrinsic term."""
        mat_ok, mask_ok = build(0.05, n_epi=200, n_other=2000)
        mat_bad, mask_bad = _misassigned(n_epi=200, n_other=2000)
        matrix = np.vstack([mat_ok, mat_bad])
        mask = np.concatenate([mask_ok, mask_bad])
        sample_id = ["s_ok"] * len(mask_ok) + ["s_degenerate"] * len(mask_bad)

        report = contamination_by_sample(
            matrix, GENES, sample_id=sample_id, cell_mask=mask
        ).set_index("sample_id")
        assert np.isnan(report.loc["s_degenerate", "contamination"])
        assert not np.isnan(report.loc["s_ok", "contamination"])
        assert report.loc["s_ok", "mask_share"] == pytest.approx(200 / 2200, abs=0.01)

    def test_one_row_per_sample_with_differing_rates(self):
        matrix, mask, sample_id = self._two_samples(0.02, 0.30)
        report = contamination_by_sample(
            matrix, GENES, sample_id=sample_id, cell_mask=mask
        )
        assert list(report["sample_id"]) == ["s_clean", "s_dirty"]
        assert (report["n_cells"] == 200).all()
        rates = report.set_index("sample_id")["contamination"]
        assert rates["s_dirty"] > rates["s_clean"]
        assert rates["s_clean"] < 0.05
        assert rates["s_dirty"] == pytest.approx(0.30, abs=0.03)

    def test_sample_with_no_masked_cells_reports_nan(self):
        matrix, mask, sample_id = self._two_samples(0.1, 0.1)
        mask[len(mask) // 2 :] = False
        report = contamination_by_sample(
            matrix, GENES, sample_id=sample_id, cell_mask=mask
        )
        assert np.isnan(report.set_index("sample_id").loc["s_dirty", "contamination"])

    def test_length_mismatch_raises(self):
        matrix, mask = build(0.1)
        with pytest.raises(AmbientError, match="entries for"):
            contamination_by_sample(matrix, GENES, sample_id=["a"], cell_mask=mask)


@pytest.mark.parametrize("fn,match", [(run_soupx, "soup profile"), (run_decontx, "celda")])
def test_correction_entry_points_are_explicit_todos(fn, match):
    with pytest.raises(NotImplementedError, match=match):
        fn()
