"""Ambient contamination measured without empty droplets. W1.

GSE178341 ships no empty droplets (docs/open_decisions.md #11), so the estimator
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
    compare_retention,
    contamination_by_sample,
    contamination_fraction,
    retention_agreement,
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


class TestDecontX:
    """The second method, and it exists because CellBender cannot run on a
    deposit with no empty droplets (#8)."""

    def _inputs(self, n_cells=60, n_genes=30, clusters=2):
        rng = np.random.default_rng(1)
        matrix = rng.poisson(3, size=(n_cells, n_genes)).astype(np.int64)
        return (
            matrix,
            [f"G{i}" for i in range(n_genes)],
            [f"c{i}" for i in range(n_cells)],
            [f"k{i % clusters}" for i in range(n_cells)],
        )

    def test_dry_run_writes_its_inputs(self, tmp_path):
        matrix, genes, barcodes, labels = self._inputs()
        out = run_decontx(matrix, genes, barcodes=barcodes, clusters=labels,
                          out_dir=tmp_path, dry_run=True)
        assert out["ran"] is False
        for name in ("matrix.mtx", "genes.tsv", "barcodes.tsv",
                     "clusters.tsv", "run_decontx.R"):
            assert (tmp_path / name).exists(), name

    def test_it_asks_for_no_empty_droplets(self, tmp_path):
        """The entire reason it replaces CellBender."""
        matrix, genes, barcodes, labels = self._inputs()
        out = run_decontx(matrix, genes, barcodes=barcodes, clusters=labels,
                          out_dir=tmp_path, dry_run=True)
        script = out["script"].read_text()
        assert "decontX(" in script
        # The ARGUMENT, not the word — the template's comment explains why it
        # is absent, and a bare substring check catches its own explanation.
        assert "background =" not in script
        assert "background=" not in script

    def test_it_emits_retention_and_per_cell_contamination(self, tmp_path):
        matrix, genes, barcodes, labels = self._inputs()
        script = run_decontx(
            matrix, genes, barcodes=barcodes, clusters=labels,
            out_dir=tmp_path, dry_run=True,
        )["script"].read_text()
        assert "decontx_retention.csv" in script
        assert "decontx_contamination.csv" in script

    def test_one_cluster_refuses(self, tmp_path):
        """Contamination is defined as counts resembling OTHER clusters, so one
        cluster leaves nothing to compare against."""
        matrix, genes, barcodes, _ = self._inputs()
        with pytest.raises(AmbientError, match="other clusters|OTHER clusters"):
            run_decontx(matrix, genes, barcodes=barcodes,
                        clusters=["k0"] * matrix.shape[0], out_dir=tmp_path,
                        dry_run=True)

    def test_cluster_length_is_checked(self, tmp_path):
        matrix, genes, barcodes, _ = self._inputs()
        with pytest.raises(AmbientError, match="entries for"):
            run_decontx(matrix, genes, barcodes=barcodes, clusters=["k0"] * 3,
                        out_dir=tmp_path, dry_run=True)

    def test_it_is_seeded(self, tmp_path):
        """Invariant 10 — every result carries a fixed seed."""
        matrix, genes, barcodes, labels = self._inputs()
        script = run_decontx(
            matrix, genes, barcodes=barcodes, clusters=labels,
            out_dir=tmp_path, seed=4242, dry_run=True,
        )["script"].read_text()
        assert "set.seed(4242)" in script


class TestSoupX:
    """Degraded mode: no empty droplets exist for this deposit (#8), so the
    profile comes from cells via setSoupProfile()."""

    def _inputs(self, n_cells=60, n_genes=30, clusters=2):
        rng = np.random.default_rng(0)
        matrix = rng.poisson(3, size=(n_cells, n_genes)).astype(np.int64)
        genes = [f"G{i}" for i in range(n_genes)]
        barcodes = [f"c{i}" for i in range(n_cells)]
        labels = [f"k{i % clusters}" for i in range(n_cells)]
        profile = pd.Series(1.0 / n_genes, index=genes)
        return matrix, genes, barcodes, labels, profile

    def test_dry_run_writes_every_input_sparse(self, tmp_path):
        matrix, genes, barcodes, labels, profile = self._inputs()
        out = run_soupx(
            matrix, genes, barcodes=barcodes, clusters=labels,
            soup_profile=profile, out_dir=tmp_path, dry_run=True,
        )
        assert out["ran"] is False
        for name in ("matrix.mtx", "genes.tsv", "barcodes.tsv",
                     "clusters.tsv", "soup_profile.csv", "run_soupx.R"):
            assert (tmp_path / name).exists(), name

    def test_the_R_uses_degraded_mode(self, tmp_path):
        """calcSoupProfile = FALSE plus setSoupProfile is the only route
        without empty droplets."""
        matrix, genes, barcodes, labels, profile = self._inputs()
        out = run_soupx(
            matrix, genes, barcodes=barcodes, clusters=labels,
            soup_profile=profile, out_dir=tmp_path, dry_run=True,
        )
        script = out["script"].read_text()
        assert "calcSoupProfile = FALSE" in script
        assert "setSoupProfile" in script
        assert "estimateSoup" not in script

    def test_it_writes_retention_not_a_corrected_matrix(self, tmp_path):
        """Decision #16 is to measure, not correct — and 62 corrected matrices
        would not fit on the project filesystem."""
        matrix, genes, barcodes, labels, profile = self._inputs()
        out = run_soupx(
            matrix, genes, barcodes=barcodes, clusters=labels,
            soup_profile=profile, out_dir=tmp_path, dry_run=True,
        )
        script = out["script"].read_text()
        assert "soupx_retention.csv" in script
        assert "adjustCounts" in script

    def test_one_cluster_refuses_unless_rho_is_fixed(self, tmp_path):
        """autoEstCont needs marker genes, which needs more than one cluster."""
        matrix, genes, barcodes, _labels, profile = self._inputs(clusters=1)
        labels = ["k0"] * matrix.shape[0]
        with pytest.raises(AmbientError, match="more than one cluster"):
            run_soupx(matrix, genes, barcodes=barcodes, clusters=labels,
                      soup_profile=profile, out_dir=tmp_path, dry_run=True)
        out = run_soupx(
            matrix, genes, barcodes=barcodes, clusters=labels,
            soup_profile=profile, out_dir=tmp_path, contamination=0.02,
            dry_run=True,
        )
        assert "0.02" in out["script"].read_text()

    def test_duplicate_symbols_in_the_profile_are_summed(self, tmp_path):
        """soup_profile_from_cells is indexed by gene SYMBOL, and this deposit
        maps several Ensembl IDs onto one symbol. reindex() refuses a duplicated
        index outright — this failed on the first real sample."""
        matrix, genes, barcodes, labels = self._inputs()[:4]
        profile = pd.Series(
            [0.02] * (len(genes) + 3),
            index=list(genes) + genes[:3],   # three symbols duplicated
        )
        run_soupx(matrix, genes, barcodes=barcodes, clusters=labels,
                  soup_profile=profile, out_dir=tmp_path, dry_run=True)
        written = pd.read_csv(tmp_path / "soup_profile.csv", index_col=0)
        assert not written.index.has_duplicates
        # Summed, not dropped — the profile is a share of one soup.
        assert written.loc[genes[0], "est"] == pytest.approx(0.04)

    def test_an_empty_profile_refuses(self, tmp_path):
        matrix, genes, barcodes, labels, _profile = self._inputs()
        with pytest.raises(AmbientError, match="soup profile is empty"):
            run_soupx(
                matrix, genes, barcodes=barcodes, clusters=labels,
                soup_profile=pd.Series(dtype=float), out_dir=tmp_path,
                dry_run=True,
            )

    def test_cluster_length_is_checked(self, tmp_path):
        matrix, genes, barcodes, _labels, profile = self._inputs()
        with pytest.raises(AmbientError, match="entries for"):
            run_soupx(matrix, genes, barcodes=barcodes, clusters=["k0"] * 3,
                      soup_profile=profile, out_dir=tmp_path, dry_run=True)


class TestRetentionComparison:
    """The week-2 deliverable is a comparison, not a winner."""

    def _pair(self, decontx_scale=2.0, n=40, shuffle=False):
        genes = [f"G{i}" for i in range(n)]
        soupx = pd.DataFrame({
            "gene": genes,
            "retention": [1 - 0.01 * i for i in range(n)],
        })
        order = list(range(n))[::-1] if shuffle else list(range(n))
        decontx = pd.DataFrame({
            "gene": genes,
            "retention": [1 - 0.01 * decontx_scale * order[i] for i in range(n)],
        })
        return soupx, decontx

    def test_agreement_on_which_genes_is_ranked_not_absolute(self):
        """DecontX strips twice as hard here, but ranks genes identically. That
        is a different kind of agreement from ranking them differently, and
        Spearman is what separates the two."""
        soupx, decontx = self._pair(decontx_scale=2.0)
        out = retention_agreement(compare_retention(soupx, decontx, sample_id="S"))
        assert out["spearman"] > 0.99
        assert out["agree"]
        # ...and the magnitude difference survives as its own number.
        assert out["median_difference"] > 0

    def test_opposite_rankings_do_not_agree(self):
        soupx, decontx = self._pair(shuffle=True)
        out = retention_agreement(compare_retention(soupx, decontx, sample_id="S"))
        assert out["spearman"] < 0
        assert not out["agree"]

    def test_only_shared_genes_are_compared(self):
        soupx, decontx = self._pair()
        decontx = decontx.iloc[:10]
        out = compare_retention(soupx, decontx, sample_id="S")
        assert len(out) == 10

    def test_no_shared_genes_refuses(self):
        soupx, decontx = self._pair()
        decontx = decontx.assign(gene=[f"X{i}" for i in range(len(decontx))])
        with pytest.raises(AmbientError, match="no genes in common"):
            compare_retention(soupx, decontx)

    def test_missing_columns_refuse(self):
        soupx, decontx = self._pair()
        with pytest.raises(AmbientError, match="missing column"):
            compare_retention(soupx.drop(columns=["retention"]), decontx)

    def test_too_few_genes_to_correlate_refuses(self):
        soupx, decontx = self._pair(n=2)
        with pytest.raises(AmbientError, match="too few"):
            retention_agreement(compare_retention(soupx, decontx))

    def test_the_sample_is_carried_through(self):
        """Per sample, never pooled — the soup belongs to one dissociation."""
        soupx, decontx = self._pair()
        out = compare_retention(soupx, decontx, sample_id="C122_N_1_1_0_c1_v2")
        assert set(out["sample_id"]) == {"C122_N_1_1_0_c1_v2"}
