"""Week-1 QC. W1.

The deliverable under test is the *thresholds table*, not a filtered matrix:
execution_plan.md §4 asks for "QC thresholds documented with rationale". So the
tests that matter are the ones asserting thresholds are computed per batch, that
the table records why, and that a shallow batch is not filtered against a deep
batch's distribution.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from src.reference.qc import (
    ABSOLUTE_MIN_COUNTS,
    ABSOLUTE_MIN_GENES,
    DEFAULT_MAX_PCT_MITO,
    DEFAULT_N_MADS,
    METRIC_COLUMNS,
    QCError,
    apply_qc,
    cell_qc_metrics,
    flag_doublets,
    qc_summary,
    qc_thresholds,
)

RNG = np.random.default_rng(20260816)

#: Small panel for the targeted metric tests, where exact percentages matter.
GENES = ["EPCAM", "PTPRC", "MT-CO1", "MT-ND1", "RPS6", "RPL13", "LGR5", "MUC2"]

#: Realistic panel for the threshold tests. The absolute floors (200 genes,
#: 500 counts) are calibrated for real cells, so a toy 8-gene fixture fails
#: every one of them — which is the floors working, not a bug.
BIG_GENES = ["MT-CO1", "MT-ND1", "RPS6", "RPL13"] + [f"G{i:04d}" for i in range(700)]


def counts(n_cells=100, scale=1.0):
    return (RNG.poisson(20.0 * scale, size=(n_cells, len(GENES)))).astype(np.int64)


def realistic_counts(n_cells=60, scale=1.0):
    """~600 genes detected, ~1400 UMI per cell — comfortably above the floors."""
    return RNG.poisson(2.0 * scale, size=(n_cells, len(BIG_GENES))).astype(np.int64)


def metrics_frame(n_per_batch=60, batches=("s1", "s2")):
    frames = []
    for i, name in enumerate(batches):
        matrix = realistic_counts(n_per_batch, scale=1.0 + i)
        frames.append(
            cell_qc_metrics(
                matrix,
                BIG_GENES,
                batch=[name] * n_per_batch,
                patient_id=[f"P{i + 1}"] * n_per_batch,
                tissue=(["tumour"] * (n_per_batch // 2) + ["normal"] * (n_per_batch // 2)),
            )
        )
    return pd.concat(frames, ignore_index=True)


class TestMetrics:
    def test_columns_and_length(self):
        m = cell_qc_metrics(counts(30), GENES, batch=["s1"] * 30)
        assert set(METRIC_COLUMNS) <= set(m.columns)
        assert len(m) == 30

    def test_mito_percentage_is_computed_from_mt_prefixed_genes(self):
        matrix = np.zeros((2, len(GENES)), dtype=np.int64)
        matrix[0, GENES.index("EPCAM")] = 75
        matrix[0, GENES.index("MT-CO1")] = 25
        matrix[1, GENES.index("MT-CO1")] = 10
        matrix[1, GENES.index("MT-ND1")] = 10
        m = cell_qc_metrics(matrix, GENES, batch=["s1", "s1"])
        assert m.loc[0, "pct_mito"] == pytest.approx(25.0)
        assert m.loc[1, "pct_mito"] == pytest.approx(100.0)

    def test_ribo_percentage_is_reported(self):
        matrix = np.zeros((1, len(GENES)), dtype=np.int64)
        matrix[0, GENES.index("RPS6")] = 30
        matrix[0, GENES.index("EPCAM")] = 70
        m = cell_qc_metrics(matrix, GENES, batch=["s1"])
        assert m.loc[0, "pct_ribo"] == pytest.approx(30.0)

    def test_genes_detected_counts_nonzero(self):
        matrix = np.zeros((1, len(GENES)), dtype=np.int64)
        matrix[0, :3] = [5, 1, 2]
        assert cell_qc_metrics(matrix, GENES, batch=["s1"]).loc[0, "n_genes"] == 3

    def test_sparse_and_dense_agree(self):
        matrix = counts(40)
        dense = cell_qc_metrics(matrix, GENES, batch=["s1"] * 40)
        spars = cell_qc_metrics(sparse.csr_matrix(matrix), GENES, batch=["s1"] * 40)
        pd.testing.assert_frame_equal(dense, spars)

    def test_zero_count_cell_does_not_divide_by_zero(self):
        matrix = np.zeros((1, len(GENES)), dtype=np.int64)
        m = cell_qc_metrics(matrix, GENES, batch=["s1"])
        assert m.loc[0, "pct_mito"] == 0.0

    def test_gene_name_length_mismatch_raises(self):
        with pytest.raises(QCError, match="cells x genes"):
            cell_qc_metrics(counts(10), GENES[:3], batch=["s1"] * 10)


class TestThresholds:
    """The week-1 deliverable: thresholds, per batch, with a written rationale."""

    def test_one_row_per_batch_and_metric(self):
        t = qc_thresholds(metrics_frame())
        assert set(t["batch"]) == {"s1", "s2"}
        assert set(t["metric"]) == {"n_counts", "n_genes", "pct_mito"}
        assert len(t) == 6

    def test_every_row_carries_a_rationale_and_a_method(self):
        t = qc_thresholds(metrics_frame())
        assert t["rationale"].str.len().gt(40).all()
        assert t["method"].str.len().gt(5).all()

    def test_thresholds_differ_between_batches_of_different_depth(self):
        """The whole point of per-batch: a deep batch gets a higher bound."""
        t = qc_thresholds(metrics_frame(n_per_batch=200))
        counts_rule = t[t["metric"] == "n_counts"].set_index("batch")
        assert counts_rule.loc["s2", "upper"] > counts_rule.loc["s1", "upper"]

    def test_absolute_floor_wins_over_a_permissive_mad_bound(self):
        """A uniformly poor batch must not keep near-empty barcodes."""
        poor = np.full((80, len(GENES)), 1, dtype=np.int64)
        m = cell_qc_metrics(poor, GENES, batch=["bad"] * 80)
        t = qc_thresholds(m)
        assert t.loc[t["metric"] == "n_genes", "lower"].iloc[0] == ABSOLUTE_MIN_GENES
        assert t.loc[t["metric"] == "n_counts", "lower"].iloc[0] == ABSOLUTE_MIN_COUNTS

    def test_mito_cap_is_hard_and_shared_across_batches(self):
        t = qc_thresholds(metrics_frame())
        mito = t[t["metric"] == "pct_mito"]
        assert (mito["upper"] == DEFAULT_MAX_PCT_MITO).all()

    def test_failure_counts_are_recorded(self):
        t = qc_thresholds(metrics_frame())
        assert (t["n_failed"] <= t["n_cells"]).all()
        assert t["n_cells"].sum() > 0

    def test_degenerate_batch_cuts_nobody_on_the_mad_rule(self):
        """Zero MAD means the rule cannot discriminate — it must not cut everything."""
        identical = np.full((50, len(GENES)), 100, dtype=np.int64)
        m = cell_qc_metrics(identical, GENES, batch=["flat"] * 50)
        t = qc_thresholds(m)
        assert np.isinf(t.loc[t["metric"] == "n_counts", "upper"].iloc[0])

    def test_missing_column_raises(self):
        with pytest.raises(QCError, match="missing column"):
            qc_thresholds(pd.DataFrame({"batch": ["s1"], "n_counts": [10]}))

    def test_defaults_match_w4(self):
        """Comparability at the gate: diverging silently would read as biology."""
        assert DEFAULT_N_MADS == 5.0
        assert DEFAULT_MAX_PCT_MITO == 20.0


class TestApplyQC:
    def test_high_mito_cells_fail(self):
        m = metrics_frame()
        m.loc[0, "pct_mito"] = 95.0
        passes = apply_qc(m, qc_thresholds(m))
        assert not passes.iloc[0]

    def test_typical_cells_pass(self):
        m = metrics_frame()
        passes = apply_qc(m, qc_thresholds(m))
        assert passes.mean() > 0.8

    def test_shallow_cell_fails_the_absolute_floor(self):
        m = metrics_frame()
        m.loc[1, ["n_counts", "n_genes"]] = [3, 2]
        passes = apply_qc(m, qc_thresholds(m))
        assert not passes.iloc[1]

    def test_a_batch_absent_from_the_table_is_left_alone(self):
        m = metrics_frame()
        thresholds = qc_thresholds(m)
        passes = apply_qc(m, thresholds[thresholds["batch"] == "s1"])
        assert passes[m["batch"] == "s2"].all()

    def test_index_is_preserved(self):
        m = metrics_frame().set_index(pd.Index(range(100, 220), name="cell"))
        passes = apply_qc(m, qc_thresholds(m))
        assert list(passes.index) == list(m.index)


class TestSummary:
    def test_counts_by_patient_and_tissue(self):
        m = metrics_frame()
        summary = qc_summary(m, apply_qc(m, qc_thresholds(m)))
        assert {"patient_id", "tissue", "n_cells", "n_passed", "retained"} <= set(summary.columns)
        assert summary["n_cells"].sum() == len(m)

    def test_retained_fraction_is_a_proportion(self):
        m = metrics_frame()
        summary = qc_summary(m, apply_qc(m, qc_thresholds(m)))
        assert summary["retained"].between(0, 1).all()

    def test_failed_plus_passed_is_total(self):
        m = metrics_frame()
        summary = qc_summary(m, apply_qc(m, qc_thresholds(m)))
        assert (summary["n_passed"] + summary["n_failed"] == summary["n_cells"]).all()

    def test_length_mismatch_raises(self):
        m = metrics_frame()
        with pytest.raises(QCError, match="entries for"):
            qc_summary(m, pd.Series([True, False]))


def test_doublet_flagging_is_an_explicit_todo():
    """Needs real matrices; the cutoff is a judgement call, not a formula."""
    with pytest.raises(NotImplementedError, match="per sample"):
        flag_doublets()
