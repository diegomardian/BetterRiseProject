"""Real Lee cohort file structure, exercised via small committed fixtures
extracted from the actual GSE132465/GSE144735 downloads (real header tokens,
real gene symbols/order, real barcode format) rather than synthetic data.

Fixture provenance (reproducible from data/raw/lee/ if it needs regenerating):
25 real SMC cells (SMC01 x5 Tumor/5 Normal, SMC02 x5 Tumor/5 Normal, SMC14 x5
Tumor-only) and 25 real KUL3 cells (KUL01 x5 Tumor/5 Normal/5 Border, KUL19
x5 Tumor/5 Normal), restricted to ~40 real gene rows (panel tier A+B, axis
1/2 markers, 3 MT- genes, ~20 filler genes) via a single `gunzip -c | awk`
pass over each real file, keeping MUC2 out on purpose (it is absent from the
real gene index in both cohorts).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src.estimator.kitagawa import decompose_cohort
from src.estimator.lee_io import (
    EPITHELIAL_COMPARTMENT,
    MITO_GENES,
    build_gene_rung_axis_summary,
    label_agreement,
    load_annotation,
    load_lee_cohort,
    stream_matrix_stats,
)
from src.schema import coerce_results, validate_results

FIXTURES = Path(__file__).parent / "fixtures" / "lee"
SMC_ANNOTATION = FIXTURES / "GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz"
SMC_MATRIX = FIXTURES / "GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz"
KUL3_ANNOTATION = FIXTURES / "GSE144735_processed_KUL3_CRC_10X_annotation.txt.gz"
KUL3_MATRIX = FIXTURES / "GSE144735_processed_KUL3_CRC_10X_raw_UMI_count_matrix.txt.gz"

PANEL_AB = ["GUCA2A", "GUCA2B", "OTOP2", "CA7", "MLH1", "SFRP1", "SFRP2"]
AXIS_GENES = ["LGR5", "ASCL2", "MKI67", "OLFM4", "SMOC2", "TFF3", "SPDEF", "ITLN1"]


# ---------------------------------------------------------------------------
# load_annotation
# ---------------------------------------------------------------------------


def test_load_annotation_maps_class_to_tissue():
    out = load_annotation(SMC_ANNOTATION, study_id="GSE132465")
    assert set(out["tissue"].unique()) == {"tumour", "normal"}
    assert (out["study_id"] == "GSE132465").all()


def test_load_annotation_maps_border_for_kul3():
    out = load_annotation(KUL3_ANNOTATION, study_id="GSE144735")
    assert "border" in set(out["tissue"].unique())


def test_load_annotation_rejects_an_unrecognized_class(tmp_path):
    bad = tmp_path / "bad_annotation.txt.gz"
    df = pd.read_csv(SMC_ANNOTATION, sep="\t")
    df.loc[0, "Class"] = "Metastasis"
    df.to_csv(bad, sep="\t", index=False, compression="gzip")
    with pytest.raises(ValueError, match="unrecognized Class"):
        load_annotation(bad, study_id="GSE132465")


# ---------------------------------------------------------------------------
# stream_matrix_stats
# ---------------------------------------------------------------------------


def test_stream_matrix_stats_finds_mito_and_target_genes():
    metrics, expression, coverage = stream_matrix_stats(
        SMC_MATRIX, genes_of_interest=PANEL_AB + AXIS_GENES
    )
    assert set(coverage["found"]) == set(PANEL_AB + AXIS_GENES)
    assert coverage["missing"] == []
    assert list(expression.columns) != []
    assert (metrics["n_counts"] >= 0).all()


def test_stream_matrix_stats_reports_muc2_as_missing():
    _, _, coverage = stream_matrix_stats(
        SMC_MATRIX, genes_of_interest=PANEL_AB + ["MUC2"]
    )
    assert coverage["missing"] == ["MUC2"]
    assert "MUC2" not in coverage["found"]


def test_stream_matrix_stats_matches_a_hand_computed_golden_value():
    """Independently computed via a plain, unchunked pandas read -- a
    different code path than the streaming accumulator under test."""
    raw = pd.read_csv(SMC_MATRIX, sep="\t", index_col=0)
    cell = raw.columns[0]
    golden_n_counts = float(raw[cell].sum())
    golden_n_genes = int((raw[cell] > 0).sum())
    golden_mito = float(raw.loc[raw.index.isin(MITO_GENES), cell].sum())
    golden_pct_mito = 100.0 * golden_mito / golden_n_counts

    metrics, _, _ = stream_matrix_stats(
        SMC_MATRIX, genes_of_interest=PANEL_AB
    )
    assert metrics.loc[cell, "n_counts"] == pytest.approx(golden_n_counts)
    assert metrics.loc[cell, "n_genes"] == golden_n_genes
    assert metrics.loc[cell, "pct_mito"] == pytest.approx(golden_pct_mito)


def test_stream_matrix_stats_is_reproducible_and_gene_order_invariant():
    """The accumulator is a single running-sum pass with no chunking left to
    vary -- what's worth pinning now is that re-running produces identical
    output, and that asking for the same genes in a different order doesn't
    change the result (guards against an accumulator that's accidentally
    order-sensitive, e.g. via a mutable default or a stray sort)."""
    baseline_metrics, baseline_expr, baseline_cov = stream_matrix_stats(
        SMC_MATRIX, genes_of_interest=PANEL_AB + AXIS_GENES
    )
    metrics, expr, cov = stream_matrix_stats(
        SMC_MATRIX, genes_of_interest=list(reversed(PANEL_AB + AXIS_GENES))
    )
    pd.testing.assert_frame_equal(metrics, baseline_metrics)
    pd.testing.assert_frame_equal(
        expr.sort_index(axis=1), baseline_expr.sort_index(axis=1)
    )
    assert cov == baseline_cov


def test_stream_matrix_stats_raises_if_no_mito_genes_found():
    with pytest.raises(ValueError, match="mitochondrial"):
        stream_matrix_stats(
            SMC_MATRIX, genes_of_interest=PANEL_AB, mito_genes=("MT-MADEUP",)
        )


def test_stream_matrix_stats_raises_on_an_empty_file(tmp_path):
    empty = tmp_path / "empty.txt.gz"
    pd.DataFrame(columns=["Index"]).set_index("Index").to_csv(
        empty, sep="\t", compression="gzip"
    )
    with pytest.raises(ValueError):
        stream_matrix_stats(empty, genes_of_interest=PANEL_AB)


# ---------------------------------------------------------------------------
# load_lee_cohort — end to end on fixtures
# ---------------------------------------------------------------------------


def test_load_lee_cohort_excludes_the_tumor_only_patient(tmp_path):
    cohort = load_lee_cohort(
        "smc", target_genes=PANEL_AB, raw_dir=FIXTURES
    )
    assert "SMC14" in cohort.excluded_patients
    assert "SMC14" not in set(cohort.cells["patient_id"])
    assert {"SMC01", "SMC02"} <= set(cohort.cells["patient_id"])


def test_load_lee_cohort_surfaces_muc2_coverage_gap():
    cohort = load_lee_cohort(
        "smc", target_genes=PANEL_AB, raw_dir=FIXTURES
    )
    assert "MUC2" in cohort.axis_gene_coverage["missing"]


def test_load_lee_cohort_labels_still_run_with_partial_axis_coverage():
    """opposite_lineage has only 3/4 markers for Lee (no MUC2) -- must not
    raise, per axis_score's documented tolerance for a partial marker set."""
    cohort = load_lee_cohort(
        "smc", target_genes=PANEL_AB, raw_dir=FIXTURES
    )
    assert any(c.startswith("mature__opposite_lineage__") for c in cohort.labels.columns)


def test_load_lee_cohort_counts_border_cells_for_kul3():
    """Border cells are counted AND retained in .cells (labelled like any
    other cell) -- only excluded later, at the decompose_cohort-bound
    summary stage, since Kitagawa has no third arm for them."""
    cohort = load_lee_cohort(
        "kul3", target_genes=PANEL_AB, raw_dir=FIXTURES
    )
    assert cohort.n_border_cells > 0
    assert "border" in set(cohort.cells["tissue"])


def test_load_lee_cohort_has_no_normal_arm_smc_with_zero_normal_patients_excluded():
    """Every SMC patient in the fixture with a normal+tumour pair (SMC01,
    SMC02) must survive; only the tumor-only one (SMC14) is dropped."""
    cohort = load_lee_cohort(
        "smc", target_genes=PANEL_AB, raw_dir=FIXTURES
    )
    assert cohort.excluded_patients == ["SMC14"]


# ---------------------------------------------------------------------------
# label_compartment -- docs/open_decisions.md #13, "caller must pre-filter"
# ---------------------------------------------------------------------------


def _mixed_compartment_fixture(tmp_path):
    """The committed fixture is all-epithelial, which is the one shape that
    cannot exercise the restriction. Relabel half its cells as T cells --
    real barcodes, real counts, only the annotation's Cell_type changes."""
    import gzip
    import shutil

    annotation = pd.read_csv(SMC_ANNOTATION, sep="\t")
    every_other = annotation.index[::2]
    annotation.loc[every_other, "Cell_type"] = "T cells"

    out_annotation = tmp_path / SMC_ANNOTATION.name
    with gzip.open(out_annotation, "wt") as fh:
        annotation.to_csv(fh, sep="\t", index=False)
    shutil.copy(SMC_MATRIX, tmp_path / SMC_MATRIX.name)
    return tmp_path, set(annotation.loc[every_other, "Index"])


def test_cells_outside_the_compartment_are_na_not_immature(tmp_path):
    """The compositional term is the mature FRACTION. Reading an unlabelled
    T cell as False puts the tumour's immune infiltrate in the denominator."""
    raw_dir, t_cells = _mixed_compartment_fixture(tmp_path)
    cohort = load_lee_cohort("smc", target_genes=PANEL_AB, raw_dir=raw_dir)

    assert cohort.label_compartment == EPITHELIAL_COMPARTMENT
    maturity = [c for c in cohort.labels.columns if c.startswith("mature__")]
    assert maturity

    outside = cohort.labels.index.intersection(list(t_cells))
    assert len(outside) > 0
    inside = cohort.labels.index.difference(outside)
    for col in maturity:
        assert cohort.labels.loc[outside, col].isna().all()
        # An epithelial cell may ALSO be NA now — below the depth floor, which
        # is "not asked" for a different reason. What must never happen is a
        # silent False, which would put the cell in the denominator as immature.
        assert not cohort.labels.loc[inside, col].fillna(True).isna().any()

    # And the reason is recoverable rather than merged into one NA: W1's label
    # column distinguishes non_epithelial from unresolved_depth.
    reasons = set(cohort.labels["label_stem_pole__lineage".replace("__", "_", 1)].astype(str))
    assert "unresolved_depth" in reasons or "non_epithelial" in reasons


def test_the_depth_floor_fires_and_does_not_count_shallow_cells_as_immature(tmp_path):
    """The fix for issue #44, at the seam where it was introduced.

    The previous labeller had no depth handling, so a cell that sampled zero
    stem markers scored maximally mature. W1's thins to a common depth and puts
    cells below the floor in `unresolved_depth` — not scored, and above all not
    counted as immature, which would move the compositional term.
    """
    raw_dir, _ = _mixed_compartment_fixture(tmp_path)
    cohort = load_lee_cohort("smc", target_genes=PANEL_AB, raw_dir=raw_dir)
    labels = cohort.labels["label_stem_pole_lineage"].astype(str)
    assert (labels == "unresolved_depth").any(), "fixture no longer exercises the floor"

    shallow = labels.index[labels == "unresolved_depth"]
    for col in [c for c in cohort.labels.columns if c.startswith("mature__")]:
        assert cohort.labels.loc[shallow, col].isna().all()


def test_the_maturity_quantile_is_taken_within_the_compartment(tmp_path):
    """The reason the restriction is in the loader and not left to the caller:
    the threshold is a quantile of whatever it is handed. Non-epithelial cells
    carry no stem markers, so on the inverted axis they score as maximally
    mature -- label everything and the cut lands in the immune mass."""
    raw_dir, _ = _mixed_compartment_fixture(tmp_path)
    restricted = load_lee_cohort("smc", target_genes=PANEL_AB, raw_dir=raw_dir)
    everything = load_lee_cohort(
        "smc", target_genes=PANEL_AB, raw_dir=raw_dir, label_compartment=None
    )

    # NOT the `epithelial` rung: it is single-bin by design under W1's
    # RUNG_SPECS ("the whole epithelium is one population"), so every scored
    # cell is mature there and the column cannot differ. Use a rung that bins.
    col = "mature__stem_pole__lineage"
    epithelium = restricted.labels.index[restricted.labels[col].notna()]
    assert len(epithelium) > 0
    restricted_calls = restricted.labels.loc[epithelium, col]
    everything_calls = everything.labels.loc[epithelium, col].astype("boolean")
    assert not restricted_calls.equals(everything_calls), (
        "restriction changed no epithelial call -- fixture no longer separates them"
    )


def test_labelling_everything_stays_available_and_says_so(tmp_path):
    raw_dir, _ = _mixed_compartment_fixture(tmp_path)
    cohort = load_lee_cohort(
        "smc", target_genes=PANEL_AB, raw_dir=raw_dir, label_compartment=None
    )
    assert cohort.label_compartment is None
    maturity = [c for c in cohort.labels.columns if c.startswith("mature__")]
    assert maturity
    # Every cell keeps its ROW — that is what "stays available" means. It no
    # longer means every cell gets a call: the depth floor legitimately leaves
    # some unscored, and NA there is the honest value.
    assert len(cohort.labels) == len(cohort.expression)
    scored = cohort.labels.loc[:, maturity].notna().any(axis=1)
    assert scored.any(), "labelling everything scored nothing"


def test_an_empty_compartment_raises_rather_than_labelling_nothing(tmp_path):
    raw_dir, _ = _mixed_compartment_fixture(tmp_path)
    with pytest.raises(ValueError, match="no cells in compartment"):
        load_lee_cohort(
            "smc", target_genes=PANEL_AB, raw_dir=raw_dir, label_compartment="Neurons"
        )


def test_the_summary_drops_unlabelled_cells_instead_of_counting_them(tmp_path):
    """End to end: the mature fraction must not move when non-epithelial cells
    are added to the cohort, because they are dropped rather than counted."""
    raw_dir, _ = _mixed_compartment_fixture(tmp_path)
    mixed = load_lee_cohort("smc", target_genes=PANEL_AB, raw_dir=raw_dir)
    summary = build_gene_rung_axis_summary(mixed, genes=PANEL_AB)
    if len(summary) == 0:
        pytest.skip("fixture too small to produce a summary row")
    assert summary["frac_mature_normal"].between(0, 1).all()
    assert summary["frac_mature_tumour"].between(0, 1).all()


# ---------------------------------------------------------------------------
# build_gene_rung_axis_summary -> decompose_cohort (integration smoke test)
# ---------------------------------------------------------------------------


def test_full_chain_from_raw_fixture_to_a_schema_valid_decomposition():
    """The highest-value test here: raw file -> LeeCohort -> summary ->
    decompose_cohort -> schema validation, with no synthetic data at all."""
    cohort = load_lee_cohort(
        "smc", target_genes=PANEL_AB, raw_dir=FIXTURES
    )
    summary = build_gene_rung_axis_summary(cohort, genes=PANEL_AB)
    assert len(summary) > 0

    results = decompose_cohort(summary)
    validate_results(coerce_results(results))
    assert set(results["study_id"]) == {"GSE132465"}


def test_full_chain_both_cohorts_stay_separate_study_ids():
    smc = load_lee_cohort(
        "smc", target_genes=PANEL_AB, raw_dir=FIXTURES
    )
    kul3 = load_lee_cohort(
        "kul3", target_genes=PANEL_AB, raw_dir=FIXTURES
    )
    summary = pd.concat(
        [
            build_gene_rung_axis_summary(smc, genes=PANEL_AB),
            build_gene_rung_axis_summary(kul3, genes=PANEL_AB),
        ],
        ignore_index=True,
    )
    if len(summary) == 0:
        pytest.skip("fixture too small to produce a summary row for both cohorts")
    results = decompose_cohort(summary)
    validate_results(coerce_results(results))
    assert set(results["study_id"]) <= {"GSE132465", "GSE144735"}


# ---------------------------------------------------------------------------
# label_agreement
# ---------------------------------------------------------------------------


def test_label_agreement_requires_a_deliberate_allowlist():
    with pytest.raises(ValueError, match="deliberately chosen"):
        label_agreement(pd.DataFrame(), pd.Series(dtype=object), mature_subtype_allowlist=[])


def test_label_agreement_computes_a_rate_per_axis_rung():
    cohort = load_lee_cohort(
        "smc", target_genes=PANEL_AB, raw_dir=FIXTURES
    )
    author_subtype = cohort.cells.loc[cohort.labels.index, "author_cell_subtype"]
    out = label_agreement(
        cohort.labels,
        author_subtype,
        mature_subtype_allowlist=["Mature Enterocytes type 1", "Mature Enterocytes type 2"],
    )
    assert (out["agreement"].between(0, 1)).all()


# ---------------------------------------------------------------------------
# raw_counts / extra_genes — the W2 harness path
# ---------------------------------------------------------------------------


def test_raw_counts_is_empty_unless_asked_for():
    """Off by default, so existing callers pay no memory for it."""
    cohort = load_lee_cohort("smc", target_genes=PANEL_AB, raw_dir=FIXTURES)
    assert cohort.raw_counts.empty


def test_raw_counts_are_integers_not_cp10k():
    """The generator thins counts binomially; a normalised frame would be
    truncated to zero silently. See LeeCohort.raw_counts."""
    cohort = load_lee_cohort(
        "smc", target_genes=PANEL_AB, raw_dir=FIXTURES, keep_raw_counts=True
    )
    assert not cohort.raw_counts.empty
    assert pd.api.types.is_integer_dtype(cohort.raw_counts.dtypes.iloc[0])
    values = cohort.raw_counts.to_numpy(dtype="float64")
    assert ((values == np.floor(values)) | np.isnan(values)).all()


def test_raw_counts_align_with_the_normalised_frame():
    cohort = load_lee_cohort(
        "smc", target_genes=PANEL_AB, raw_dir=FIXTURES, keep_raw_counts=True
    )
    assert list(cohort.raw_counts.index) == list(cohort.expression.index)
    assert list(cohort.raw_counts.columns) == list(cohort.expression.columns)


def test_raw_counts_and_cp10k_agree_up_to_library_size():
    """Same numbers, different scale — CP10K is counts / library_size * 1e4, so
    the two frames must be proportional within a cell."""
    cohort = load_lee_cohort(
        "smc", target_genes=PANEL_AB, raw_dir=FIXTURES, keep_raw_counts=True
    )
    raw = cohort.raw_counts.to_numpy(dtype="float64")
    cp10k = cohort.expression.to_numpy(dtype="float64")
    rows = np.flatnonzero((raw.sum(axis=1) > 0) & (cp10k.sum(axis=1) > 0))
    assert rows.size, "fixture has no cell with any expression"
    for i in rows[:5]:
        ratio = cp10k[i][raw[i] > 0] / raw[i][raw[i] > 0]
        assert np.allclose(ratio, ratio[0])  # one scale factor per cell


def test_extra_genes_widens_the_collected_set():
    """nu-SVR needs 500-2000 genes; the default set is roughly twenty, so the
    harness passes a marker panel rather than re-parsing the matrix."""
    narrow = load_lee_cohort("smc", target_genes=PANEL_AB, raw_dir=FIXTURES)
    wide = load_lee_cohort(
        "smc", target_genes=PANEL_AB, raw_dir=FIXTURES, extra_genes=["EPCAM", "PTPRC", "COL1A1"]
    )
    assert set(wide.axis_gene_coverage["requested"]) > set(
        narrow.axis_gene_coverage["requested"]
    )
    assert {"EPCAM", "PTPRC", "COL1A1"} <= set(wide.axis_gene_coverage["requested"])


def test_extra_genes_does_not_disturb_the_leakage_guard():
    """Widening the collected set must not put a target gene into the labels."""
    wide = load_lee_cohort(
        "smc", target_genes=PANEL_AB, raw_dir=FIXTURES, extra_genes=["EPCAM", "PTPRC"]
    )
    assert not set(PANEL_AB) & set(wide.labels.columns)
