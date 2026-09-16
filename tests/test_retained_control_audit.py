"""The retained-control audit must be able to convict the labelling axis.

An audit that can only exonerate is not an audit. Every test here that asserts
REFUTED has a sibling feeding data where the same check FIRES.
"""

import numpy as np
import pandas as pd
import pytest

from src.reference.jobs.retained_control_audit import (
    AXIS_CONTROL_GENE,
    RETAINED_GENE,
    abundance_inversions,
    abundance_rank_correlation,
    adjudicate,
    relative_change,
    stratum_table,
    verdict,
)


def _rows(gene, rung, axis, *, base, tip, frac_n=1.0, frac_t=1.0, n_mature=500,
          n_patients=4):
    return [
        {
            "patient_id": f"P{i}",
            "study_id": "TEST",
            "gene": gene,
            "granularity_rung": rung,
            "labeling_axis": axis,
            "frac_mature_normal": frac_n,
            "frac_mature_tumour": frac_t,
            "mean_normal": base,
            "mean_tumour": tip,
            "n_cells_mature": n_mature,
        }
        for i in range(n_patients)
    ]


def _frames(retained_tip, *, control_tip=0.9, frac_t=1.0, n_mature=500,
            matched_tip=None):
    """A minimal two-source grid: one unmatched read, one depth-matched read."""
    def build(tip):
        rows = []
        for rung in ("epithelial", "lineage"):
            rows += _rows(RETAINED_GENE, rung, "stem_pole", base=2.5, tip=tip,
                          frac_t=frac_t, n_mature=n_mature)
            rows += _rows(AXIS_CONTROL_GENE, rung, "stem_pole", base=1.3,
                          tip=control_tip, frac_t=frac_t, n_mature=n_mature)
            rows += _rows("GUCA2A", rung, "stem_pole", base=25.0, tip=0.5,
                          frac_t=frac_t, n_mature=n_mature)
        return pd.DataFrame(rows)

    return {
        ("TEST", "unmatched"): build(retained_tip),
        ("TEST", "depth_matched"): build(
            retained_tip if matched_tip is None else matched_tip
        ),
    }


def _bulk(retained_log2=-8.176, control_log2=-0.237):
    return pd.DataFrame([
        {"gene": "GUCA2A", "log2_fold_change": -5.858, "fold_change": 0.0172},
        {"gene": AXIS_CONTROL_GENE, "log2_fold_change": control_log2,
         "fold_change": 2.0 ** control_log2},
        {"gene": RETAINED_GENE, "log2_fold_change": retained_log2,
         "fold_change": 2.0 ** retained_log2},
    ])


# --------------------------------------------------------------------------
# relative_change — the quantity both papers report
# --------------------------------------------------------------------------

def test_relative_change_is_floored_at_minus_one_only_by_a_zero_numerator():
    df = pd.DataFrame({"mean_normal": [2.0, 2.0], "mean_tumour": [0.0, 1.0]})
    assert relative_change(df).tolist() == [-1.0, -0.5]


def test_a_zero_baseline_is_not_silently_called_minus_one():
    """The confusion this whole audit is about. A gene with no normal-arm
    signal has no baseline to fall from; that is NaN, not a 100% loss."""
    df = pd.DataFrame({"mean_normal": [0.0], "mean_tumour": [0.0]})
    assert np.isnan(relative_change(df)[0])


# --------------------------------------------------------------------------
# stratum_table — label_selects_nothing must be measured, not assumed
# --------------------------------------------------------------------------

def test_label_selects_nothing_is_true_only_when_both_fractions_are_one():
    grid = stratum_table(_frames(0.08))
    assert grid["label_selects_nothing"].all()


def test_label_selects_nothing_is_false_when_the_axis_partitions():
    grid = stratum_table(_frames(0.08, frac_t=0.3))
    assert not grid["label_selects_nothing"].any()


def test_a_tiny_tumour_compartment_is_flagged_degenerate():
    grid = stratum_table(_frames(0.08, n_mature=2))
    assert grid["degenerate_stratum"].all()
    assert not stratum_table(_frames(0.08))["degenerate_stratum"].any()


# --------------------------------------------------------------------------
# hypothesis (d) — abundance ordering
#
# Sign convention, because the first draft of this job had it backwards:
# dropout eats the LEAST abundant genes hardest, so (d) predicts a NEGATIVE
# Spearman rho between log10 baseline and drop magnitude.
# --------------------------------------------------------------------------

def test_rank_correlation_is_negative_when_dropout_really_does_drive_it():
    """The least abundant gene falls hardest: rho negative, (d)'s own shape."""
    rows = []
    for gene, base, tip in (("LOW", 0.5, 0.005), ("MID", 5.0, 2.5),
                            ("HIGH", 50.0, 45.0)):
        rows += _rows(gene, "epithelial", "stem_pole", base=base, tip=tip)
    grid = stratum_table({("TEST", "unmatched"): pd.DataFrame(rows)})
    assert abundance_rank_correlation(grid) == pytest.approx(-1.0)


def test_rank_correlation_is_not_what_the_verdict_rests_on():
    """On the committed shape rho runs the same way dropout would predict, and
    the verdict must NOT follow it — the inversion test is what decides."""
    grid = stratum_table(_frames(0.08))
    ledger = adjudicate(grid, _bulk())
    row = ledger[ledger["hypothesis"] == "d_detection_floor"].iloc[0]
    assert row["verdict"] == "REFUTED"
    assert len(abundance_inversions(grid)) >= 1


def test_degenerate_strata_cannot_vote_in_the_correlation():
    ok = abundance_rank_correlation(stratum_table(_frames(0.08)))
    with_junk = abundance_rank_correlation(
        stratum_table(_frames(0.08, n_mature=1))
    )
    assert np.isnan(with_junk)
    assert not np.isnan(ok)


def test_an_inversion_is_a_less_abundant_gene_that_drops_less():
    """CDX2's real shape: half MS4A12's baseline, a third of its drop."""
    inv = abundance_inversions(stratum_table(_frames(0.08)))
    assert AXIS_CONTROL_GENE in inv.index
    assert inv.loc[AXIS_CONTROL_GENE, "baseline"] < 2.5


def test_no_inversion_when_abundance_does_order_the_drop():
    """Make every less-abundant gene fall harder and (d) must survive."""
    rows = []
    for rung in ("epithelial", "lineage"):
        rows += _rows(RETAINED_GENE, rung, "stem_pole", base=2.5, tip=0.08)
        rows += _rows(AXIS_CONTROL_GENE, rung, "stem_pole", base=1.3, tip=0.001)
        rows += _rows("GUCA2A", rung, "stem_pole", base=25.0, tip=0.5)
    grid = stratum_table({("TEST", "unmatched"): pd.DataFrame(rows)})
    assert len(abundance_inversions(grid)) == 0
    ledger = adjudicate(grid, _bulk())
    assert ledger[ledger["hypothesis"] == "d_detection_floor"].iloc[0]["verdict"] == "NOT REFUTED"


# --------------------------------------------------------------------------
# adjudicate / verdict — the audit must be able to convict
# --------------------------------------------------------------------------

def test_the_committed_shape_exonerates_the_axis():
    ledger = adjudicate(stratum_table(_frames(0.08)), _bulk())
    axis = ledger[ledger["hypothesis"].str.startswith("a_")]
    assert (axis["verdict"] == "REFUTED").all()
    assert verdict(ledger)["verdict"] == "THE LABELLING AXIS IS NOT DRIVING THE CONTROL"


def test_the_axis_is_convicted_when_the_control_gene_collapses_too():
    """If CDX2 also read near zero on the same label-free cells, hypothesis (a)
    would be the live one and the verdict must say so."""
    ledger = adjudicate(
        stratum_table(_frames(0.08, control_tip=0.005)), _bulk()
    )
    assert verdict(ledger)["verdict"] != "THE LABELLING AXIS IS NOT DRIVING THE CONTROL"


def test_the_depth_leg_is_refuted_only_when_matching_leaves_the_drop_intact():
    ledger = adjudicate(stratum_table(_frames(0.08)), _bulk())
    row = ledger[ledger["hypothesis"] == "a_axis_driven_via_depth"].iloc[0]
    assert row["verdict"] == "REFUTED"


def test_the_axis_is_convicted_when_depth_matching_removes_the_drop():
    """Depth-match the arms and the drop mostly goes away: depth was doing the
    work, and the verdict must flip. This is the test that caught the first
    draft of this job asserting REFUTED in the source."""
    ledger = adjudicate(
        stratum_table(_frames(0.08, matched_tip=2.4)), _bulk()
    )
    row = ledger[ledger["hypothesis"] == "a_axis_driven_via_depth"].iloc[0]
    assert row["verdict"] == "NOT REFUTED"
    assert verdict(ledger)["verdict"] == "THE LABELLING AXIS MAY BE DRIVING THE CONTROL"


def test_the_denominator_leg_is_convicted_when_the_baseline_is_at_the_floor():
    rows = []
    for rung in ("epithelial", "lineage"):
        rows += _rows(RETAINED_GENE, rung, "stem_pole", base=0.02, tip=0.0)
        rows += _rows(AXIS_CONTROL_GENE, rung, "stem_pole", base=1.3, tip=0.9)
        rows += _rows("GUCA2A", rung, "stem_pole", base=25.0, tip=0.5)
    grid = stratum_table({("TEST", "unmatched"): pd.DataFrame(rows)})
    ledger = adjudicate(grid, _bulk())
    row = ledger[ledger["hypothesis"] == "c_denominator_collapse"].iloc[0]
    assert row["verdict"] == "NOT REFUTED"


def test_bulk_replication_is_what_supports_the_not_retained_reading():
    """Strip the label-free bulk fall and (b) must stop being SUPPORTED."""
    ledger = adjudicate(stratum_table(_frames(0.08)), _bulk(retained_log2=-0.05))
    row = ledger[ledger["hypothesis"] == "b_gene_not_retained"].iloc[0]
    assert row["verdict"] == "NOT SUPPORTED"
    assert verdict(ledger)["verdict"] == "UNRESOLVED"


def test_the_endpoint_floor_is_reported_separately_and_not_as_a_refutation():
    ledger = adjudicate(stratum_table(_frames(0.0, n_mature=1)), _bulk())
    row = ledger[ledger["hypothesis"] == "c_degenerate_stratum_floor"].iloc[0]
    assert row["verdict"] == "UPHELD_FOR_THE_ENDPOINT_ONLY"


def test_every_ledger_row_states_its_prediction_before_its_observation():
    ledger = adjudicate(stratum_table(_frames(0.08)), _bulk())
    for col in ("statement", "prediction", "observed", "verdict"):
        assert ledger[col].map(bool).all(), f"{col} is empty somewhere"
