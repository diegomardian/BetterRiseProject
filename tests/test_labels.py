"""Per-cell maturity labelling. CLAUDE.md invariant 2 lives here too: axis
markers must never overlap the target genes under test for a run.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.common.panel import granularity_rungs
from src.estimator.labels import (
    PROVISIONAL_RUNG_QUANTILES,
    axis_score,
    classify_maturity,
    label_cohort,
)
from src.reference.signature import LeakageError


def _synthetic_expression(n_cells=400, seed=0):
    """Cells cleanly separated along the stem-pole axis genes, so a maturity
    call has an unambiguous right answer to check against."""
    rng = np.random.default_rng(seed)
    stem_genes = ["LGR5", "ASCL2", "MKI67", "OLFM4", "SMOC2"]
    goblet_genes = ["MUC2", "TFF3", "SPDEF", "ITLN1"]
    other_genes = [f"GENE{i:04d}" for i in range(50)]

    # First half: near the stem pole (high stem score). Second half: mature
    # (low stem score, low goblet score too).
    is_stem_like = np.arange(n_cells) < n_cells // 2
    stem_expr = np.where(is_stem_like[:, None], rng.normal(8, 0.5, (n_cells, len(stem_genes))),
                          rng.normal(1, 0.5, (n_cells, len(stem_genes))))
    goblet_expr = rng.normal(1, 0.5, (n_cells, len(goblet_genes)))
    other_expr = rng.normal(3, 1, (n_cells, len(other_genes)))

    columns = stem_genes + goblet_genes + other_genes
    data = np.hstack([stem_expr, goblet_expr, other_expr])
    return pd.DataFrame(data, columns=columns), is_stem_like


def test_axis_score_averages_the_axis_marker_genes():
    expression, is_stem_like = _synthetic_expression()
    score = axis_score(expression, "stem_pole", target_genes=["GUCA2A"])
    assert score[is_stem_like].mean() > score[~is_stem_like].mean()


def test_axis_score_guards_against_target_leakage():
    expression, _ = _synthetic_expression()
    with pytest.raises(LeakageError):
        axis_score(expression, "stem_pole", target_genes=["LGR5"])


def test_axis_score_raises_for_an_axis_with_no_transcript_markers():
    """Chromatin and spatial (axis 3) carry no marker genes yet -- week 13+.
    This must raise, not silently score everything as zero."""
    expression, _ = _synthetic_expression()
    with pytest.raises(ValueError, match="chromatin and spatial"):
        axis_score(expression, "chromatin", target_genes=["GUCA2A"])


def test_classify_maturity_calls_stem_cells_immature():
    """Cells near the stem pole should not be classified as mature at any
    rung -- the axis is inverted before thresholding."""
    expression, is_stem_like = _synthetic_expression()
    mature = classify_maturity(expression, "stem_pole", "lineage", target_genes=["GUCA2A"])
    # Every stem-like cell should score below the mature threshold.
    assert not mature[is_stem_like].any()
    assert mature[~is_stem_like].mean() > 0.5


def test_classify_maturity_rejects_an_unknown_rung():
    expression, _ = _synthetic_expression()
    with pytest.raises(KeyError, match="unknown rung"):
        classify_maturity(expression, "stem_pole", "single_cell", target_genes=["GUCA2A"])


def test_finer_rungs_call_a_narrower_mature_population():
    """best4 (0.95) should flag fewer cells mature than epithelial (0.50) --
    README design decision 3: granularity narrows the compartment."""
    expression, _ = _synthetic_expression(n_cells=1000)
    epithelial = classify_maturity(expression, "stem_pole", "epithelial", target_genes=["GUCA2A"])
    best4 = classify_maturity(expression, "stem_pole", "best4", target_genes=["GUCA2A"])
    assert best4.sum() < epithelial.sum()


def test_rung_quantiles_are_monotonically_narrowing():
    ordered = [PROVISIONAL_RUNG_QUANTILES[r] for r in granularity_rungs()]
    assert ordered == sorted(ordered)


def test_label_cohort_produces_one_column_per_axis_and_rung_never_overwriting():
    expression, _ = _synthetic_expression(n_cells=100)
    out = label_cohort(
        expression,
        patient_id=[f"P{i % 5}" for i in range(100)],
        tissue=["normal"] * 50 + ["tumour"] * 50,
        target_genes=["GUCA2A"],
        axes=("stem_pole",),
        rungs=("epithelial", "best4"),
    )
    assert "mature__stem_pole__epithelial" in out.columns
    assert "mature__stem_pole__best4" in out.columns
    assert not out["mature__stem_pole__epithelial"].equals(out["mature__stem_pole__best4"])


def test_label_cohort_rejects_mismatched_metadata_length():
    expression, _ = _synthetic_expression(n_cells=100)
    with pytest.raises(ValueError, match="one entry per cell"):
        label_cohort(
            expression,
            patient_id=["P1"],
            tissue=["normal"] * 100,
            target_genes=["GUCA2A"],
        )


def test_a_run_testing_muc2_must_not_use_opposite_lineage_axis():
    """docs/open_decisions.md #1: MUC2/TFF3 are both tier-E targets and
    axis-2 (opposite_lineage) labels."""
    expression, _ = _synthetic_expression()
    with pytest.raises(LeakageError):
        classify_maturity(expression, "opposite_lineage", "lineage", target_genes=["MUC2"])
