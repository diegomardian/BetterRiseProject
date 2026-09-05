"""A1's machinery, against the inputs that would make its answer meaningless.

The experiment reports a difference between two references. The ways that
difference can be fake are: the two references not actually differing, the
ground truth not varying, the deconvolution seeing no genes, or the truth being
scored against the requested draw rather than the realised one. Each has a test.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.harness.pseudobulk_recovery import (
    EPITHELIAL_LABELS,
    RECIPES,
    RecoveryError,
    build_reference,
    make_pseudobulk,
    non_epithelial,
    recipe_gap,
    recover,
)

COMPARTMENTS = ["Epithelial cells", "T cells", "Stromal cells", "Myeloids"]


@pytest.fixture
def cells():
    """600 cells, four compartments, each with genes it dominates.

    Deliberately over-dispersed within compartment: the two recipes differ by
    Jensen's gap, which is zero when every cell of a type is identical. A
    fixture with no within-type spread would make the experiment measure
    nothing and every test pass.
    """
    rng = np.random.default_rng(11)
    n_genes, per = 300, 150
    rows, labels = [], []
    for j, compartment in enumerate(COMPARTMENTS):
        base = rng.lognormal(0.5, 0.4, n_genes)
        base[j * 60:(j + 1) * 60] *= 25.0
        for _ in range(per):
            rows.append(base * rng.lognormal(0.0, 0.9, n_genes))
            labels.append(compartment)
    expression = pd.DataFrame(
        rows, columns=[f"G{i:04d}" for i in range(n_genes)],
        index=[f"c{i:05d}" for i in range(len(rows))],
    )
    return expression, pd.Series(labels, index=expression.index)


# ---------------------------------------------------------------------------
# The two recipes must actually differ


def test_the_two_recipes_give_different_references(cells):
    """If they agree, the experiment has no signal and its gap is noise."""
    expression, labels = cells
    linear = build_reference(expression, labels, recipe="linear")
    log1p = build_reference(expression, labels, recipe="log1p")
    assert list(linear.columns) == list(log1p.columns)
    relative = ((linear - log1p).abs() / linear.replace(0, np.nan)).stack().median()
    assert relative > 0.05, (
        f"the recipes differ by a median {relative:.3%}. With no within-type "
        f"dispersion Jensen's gap vanishes and this experiment measures nothing."
    )


def test_the_geometric_mean_is_biased_low_as_jensen_says(cells):
    """Direction, not just magnitude -- it is the whole mechanism."""
    expression, labels = cells
    linear = build_reference(expression, labels, recipe="linear")
    log1p = build_reference(expression, labels, recipe="log1p")
    assert (log1p <= linear + 1e-9).to_numpy().mean() > 0.95


def test_an_unknown_recipe_is_refused(cells):
    expression, labels = cells
    with pytest.raises(RecoveryError, match="recipe must be"):
        build_reference(expression, labels, recipe="sqrt")


def test_a_cell_without_a_label_is_refused(cells):
    expression, labels = cells
    labels = labels.copy()
    labels.iloc[0] = np.nan
    with pytest.raises(RecoveryError, match="no compartment label"):
        build_reference(expression, labels, recipe="linear")


# ---------------------------------------------------------------------------
# The ground truth


def test_pseudobulk_truth_varies_and_sums_to_one(cells):
    expression, labels = cells
    bulk, truth = make_pseudobulk(expression, labels, n_samples=40, n_cells=200, seed=3)
    assert bulk.shape == (40, expression.shape[1])
    assert np.allclose(truth.sum(axis=1), 1.0)
    assert truth.std().min() > 0.01, (
        "a constant truth makes every correlation undefined, which is the "
        "degenerate case this project keeps finding"
    )


def test_truth_is_the_realised_draw_not_the_requested_one(cells):
    """The defect the WMHS paper documents at length, in a new generator.

    The multinomial does not land on the Dirichlet exactly. Scoring against the
    request rather than the draw makes the statistic a function of the
    generator, and the recovered fractions can then never disagree with it.
    """
    expression, labels = cells
    _, truth = make_pseudobulk(expression, labels, n_samples=200, n_cells=50, seed=5)
    # With only 50 cells the realised fractions land on multiples of 1/50.
    scaled = truth.to_numpy() * 50
    residual = np.minimum(scaled % 1, 1 - (scaled % 1))   # nearest integer, either side
    assert np.allclose(residual, 0, atol=1e-6), (
        "the truth is not a realised cell count, so it is the requested "
        "parameter rather than the draw"
    )


def test_a_compartment_too_thin_to_draw_is_refused(cells):
    expression, labels = cells
    labels = labels.copy()
    labels[labels == "Myeloids"] = "T cells"
    labels.iloc[:5] = "Myeloids"
    with pytest.raises(RecoveryError, match="under 20 cells"):
        make_pseudobulk(expression, labels, n_samples=10, n_cells=100, seed=1)


def test_one_compartment_is_refused(cells):
    expression, labels = cells
    with pytest.raises(RecoveryError, match="at least two compartments"):
        make_pseudobulk(expression, pd.Series("Epithelial cells", index=labels.index),
                        n_samples=5, n_cells=50, seed=1)


def test_non_epithelial_is_everything_that_is_not_epithelium(cells):
    frame = pd.DataFrame([[0.4, 0.3, 0.2, 0.1]], columns=COMPARTMENTS)
    assert non_epithelial(frame).iloc[0] == pytest.approx(0.6)
    assert EPITHELIAL_LABELS[0] in COMPARTMENTS
    with pytest.raises(RecoveryError, match="no non-epithelial"):
        non_epithelial(pd.DataFrame([[1.0]], columns=[EPITHELIAL_LABELS[0]]))


# ---------------------------------------------------------------------------
# Recovery, end to end


def test_a_correct_reference_recovers_the_fractions_it_was_given(cells):
    """The experiment must be able to succeed, or its failure means nothing."""
    expression, labels = cells
    bulk, truth = make_pseudobulk(expression, labels, n_samples=60, n_cells=400, seed=7)
    reference = build_reference(expression, labels, recipe="linear")
    results = recover(bulk, truth, reference, leg="self", recipe="linear")
    assert results
    for r in results:
        assert r.r_non_epithelial > 0.8, (
            f"{r.method} recovers the non-epithelial fraction at only "
            f"r={r.r_non_epithelial:.3f} from the reference that generated it"
        )


def test_the_gap_table_carries_both_recipes_and_their_difference(cells):
    expression, labels = cells
    bulk, truth = make_pseudobulk(expression, labels, n_samples=40, n_cells=300, seed=9)
    results = []
    for recipe in RECIPES:
        reference = build_reference(expression, labels, recipe=recipe)
        results += recover(bulk, truth, reference, leg="self", recipe=recipe)
    gaps = recipe_gap(results)
    assert {"leg", "method"} <= set(gaps.columns)
    assert "r_non_epithelial_gap" in gaps.columns
    for _, row in gaps.iterrows():
        assert row["r_non_epithelial_gap"] == pytest.approx(
            row["r_non_epithelial_linear"] - row["r_non_epithelial_log1p"]
        )


# ---------------------------------------------------------------------------
# The driver's one silent-failure route


def test_the_signature_resolves_to_symbols_or_raises():
    """Ensembl-keyed S matrix against a symbol-keyed cohort is an EMPTY join.

    Not an error -- an empty gene set, a deconvolution of nothing, and a
    reported number. The same identifier-space class the WMHS appendix records
    twice.
    """
    from src.common.paths import RESULTS_DIR
    from src.harness.jobs.run_pseudobulk_recovery import signature_symbols

    path = RESULTS_DIR / "2026-08-26_63ead2e" / "S_matrix_lineage_1.0.0.parquet"
    if not path.exists():
        pytest.skip("S matrix not committed")
    symbols = signature_symbols(path)
    assert len(symbols) > 700
    assert all(not s.startswith("ENSG") for s in symbols)


def test_a_signature_that_resolves_to_nothing_raises(tmp_path):
    from src.harness.jobs.run_pseudobulk_recovery import signature_symbols

    frame = pd.DataFrame({"gene": ["NOT_AN_ENSEMBL_ID"], "differentiated": [1.0]})
    path = tmp_path / "s.parquet"
    frame.to_parquet(path)
    with pytest.raises(RecoveryError, match="empty gene set"):
        signature_symbols(path)
