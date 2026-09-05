"""Fractions from bulk, against the inputs that must defeat them.

The headline case is the one the pre-committed instrument gate cannot see: a
log-scale reference against linear bulk returns a mature fraction of exactly
0.0 on every sample, while the non-epithelial aggregate the gate reads comes
back at r ~ 0.88 and passes. Every guard here is tested against the input that
would have made it pass wrongly, not only against clean data.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src.bulk.deconvolution import (
    LINEAR_CP10K,
    LOG1P_CP10K,
    MIN_FRACTION_SD,
    NON_EPITHELIAL,
    DeconvolutionError,
    Reference,
    align,
    assert_scales_agree,
    check_predictor,
    deconvolve_cohort,
    default_methods,
    linearise,
    load_reference,
    mature_column,
    require_usable_predictor,
    summarise_fractions,
)
from src.common.paths import RESULTS_DIR

COMMITTED = RESULTS_DIR / "2026-08-26_63ead2e"
TYPES = ["differentiated", "endothelial", "epithelial_unscored", "immune", "stem_like", "stromal"]


# ---------------------------------------------------------------------------
# Fixtures: a small synthetic world whose truth is known in closed form


@pytest.fixture
def truth() -> pd.DataFrame:
    """Known fractions for 60 samples, mature-rich, none degenerate."""
    rng = np.random.default_rng(20260905)
    alpha = [3.0 if t == "differentiated" else 1.5 for t in TYPES]
    return pd.DataFrame(
        rng.dirichlet(alpha, size=60),
        columns=TYPES,
        index=[f"S{i:03d}" for i in range(60)],
    )


@pytest.fixture
def linear_reference() -> Reference:
    """A linear-scale signature: each type has genes it dominates."""
    rng = np.random.default_rng(7)
    n_genes = 400
    base = rng.lognormal(1.0, 0.6, size=(n_genes, len(TYPES)))
    for j in range(len(TYPES)):  # a private marker block per type
        base[j * 50:(j + 1) * 50, j] *= 30.0
    matrix = pd.DataFrame(
        base, columns=TYPES, index=[f"ENSG{i:08d}" for i in range(n_genes)]
    )
    return Reference(matrix=matrix, rung="lineage", scale=LINEAR_CP10K,
                     source=Path("synthetic.parquet"))


def _bulk_from(truth: pd.DataFrame, reference: Reference, *, noise: float = 0.0) -> pd.DataFrame:
    """Linear mixture, which is what TPM/CPM actually is."""
    values = truth[reference.matrix.columns].to_numpy() @ reference.matrix.to_numpy().T
    if noise:
        rng = np.random.default_rng(11)
        values = values * rng.lognormal(0.0, noise, size=values.shape)
    return pd.DataFrame(values, index=truth.index, columns=reference.matrix.index)


# ---------------------------------------------------------------------------
# It recovers a known truth when the scales agree


def test_known_fractions_are_recovered_when_the_scales_agree(truth, linear_reference):
    bulk = _bulk_from(truth, linear_reference)
    long, skipped = deconvolve_cohort(bulk, linear_reference, bulk_scale=LINEAR_CP10K)
    assert not skipped, skipped
    summary = summarise_fractions(long, "lineage")
    for method in summary["method"].unique():
        got = summary[summary["method"] == method].set_index("sample_id")
        r = np.corrcoef(
            got.loc[truth.index, "mature_colonocyte_fraction"], truth["differentiated"]
        )[0, 1]
        assert r > 0.9, f"{method} recovers the mature fraction at only r={r:.3f}"


def test_recovery_survives_multiplicative_noise(truth, linear_reference):
    bulk = _bulk_from(truth, linear_reference, noise=0.2)
    long, _ = deconvolve_cohort(bulk, linear_reference, bulk_scale=LINEAR_CP10K)
    summary = summarise_fractions(long, "lineage")
    got = summary[summary["method"] == "nnls"].set_index("sample_id")
    r = np.corrcoef(
        got.loc[truth.index, "mature_colonocyte_fraction"], truth["differentiated"]
    )[0, 1]
    assert r > 0.8


# ---------------------------------------------------------------------------
# THE FINDING: the scale mismatch, and the gate that does not see it


def test_a_log_reference_against_linear_bulk_is_refused_before_it_runs(truth, linear_reference):
    log_reference = Reference(
        matrix=np.log1p(linear_reference.matrix), rung="lineage",
        scale=LOG1P_CP10K, source=Path("synthetic.parquet"),
    )
    bulk = _bulk_from(truth, linear_reference)
    with pytest.raises(DeconvolutionError, match="exactly 0.0 on every sample"):
        deconvolve_cohort(bulk, log_reference, bulk_scale=LINEAR_CP10K)


def test_the_mismatch_zeroes_the_predictor_on_the_COMMITTED_matrix(truth):
    """The input that justifies the refusal above, on the real artifact.

    Deliberately mislabels the committed log matrix as linear to get past the
    guard, then asserts BOTH halves of the finding: the mature fraction is
    constant zero, AND the non-epithelial aggregate -- the quantity the
    pre-committed instrument gate correlates against ABSOLUTE purity -- is
    recovered well enough to pass its own r >= 0.5 threshold.

    Without this the guard would be a rule with an anecdote behind it.
    """
    path = COMMITTED / "S_matrix_lineage_1.0.0.parquet"
    if not path.exists():
        pytest.skip("S matrix not committed")
    committed = load_reference(path, rung="lineage", scale=LOG1P_CP10K)
    linear = linearise(committed)

    # Bulk is a LINEAR mixture, which is what TPM is.
    frame = truth[linear.matrix.columns]
    bulk = pd.DataFrame(
        frame.to_numpy() @ linear.matrix.to_numpy().T,
        index=truth.index, columns=linear.matrix.index,
    )
    mislabelled = Reference(
        matrix=committed.matrix, rung="lineage",
        scale=LINEAR_CP10K, source=committed.source,
    )
    long, _ = deconvolve_cohort(
        bulk, mislabelled, bulk_scale=LINEAR_CP10K, methods=[default_methods()[0]]
    )
    summary = summarise_fractions(long, "lineage").set_index("sample_id")

    mature = summary.loc[truth.index, "mature_colonocyte_fraction"]
    assert (mature == 0.0).all(), (
        f"the mismatch no longer zeroes the mature fraction on the committed "
        f"matrix ({int((mature == 0.0).sum())}/{len(mature)} zero). The guard's "
        f"stated justification no longer holds -- re-derive it before relaxing it."
    )

    non_epi = summary.loc[truth.index, "non_epithelial_fraction"]
    truth_non_epi = truth[list(NON_EPITHELIAL)].sum(axis=1)
    r = float(np.corrcoef(non_epi, truth_non_epi)[0, 1])
    assert r >= 0.5, (
        f"the gate's aggregate came back at r={r:.3f}, below its own threshold. "
        f"If the gate CAN fail here, the argument for a separate predictor check "
        f"is weaker and this module's docstring must be corrected."
    )

    # And the documented repair works on the same bulk.
    fixed, _ = deconvolve_cohort(
        bulk, linear, bulk_scale=LINEAR_CP10K, methods=[default_methods()[0]]
    )
    got = summarise_fractions(fixed, "lineage").set_index("sample_id")
    assert got.loc[truth.index, "mature_colonocyte_fraction"].std() > MIN_FRACTION_SD


def test_a_clean_reference_survives_the_same_mismatch(truth, linear_reference):
    """The narrower claim, pinned so it cannot quietly widen.

    "A log reference zeroes the predictor" is FALSE in general. Where each cell
    type owns a private block of marker genes, the same mismatch still recovers
    the mature fraction. It takes the mismatch AND a reference whose columns sit
    close together -- which the committed one does and this one does not. If
    this test ever starts failing, the module docstring is overclaiming and must
    be widened deliberately rather than by drift.
    """
    log_reference = Reference(
        matrix=np.log1p(linear_reference.matrix), rung="lineage",
        scale=LINEAR_CP10K, source=Path("synthetic.parquet"),
    )
    bulk = _bulk_from(truth, linear_reference)
    long, _ = deconvolve_cohort(
        bulk, log_reference, bulk_scale=LINEAR_CP10K, methods=[default_methods()[0]]
    )
    summary = summarise_fractions(long, "lineage").set_index("sample_id")
    mature = summary.loc[truth.index, "mature_colonocyte_fraction"]
    assert not (mature == 0.0).all()
    assert mature.std() > MIN_FRACTION_SD


def test_the_committed_matrix_is_near_collinear_where_the_analysis_reads_it():
    """Why the committed matrix is the fragile one: 0.982 against 0.929.

    Log compression pulls `differentiated` almost parallel to
    `epithelial_unscored`, and that pair is exactly the epithelial-internal
    split Stage 4's predictor comes from. The instrument gate reads across the
    epithelial boundary instead, which is why it survives.
    """
    path = COMMITTED / "S_matrix_lineage_1.0.0.parquet"
    if not path.exists():
        pytest.skip("S matrix not committed")
    committed = load_reference(path, rung="lineage", scale=LOG1P_CP10K)

    def cosine(matrix: pd.DataFrame, a: str, b: str) -> float:
        x, y = matrix[a].to_numpy(float), matrix[b].to_numpy(float)
        return float(x @ y / (np.linalg.norm(x) * np.linalg.norm(y)))

    log_cos = cosine(committed.matrix, "differentiated", "epithelial_unscored")
    lin_cos = cosine(linearise(committed).matrix, "differentiated", "epithelial_unscored")
    assert round(log_cos, 3) == 0.982
    assert round(lin_cos, 3) == 0.929
    assert log_cos > lin_cos, "log compression no longer worsens this pair"


def test_check_predictor_refuses_a_constant_zero_column():
    summary = pd.DataFrame({
        "sample_id": [f"S{i}" for i in range(20)],
        "method": "nnls",
        "mature_colonocyte_fraction": 0.0,
    })
    check = check_predictor(summary, rung="lineage", method="nnls")
    assert check.verdict == "refused"
    assert not check.usable
    assert check.is_constant and check.n_exact_zero == 20
    assert "invariant 1" in check.detail


def test_check_predictor_refuses_a_constant_nonzero_column():
    """A constant 0.4 is as unusable as a constant 0.0 and less obvious."""
    summary = pd.DataFrame({
        "sample_id": [f"S{i}" for i in range(20)],
        "method": "nnls",
        "mature_colonocyte_fraction": 0.4,
    })
    assert check_predictor(summary, rung="lineage", method="nnls").verdict == "refused"


def test_check_predictor_flags_a_near_constant_column_without_refusing_it():
    rng = np.random.default_rng(3)
    summary = pd.DataFrame({
        "sample_id": [f"S{i}" for i in range(50)],
        "method": "nnls",
        "mature_colonocyte_fraction": 0.3 + rng.normal(0, 0.001, 50),
    })
    check = check_predictor(summary, rung="lineage", method="nnls")
    assert check.verdict == "degenerate"
    assert not check.usable


def test_check_predictor_passes_a_real_one(truth, linear_reference):
    bulk = _bulk_from(truth, linear_reference)
    long, _ = deconvolve_cohort(bulk, linear_reference, bulk_scale=LINEAR_CP10K)
    summary = summarise_fractions(long, "lineage")
    checks = [check_predictor(summary, rung="lineage", method=m)
              for m in summary["method"].unique()]
    assert all(c.usable for c in checks), [c.detail for c in checks]
    require_usable_predictor(checks)  # does not raise


def test_one_surviving_method_is_enough_but_none_is_not():
    good = pd.DataFrame({"sample_id": [f"S{i}" for i in range(20)], "method": "nusvr",
                         "mature_colonocyte_fraction": np.linspace(0.1, 0.6, 20)})
    bad = pd.DataFrame({"sample_id": [f"S{i}" for i in range(20)], "method": "nnls",
                        "mature_colonocyte_fraction": 0.0})
    ok = check_predictor(good, rung="lineage", method="nusvr")
    dead = check_predictor(bad, rung="lineage", method="nnls")
    require_usable_predictor([ok, dead])  # method-specific failure is reportable
    with pytest.raises(DeconvolutionError, match="no method produced"):
        require_usable_predictor([dead])


# ---------------------------------------------------------------------------
# The scale guard itself


def test_assert_scales_agree_catches_both_directions(linear_reference):
    assert_scales_agree(linear_reference, LINEAR_CP10K)  # does not raise
    with pytest.raises(DeconvolutionError, match="linear mixture model"):
        assert_scales_agree(linear_reference, LOG1P_CP10K)
    log_reference = Reference(
        matrix=np.log1p(linear_reference.matrix), rung="lineage",
        scale=LOG1P_CP10K, source=Path("s.parquet"),
    )
    with pytest.raises(DeconvolutionError, match="linear mixture model"):
        assert_scales_agree(log_reference, LINEAR_CP10K)
    assert_scales_agree(log_reference, LOG1P_CP10K)


def test_linearise_is_the_exact_inverse_and_says_where_it_came_from(linear_reference):
    log_reference = Reference(
        matrix=np.log1p(linear_reference.matrix), rung="lineage",
        scale=LOG1P_CP10K, source=Path("s.parquet"),
    )
    back = linearise(log_reference)
    assert back.scale == LINEAR_CP10K
    assert back.derived_from_scale == LOG1P_CP10K
    pd.testing.assert_frame_equal(back.matrix, linear_reference.matrix)
    assert_scales_agree(back, LINEAR_CP10K)


def test_linearising_an_already_linear_reference_is_a_no_op(linear_reference):
    assert linearise(linear_reference) is linear_reference


# ---------------------------------------------------------------------------
# The committed matrices


@pytest.mark.parametrize("rung", ["epithelial", "lineage", "crypt_position", "best4"])
def test_every_committed_s_matrix_loads_and_is_log_scale(rung):
    path = COMMITTED / f"S_matrix_{rung}_1.0.0.parquet"
    if not path.exists():
        pytest.skip(f"{path.name} not committed")
    ref = load_reference(path, rung=rung, scale=LOG1P_CP10K)
    assert ref.matrix.shape[0] == 800
    assert ref.matrix.index.name == "gene"
    # It is the mean of log1p(CP10K), so it cannot exceed log1p(1e4).
    assert ref.matrix.to_numpy().max() <= np.log1p(1e4)
    mature = mature_column(rung)
    assert mature is None or mature in ref.matrix.columns


def test_the_committed_matrix_is_rejected_against_linear_bulk():
    """The whole reason this module exists, asserted on the real artifact."""
    path = COMMITTED / "S_matrix_lineage_1.0.0.parquet"
    if not path.exists():
        pytest.skip("S matrix not committed")
    ref = load_reference(path, rung="lineage", scale=LOG1P_CP10K)
    with pytest.raises(DeconvolutionError, match="linear mixture model"):
        assert_scales_agree(ref, LINEAR_CP10K)
    assert_scales_agree(linearise(ref), LINEAR_CP10K)  # the documented repair


def test_a_positional_index_is_refused_rather_than_aligned_by_order(tmp_path):
    """23b1d83's leakage class: aligning by position instead of identity."""
    frame = pd.DataFrame(np.ones((10, 2)), columns=["differentiated", "immune"])
    path = tmp_path / "S_matrix_lineage_1.0.0.parquet"
    frame.to_parquet(path)
    with pytest.raises(DeconvolutionError, match="align by ORDER"):
        load_reference(path, rung="lineage")


def test_a_target_gene_in_the_reference_is_refused(tmp_path):
    frame = pd.DataFrame({
        "gene": ["ENSG1", "GUCA2A"],
        "differentiated": [1.0, 2.0], "immune": [0.5, 0.1],
    })
    path = tmp_path / "s.parquet"
    frame.to_parquet(path)
    with pytest.raises(DeconvolutionError, match="invariant 2"):
        load_reference(path, rung="lineage", targets=["GUCA2A"])


def test_duplicate_genes_are_refused(tmp_path):
    frame = pd.DataFrame({
        "gene": ["ENSG1", "ENSG1"], "differentiated": [1.0, 2.0], "immune": [0.5, 0.1],
    })
    path = tmp_path / "s.parquet"
    frame.to_parquet(path)
    with pytest.raises(DeconvolutionError, match="duplicate genes"):
        load_reference(path, rung="lineage")


# ---------------------------------------------------------------------------
# The join


def test_align_uses_identity_not_position(linear_reference, truth):
    """Shuffled bulk columns must give identical fractions, not different ones."""
    bulk = _bulk_from(truth, linear_reference)
    shuffled = bulk.sample(axis=1, frac=1.0, random_state=5)
    a, _ = deconvolve_cohort(bulk, linear_reference, bulk_scale=LINEAR_CP10K,
                             methods=[default_methods()[0]])
    b, _ = deconvolve_cohort(shuffled, linear_reference, bulk_scale=LINEAR_CP10K,
                             methods=[default_methods()[0]])
    pd.testing.assert_frame_equal(
        a.sort_values(["sample_id", "cell_type"]).reset_index(drop=True),
        b.sort_values(["sample_id", "cell_type"]).reset_index(drop=True),
    )


def test_a_disjoint_gene_space_is_refused_not_silently_empty(linear_reference, truth):
    bulk = _bulk_from(truth, linear_reference)
    bulk.columns = [f"SYMBOL{i}" for i in range(bulk.shape[1])]
    with pytest.raises(DeconvolutionError, match="share no genes"):
        align(bulk, linear_reference)


def test_too_few_shared_genes_is_refused(linear_reference, truth):
    bulk = _bulk_from(truth, linear_reference).iloc[:, :40]
    with pytest.raises(DeconvolutionError, match="high-dimensionality"):
        align(bulk, linear_reference)


# ---------------------------------------------------------------------------
# Invariant 1 at the rung with no maturity call


def test_the_epithelial_rung_yields_none_not_zero(truth, linear_reference):
    """One bin means no maturity call, and that is `None` with a reason."""
    assert mature_column("epithelial") is None
    bulk = _bulk_from(truth, linear_reference)
    long, _ = deconvolve_cohort(bulk, linear_reference, bulk_scale=LINEAR_CP10K,
                                methods=[default_methods()[0]])
    summary = summarise_fractions(long, "epithelial")
    assert summary["mature_colonocyte_fraction"].isna().all()
    assert (summary["estimability"] == "not_estimable").all()
    assert (summary["mature_colonocyte_fraction"] == 0.0).sum() == 0
    check = check_predictor(summary, rung="epithelial", method="nnls")
    assert check.verdict == "not_estimable" and not check.usable


# ---------------------------------------------------------------------------
# The driver, end to end. Phase-5 plumbing: S matrix -> synthetic bulk -> both
# deconvolvers -> versioned table, without touching the cluster.


def _write_synthetic_bulk(path: Path, truth: pd.DataFrame, reference: Reference) -> None:
    """A linear mixture on the reference's own gene index, scaled to TPM.

    ``assert_linear_scale`` reads the driver's input, so the fixture has to be a
    plausible linear matrix rather than any positive array.
    """
    values = truth[reference.matrix.columns].to_numpy() @ reference.matrix.to_numpy().T
    frame = pd.DataFrame(values, index=truth.index, columns=reference.matrix.index)
    frame = frame.div(frame.sum(axis=1), axis=0) * 1e6
    frame.to_parquet(path)


def test_the_driver_runs_end_to_end_and_stamps_a_table(truth, tmp_path):
    """S matrix -> bulk -> NNLS + nu-SVR -> versioned parquet, both rungs."""
    from src.bulk.run_deconvolution import main

    s_path = COMMITTED / "S_matrix_lineage_1.0.0.parquet"
    if not s_path.exists():
        pytest.skip("S matrix not committed")
    reference = linearise(load_reference(s_path, rung="lineage", scale=LOG1P_CP10K))

    bulk_path = tmp_path / "tcga_tpm_1.0.0.parquet"
    _write_synthetic_bulk(bulk_path, truth, reference)
    results = tmp_path / "results"

    code = main([
        "--rung", "lineage", "--bulk", str(bulk_path),
        "--s-matrix-dir", str(COMMITTED), "--linearise-reference",
        "--results-dir", str(results), "--allow-dirty",
    ])
    assert code == 0, "the driver refused a run it should have completed"

    written = sorted(results.glob("*/stage4_fractions.parquet"))
    assert written, f"no table under {results}"
    table = pd.read_parquet(written[0])
    assert set(table["method"]) == {"nnls", "nusvr"}
    assert (table["granularity_rung"] == "lineage").all()
    assert (table["mature_cell_type"] == "differentiated").all()
    assert table["mature_colonocyte_fraction"].between(0, 1).all()
    assert (table["estimability"] == "estimated").all()

    checks = pd.read_parquet(sorted(results.glob("*/stage4_predictor_checks.parquet"))[0])
    assert set(checks["method"]) == {"nnls", "nusvr"}
    assert (checks["verdict"] == "usable").all(), checks[["method", "detail"]].to_dict()


def test_the_driver_refuses_the_committed_scale_and_says_why(truth, tmp_path):
    """Without --linearise-reference the run stops before producing a number."""
    from src.bulk.run_deconvolution import main

    s_path = COMMITTED / "S_matrix_lineage_1.0.0.parquet"
    if not s_path.exists():
        pytest.skip("S matrix not committed")
    reference = linearise(load_reference(s_path, rung="lineage", scale=LOG1P_CP10K))
    bulk_path = tmp_path / "tcga_tpm_1.0.0.parquet"
    _write_synthetic_bulk(bulk_path, truth, reference)

    code = main([
        "--rung", "lineage", "--bulk", str(bulk_path),
        "--s-matrix-dir", str(COMMITTED),
        "--results-dir", str(tmp_path / "results"), "--allow-dirty",
    ])
    assert code == 2, "a scale mismatch must refuse, not warn"
    assert not list((tmp_path / "results").glob("*/stage4_fractions.parquet"))


def test_the_driver_records_a_refusal_rather_than_only_raising(truth, tmp_path):
    """A refused predictor is a Stage 4 result and gets a stamped table.

    The epithelial rung has no maturity call, so every method's predictor is
    not-estimable and `require_usable_predictor` fires. The run must still leave
    both tables behind: a traceback is not a record, and 'we ran it and got
    nothing' has to be as auditable as any other outcome.
    """
    from src.bulk.run_deconvolution import main

    s_path = COMMITTED / "S_matrix_epithelial_1.0.0.parquet"
    if not s_path.exists():
        pytest.skip("S matrix not committed")
    reference = linearise(load_reference(s_path, rung="epithelial", scale=LOG1P_CP10K))
    bulk_path = tmp_path / "tcga_tpm_1.0.0.parquet"
    epithelial_truth = pd.DataFrame(
        np.repeat([[0.4, 0.1, 0.2, 0.2, 0.1]], len(truth), axis=0),
        index=truth.index, columns=reference.matrix.columns,
    )
    _write_synthetic_bulk(bulk_path, epithelial_truth, reference)
    results = tmp_path / "results"

    code = main([
        "--rung", "epithelial", "--bulk", str(bulk_path),
        "--s-matrix-dir", str(COMMITTED), "--linearise-reference",
        "--results-dir", str(results), "--allow-dirty",
    ])
    assert code == 3
    table = pd.read_parquet(sorted(results.glob("*/stage4_fractions.parquet"))[0])
    assert table["mature_colonocyte_fraction"].isna().all()
    assert (table["mature_colonocyte_fraction"] == 0.0).sum() == 0  # invariant 1
    checks = pd.read_parquet(sorted(results.glob("*/stage4_predictor_checks.parquet"))[0])
    assert (checks["verdict"] == "not_estimable").all()
