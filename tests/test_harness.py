"""The harness itself, before it is trusted to judge the estimator.

The harness makes claims about whether the estimator works. Nothing it says
means anything until it has been shown to (a) produce a self-consistent truth
and (b) reject something. Test order below follows that.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.common.io import DirtyTreeError, write_versioned_table
from src.harness.controls import (
    assert_housekeeping_are_not_panel_genes,
    housekeeping_panel,
    permutation_preserved_counts,
    permute_labels_within_patient,
)
from src.harness.deconvolve import NNLSDeconvolver, available_methods
from src.harness.positivity import classify_estimability
from src.harness.pseudobulk import generate_pseudobulk, patient_holdout
from src.harness.results import (
    HarnessTableError,
    empty_harness_table,
    validate_harness_table,
    write_harness_table,
)
from src.harness.truth import (
    GroundTruth,
    analytic_terms,
    assert_identity_closes,
    identity_residual,
)

MATURE = "mature_colonocyte"
TARGET = "GUCA2A"
GENES = [TARGET, "MLH1", "ACTB", "FILLER1", "FILLER2"]
TYPES = [MATURE, "stem", "stromal", "immune"]


# ---------------------------------------------------------------------------
# Fixtures — a small synthetic cohort. Real cells arrive in week 2 (KUL3).
# ---------------------------------------------------------------------------


@pytest.fixture
def cohort():
    """8 patients x 4 cell types x 40 cells, with GUCA2A high in mature cells."""
    rng = np.random.default_rng(7)
    rows, ctypes, patients = [], [], []
    for p in range(8):
        for t in TYPES:
            for _ in range(40):
                base = np.array([2.0, 8.0, 50.0, 5.0, 5.0])
                if t == MATURE:
                    base[0] = 60.0  # GUCA2A is a mature-cell marker
                rows.append(rng.poisson(base))
                ctypes.append(t)
                patients.append(f"P{p:02d}")
    return np.array(rows), np.array(ctypes), np.array(patients)


def _comp(mature_frac: float) -> dict[str, float]:
    rest = (1.0 - mature_frac) / 3.0
    return {MATURE: mature_frac, "stem": rest, "stromal": rest, "immune": rest}


# ---------------------------------------------------------------------------
# 1 · The truth must be self-consistent
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("weighting", ["normal", "tumour"])
@pytest.mark.parametrize("f_n", [0.05, 0.4])
@pytest.mark.parametrize("f_t", [0.0, 0.01, 0.4])
@pytest.mark.parametrize("shift", [0.25, 0.5, 0.8, 1.0, 1.5])
def test_identity_closes_across_the_grid(weighting, f_n, f_t, shift):
    terms = analytic_terms(f_n, f_t, mean_normal=60.0, shift=shift, weighting=weighting)
    assert identity_residual(terms) == pytest.approx(0.0, abs=1e-9)


@pytest.mark.parametrize("weighting", ["normal", "tumour"])
def test_null_shift_gives_exactly_zero_intrinsic(weighting):
    """Exactly zero, not small. The null must be a genuine no-op."""
    terms = analytic_terms(0.4, 0.1, mean_normal=60.0, shift=1.0, weighting=weighting)
    assert terms["intrinsic"] == 0.0
    assert terms["interaction"] == 0.0
    assert terms["compositional"] != 0.0  # cells still left


@pytest.mark.parametrize("weighting", ["normal", "tumour"])
def test_pure_intrinsic_case_has_zero_compositional(weighting):
    """Cells present, survivors silenced. Tier B's expectation."""
    terms = analytic_terms(0.4, 0.4, mean_normal=60.0, shift=0.5, weighting=weighting)
    assert terms["compositional"] == pytest.approx(0.0)
    assert terms["intrinsic"] < 0


def test_the_two_weightings_disagree():
    """The split is not unique. Both get reported; neither is 'the' answer."""
    kw = dict(frac_mature_normal=0.4, frac_mature_tumour=0.1, mean_normal=60.0, shift=0.5)
    n = analytic_terms(**kw, weighting="normal")
    t = analytic_terms(**kw, weighting="tumour")
    assert n["intrinsic"] != pytest.approx(t["intrinsic"])
    assert n["total"] == pytest.approx(t["total"])


def test_assert_identity_closes_actually_fires():
    with pytest.raises(AssertionError, match="does not close"):
        assert_identity_closes(
            {"compositional": 1.0, "intrinsic": 1.0, "interaction": 1.0, "total": 99.0}
        )


@pytest.mark.parametrize("bad", [-0.1, -1.0])
def test_negative_shift_is_rejected(bad):
    with pytest.raises(ValueError, match="negative"):
        analytic_terms(0.4, 0.2, 60.0, bad)


def test_ground_truth_rejects_a_composition_that_does_not_sum_to_one():
    with pytest.raises(ValueError, match="sums to"):
        GroundTruth(
            composition_normal={MATURE: 0.5, "stem": 0.2},
            composition_tumour=_comp(0.1),
            shift={TARGET: 0.5},
            n_cells_mature=10,
            patient_ids=("P00",),
            seed=1,
        )


def test_ground_truth_rejects_an_unknown_mature_label():
    with pytest.raises(ValueError, match="not a cell type"):
        GroundTruth(
            composition_normal=_comp(0.4),
            composition_tumour=_comp(0.1),
            shift={TARGET: 0.5},
            n_cells_mature=10,
            patient_ids=("P00",),
            seed=1,
            mature_label="goblet",
        )


# ---------------------------------------------------------------------------
# 2 · The generator
# ---------------------------------------------------------------------------


def test_patient_holdout_splits_patients_not_cells(cohort):
    _, _, patients = cohort
    train, held = patient_holdout(patients, n_held_out=3, seed=1)
    assert len(held) == 3
    assert not set(train) & set(held)
    assert set(train) | set(held) == set(patients.tolist())


def test_patient_holdout_is_deterministic(cohort):
    _, _, patients = cohort
    assert patient_holdout(patients, n_held_out=3, seed=1) == patient_holdout(
        patients, n_held_out=3, seed=1
    )


def test_generator_draws_only_from_held_out_patients(cohort):
    counts, ctypes, patients = cohort
    with pytest.raises(ValueError, match="not present in the data"):
        generate_pseudobulk(
            counts, ctypes, patients, GENES,
            composition_normal=_comp(0.4), composition_tumour=_comp(0.1),
            shift={TARGET: 0.5}, held_out_patients=["P99"], n_cells=200, seed=1,
        )


def test_generator_is_deterministic_given_a_seed(cohort):
    counts, ctypes, patients = cohort
    kw = dict(
        composition_normal=_comp(0.4), composition_tumour=_comp(0.1),
        shift={TARGET: 0.5}, held_out_patients=["P00", "P01"], n_cells=400, seed=42,
    )
    a = generate_pseudobulk(counts, ctypes, patients, GENES, **kw)
    b = generate_pseudobulk(counts, ctypes, patients, GENES, **kw)
    assert np.array_equal(a.bulk_tumour, b.bulk_tumour)
    assert np.array_equal(a.bulk_normal, b.bulk_normal)
    assert a.truth.n_cells_mature == b.truth.n_cells_mature


def test_null_shift_leaves_the_tumour_counts_untouched(cohort):
    """s = 1.0 is a no-op, so the realised intrinsic term is exactly zero."""
    counts, ctypes, patients = cohort
    sample = generate_pseudobulk(
        counts, ctypes, patients, GENES,
        composition_normal=_comp(0.4), composition_tumour=_comp(0.4),
        shift={TARGET: 1.0}, held_out_patients=["P00", "P01"], n_cells=400, seed=3,
    )
    assert sample.truth.is_null
    for w in ("normal", "tumour"):
        assert sample.truth.parametric[TARGET][w]["intrinsic"] == 0.0


def test_shift_reduces_the_target_and_leaves_other_genes_alone(cohort):
    counts, ctypes, patients = cohort
    kw = dict(
        composition_normal=_comp(0.4), composition_tumour=_comp(0.4),
        held_out_patients=["P00", "P01", "P02"], n_cells=800, seed=5,
    )
    null = generate_pseudobulk(counts, ctypes, patients, GENES, shift={TARGET: 1.0}, **kw)
    hit = generate_pseudobulk(counts, ctypes, patients, GENES, shift={TARGET: 0.25}, **kw)
    j = GENES.index(TARGET)
    k = GENES.index("ACTB")
    assert hit.bulk_tumour[j] < null.bulk_tumour[j] * 0.6
    assert hit.bulk_tumour[k] == null.bulk_tumour[k]  # untouched gene is identical


def test_realised_truth_is_recorded_alongside_parametric(cohort):
    counts, ctypes, patients = cohort
    s = generate_pseudobulk(
        counts, ctypes, patients, GENES,
        composition_normal=_comp(0.4), composition_tumour=_comp(0.1),
        shift={TARGET: 0.5}, held_out_patients=["P00", "P01"], n_cells=400, seed=9,
    )
    assert set(s.truth.parametric) == set(s.truth.realised) == {TARGET}
    for w in ("normal", "tumour"):
        assert_identity_closes(s.truth.realised[TARGET][w])
        # They agree in direction but not to the digit — that gap is the point.
        assert np.sign(s.truth.realised[TARGET][w]["intrinsic"]) == np.sign(
            s.truth.parametric[TARGET][w]["intrinsic"]
        )


def test_zero_mature_cells_is_not_estimable_and_does_not_crash(cohort):
    """The third segment's edge case. Compositional is still reported."""
    counts, ctypes, patients = cohort
    s = generate_pseudobulk(
        counts, ctypes, patients, GENES,
        composition_normal=_comp(0.4), composition_tumour=_comp(0.0),
        shift={TARGET: 0.5}, held_out_patients=["P00", "P01"], n_cells=400, seed=11,
    )
    assert s.truth.n_cells_mature == 0
    assert s.truth.expected_estimability() == "not_estimable"
    assert s.truth.parametric[TARGET]["normal"]["compositional"] != 0.0


def test_shift_naming_a_missing_gene_is_an_error(cohort):
    counts, ctypes, patients = cohort
    with pytest.raises(KeyError, match="not in the matrix"):
        generate_pseudobulk(
            counts, ctypes, patients, GENES,
            composition_normal=_comp(0.4), composition_tumour=_comp(0.1),
            shift={"NOTAGENE": 0.5}, held_out_patients=["P00"], n_cells=100, seed=1,
        )


def test_mature_count_falls_as_the_swept_fraction_falls(cohort):
    """The sweep's x-axis has to actually move the estimability verdict."""
    counts, ctypes, patients = cohort
    seen = []
    for f in (0.4, 0.1, 0.02, 0.0):
        s = generate_pseudobulk(
            counts, ctypes, patients, GENES,
            composition_normal=_comp(0.4), composition_tumour=_comp(f),
            shift={TARGET: 0.5}, held_out_patients=["P00", "P01"], n_cells=500, seed=13,
        )
        seen.append(s.truth.n_cells_mature)
    assert seen == sorted(seen, reverse=True)
    assert classify_estimability(seen[0]) == "ok"
    assert classify_estimability(seen[-1]) == "not_estimable"


# ---------------------------------------------------------------------------
# 3 · Can the harness reject anything?
# ---------------------------------------------------------------------------


def _broken_estimator(sample) -> float:
    """Always says 'no intrinsic effect'. Must not survive G3."""
    return 0.0


def _honest_estimator(sample) -> float:
    """Reads the realised truth. Stands in for a working estimator."""
    return sample.truth.realised[TARGET]["normal"]["intrinsic"]


def _g3_verdict(samples, estimator, *, tol=0.15) -> bool:
    """Passes only if the estimator tracks truth on samples where truth is nonzero."""
    truths, hats = [], []
    for s in samples:
        if s.truth.expected_estimability() != "ok":
            continue
        truths.append(s.truth.realised[TARGET]["normal"]["intrinsic"])
        hats.append(estimator(s))
    if not truths:
        raise AssertionError("no estimable samples — the G3 check would be vacuous")
    truths, hats = np.array(truths), np.array(hats)
    scale = np.abs(truths).mean()
    return bool(np.abs(hats - truths).mean() <= tol * scale)


@pytest.fixture
def estimable_samples(cohort):
    counts, ctypes, patients = cohort
    return [
        generate_pseudobulk(
            counts, ctypes, patients, GENES,
            composition_normal=_comp(0.4), composition_tumour=_comp(0.35),
            shift={TARGET: 0.5}, held_out_patients=["P00", "P01"], n_cells=600, seed=100 + i,
        )
        for i in range(8)
    ]


def test_harness_passes_an_honest_estimator(estimable_samples):
    assert _g3_verdict(estimable_samples, _honest_estimator)


def test_harness_rejects_a_broken_estimator(estimable_samples):
    """THE test. A harness that has never rejected anything is not evidence."""
    assert not _g3_verdict(estimable_samples, _broken_estimator)


# ---------------------------------------------------------------------------
# 4 · Negative controls
# ---------------------------------------------------------------------------


def test_permutation_preserves_per_patient_composition(cohort):
    _, ctypes, patients = cohort
    permuted = permute_labels_within_patient(ctypes, patients, seed=1)
    assert permutation_preserved_counts(ctypes, permuted, patients)


def test_permutation_actually_moves_labels(cohort):
    _, ctypes, patients = cohort
    permuted = permute_labels_within_patient(ctypes, patients, seed=1)
    assert (permuted != ctypes).mean() > 0.5


def test_permutation_is_deterministic(cohort):
    _, ctypes, patients = cohort
    a = permute_labels_within_patient(ctypes, patients, seed=4)
    b = permute_labels_within_patient(ctypes, patients, seed=4)
    assert np.array_equal(a, b)


def test_permutation_destroys_the_marker_signal(cohort):
    """After shuffling, mature cells are no longer GUCA2A-high."""
    counts, ctypes, patients = cohort
    j = GENES.index(TARGET)
    real_gap = counts[ctypes == MATURE, j].mean() - counts[ctypes != MATURE, j].mean()
    permuted = permute_labels_within_patient(ctypes, patients, seed=2)
    perm_gap = counts[permuted == MATURE, j].mean() - counts[permuted != MATURE, j].mean()
    assert real_gap > 40
    assert abs(perm_gap) < real_gap * 0.1


def test_housekeeping_controls_are_not_panel_or_label_genes():
    assert_housekeeping_are_not_panel_genes()
    assert len(housekeeping_panel()) >= 8
    assert "ACTB" not in housekeeping_panel(exclude=["ACTB"])


# ---------------------------------------------------------------------------
# 5 · Deconvolution adapters
# ---------------------------------------------------------------------------


def test_nnls_recovers_fractions_from_a_clean_mixture():
    rng = np.random.default_rng(0)
    sig = pd.DataFrame(
        rng.gamma(2.0, 5.0, size=(600, 4)),
        columns=TYPES,
        index=[f"G{i:04d}" for i in range(600)],
    )
    true = np.array([0.4, 0.3, 0.2, 0.1])
    # Mix on the column-normalised signature, which is what the adapter fits.
    normalised = sig.to_numpy() / sig.sum(axis=0).to_numpy()
    bulk = normalised @ (true * 1e6)
    hat = NNLSDeconvolver().fit_predict(bulk, sig)
    assert hat.sum() == pytest.approx(1.0)
    np.testing.assert_allclose(hat.loc[TYPES].to_numpy(), true, atol=0.02)


def test_deconvolver_rejects_a_gene_index_mismatch():
    sig = pd.DataFrame(np.ones((10, 2)), columns=["a", "b"])
    with pytest.raises(ValueError, match="shared gene index"):
        NNLSDeconvolver().fit_predict(np.ones(9), sig)


def test_available_methods_reports_skips_by_name():
    class Missing(NNLSDeconvolver):
        name = "cibersortx"

        def is_available(self):
            return False, "no token configured"

    usable, skipped = available_methods([NNLSDeconvolver(), Missing()])
    assert [m.name for m in usable] == ["nnls"]
    assert skipped == {"cibersortx": "no token configured"}


# ---------------------------------------------------------------------------
# 6 · Persistence, without touching the frozen schema
# ---------------------------------------------------------------------------


def test_harness_table_shape_is_validated():
    df = empty_harness_table("calibration")
    validate_harness_table(df, "calibration")
    with pytest.raises(HarnessTableError, match="unexpected columns"):
        validate_harness_table(df.assign(pvalue=[]), "calibration")
    with pytest.raises(HarnessTableError, match="missing columns"):
        validate_harness_table(df.drop(columns=["verdict"]), "calibration")


def test_harness_table_writes_with_provenance(tmp_path):
    df = pd.DataFrame(
        [
            {
                "n_cells_mature": 50, "shift": 0.5, "coverage": 0.93,
                "discrimination": 0.82, "median_ci_width": 0.4,
                "n_replicates": 100, "verdict": "ok",
            }
        ]
    )
    path = write_harness_table(
        df, "calibration", seed=20260815, results_dir=tmp_path,
        extra_meta={"skipped_methods": {"cibersortx": "no token"}}, allow_dirty=True,
    )
    import json

    meta = json.loads((path.parent / "harness_calibration.meta.json").read_text())
    assert meta["seed"] == 20260815
    assert meta["harness_table_kind"] == "calibration"
    assert meta["skipped_methods"] == {"cibersortx": "no token"}
    assert meta["git_sha"]


def test_versioned_writer_refuses_a_dirty_tree(tmp_path, monkeypatch):
    monkeypatch.setattr("src.common.io.provenance_record", lambda **kw: {
        "date": "2026-08-15", "git_sha_short": "abc1234", "git_dirty": True,
    })
    with pytest.raises(DirtyTreeError, match="dirty"):
        write_versioned_table(pd.DataFrame({"a": [1]}), "x", seed=1, results_dir=tmp_path)


def test_harness_tables_do_not_go_through_the_frozen_schema():
    """The frozen contract is for decomposition rows only."""
    from src.harness.results import ATTENUATION_COLUMNS
    from src.schema import REQUIRED_COLUMNS

    assert set(ATTENUATION_COLUMNS) - set(REQUIRED_COLUMNS)
