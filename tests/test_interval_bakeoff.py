"""The per-patient interval, the bake-off, and the negative-control runner."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.harness.bakeoff import (
    marker_ranked_genes,
    rank_methods,
    run_bakeoff,
    signature_width_comparison,
)
from src.harness.controls import run_negative_controls, summarise_negative_controls
from src.harness.deconvolve import NNLSDeconvolver, NuSVRDeconvolver
from src.harness.interval import ci_width, within_patient_intrinsic_ci
from src.harness.pseudobulk import generate_pseudobulk
from src.harness.results import validate_harness_table
from src.reference.signature import LeakageError

MATURE = "mature_colonocyte"
TARGET = "GUCA2A"
TYPES = [MATURE, "stem", "stromal", "immune", "endothelial"]


# ---------------------------------------------------------------------------
# within_patient_intrinsic_ci
# ---------------------------------------------------------------------------


def _mature(n, mean, seed=0):
    return np.random.default_rng(seed).poisson(mean, size=n).astype(float)


def test_interval_brackets_the_point_estimate():
    lo, hi = within_patient_intrinsic_ci(
        _mature(200, 60.0), _mature(200, 30.0, seed=1),
        frac_mature_normal=0.4, frac_mature_tumour=0.4, seed=1,
    )
    point = 0.4 * (_mature(200, 30.0, seed=1).mean() - _mature(200, 60.0).mean())
    assert lo < point < hi


def test_interval_narrows_as_mature_cells_accumulate():
    """The dependence the cutpoints are calibrated on. Without it the whole
    coverage-vs-n curve is flat and no cutpoint exists."""
    widths = []
    for n in (10, 40, 160, 640):
        lo, hi = within_patient_intrinsic_ci(
            _mature(n, 60.0), _mature(n, 30.0, seed=1),
            frac_mature_normal=0.4, frac_mature_tumour=0.4, seed=2, n_boot=300,
        )
        widths.append(ci_width(lo, hi))
    assert widths == sorted(widths, reverse=True)
    # Roughly 1/sqrt(n): 64x more cells should be about 8x tighter.
    assert 4.0 < widths[0] / widths[-1] < 16.0


def test_interval_is_undefined_with_no_mature_cells():
    """None, not a zero-width interval. Invariant 1."""
    lo, hi = within_patient_intrinsic_ci(
        _mature(50, 60.0), np.array([]),
        frac_mature_normal=0.4, frac_mature_tumour=0.0, seed=1,
    )
    assert lo is None and hi is None
    assert np.isnan(ci_width(lo, hi))


def test_interval_is_reproducible_given_a_seed():
    kw = dict(frac_mature_normal=0.4, frac_mature_tumour=0.3, seed=7)
    a = within_patient_intrinsic_ci(_mature(80, 60.0), _mature(80, 30.0, seed=1), **kw)
    b = within_patient_intrinsic_ci(_mature(80, 60.0), _mature(80, 30.0, seed=1), **kw)
    assert a == b


def test_interval_rejects_nonsense_arguments():
    with pytest.raises(ValueError, match="n_boot"):
        within_patient_intrinsic_ci(
            _mature(10, 5.0), _mature(10, 5.0),
            frac_mature_normal=0.4, frac_mature_tumour=0.4, seed=1, n_boot=0,
        )
    with pytest.raises(ValueError, match="alpha"):
        within_patient_intrinsic_ci(
            _mature(10, 5.0), _mature(10, 5.0),
            frac_mature_normal=0.4, frac_mature_tumour=0.4, seed=1, alpha=1.5,
        )


# ---------------------------------------------------------------------------
# bake-off
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def cohort():
    rng = np.random.default_rng(5)
    n_bg = 300
    genes = [TARGET] + [f"BG{i:04d}" for i in range(n_bg)]
    base = rng.gamma(2.0, 6.0, size=(n_bg, len(TYPES)))
    for k in range(len(TYPES)):
        base[k * 40 : (k + 1) * 40, k] *= 8.0
    profiles = {
        t: np.concatenate([[60.0 if t == MATURE else 1.5], base[:, k]])
        for k, t in enumerate(TYPES)
    }
    rows, ctypes, patients = [], [], []
    for p in range(8):
        for t in TYPES:
            for _ in range(40):
                rows.append(rng.poisson(profiles[t]))
                ctypes.append(t)
                patients.append(f"P{p:02d}")
    return np.array(rows), np.array(ctypes), np.array(patients), genes


def _comp(mature_frac):
    rest = (1.0 - mature_frac) / (len(TYPES) - 1)
    return {MATURE: mature_frac} | {t: rest for t in TYPES if t != MATURE}


@pytest.fixture(scope="module")
def samples_and_signature(cohort):
    counts, ctypes, patients, genes = cohort
    from src.harness.bulk_recovery import reference_profiles

    train = [f"P{p:02d}" for p in range(6)]
    held = [f"P{p:02d}" for p in range(6, 8)]
    train_rows = np.isin(patients, train)
    signature = reference_profiles(
        counts[train_rows], ctypes[train_rows], genes, exclude_genes=[TARGET]
    )
    samples = [
        generate_pseudobulk(
            counts, ctypes, patients, genes,
            composition_normal=_comp(0.4), composition_tumour=_comp(f),
            shift={TARGET: 0.5}, held_out_patients=held, n_cells=1500, seed=100 + i,
        )
        for i, f in enumerate((0.35, 0.20, 0.10, 0.05))
    ]
    return samples, signature


def test_bakeoff_table_is_shaped_correctly(samples_and_signature):
    samples, signature = samples_and_signature
    table, skipped = run_bakeoff(
        samples, signature, [NNLSDeconvolver()], seed=1, target_genes=[TARGET]
    )
    validate_harness_table(table, "bakeoff")
    assert skipped == {}
    assert len(table) == len(samples) * len(TYPES)


def test_bakeoff_refuses_a_signature_containing_a_target(samples_and_signature):
    """Invariant 2 — the deconvolution must not be informed by the gene under test."""
    samples, signature = samples_and_signature
    leaky = pd.concat([signature, signature.iloc[[0]].rename(index={signature.index[0]: TARGET})])
    with pytest.raises(LeakageError, match="invariant 2"):
        run_bakeoff(samples, leaky, [NNLSDeconvolver()], seed=1, target_genes=[TARGET])


def test_bakeoff_ranks_methods_and_reports_the_mature_error(samples_and_signature):
    samples, signature = samples_and_signature
    table, _ = run_bakeoff(
        samples, signature, [NNLSDeconvolver(), NuSVRDeconvolver()],
        seed=1, target_genes=[TARGET],
    )
    ranked = rank_methods(table, mature_label=MATURE)
    assert set(ranked["method"]) == {"nnls", "nusvr"}
    assert ranked["rmse"].is_monotonic_increasing  # best first
    assert ranked["rmse_mature"].notna().all()


def test_bakeoff_reports_a_skipped_method_rather_than_dropping_it(samples_and_signature):
    samples, signature = samples_and_signature

    class Absent(NNLSDeconvolver):
        name = "cibersortx"

        def is_available(self):
            return False, "no token configured"

    _, skipped = run_bakeoff(
        samples, signature, [NNLSDeconvolver(), Absent()], seed=1, target_genes=[TARGET]
    )
    assert skipped == {"cibersortx": "no token configured"}


def test_bakeoff_raises_when_nothing_is_available(samples_and_signature):
    samples, signature = samples_and_signature

    class Absent(NNLSDeconvolver):
        name = "x"

        def is_available(self):
            return False, "nope"

    with pytest.raises(ValueError, match="no deconvolution method"):
        run_bakeoff(samples, signature, [Absent()], seed=1)


# ---------------------------------------------------------------------------
# signature width — §2.1 error #4
# ---------------------------------------------------------------------------


def test_marker_ranking_puts_specific_genes_first(samples_and_signature):
    _, signature = samples_and_signature
    ranked = marker_ranked_genes(signature)
    values = signature.to_numpy(dtype=float)
    spec = values.max(axis=1) / values.mean(axis=1)
    top = spec[[signature.index.get_loc(g) for g in ranked[:20]]]
    bottom = spec[[signature.index.get_loc(g) for g in ranked[-20:]]]
    assert top.min() >= bottom.max()


def test_signature_width_comparison_covers_every_requested_width(samples_and_signature):
    samples, signature = samples_and_signature
    widths = (11, 50, 200)
    out = signature_width_comparison(
        samples, signature, NNLSDeconvolver(), widths, seed=1, target_genes=[TARGET]
    )
    validate_harness_table(out, "bakeoff")
    assert set(out["n_signature_genes"]) == set(widths)


def test_signature_width_comparison_rejects_a_width_it_cannot_serve(samples_and_signature):
    samples, signature = samples_and_signature
    with pytest.raises(ValueError, match="exceeds"):
        signature_width_comparison(
            samples, signature, NNLSDeconvolver(), (99_999,), seed=1
        )


# ---------------------------------------------------------------------------
# negative controls, run end to end
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def controls_table(cohort):
    counts, ctypes, patients, genes = cohort
    # Two flat genes to stand in for housekeeping in the synthetic cohort.
    return run_negative_controls(
        counts, ctypes, patients, genes,
        target_gene=TARGET, housekeeping=["BG0290", "BG0291"],
        composition_normal=_comp(0.4), composition_tumour=_comp(0.15),
        n_cells=1200, n_replicates=6, seed=11,
    )


def test_controls_table_is_shaped_correctly(controls_table):
    validate_harness_table(controls_table, "controls")
    assert set(controls_table["control"]) == {"target", "housekeeping", "permuted"}


def test_permutation_collapses_the_intrinsic_term(controls_table):
    """If it survives, the estimator is reading something other than the labels."""
    summary = summarise_negative_controls(controls_table)
    permuted = summary[
        (summary["control"] == "permuted") & (summary["term"] == "intrinsic")
    ]
    assert permuted["relative_to_target"].iloc[0] < 0.30


def test_housekeeping_genes_show_neither_term(controls_table):
    summary = summarise_negative_controls(controls_table)
    hk = summary[summary["control"] == "housekeeping"]
    assert (hk["relative_to_target"] < 0.10).all()


def test_target_arm_actually_shows_an_intrinsic_term(controls_table):
    """The positive reference. Without it the other two arms prove nothing."""
    summary = summarise_negative_controls(controls_table)
    target = summary[(summary["control"] == "target") & (summary["term"] == "intrinsic")]
    assert target["median_abs"].iloc[0] > 1.0
