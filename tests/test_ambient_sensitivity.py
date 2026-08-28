"""Ambient sensitivity — turning "degraded correction" into a number. W2, §5.7.

The finding these pin is an asymmetry: contamination manufactures **intrinsic**
signal out of **compositional** truth, and never the other way round. If that
ever stops being true, the gate memo's §10 is wrong and so is what W2 told the
gate about G1.
"""

from __future__ import annotations

import numpy as np
import pytest

from src.harness.ambient_sensitivity import (
    W1_EXCLUSION_THRESHOLD,
    AmbientRegime,
    ambient_sweep,
    contaminate,
    cost_at_threshold,
    default_regimes,
    soup_rate,
    summarise_ambient,
)

MATURE = "mature_colonocyte"
TARGET = "GUCA2A"
GENES = [TARGET, "MLH1", "ACTB"]
TYPES = [MATURE, "stem", "stromal", "immune"]

#: GUCA2A is mature-restricted — that is what makes the two arms' soups differ,
#: and the whole mechanism depends on it.
PROFILE = {
    MATURE: [60.0, 8.0, 50.0],
    "stem": [2.0, 8.0, 50.0],
    "stromal": [1.0, 6.0, 55.0],
    "immune": [1.0, 7.0, 45.0],
}
GRID = (0.0, 0.10, 0.30)


@pytest.fixture(scope="module")
def cohort():
    rng = np.random.default_rng(11)
    rows, ctypes, patients = [], [], []
    for p in range(6):
        for t in TYPES:
            for _ in range(60):
                rows.append(rng.poisson(PROFILE[t]))
                ctypes.append(t)
                patients.append(f"P{p:02d}")
    return np.array(rows), np.array(ctypes), np.array(patients)


@pytest.fixture(scope="module")
def summary(cohort):
    counts, ctypes, patients = cohort
    sweep = ambient_sweep(
        counts,
        ctypes,
        patients,
        GENES,
        target_gene=TARGET,
        seed=20260827,
        held_out_patients=["P04", "P05"],
        grid=GRID,
        n_cells=1200,
        n_replicates=8,
    )
    return summarise_ambient(sweep)


def _row(summary, regime, term, fraction):
    hit = summary[
        (summary["regime"] == regime)
        & (summary["term"] == term)
        & np.isclose(summary["ambient_fraction"], fraction)
    ]
    assert len(hit) == 1, (regime, term, fraction)
    return hit.iloc[0]


# ---------------------------------------------------------------------------
# The finding
# ---------------------------------------------------------------------------


def test_contamination_manufactures_intrinsic_from_compositional_truth(summary):
    """The mechanism the sweep exists to measure.

    In a world where mature cells are gone and the survivors are untouched, the
    true intrinsic term is exactly zero. Each arm's mature cells are pulled
    toward *their own arm's* soup, the normal arm's soup is richer in the gene,
    and the means separate where the truth says they must not.
    """
    at_zero = _row(summary, "compositional_only", "intrinsic", 0.0)
    at_cap = _row(summary, "compositional_only", "intrinsic", W1_EXCLUSION_THRESHOLD)
    assert at_zero["truth_is_zero"]
    assert at_cap["median_manufactured"] > at_zero["median_manufactured"]
    assert abs(at_cap["median_contaminated"]) > abs(at_zero["median_contaminated"])


def test_contamination_never_manufactures_compositional_from_intrinsic_truth(summary):
    """The asymmetry, and it is structural rather than lucky.

    The compositional term is a function of the mature-fraction *difference*.
    Contamination in this model moves means, not fractions, so where the two
    arms have the same mature fraction the compositional term stays identically
    zero however much soup is added. Ambient can invent silencing; it cannot
    invent depletion.
    """
    for fraction in GRID:
        row = _row(summary, "intrinsic_only", "compositional", fraction)
        assert row["truth_is_zero"]
        assert row["median_contaminated"] == pytest.approx(0.0, abs=1e-9)
        assert row["median_manufactured"] == pytest.approx(0.0, abs=1e-9)


def test_manufacture_grows_with_contamination(summary):
    made = [
        _row(summary, "compositional_only", "intrinsic", f)["median_manufactured"]
        for f in GRID
    ]
    assert made == sorted(made)


def test_real_terms_are_attenuated_but_not_destroyed(summary):
    """Contamination shrinks a true term. It should not annihilate or flip it."""
    for regime, term in (
        ("compositional_only", "compositional"),
        ("intrinsic_only", "intrinsic"),
        ("both", "intrinsic"),
        ("both", "compositional"),
    ):
        retention = [_row(summary, regime, term, f)["median_retention"] for f in GRID]
        assert retention[0] == pytest.approx(1.0)
        assert retention == sorted(retention, reverse=True), (regime, term)
        assert retention[-1] > 0.5, (regime, term)


def test_the_cost_at_the_pre_committed_exclusion_cap_is_reported(summary):
    """Decision #16 keeps a sample at exactly 10%, so 10% is a real worst case."""
    at_cap = cost_at_threshold(summary)
    assert not at_cap.empty
    assert set(at_cap["regime"]) == {"compositional_only", "intrinsic_only", "both"}
    made = _row(summary, "compositional_only", "intrinsic", W1_EXCLUSION_THRESHOLD)
    # Small enough not to dominate a real signal, large enough to be worth saying.
    assert 0.0 < made["median_manufactured"] < 0.25


# ---------------------------------------------------------------------------
# The bug this module nearly shipped
# ---------------------------------------------------------------------------


def test_truth_is_zero_comes_from_the_design_not_from_the_numbers():
    """The realised intrinsic term in a compositional-only world is noise, not 0.

    Inferring "the truth is zero here" from the realised value calls that noise
    non-zero, forms a retention ratio against a ~0.05 denominator, and reports a
    confident 2.4x where the honest statement is "the truth is zero and 0.79
    appeared". Same mistake as scoring coverage against realised truth.
    """
    compositional, intrinsic, both = default_regimes(TARGET, mature_label=MATURE)
    assert compositional.parametric_zero_terms(TARGET, MATURE) == {
        "intrinsic",
        "interaction",
    }
    assert intrinsic.parametric_zero_terms(TARGET, MATURE) == {
        "compositional",
        "interaction",
    }
    assert both.parametric_zero_terms(TARGET, MATURE) == set()


def test_summarise_refuses_a_sweep_that_lost_the_design_flag(summary, cohort):
    counts, ctypes, patients = cohort
    sweep = ambient_sweep(
        counts, ctypes, patients, GENES, target_gene=TARGET, seed=1,
        held_out_patients=["P04", "P05"], grid=(0.0,), n_cells=400, n_replicates=1,
    )
    with pytest.raises(ValueError, match="truth_is_zero"):
        summarise_ambient(sweep.drop(columns=["truth_is_zero"]))


# ---------------------------------------------------------------------------
# The contamination model
# ---------------------------------------------------------------------------


def test_zero_contamination_is_the_identity():
    values = np.array([0, 3, 17, 250])
    out = contaminate(values, fraction=0.0, rate=99.0, rng=np.random.default_rng(0))
    assert np.array_equal(out, values)


def test_contamination_keeps_counts_integral_and_non_negative():
    rng = np.random.default_rng(2)
    out = contaminate(np.array([0, 1, 5, 400]), fraction=0.2, rate=10.0, rng=rng)
    assert out.dtype.kind in "iu"
    assert (out >= 0).all()


def test_contamination_lifts_a_silenced_cell_off_the_floor():
    """Why a silenced mature cell stops reading as silenced. The whole point."""
    rng = np.random.default_rng(5)
    silenced = np.zeros(2000, dtype=int)
    out = contaminate(silenced, fraction=0.10, rate=40.0, rng=rng)
    assert out.mean() == pytest.approx(4.0, rel=0.15)


@pytest.mark.parametrize(
    ("fraction", "rate"), [(-0.1, 1.0), (1.5, 1.0), (0.1, -1.0)]
)
def test_contaminate_rejects_impossible_inputs(fraction, rate):
    with pytest.raises(ValueError):
        contaminate(np.array([1, 2]), fraction=fraction, rate=rate, rng=np.random.default_rng(0))


def test_soup_is_drawn_from_every_cell_not_only_the_mature_ones(cohort):
    """The single line the finding rests on.

    If the soup were the mature-cell mean, both arms would share it and nothing
    would separate. It is the sample average, so a depleted tumour has a poorer
    soup than its matched normal — and that difference is the artefact.
    """
    from src.harness.pseudobulk import generate_pseudobulk

    counts, ctypes, patients = cohort
    regime = AmbientRegime(
        name="probe",
        composition_normal={MATURE: 0.40, "stem": 0.20, "stromal": 0.20, "immune": 0.20},
        composition_tumour={MATURE: 0.05, "stem": 0.317, "stromal": 0.317, "immune": 0.316},
        shift={TARGET: 1.0},
        why="",
    )
    sample = generate_pseudobulk(
        counts, ctypes, patients, GENES,
        composition_normal=dict(regime.composition_normal),
        composition_tumour=dict(regime.composition_tumour),
        shift=dict(regime.shift), held_out_patients=["P04", "P05"],
        n_cells=1500, seed=7, mature_label=MATURE,
    )
    normal_soup = soup_rate(sample, TARGET, "normal")
    tumour_soup = soup_rate(sample, TARGET, "tumour")
    mature_mean = float(
        sample.drawn_expression[TARGET]["normal"][sample.drawn_is_mature["normal"]].mean()
    )
    assert normal_soup > tumour_soup           # the depleted arm has a poorer soup
    assert normal_soup < mature_mean           # and the soup is not the mature mean
