"""Kitagawa identity and the positivity rule — analytically known answers.

W4 week 2-3 requires the estimator be unit-tested on synthetic data with
analytically known answers. This is the scalar floor of that; W4 extends it.
"""

from __future__ import annotations

import pytest

from src.estimator.kitagawa import decompose
from src.harness.positivity import (
    CUTPOINTS,
    classify_estimability,
    gate_g4_verdict,
)

# A tumour that lost half its mature cells AND silenced the survivors.
CASE = dict(
    frac_mature_normal=0.40,
    frac_mature_tumour=0.20,
    mean_normal=10.0,
    mean_tumour=6.0,
)


def test_kitagawa_identity_holds_under_normal_weighting():
    d = decompose(**CASE, n_cells_mature=200, weighting="normal")
    total = (
        CASE["frac_mature_tumour"] * CASE["mean_tumour"]
        - CASE["frac_mature_normal"] * CASE["mean_normal"]
    )
    assert d.compositional + d.intrinsic + d.interaction == pytest.approx(total)


def test_kitagawa_identity_holds_under_tumour_weighting():
    d = decompose(**CASE, n_cells_mature=200, weighting="tumour")
    total = (
        CASE["frac_mature_tumour"] * CASE["mean_tumour"]
        - CASE["frac_mature_normal"] * CASE["mean_normal"]
    )
    assert d.compositional + d.intrinsic + d.interaction == pytest.approx(total)


def test_the_two_weightings_disagree_and_the_difference_is_the_interaction():
    """The split is not unique. This is why both are reported."""
    n = decompose(**CASE, n_cells_mature=200, weighting="normal")
    t = decompose(**CASE, n_cells_mature=200, weighting="tumour")
    assert n.compositional != pytest.approx(t.compositional)
    assert n.intrinsic != pytest.approx(t.intrinsic)


def test_pure_compositional_case_has_zero_intrinsic():
    """Cells left, survivors unchanged. Tier A's expectation."""
    d = decompose(0.40, 0.10, 10.0, 10.0, n_cells_mature=200)
    assert d.intrinsic == pytest.approx(0.0)
    assert d.interaction == pytest.approx(0.0)
    assert d.compositional < 0


def test_pure_intrinsic_case_has_zero_compositional():
    """Cells present, silenced. Tier B's expectation — MLH1."""
    d = decompose(0.40, 0.40, 10.0, 4.0, n_cells_mature=200)
    assert d.compositional == pytest.approx(0.0)
    assert d.interaction == pytest.approx(0.0)
    assert d.intrinsic < 0


def test_neither_case_is_flat():
    """Tier D — MS4A12, colonocyte-restricted yet maintained."""
    d = decompose(0.40, 0.40, 10.0, 10.0, n_cells_mature=200)
    assert (d.compositional, d.intrinsic, d.interaction) == pytest.approx((0.0, 0.0, 0.0))


def test_doubly_robust_is_not_silently_the_plain_version():
    with pytest.raises(NotImplementedError):
        decompose(**CASE, n_cells_mature=200, weighting="doubly_robust")


# ---------------------------------------------------------------------------
# Positivity — where the third segment comes from
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "n,expected",
    [
        (10_000, "ok"),
        (50, "ok"),
        (49, "wide_interval"),
        (20, "wide_interval"),
        (19, "not_estimable"),
        (0, "not_estimable"),
    ],
)
def test_provisional_cutpoints(n, expected):
    assert classify_estimability(n) == expected


def test_cutpoints_are_ordered():
    assert CUTPOINTS.ok > CUTPOINTS.wide > 0


def test_negative_cell_count_is_an_error_not_a_verdict():
    with pytest.raises(ValueError):
        classify_estimability(-1)


def test_gate_g4_passes_when_most_patients_are_estimable():
    v = gate_g4_verdict([200, 150, 80, 60, 5])
    assert v["passes"] and v["fraction_below"] == pytest.approx(0.2)


def test_gate_g4_fails_and_names_the_pre_committed_consequence():
    v = gate_g4_verdict([2, 3, 4, 200])
    assert not v["passes"]
    assert "headline" in v["consequence"]


def test_gate_g4_needs_patients():
    with pytest.raises(ValueError):
        gate_g4_verdict([])
