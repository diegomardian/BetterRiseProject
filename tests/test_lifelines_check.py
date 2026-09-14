from __future__ import annotations

from src.harness.lifelines_check import FIXED_TOLERANCE, compare


def test_independent_lifelines_comparison_is_valid_and_within_fixed_tolerance():
    table = compare(n_patients=1500, seed=5)
    assert len(table) == 8
    assert table["valid_pair"].all()
    assert table["within_fixed_tolerance"].all()
    assert (table["fixed_tolerance"] == FIXED_TOLERANCE).all()
    assert table["absolute_difference"].max() <= FIXED_TOLERANCE
