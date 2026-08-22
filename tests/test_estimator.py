"""Kitagawa identity and the positivity rule — analytically known answers.

W4 week 2-3 requires the estimator be unit-tested on synthetic data with
analytically known answers. This is the scalar floor of that; W4 extends it.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.estimator.kitagawa import (
    ADDITIVE_WEIGHTINGS,
    attach_intrinsic_ci,
    bootstrap_over_patients,
    decompose,
    decompose_cohort,
    identity_residual,
)
from src.harness.positivity import (
    CUTPOINTS,
    classify_estimability,
    gate_g4_verdict,
)
from src.schema import WEIGHTINGS, coerce_results, validate_results

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


def test_pooled_reference_closes_on_two_terms_not_three():
    """docs/open_decisions.md #9. The pooled-reference split ("doubly_robust")
    puts half the cross term into each arm, so compositional + intrinsic is
    already the total and the third column is disclosure, not a summand."""
    d = decompose(**CASE, n_cells_mature=200, weighting="doubly_robust")
    total = (
        CASE["frac_mature_tumour"] * CASE["mean_tumour"]
        - CASE["frac_mature_normal"] * CASE["mean_normal"]
    )
    assert d.compositional + d.intrinsic == pytest.approx(total)
    assert d.compositional + d.intrinsic + d.interaction != pytest.approx(total)
    assert "doubly_robust" not in ADDITIVE_WEIGHTINGS


def test_pooled_reference_reports_the_real_cross_term_not_zero():
    """The bug W2 caught: `interaction = 0.0` read as "no interaction here"
    when the honest statement is "it is in the two arms, half each"."""
    dr = decompose(**CASE, n_cells_mature=200, weighting="doubly_robust")
    n = decompose(**CASE, n_cells_mature=200, weighting="normal")

    assert dr.interaction != 0.0
    assert dr.interaction == pytest.approx(n.interaction)


def test_the_pooled_arms_are_the_normal_arms_plus_half_the_cross_term():
    """Makes the folding legible: this is the algebra decision #9 turns on, so
    it is asserted rather than left in a docstring."""
    n = decompose(**CASE, n_cells_mature=200, weighting="normal")
    dr = decompose(**CASE, n_cells_mature=200, weighting="doubly_robust")
    half = n.interaction / 2.0

    assert dr.compositional == pytest.approx(n.compositional + half)
    assert dr.intrinsic == pytest.approx(n.intrinsic + half)


def test_identity_residual_uses_each_weightings_own_identity():
    """The helper exists so callers stop hard-coding a three-column sum."""
    total = (
        CASE["frac_mature_tumour"] * CASE["mean_tumour"]
        - CASE["frac_mature_normal"] * CASE["mean_normal"]
    )
    for weighting in WEIGHTINGS:
        d = decompose(**CASE, n_cells_mature=200, weighting=weighting)
        assert identity_residual(d, total) == pytest.approx(0.0, abs=1e-12)


def test_identity_residual_refuses_a_none_arm_rather_than_answering_zero():
    """Invariant 1 again: an unestimable term has no residual, not a zero one."""
    d = decompose(**CASE, n_cells_mature=3, weighting="normal")
    unestimable = type(d)(
        compositional=d.compositional,
        intrinsic=None,
        interaction=d.interaction,
        weighting="normal",
        n_cells_mature=3,
    )
    with pytest.raises(ValueError, match="not a residual of zero"):
        identity_residual(unestimable, -3.6)


def test_pooled_reference_is_not_silently_either_plain_version():
    """The whole point of the third weighting: it disagrees with both."""
    n = decompose(**CASE, n_cells_mature=200, weighting="normal")
    t = decompose(**CASE, n_cells_mature=200, weighting="tumour")
    dr = decompose(**CASE, n_cells_mature=200, weighting="doubly_robust")
    assert dr.compositional != pytest.approx(n.compositional)
    assert dr.compositional != pytest.approx(t.compositional)
    assert dr.intrinsic != pytest.approx(n.intrinsic)
    assert dr.intrinsic != pytest.approx(t.intrinsic)


def test_the_weightings_agree_in_degenerate_cases():
    """Tier A/B-style cases (only one side moves) collapse the reference
    choice, so all three weightings should agree. The cross term is genuinely
    zero here -- one of its two factors is -- so the pooled split's disclosed
    interaction is zero for a real reason, not by construction."""
    pure_compositional = dict(
        frac_mature_normal=0.40, frac_mature_tumour=0.10, mean_normal=10.0, mean_tumour=10.0
    )
    for weighting in WEIGHTINGS:
        d = decompose(**pure_compositional, n_cells_mature=200, weighting=weighting)
        assert d.intrinsic == pytest.approx(0.0)
        assert d.interaction == pytest.approx(0.0)

    pure_intrinsic = dict(
        frac_mature_normal=0.40, frac_mature_tumour=0.40, mean_normal=10.0, mean_tumour=4.0
    )
    for weighting in WEIGHTINGS:
        d = decompose(**pure_intrinsic, n_cells_mature=200, weighting=weighting)
        assert d.compositional == pytest.approx(0.0)
        assert d.interaction == pytest.approx(0.0)


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


# ---------------------------------------------------------------------------
# decompose_cohort — per-patient wrapper, schema-shaped
# ---------------------------------------------------------------------------


def _summary_frame():
    """Three synthetic patients spanning the positivity rule: comfortably
    estimable, wide-interval, and below the not-estimable cutpoint."""
    return pd.DataFrame(
        [
            dict(
                patient_id="P1", study_id="LEE_SMC", gene="GUCA2A",
                granularity_rung="lineage", labeling_axis="stem_pole",
                frac_mature_normal=0.40, frac_mature_tumour=0.10,
                mean_normal=10.0, mean_tumour=10.0, n_cells_mature=200,
            ),
            dict(
                patient_id="P2", study_id="LEE_SMC", gene="GUCA2A",
                granularity_rung="lineage", labeling_axis="stem_pole",
                frac_mature_normal=0.40, frac_mature_tumour=0.15,
                mean_normal=9.0, mean_tumour=9.0, n_cells_mature=30,
            ),
            dict(
                patient_id="P3", study_id="LEE_SMC", gene="GUCA2A",
                granularity_rung="lineage", labeling_axis="stem_pole",
                frac_mature_normal=0.40, frac_mature_tumour=0.05,
                mean_normal=11.0, mean_tumour=11.0, n_cells_mature=5,
            ),
        ]
    )  # fmt: skip


def test_decompose_cohort_fans_out_three_weightings_per_patient():
    out = decompose_cohort(_summary_frame())
    assert set(out["weighting"]) == set(WEIGHTINGS)
    assert len(out) == 3 * len(WEIGHTINGS)


def test_decompose_cohort_nulls_intrinsic_below_the_wide_cutpoint():
    """P3 has 5 mature cells -- below the wide-interval floor. Intrinsic is
    None, not zero; compositional stays estimable (README §Positivity)."""
    out = decompose_cohort(_summary_frame())
    p3 = out[out["patient_id"] == "P3"]
    assert (p3["estimability"] == "not_estimable").all()
    assert p3["intrinsic"].isna().all()
    assert p3["compositional"].notna().all()


def test_decompose_cohort_flags_wide_interval_without_dropping_the_estimate():
    out = decompose_cohort(_summary_frame())
    p2 = out[out["patient_id"] == "P2"]
    assert (p2["estimability"] == "wide_interval").all()
    assert p2["intrinsic"].notna().all()


def test_decompose_cohort_estimability_does_not_vary_by_weighting():
    """The mature-cell count, and therefore the verdict, doesn't depend on
    which weighting is used to report the split."""
    out = decompose_cohort(_summary_frame())
    for _patient_id, group in out.groupby("patient_id"):
        assert group["estimability"].nunique() == 1


def test_decompose_cohort_output_validates_against_the_frozen_schema():
    out = decompose_cohort(_summary_frame())
    validate_results(coerce_results(out))


def test_decompose_cohort_rejects_a_summary_missing_columns():
    bad = _summary_frame().drop(columns=["mean_tumour"])
    with pytest.raises(ValueError, match="missing column"):
        decompose_cohort(bad)


# ---------------------------------------------------------------------------
# bootstrap_over_patients / attach_intrinsic_ci
# ---------------------------------------------------------------------------


def _bootstrap_cohort(n_patients=8, seed=0):
    """Patients with heterogeneous but comfortably-estimable mature counts, so
    the bootstrap has something to resample without hitting NaNs."""
    rng = np.random.default_rng(seed)
    rows = [
        dict(
            patient_id=f"P{i}", study_id="LEE_SMC", gene="MLH1",
            granularity_rung="lineage", labeling_axis="stem_pole",
            frac_mature_normal=0.40, frac_mature_tumour=0.40,
            mean_normal=10.0, mean_tumour=float(rng.normal(4.0, 0.5)),
            n_cells_mature=int(rng.integers(80, 200)),
        )
        for i in range(n_patients)
    ]  # fmt: skip
    return pd.DataFrame(rows)


def test_bootstrap_ci_is_bounded_by_the_observed_patients_estimates():
    """Resampling patients can only reweight the observed set, never invent a
    value outside it -- the defining difference from resampling cells."""
    summary = _bootstrap_cohort()
    point = decompose_cohort(summary)
    ci = bootstrap_over_patients(summary, n_boot=200, seed=0)

    observed = point.loc[point["weighting"] == "normal", "intrinsic"]
    band = ci[(ci["weighting"] == "normal") & (ci["term"] == "intrinsic")].iloc[0]
    assert observed.min() - 1e-9 <= band["ci_low"] <= band["ci_high"] <= observed.max() + 1e-9


def test_bootstrap_requires_an_explicit_seed():
    with pytest.raises(TypeError):
        bootstrap_over_patients(_bootstrap_cohort(), n_boot=10)  # type: ignore[call-arg]


def test_bootstrap_rejects_non_positive_n_boot():
    with pytest.raises(ValueError, match="n_boot"):
        bootstrap_over_patients(_bootstrap_cohort(), n_boot=0, seed=0)


def test_bootstrap_is_reproducible_given_the_same_seed():
    summary = _bootstrap_cohort()
    a = bootstrap_over_patients(summary, n_boot=50, seed=42)
    b = bootstrap_over_patients(summary, n_boot=50, seed=42)
    pd.testing.assert_frame_equal(a, b)


def test_bootstrap_returns_one_row_per_group_weighting_term():
    summary = _bootstrap_cohort()  # one (study, gene, rung, axis) key
    ci = bootstrap_over_patients(summary, n_boot=20, seed=0)
    assert len(ci) == len(WEIGHTINGS) * 3


def test_attach_intrinsic_ci_broadcasts_the_cohort_band_to_every_patient():
    summary = _bootstrap_cohort()
    point = decompose_cohort(summary)
    ci = bootstrap_over_patients(summary, n_boot=50, seed=0)
    merged = attach_intrinsic_ci(point, ci)

    band = ci[(ci["weighting"] == "normal") & (ci["term"] == "intrinsic")].iloc[0]
    normal_rows = merged[merged["weighting"] == "normal"]
    assert np.allclose(normal_rows["ci_low"].to_numpy(dtype=float), float(band["ci_low"]))
    assert np.allclose(normal_rows["ci_high"].to_numpy(dtype=float), float(band["ci_high"]))


def test_attach_intrinsic_ci_output_validates_against_the_frozen_schema():
    summary = _bootstrap_cohort()
    point = decompose_cohort(summary)
    ci = bootstrap_over_patients(summary, n_boot=30, seed=0)
    merged = attach_intrinsic_ci(point, ci)
    validate_results(coerce_results(merged))
