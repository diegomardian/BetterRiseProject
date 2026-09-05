"""Stage 4's variance arm, against inputs constructed to defeat each claim.

The arm's whole design is a response to issue #54: a raw R-squared measures the
assay floor, so the statistic is a percentile within an abundance-matched null.
The tests that matter most are the ones that reproduce #54's simulation and
check that the percentile is NOT fooled where the raw comparison was.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.bulk.prespec import load_prespec
from src.bulk.variance_arm import (
    Attrition,
    GeneFit,
    VarianceArmError,
    benjamini_hochberg,
    build_design,
    compare_to_null,
    gene_r_squared,
    negative_control_verdict,
    primary_verdict,
    resolve_r_squared_kinds,
    secondary_verdict,
)

RNG = np.random.default_rng(20260905)


def _design(n: int = 300, *, with_covariates: bool = False) -> pd.DataFrame:
    frame = pd.DataFrame({
        "sample_id": [f"S{i:04d}" for i in range(n)],
        "mature_colonocyte_fraction": RNG.uniform(0.05, 0.7, n),
    })
    if with_covariates:
        frame["purity"] = RNG.uniform(0.3, 0.9, n)
        frame["plate"] = RNG.choice([f"P{i}" for i in range(8)], n)
    return frame


def _gene(design: pd.DataFrame, *, share: float, floor: float = 0.0,
          name: str = "G") -> pd.Series:
    """Expression with a known share of variance from fraction, then floored.

    `floor` is the assay floor issue #54 is about: values below it read as
    noise, which destroys R-squared for a low-abundance gene whatever its
    biology.
    """
    f = design["mature_colonocyte_fraction"].to_numpy()
    z = (f - f.mean()) / f.std()
    signal = np.sqrt(share) * z
    noise = np.sqrt(max(0.0, 1 - share)) * RNG.normal(0, 1, len(f))
    values = 10.0 + signal + noise
    if floor:
        values = np.where(values < floor, floor + RNG.normal(0, 0.01, len(f)), values)
    return pd.Series(values, index=design.index, name=name)


# ---------------------------------------------------------------------------
# The R-squared itself


def test_marginal_r_squared_recovers_a_known_share():
    design = _design(600)
    fit = gene_r_squared(_gene(design, share=0.4), design, [])
    assert 0.33 < fit.marginal_r2 < 0.47


def test_a_gene_unrelated_to_fraction_returns_near_zero():
    design = _design(600)
    unrelated = pd.Series(RNG.normal(0, 1, len(design)), name="U")
    assert gene_r_squared(unrelated, design, []).marginal_r2 < 0.05


def test_partial_r_squared_removes_a_covariate_that_drives_both():
    """Purity drives fraction and expression; marginal is inflated, partial is not."""
    n = 500
    purity = RNG.uniform(0.3, 0.9, n)
    design = pd.DataFrame({
        "sample_id": [f"S{i}" for i in range(n)],
        "mature_colonocyte_fraction": 0.6 * (1 - purity) + RNG.normal(0, 0.03, n),
        "purity": purity,
    })
    expression = pd.Series(5.0 + 4.0 * (1 - purity) + RNG.normal(0, 0.5, n), name="G")
    fit = gene_r_squared(expression, design, ["purity"])
    assert fit.marginal_r2 > 0.5, "the confounded marginal should be large"
    assert fit.partial_r2 < 0.2, (
        f"partial R2 is {fit.partial_r2:.3f}; the covariate should absorb it"
    )


def test_the_same_function_serves_target_and_null():
    """Not a behavioural test -- an identity one. Two code paths that agree by
    review is exactly how a target ends up adjusted and its null not."""
    import inspect

    from src.bulk import variance_arm

    source = inspect.getsource(variance_arm.compare_to_null)
    assert "f.value(kind)" in source
    assert "lstsq" not in source, (
        "compare_to_null fits its own model. The null must go through "
        "gene_r_squared, the same function the target uses."
    )


# ---------------------------------------------------------------------------
# Issue #54: the percentile must not be fooled where the raw value was


def test_the_assay_floor_destroys_raw_r_squared_for_a_low_abundance_gene():
    """#54's simulation, reproduced. This is the failure the null exists for."""
    design = _design(600)
    unfloored = gene_r_squared(_gene(design, share=0.6), design, []).marginal_r2
    floored = gene_r_squared(_gene(design, share=0.6, floor=10.5), design, []).marginal_r2
    assert unfloored > 0.5
    assert floored < unfloored / 2, (
        "the floor no longer destroys R-squared, so the matched null's "
        "justification needs re-deriving before it is relaxed"
    )


def test_the_matched_null_percentile_survives_the_floor_that_broke_the_raw_value():
    """A floored gene with real signal still beats floored genes without it.

    The raw R-squared collapses; the percentile does not, because every null
    gene is floored the same way. This is the whole reason the statistic was
    changed, so it is asserted rather than assumed.
    """
    design = _design(600)
    target = gene_r_squared(_gene(design, share=0.6, floor=10.5, name="T"), design, [])
    nulls = [
        gene_r_squared(_gene(design, share=0.0, floor=10.5, name=f"N{i}"), design, [])
        for i in range(60)
    ]
    comparison = compare_to_null(target, nulls, "marginal")
    assert comparison.r2 < 0.35, "the target's raw R-squared should be depressed"
    assert comparison.exceeds_null, (
        f"percentile {comparison.percentile:.3f}: the floored target with real "
        f"signal must still exceed floored nulls without it"
    )


def test_a_gene_with_no_signal_does_not_exceed_its_null():
    design = _design(600)
    target = gene_r_squared(_gene(design, share=0.0, name="T"), design, [])
    nulls = [gene_r_squared(_gene(design, share=0.0, name=f"N{i}"), design, [])
             for i in range(60)]
    assert not compare_to_null(target, nulls, "marginal").exceeds_null


def test_a_null_too_small_to_be_a_percentile_is_refused():
    design = _design(200)
    target = gene_r_squared(_gene(design, share=0.3), design, [])
    nulls = [gene_r_squared(_gene(design, share=0.0, name=f"N{i}"), design, [])
             for i in range(5)]
    with pytest.raises(VarianceArmError, match="not one"):
        compare_to_null(target, nulls, "marginal")


# ---------------------------------------------------------------------------
# The pre-registered verdicts, applied verbatim


def _comparison(percentile: float, excess: float = 0.0, gene: str = "G",
                kind: str = "marginal"):
    from src.bulk.variance_arm import NullComparison

    return NullComparison(
        gene=gene, kind=kind, r2=0.1, n_null=200, null_median=0.05,
        null_p05=0.0, null_p95=0.2, percentile=percentile, excess=excess,
        exceeds_null=percentile > 0.95,
    )


def test_the_primary_verdict_matches_the_locked_spec_in_all_four_cases():
    within, above = _comparison(0.50, gene="GUCA2A"), _comparison(0.99, gene="CDX2")
    assert primary_verdict(within, above)[0] == "confirmed"
    assert primary_verdict(_comparison(0.99), _comparison(0.50))[0] == "disconfirmed"
    assert primary_verdict(_comparison(0.99), _comparison(0.99))[0] == "indeterminate"
    assert primary_verdict(_comparison(0.20), _comparison(0.20))[0] == "indeterminate"


def test_a_dead_instrument_reads_indeterminate_and_not_confirmed():
    """A constant predictor gives both genes R-squared 0 and neither exceeds.

    The locked spec's `disconfirmed_if` clause is what stops total instrument
    failure reading as the predicted result -- worth pinning, because the
    prediction's expected direction for GUCA2A is exactly what a broken
    instrument produces.
    """
    dead = _comparison(0.0, gene="GUCA2A"), _comparison(0.0, gene="CDX2")
    verdict, detail = primary_verdict(*dead)
    assert verdict == "indeterminate"
    assert "neither" in detail


def test_verdicts_refuse_to_mix_a_partial_with_a_marginal():
    with pytest.raises(VarianceArmError, match="cannot compare"):
        primary_verdict(_comparison(0.5, kind="marginal"),
                        _comparison(0.99, kind="partial"))


def test_the_secondary_verdict_is_an_excess_comparison():
    assert secondary_verdict(_comparison(0.5, excess=0.01),
                             _comparison(0.99, excess=0.20))[0] == "confirmed"
    assert secondary_verdict(_comparison(0.5, excess=0.30),
                             _comparison(0.99, excess=0.20))[0] == "disconfirmed"


def test_both_negative_controls_are_reported_and_neither_is_silent():
    clean = [GeneFit("ACTB", 300, 0.02, 0.01, 0.3), GeneFit("GAPDH", 300, 0.03, 0.02, 0.3)]
    assert negative_control_verdict(clean, 0.04, "marginal")[0] == "clean"

    verdict, detail = negative_control_verdict(
        [GeneFit("ACTB", 300, 0.31, 0.2, 0.3)], 0.04, "marginal")
    assert verdict == "breached" and "UPPER BOUND" in detail

    verdict, detail = negative_control_verdict(clean, 0.22, "marginal")
    assert verdict == "breached" and "INDETERMINATE" in detail


# ---------------------------------------------------------------------------
# Invariant 1 and the attrition table


def test_a_not_estimable_fraction_is_dropped_and_never_entered_as_zero():
    design = _design(100)
    design.loc[:19, "mature_colonocyte_fraction"] = np.nan
    attrition = Attrition()
    built = build_design(design, None, covariate_names=[], attrition=attrition)
    assert len(built) == 80
    assert (built["mature_colonocyte_fraction"] == 0.0).sum() == 0
    row = attrition.to_frame().iloc[0]
    assert row["n_dropped"] == 20 and "Invariant 1" in row["reason"]


def test_entering_them_as_zero_would_manufacture_the_signal():
    """The input the guard exists for: zeros correlate with everything.

    Not a test of our code path -- a demonstration that the dropped alternative
    is not neutral. Coding the unestimable as 0.0 creates a low-fraction group
    whose expression differs, which is the compositional signal the analysis is
    supposed to be measuring.
    """
    design = _design(300)
    honest = gene_r_squared(_gene(design, share=0.0, name="G"), design, [])
    coerced = design.copy()
    coerced.loc[:99, "mature_colonocyte_fraction"] = 0.0
    expression = _gene(design, share=0.0, name="G").copy()
    expression.iloc[:100] -= 2.0          # those samples happen to be lower
    manufactured = gene_r_squared(expression, coerced, [])
    assert honest.marginal_r2 < 0.05
    # Observed 0.237 against an honest <0.05: a five-fold inflation out of a
    # gene with no relationship to fraction at all. The ratio is the claim; the
    # absolute value depends on how far the coerced group happens to sit.
    assert manufactured.marginal_r2 > 5 * honest.marginal_r2, (
        f"coercion inflated R-squared from {honest.marginal_r2:.4f} to only "
        f"{manufactured.marginal_r2:.4f}. If zeros no longer manufacture a "
        f"signal here, re-derive invariant 1's justification before relying on it."
    )
    assert manufactured.marginal_r2 > 0.15


def test_the_design_refuses_a_missing_prespecified_covariate():
    design = _design(100, with_covariates=True)
    covariates = pd.DataFrame({"sample_id": design["sample_id"], "purity": 0.5})
    with pytest.raises(VarianceArmError, match="missing"):
        build_design(design[["sample_id", "mature_colonocyte_fraction"]], covariates,
                     covariate_names=["purity", "plate"], attrition=Attrition())


def test_too_few_surviving_samples_is_refused():
    design = _design(25)
    design.loc[:14, "mature_colonocyte_fraction"] = np.nan
    with pytest.raises(VarianceArmError, match="not interpretable"):
        build_design(design, None, covariate_names=[], attrition=Attrition())


# ---------------------------------------------------------------------------
# Multiplicity, and the gap in the lock


def test_benjamini_hochberg_matches_the_known_values():
    # m=3. Raw: a 0.01*3/1=0.03, c 0.03*3/2=0.045, b 0.04*3/3=0.04. The step-up
    # then enforces monotonicity from the largest down, so c is capped at b's
    # 0.04 -- an adjusted value may not exceed one from a larger raw p.
    got = benjamini_hochberg({"a": 0.01, "b": 0.04, "c": 0.03})
    assert got["a"] == pytest.approx(0.03)
    assert got["c"] == pytest.approx(0.04)
    assert got["b"] == pytest.approx(0.04)
    assert got["a"] <= got["c"] <= got["b"], "BH output must stay monotone"
    assert benjamini_hochberg({}) == {}


def test_both_r_squared_kinds_are_carried_because_the_lock_is_silent():
    """The locked spec asks for both and its arms say only "R-squared".

    Resolving that after seeing the numbers is the move pre-specification
    exists to prevent, so both are computed and both verdicts reported.
    """
    kinds = resolve_r_squared_kinds(load_prespec())
    assert set(kinds) == {"partial", "marginal"}
