"""Random-effects meta-analysis, checked against properties it must have.

There is no reference implementation in this repo to compare against, so these
test identities rather than remembered numbers: with no heterogeneity the
random-effects estimate must equal the fixed-effect one, identical studies must
pool to their common value, and a study with a huge SE must barely move the
answer. Each is a way the arithmetic can be wrong that a single worked example
would not catch.
"""

from __future__ import annotations

import numpy as np
import pytest

from src.harness.meta import (
    MAX_I_SQUARED,
    MIN_STUDIES,
    MetaError,
    meta_analyse,
    premise_verdict,
    se_from_interval,
)

# ---------------------------------------------------------------------------
# Identities the arithmetic must satisfy


def test_identical_studies_pool_to_their_common_value_with_no_heterogeneity():
    result = meta_analyse([0.4] * 6, [0.1] * 6)
    assert result.pooled == pytest.approx(0.4)
    assert result.tau_squared == 0.0
    assert result.i_squared == 0.0
    assert result.q == pytest.approx(0.0, abs=1e-9)
    # Six studies at SE 0.1 carry the precision of one at 0.1/sqrt(6).
    assert result.se == pytest.approx(0.1 / np.sqrt(6))


def test_with_no_heterogeneity_random_effects_equals_fixed_effects():
    """tau^2 = 0 makes the two weightings identical. If they disagree there,
    the weights are wrong somewhere."""
    result = meta_analyse([0.30, 0.32, 0.28, 0.31], [0.10, 0.11, 0.09, 0.10])
    assert result.tau_squared == 0.0
    assert result.pooled == pytest.approx(result.pooled_fixed)


def test_a_precise_study_outweighs_an_imprecise_one():
    result = meta_analyse([1.0, 0.0, 0.0], [0.01, 10.0, 10.0])
    assert result.pooled > 0.9, "inverse-variance weighting is not being applied"


def test_heterogeneity_pulls_the_pooled_estimate_toward_the_unweighted_mean():
    """tau^2 > 0 flattens the weights: that is what random effects DO."""
    spread = meta_analyse([-0.7, 0.1, 0.6, 0.4, -0.3], [0.05] * 5)
    assert spread.tau_squared > 0
    assert spread.i_squared > 0.5
    assert spread.pooled == pytest.approx(np.mean([-0.7, 0.1, 0.6, 0.4, -0.3]), abs=0.05)


def test_i_squared_and_tau_squared_are_never_negative():
    """Q below its expectation means LESS variation than sampling predicts,
    which is zero between-study variance, not a negative one."""
    result = meta_analyse([0.5, 0.5, 0.5, 0.5], [1.0, 1.0, 1.0, 1.0])
    assert result.tau_squared >= 0.0
    assert result.i_squared >= 0.0


def test_the_interval_narrows_as_studies_accumulate():
    few = meta_analyse([0.4] * 3, [0.2] * 3)
    many = meta_analyse([0.4] * 14, [0.2] * 14)
    assert many.se < few.se
    assert (many.ci_high - many.ci_low) < (few.ci_high - few.ci_low)


# ---------------------------------------------------------------------------
# Refusals


def test_too_few_studies_is_refused():
    with pytest.raises(MetaError, match=f"{MIN_STUDIES} minimum"):
        meta_analyse([0.1, 0.2], [0.1, 0.1])


def test_non_finite_studies_are_dropped_not_propagated():
    """A NaN in a weighted mean takes the whole pooled value with it."""
    result = meta_analyse([0.4, np.nan, 0.4, 0.4, 0.4], [0.1, 0.1, 0.1, 0.1, 0.1])
    assert result.k == 4
    assert np.isfinite(result.pooled)
    assert result.pooled == pytest.approx(0.4)


def test_a_zero_standard_error_is_dropped_rather_than_infinitely_weighted():
    result = meta_analyse([9.9, 0.4, 0.4, 0.4], [0.0, 0.1, 0.1, 0.1])
    assert result.k == 3
    assert result.pooled == pytest.approx(0.4)


def test_mismatched_input_lengths_are_refused():
    with pytest.raises(MetaError, match="against"):
        meta_analyse([0.1, 0.2, 0.3], [0.1, 0.1])


def test_se_from_interval_inverts_a_symmetric_interval():
    assert se_from_interval(0.2, 0.8) == pytest.approx(0.3 / 1.959963984540054)
    assert np.isnan(se_from_interval(np.nan, 0.5))
    with pytest.raises(MetaError, match="inverted"):
        se_from_interval(0.8, 0.2)


# ---------------------------------------------------------------------------
# The verdict, and the state that only exists at this level


def test_heterogeneous_studies_refuse_a_pooled_reading_that_would_say_HOLDS():
    """The failure this gate exists for, at the k we will actually have.

    Fourteen studies whose controls swing between -0.7 and +0.6 -- no two
    agreeing on whether the arms are comparable, several individually beyond
    any tolerance -- pool to +0.002 with an interval of [-0.304, +0.308]. That
    sits ENTIRELY inside a 0.5 tolerance, so the straddle rule alone returns
    HOLDS. The pooled number is arithmetically correct and substantively
    meaningless.

    At k = 5 random effects widens the interval enough to catch this on its
    own, which is why the gate looks redundant on a small example. It is not
    redundant at fourteen: more studies shrink the standard error even when
    tau^2 is large, so the more data you have the more confidently the naive
    reading is wrong.
    """
    values = [-0.7, 0.6, -0.5, 0.55, -0.62, 0.58, -0.48,
              0.52, -0.55, 0.61, -0.6, 0.5, -0.45, 0.57]
    result = meta_analyse(values, [0.05] * len(values))
    assert result.k == 14

    # The pooled estimate looks entirely benign.
    assert abs(result.pooled) < 0.05
    reach = max(abs(result.ci_low), abs(result.ci_high))
    assert reach < 0.5, (
        "the interval no longer sits inside the tolerance, so this example no "
        "longer demonstrates what the gate is for -- rebuild it before relying "
        "on the gate's justification"
    )

    # Which is exactly why the gate has to fire.
    assert result.i_squared > MAX_I_SQUARED
    verdict, detail = premise_verdict(result, tolerance=0.5)
    assert verdict == "UNRESOLVED"
    assert "not estimating a common quantity" in detail


def test_a_tight_homogeneous_result_inside_the_tolerance_holds():
    result = meta_analyse([0.10, 0.12, 0.09, 0.11, 0.10], [0.03] * 5)
    verdict, _ = premise_verdict(result, tolerance=0.5)
    assert verdict == "HOLDS"


def test_a_tight_homogeneous_result_beyond_the_tolerance_refuses():
    result = meta_analyse([0.90, 0.92, 0.89, 0.91, 0.90], [0.03] * 5)
    verdict, detail = premise_verdict(result, tolerance=0.5)
    assert verdict == "REFUSED" and "not comparable" in detail


def test_an_interval_straddling_the_tolerance_is_undecided():
    result = meta_analyse([0.45, 0.55, 0.50, 0.48, 0.52], [0.20] * 5)
    verdict, detail = premise_verdict(result, tolerance=0.5)
    assert verdict == "UNRESOLVED" and "straddles" in detail


def test_a_negative_interval_straddling_the_tolerance_is_also_undecided():
    """|shift| is what the tolerance reads, so signs must not flip the logic.
    This is the bug the Stage 4 straddle check had."""
    result = meta_analyse([-0.45, -0.55, -0.50, -0.48, -0.52], [0.20] * 5)
    verdict, _ = premise_verdict(result, tolerance=0.5)
    assert verdict == "UNRESOLVED"


def test_the_three_verdicts_are_the_same_three_the_per_study_check_uses():
    """One vocabulary across both levels, or a reader has to learn two.

    Checked against the per-study module's SOURCE rather than its docstrings --
    the strings it actually returns are what a results table carries.
    """
    import inspect

    from src.reference.jobs import coexpression_silencing as per_study

    source = inspect.getsource(per_study)
    for state in ("REFUSED", "UNRESOLVED"):
        assert f'"{state}' in source or f"'{state}" in source, (
            f"the per-study check no longer emits {state!r}; the two levels "
            f"have drifted apart"
        )
