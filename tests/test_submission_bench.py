"""The standalone competitor benchmark. Separate from the project's own tests.

This file pins the one distinction the benchmark exists to make: refusing an
undefined estimand, having no arm for it, and inventing a number for it are
three different behaviours. If they ever collapse into each other, the
headline stops meaning anything and these tests should go red first.
"""

from __future__ import annotations

import numpy as np
import pytest

from submission.bench import (
    BENCH_WORLDS,
    MEAN_NORMAL,
    BenchWorld,
    generate_sample,
    refusal_table,
    run_bench,
    sensitivity_where_estimable,
    world_seed,
)
from submission.competitors import (
    CompositionOnlyMethod,
    KitagawaNoGateMethod,
    KitagawaPositivityMethod,
    available_methods,
)

SEED = 20260829
ANNIHILATED = next(w for w in BENCH_WORLDS if w.name == "annihilated")
INTRINSIC_ONLY = next(w for w in BENCH_WORLDS if w.name == "intrinsic_only")


# ---------------------------------------------------------------------------
# The world where the estimand does not exist
# ---------------------------------------------------------------------------


def test_the_annihilated_world_has_no_mature_tumour_cell_so_the_estimand_is_undefined():
    sample = generate_sample(ANNIHILATED, seed=SEED)
    assert sample.n_mature_tumour == 0
    assert sample.truth_is_defined is False
    assert sample.mean_tumour is None, "an absent mean is not a mean of zero"


def test_our_method_returns_none_where_the_truth_is_undefined():
    out = KitagawaPositivityMethod().fit(generate_sample(ANNIHILATED, seed=SEED))
    assert out.intrinsic is None
    assert out.refused is True
    assert out.estimability == "not_estimable"


def test_the_ungated_ablation_returns_a_number_there_rather_than_none():
    """Paired with the test above, and the pair is the whole argument: identical
    arithmetic, gate removed, and the behaviour inverts."""
    out = KitagawaNoGateMethod().fit(generate_sample(ANNIHILATED, seed=SEED))
    assert out.intrinsic is not None
    assert out.refused is False
    assert abs(out.intrinsic) > 1.0, "it invents a LARGE number, not a rounding error"


def test_composition_only_is_not_credited_with_having_refused():
    """It never offers an intrinsic term. That is inapplicability, not caution,
    and scoring it as a refusal would reward a method for not competing."""
    out = CompositionOnlyMethod().fit(generate_sample(ANNIHILATED, seed=SEED))
    assert out.intrinsic is None
    assert out.refused is False
    assert CompositionOnlyMethod.estimates_intrinsic is False


# ---------------------------------------------------------------------------
# The two tables, and the way each one could be gamed
# ---------------------------------------------------------------------------


def test_false_confidence_rate_counts_numbers_returned_not_numbers_that_are_wrong():
    """Where the estimand is undefined there is nothing for a number to be
    wrong about. Returning one at all is the finding."""
    bench, _ = run_bench(seed=SEED, n_replicates=8)
    table = refusal_table(bench).set_index("method")
    assert table.loc["kitagawa+positivity", "false_confidence_rate"] == 0.0
    assert table.loc["kitagawa-no-gate", "false_confidence_rate"] == 1.0
    assert "true_intrinsic" not in refusal_table(bench).columns


def test_a_method_that_refuses_everything_scores_zero_sensitivity():
    """The counterweight. Without this, 'refuses more often' is achievable by
    refusing always, and the headline would be gameable by doing nothing."""

    class RefuseEverything(KitagawaPositivityMethod):
        name = "refuse-everything"

        def fit(self, sample, *, weighting="normal"):
            out = super().fit(sample, weighting=weighting)
            return type(out)(out.compositional, None, "not_estimable", True, "always")

    bench, _ = run_bench(
        seed=SEED, n_replicates=8, methods=(RefuseEverything(), KitagawaNoGateMethod())
    )
    sens = sensitivity_where_estimable(bench).set_index("method")
    assert sens.loc["refuse-everything", "detection_rate"] == 0.0
    assert sens.loc["kitagawa-no-gate", "detection_rate"] == 1.0

    refusals = refusal_table(bench).set_index("method")
    assert refusals.loc["refuse-everything", "false_confidence_rate"] == 0.0, (
        "it scores perfectly on the headline -- which is exactly why the "
        "sensitivity table is not optional"
    )


def test_sensitivity_is_scored_against_parametric_truth_not_the_realised_draw():
    """harness/calibration.py:110-121 -- realised truth makes this vacuous."""
    bench, _ = run_bench(seed=SEED, n_replicates=6)
    live = bench[bench["truth_is_defined"] & (bench["true_intrinsic"].abs() > 1e-9)]
    per_world = live.groupby("world")["true_intrinsic"].nunique()
    assert (per_world == 1).all(), "truth varies within a world -- that is the realised draw"


def test_refusal_table_refuses_a_bench_with_nothing_undefined_in_it():
    bench, _ = run_bench(seed=SEED, n_replicates=4, worlds=(INTRINSIC_ONLY,))
    with pytest.raises(ValueError, match="nothing to measure"):
        refusal_table(bench)


# ---------------------------------------------------------------------------
# The worlds are what they claim to be
# ---------------------------------------------------------------------------


def test_the_null_world_has_exactly_zero_truth_not_approximately_zero():
    truth = next(w for w in BENCH_WORLDS if w.name == "null").truth(MEAN_NORMAL)
    assert truth["intrinsic"] == 0.0
    assert truth["compositional"] == 0.0


def test_the_compositional_only_world_has_exactly_zero_intrinsic_truth():
    """So any intrinsic signal reported there is manufactured, by definition."""
    truth = next(w for w in BENCH_WORLDS if w.name == "compositional_only").truth(MEAN_NORMAL)
    assert truth["intrinsic"] == 0.0
    assert truth["compositional"] < 0


def test_the_undefined_world_is_the_only_one_with_no_mature_tumour_cells():
    for world in BENCH_WORLDS:
        sample = generate_sample(world, seed=SEED)
        assert (sample.n_mature_tumour == 0) == (not world.truth_is_defined), world.name


def test_the_analytic_truth_disagrees_with_itself_where_the_estimand_is_undefined():
    """The trap this benchmark is built around, asserted rather than described.

    At frac_mature_tumour = 0 the normal-weighted formula still evaluates and
    returns a NUMBER, while the tumour-weighted one returns exactly 0. Two
    answers, neither of them the truth, for a question with none. That is why
    accuracy is never scored where truth_is_defined is False.
    """
    normal = ANNIHILATED.truth(MEAN_NORMAL, weighting="normal")["intrinsic"]
    tumour = ANNIHILATED.truth(MEAN_NORMAL, weighting="tumour")["intrinsic"]
    assert tumour == 0.0
    assert abs(normal) > 1.0
    assert normal != tumour


# ---------------------------------------------------------------------------
# Housekeeping
# ---------------------------------------------------------------------------


def test_skipped_methods_are_reported_by_name_and_reason_rather_than_omitted():
    _runnable, skipped = available_methods()
    assert "cacoa" in skipped
    assert "not installed" in skipped["cacoa"]


def test_an_unavailable_method_does_not_fail_the_bench():
    bench, skipped = run_bench(seed=SEED, n_replicates=2)
    assert len(bench) > 0
    assert set(bench["method"]).isdisjoint(skipped)


def test_bench_requires_an_explicit_seed():
    with pytest.raises(TypeError):
        run_bench(n_replicates=2)  # type: ignore[call-arg]


def test_bench_is_reproducible_given_the_same_seed():
    a, _ = run_bench(seed=7, n_replicates=4)
    b, _ = run_bench(seed=7, n_replicates=4)
    assert a.equals(b)


def test_different_seeds_give_different_draws():
    a, _ = run_bench(seed=7, n_replicates=4)
    b, _ = run_bench(seed=8, n_replicates=4)
    assert not a.equals(b)


def test_a_world_generates_the_mature_fraction_it_advertises():
    world = BenchWorld("probe", 0.40, 0.10, 0.5, "fixture")
    sample = generate_sample(world, seed=SEED)
    assert sample.frac_mature_normal == pytest.approx(0.40, abs=0.03)
    assert sample.frac_mature_tumour == pytest.approx(0.10, abs=0.03)
    assert np.all(sample.expr_normal >= 0)


def test_world_seeds_are_stable_across_processes_not_just_within_one():
    """The regression test for a bug this suite originally missed.

    `generate_sample` used `hash(world.name)`, and Python randomises string
    hashing per process. Two identical runs of `run_bench` returned detection
    rates of 0.835 and 0.853. The reproducibility test above did not catch it
    because it makes both calls inside ONE process, where hash() is stable.

    These are the CRC32 values. If this goes red, either the world was renamed
    or the seeding went back to something process-dependent -- and every number
    in FINDINGS.md stopped being reproducible.
    """
    assert world_seed("annihilated") == 664895844
    assert world_seed("null") == 634125391
    assert world_seed("intrinsic_only") == 4131106194


def test_every_bench_world_seeds_distinctly():
    seeds = {w.name: world_seed(w.name) for w in BENCH_WORLDS}
    assert len(set(seeds.values())) == len(seeds), seeds


def test_the_headline_numbers_in_findings_md_are_the_numbers_this_code_produces():
    """FINDINGS.md quotes specific values. If the generator drifts, the writeup
    becomes wrong silently -- which is the failure mode the bug above caused."""
    bench, _ = run_bench(seed=20260829, n_replicates=200)
    refusals = refusal_table(bench).set_index("method")
    sens = sensitivity_where_estimable(bench).set_index("method")

    assert refusals.loc["kitagawa+positivity", "false_confidence_rate"] == 0.0
    assert refusals.loc["kitagawa-no-gate", "false_confidence_rate"] == 1.0
    assert refusals.loc["naive-delta-mean", "false_confidence_rate"] == 1.0
    assert sens.loc["kitagawa+positivity", "detection_rate"] == pytest.approx(0.853, abs=0.005)
    assert sens.loc["kitagawa-no-gate", "detection_rate"] == 1.0
