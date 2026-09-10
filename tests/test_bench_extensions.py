"""Every guard the three benchmark extensions add, against the input that
forces it to fail.

Same rule as ``tests/test_checks_can_fail.py``, applied to new code: **a check
unable to fail is worse than no check, because it turns an absence of evidence
into a green light.** Each guard here gets a constructed input it MUST catch,
plus a positive control on real output so that "it always raises" is not how it
passes.

Three of these are not guards over data but guards over an ARGUMENT:

* ``test_annihilated_theorem_*`` mechanises the proof in
  ``annihilated_invariance_proof``. The proof turns on ``x < 0.0`` being False
  for every ``x`` ``Generator.random`` can return, and on the width gate
  checking ``n_t < 1`` *before* it computes a moment. Either could be edited
  away while every other test still passed.
* ``test_poisson_path_is_bit_identical`` pins the refactor. The count-model
  seam had to leave the committed benchmark untouched, and "it looked fine" is
  not a measurement.
* ``test_every_model_preserves_the_marginal_mean`` pins the one design decision
  in ``expression_models``: if a misspecified arm moved the mean, the truth
  column would be wrong and the accuracy numbers would be measuring the
  generator's error while reading as the estimator's bias.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.reference.interval_calibration import expected_false_positive_rate
from submission.bench import BENCH_WORLDS, MEAN_NORMAL, generate_sample, run_bench
from submission.competitors import (
    KitagawaNoGateMethod,
    KitagawaPositivityMethod,
    VarianceGateKitagawaMethod,
)
from submission.expression_models import (
    EXPRESSION_MODELS,
    POISSON,
    ExpressionModel,
    ExpressionModelError,
)
from submission.extensions import (
    CLOSED_FORM_AGREEMENT_TOLERANCE,
    COMMITTED_DETECTABLE_SHIFT,
    DETECTABLE_SHIFT_GRID,
    ClosedFormDivergence,
    InvarianceViolation,
    SweepDoesNotBracketError,
    check_annihilated_refusal_is_model_invariant,
    check_closed_form_tracks_simulation,
    check_committed_value_is_in_grid,
    check_seeds_separate_model_from_stream,
    check_step_matches_grid,
    check_sweep_brackets_the_step,
    critical_detectable_shift,
    detectable_shift_sweep,
    interval_width_curve,
    misspecification_sweep,
    model_effect_vs_seed_noise,
    variance_gate_at,
)
from submission.run_bench import SEED, SIMULATION_REFERENCE

#: Small enough to keep the file fast, large enough that the step in
#: ``depleted_wide`` is still bracketed by the committed grid.
N_REPS = 40

ANNIHILATED = next(w for w in BENCH_WORLDS if w.frac_mature_tumour == 0.0)
DEPLETED_WIDE = next(w for w in BENCH_WORLDS if w.name == "depleted_wide")


@pytest.fixture(scope="module")
def sweep_and_steps():
    return detectable_shift_sweep(seed=SEED, n_replicates=N_REPS)


@pytest.fixture(scope="module")
def misspec():
    return misspecification_sweep(n_replicates=N_REPS)


@pytest.fixture(scope="module")
def misspec_full():
    """The sweep at the COMMITTED replicate count.

    The fast fixture above is enough to exercise the guards, which are
    structural. It is NOT enough for the two attribution tests: separating a
    model effect from seed noise is a comparison of two spreads, and at 40
    replicates both are small integers and the ratio is dominated by rounding.
    A substantive finding has to be pinned at the count it is reported at.
    """
    return misspecification_sweep(n_replicates=200)


@pytest.fixture(scope="module")
def committed_simulation():
    from src.common.paths import RESULTS_DIR

    path = RESULTS_DIR / SIMULATION_REFERENCE
    if not path.exists():
        pytest.skip(f"{path} absent — the cross-check reference is not in this tree")
    return pd.read_parquet(path)


# ---------------------------------------------------------------------------
# 1. The detectable-effect sweep
# ---------------------------------------------------------------------------


def test_sweep_bracket_guard_fails_on_a_grid_that_stops_short():
    """THE FAILING INPUT. A grid entirely below the transition.

    This is the defect the guard exists for: at the committed s_detect the
    width gate refuses 0/200 in every estimable world, so a sweep run only
    near 0.5 reports "the width criterion never binds" — which is true of the
    grid and says nothing about the criterion. Without this guard that table
    is publishable and looks like a negative result.
    """
    short_grid = (0.50, 0.60, 0.70)
    sweep, _ = detectable_shift_sweep(
        seed=SEED, n_replicates=N_REPS, grid=short_grid
    )
    # Truth-defined worlds only. `annihilated` refuses at rate 1.0 at every
    # s_detect via the finiteness guard, which is precisely why the guard under
    # test excludes it — a max taken over all worlds would be 1.0 here and
    # would say nothing about whether the grid reached the width threshold.
    live = sweep[sweep["truth_is_defined"]]
    assert live["width_gate_refusal_rate"].max() == 0.0, (
        "the constructed input is meant to contain NO refusals in any "
        "truth-defined world; if the gate now binds at s_detect <= 0.70 this "
        "test is no longer exercising the guard"
    )
    with pytest.raises(SweepDoesNotBracketError, match="never-binding"):
        check_sweep_brackets_the_step(sweep)


def test_sweep_bracket_guard_fails_when_only_annihilated_is_present():
    """A sweep of the undefined world alone brackets nothing.

    ``annihilated`` refuses at every s_detect via the finiteness guard, so its
    column is saturated at 200 regardless of the parameter. A guard reading
    only "did we see both 0 and 1" over all worlds would pass on that and
    certify a sweep that never varied anything.
    """
    sweep, _ = detectable_shift_sweep(
        seed=SEED, n_replicates=N_REPS, worlds=(ANNIHILATED,)
    )
    with pytest.raises(SweepDoesNotBracketError, match="no truth-defined world"):
        check_sweep_brackets_the_step(sweep)


def test_sweep_bracket_guard_passes_on_the_committed_grid(sweep_and_steps):
    """POSITIVE CONTROL."""
    sweep, _ = sweep_and_steps
    check_sweep_brackets_the_step(sweep)


def test_committed_value_guard_fails_when_the_shipped_value_is_absent():
    """THE FAILING INPUT. A grid that omits the value the method ships with.

    "The width gate does the job for free" is a claim about DETECTABLE_SHIFT =
    0.5 specifically. A sweep over 0.90..0.99 cannot support or refute it, and
    would be read as if it could.
    """
    sweep, _ = detectable_shift_sweep(
        seed=SEED, n_replicates=N_REPS, grid=(0.90, 0.95, 0.98)
    )
    with pytest.raises(SweepDoesNotBracketError, match="not a point on this grid"):
        check_committed_value_is_in_grid(sweep)


def test_committed_value_guard_passes_on_the_committed_grid(sweep_and_steps):
    """POSITIVE CONTROL."""
    sweep, _ = sweep_and_steps
    check_committed_value_is_in_grid(sweep)
    assert COMMITTED_DETECTABLE_SHIFT in DETECTABLE_SHIFT_GRID


def test_step_matches_grid_guard_fails_on_a_tampered_sweep(sweep_and_steps):
    """THE FAILING INPUT. One refusal count moved by one.

    The step table is derived from a closed form for the critical s_detect, not
    from the grid. That is exact and fast and it is also a second
    implementation of the gate's arithmetic — which is the standard way a
    benchmark comes to agree with itself. This guard is the only thing that
    would notice the rearrangement drifting from the shipped ``fit``.
    """
    sweep, steps = sweep_and_steps
    tampered = sweep.copy()
    target = tampered.index[
        (tampered["world"] == "depleted_wide") & (tampered["s_detect"] > 0.93)
    ][0]
    tampered.loc[target, "width_gate_refusals"] += 1
    with pytest.raises(SweepDoesNotBracketError, match="does not reproduce"):
        check_step_matches_grid(
            tampered, steps, seed=SEED, n_replicates=N_REPS
        )


def test_step_matches_grid_guard_passes_on_real_output(sweep_and_steps):
    """POSITIVE CONTROL: the closed form predicts every cell the gate produced."""
    sweep, steps = sweep_and_steps
    check_step_matches_grid(sweep, steps, seed=SEED, n_replicates=N_REPS)


def test_critical_shift_is_minus_infinity_where_no_mature_cell_survives():
    """s* has no finite value when the finiteness guard is what fires.

    Returning a finite number here — or 0.0, the invariant-1 mistake in this
    module's own coinage — would put ``annihilated`` on the step curve and let
    it be averaged into a statement about where the width threshold binds. It
    does not bind there; a different branch fires.
    """
    sample = generate_sample(ANNIHILATED, seed=SEED, replicate=0)
    assert critical_detectable_shift(sample) == float("-inf")


def test_variance_gate_subclass_does_not_leak_into_the_committed_class():
    """Sweeping a class attribute must not mutate the class being swept.

    Setting ``VarianceGateKitagawaMethod.DETECTABLE_SHIFT`` in a loop is the
    obvious way to write this sweep and it would leave the committed default
    changed for everything later in the process — including the misspecification
    run, if both are invoked by ``run_bench all``.
    """
    before = VarianceGateKitagawaMethod.DETECTABLE_SHIFT
    gate = variance_gate_at(0.97)
    assert gate.DETECTABLE_SHIFT == 0.97
    assert VarianceGateKitagawaMethod.DETECTABLE_SHIFT == before
    assert isinstance(gate, VarianceGateKitagawaMethod)
    # and the inherited fit is the committed one, not an override
    assert type(gate).fit is VarianceGateKitagawaMethod.fit


def test_the_step_is_narrow_relative_to_the_parameter_range(sweep_and_steps):
    """The finding itself, pinned: the live range is a small slice of [0, 1).

    Stated as a property rather than as the exact numbers, which move with
    n_replicates. What must not change quietly is that the committed value sits
    entirely below the range where the gate does anything.
    """
    _, steps = sweep_and_steps
    wide = steps[steps["world"] == "depleted_wide"].iloc[0]
    assert wide["s_star_min"] > COMMITTED_DETECTABLE_SHIFT + 0.3, (
        "the committed detectable effect is supposed to sit far below the "
        "range where the width gate binds at all"
    )
    assert wide["step_width"] < 0.15
    assert 0.0 < wide["headroom_below_collapse"] < 0.02, (
        "the band separating 'refuses in the wide regime' from 'refuses where "
        "the estimand plainly exists' is supposed to be very narrow — that is "
        "the result"
    )


# ---------------------------------------------------------------------------
# 2. Misspecification, and the theorem underneath it
# ---------------------------------------------------------------------------


def test_annihilated_theorem_no_model_ever_produces_a_mature_tumour_cell():
    """Step 1 of the proof, mechanised over every model and replicate.

    ``rng.random`` returns values in [0, 1), so ``x < 0.0`` is False for all of
    them. If the generator's comparison ever became ``<=``, or the world's
    fraction became a small positive number, the whole invariance argument
    would collapse and nothing else in the suite would notice.
    """
    for model in EXPRESSION_MODELS:
        for replicate in range(25):
            sample = generate_sample(
                ANNIHILATED, seed=SEED, replicate=replicate, model=model
            )
            assert sample.n_mature_tumour == 0
            # Step 2: a draw of size 0 leaves the arm exactly zero.
            assert not sample.expr_tumour.any()
            # And the absent mean is None, never 0.0 (invariant 1).
            assert sample.mean_tumour is None


def test_annihilated_theorem_every_gate_decides_without_reading_the_model():
    """Step 3 of the proof: the three gates' verdicts, per model, identical.

    Includes both zero-inflated arms, which is where the claim was least
    obvious: zero inflation adds mass at zero, and a gate reading the tumour
    arm's zeros could plausibly have behaved differently. None of them reads it.
    """
    count_gate, width_gate, no_gate = (
        KitagawaPositivityMethod(),
        VarianceGateKitagawaMethod(),
        KitagawaNoGateMethod(),
    )
    for model in EXPRESSION_MODELS:
        sample = generate_sample(ANNIHILATED, seed=SEED, replicate=3, model=model)
        assert count_gate.fit(sample).refused is True
        assert count_gate.fit(sample).intrinsic is None

        width = width_gate.fit(sample)
        assert width.refused is True
        assert width.intrinsic is None
        # It must be the FINITENESS branch, not the threshold branch — the
        # threshold is the only place a variance (and so the model) enters.
        assert "standard error of the per-cell mean is undefined" in (
            width.refusal_reason or ""
        )

        ungated = no_gate.fit(sample)
        assert ungated.refused is False
        assert ungated.intrinsic is not None


def test_invariance_guard_fails_on_a_tampered_sweep(misspec):
    """THE FAILING INPUT. One model's annihilated refusal count moved.

    This is what a real regression would look like: a width gate that computed
    its moments before checking finiteness would return NaN half-widths, and
    ``nan > threshold`` is False, so it would silently stop refusing under some
    models and not others.
    """
    tampered = misspec.copy()
    target = tampered.index[
        (tampered["world"] == "annihilated")
        & (tampered["gate"] == "width_gate")
        & (tampered["expression_model"] == "nb_disp0.5")
    ][0]
    tampered.loc[target, "n_refused"] = 0
    with pytest.raises(InvarianceViolation, match="moved with the count model"):
        check_annihilated_refusal_is_model_invariant(tampered)


def test_invariance_guard_fails_on_a_single_model_sweep():
    """A one-model sweep cannot falsify an invariance claim.

    The guard would pass trivially, which is the failure mode this whole file
    is about: it would certify invariance having compared nothing.
    """
    single = misspecification_sweep(
        n_replicates=5, models=(POISSON,), worlds=(ANNIHILATED,)
    )
    with pytest.raises(InvarianceViolation, match="at least\n?\\s*two count models|two count models"):
        check_annihilated_refusal_is_model_invariant(single)


def test_invariance_guard_fails_when_the_undefined_world_is_absent():
    """A sweep with no frac_mature_tumour == 0 world has nothing to check."""
    without = misspecification_sweep(
        n_replicates=5, worlds=(DEPLETED_WIDE,), models=EXPRESSION_MODELS[:2]
    )
    with pytest.raises(InvarianceViolation, match="nothing here"):
        check_annihilated_refusal_is_model_invariant(without)


def test_invariance_guard_passes_on_the_real_sweep(misspec):
    """POSITIVE CONTROL. The theorem's corollary, measured."""
    check_annihilated_refusal_is_model_invariant(misspec)
    block = misspec[misspec["world"] == "annihilated"]
    assert set(block.loc[block["gate"] == "count_gate", "refusal_rate"]) == {1.0}
    assert set(block.loc[block["gate"] == "width_gate", "refusal_rate"]) == {1.0}
    assert set(block.loc[block["gate"] == "no_gate", "refusal_rate"]) == {0.0}


def test_seed_guard_fails_on_a_single_seed_sweep():
    """THE FAILING INPUT. One seed, where model and stream are confounded.

    ``misspecification_sweep(seed=...)`` deliberately REPLACES the seed list
    rather than joining it, so a caller who asks for one seed gets one seed and
    trips this — rather than quietly receiving three and believing it measured
    one.
    """
    one = misspecification_sweep(
        seed=SEED, n_replicates=5, worlds=(DEPLETED_WIDE,)
    )
    assert one["seed"].nunique() == 1
    with pytest.raises(InvarianceViolation, match="perfectly confounded"):
        check_seeds_separate_model_from_stream(one)


def test_seed_guard_passes_on_the_real_sweep(misspec):
    """POSITIVE CONTROL."""
    check_seeds_separate_model_from_stream(misspec)
    assert misspec["seed"].nunique() >= 2


def test_count_gate_variation_is_not_attributed_to_the_generator(misspec_full):
    """The reading the separation table exists to prevent.

    The count gate's refusal count in ``depleted_wide`` moves with the model —
    because a different model consumes a different number of variates in the
    normal arm and so shifts the maturity draw, not because the count gate reads
    expression. It cannot: ``classify_estimability`` takes an integer.

    ``mean_n_mature_tumour`` is the corroboration: the distribution being drawn
    from, Binomial(2000, 0.01), has expectation 20 under every model.
    """
    separation = model_effect_vs_seed_noise(misspec_full)
    row = separation[
        (separation["world"] == "depleted_wide")
        & (separation["gate"] == "count_gate")
    ].iloc[0]
    assert row["attribution"] in ("SEED_CONFOUNDED", "AMBIGUOUS"), (
        "the count gate reads an integer, not an expression value; a "
        "GENERATOR attribution here would mean the sweep is measuring the RNG "
        f"stream and calling it misspecification (got {row['attribution']})"
    )
    means = misspec_full.loc[
        misspec_full["world"] == "depleted_wide", "mean_n_mature_tumour"
    ]
    assert means.between(18.0, 22.0).all()


def test_width_gate_overdispersion_effect_is_attributed_to_the_generator(misspec_full):
    """The other half: a real model effect, separated from the stream.

    Under Poisson the width gate refuses 0 times in ``depleted_wide``. Under
    NB(r=0.5) — variance 41x Poisson at mu=20 — it begins to bind. That is the
    part of ``FINDINGS.md`` limitation 6 that turns out to be true: the width
    gate's non-binding is NOT invariant to the expression model, while the
    finiteness guard's firing is.
    """
    block = misspec_full[
        (misspec_full["world"] == "depleted_wide") & (misspec_full["gate"] == "width_gate")
    ]
    poisson = block[block["expression_model"] == "poisson"]["n_refused"]
    severe = block[block["expression_model"] == "nb_disp0.5"]["n_refused"]
    assert (poisson == 0).all(), "Poisson is the published baseline: 0 refusals"
    assert (severe > 0).all(), (
        "at dispersion 0.5 the width gate is expected to bind at every seed"
    )
    separation = model_effect_vs_seed_noise(misspec_full)
    row = separation[
        (separation["world"] == "depleted_wide")
        & (separation["gate"] == "width_gate")
    ].iloc[0]
    assert row["attribution"] == "GENERATOR"


# ---------------------------------------------------------------------------
# The count models themselves
# ---------------------------------------------------------------------------


def test_poisson_path_is_bit_identical_to_the_committed_benchmark():
    """The refactor must not have moved a single number.

    Compared against the committed ``submission/results/bench_raw.parquet`` on
    every column that existed before the count-model seam was added. A refactor
    that "looks inert" and a refactor that is inert are different claims and
    only one of them is checkable.
    """
    from pathlib import Path

    committed_path = (
        Path(__file__).resolve().parents[1]
        / "submission"
        / "results"
        / "bench_raw.parquet"
    )
    if not committed_path.exists():
        pytest.skip("committed bench_raw.parquet absent")
    committed = pd.read_parquet(committed_path)
    fresh, _ = run_bench(seed=SEED, n_replicates=200)
    shared = [c for c in committed.columns if c in fresh.columns]
    assert len(shared) == len(committed.columns)
    pd.testing.assert_frame_equal(
        committed[shared].reset_index(drop=True),
        fresh[shared].reset_index(drop=True),
    )


def test_every_model_preserves_the_marginal_mean():
    """THE LOAD-BEARING DESIGN DECISION in ``expression_models``.

    ``BenchWorld.truth`` is a statement about means. A zero-inflated arm drawn
    the naive way has mean ``(1 - pi) * mu``, so the truth column would be
    wrong by 30% at pi = 0.3 — and every accuracy number computed against it
    would be measuring the generator's own error while reading as the
    estimator's bias under misspecification. That is a wrong finding, not a
    noisy one.
    """
    rng = np.random.default_rng(0)
    for model in EXPRESSION_MODELS:
        draws = model.draw(rng, MEAN_NORMAL, 200_000)
        assert abs(draws.mean() - MEAN_NORMAL) < 0.25, model.name
        # and the sampler matches the variance it claims in closed form
        claimed = model.theoretical_variance(MEAN_NORMAL)
        assert abs(draws.var() / claimed - 1.0) < 0.05, model.name


def test_overdispersion_is_ordered_and_poisson_is_the_floor():
    """The arms must actually be more dispersed, in the stated order."""
    variances = [m.theoretical_variance(MEAN_NORMAL) for m in EXPRESSION_MODELS]
    assert variances[0] == pytest.approx(MEAN_NORMAL)
    assert min(variances) == variances[0], "Poisson must be the least dispersed arm"


def test_a_draw_of_size_zero_consumes_nothing_from_the_stream():
    """Step 2 of the theorem, at the level of the sampler.

    If a model consumed variates on an empty draw, the tumour arm's absence
    would perturb everything after it and the Poisson path would not be
    bit-identical.
    """
    for model in EXPRESSION_MODELS:
        a = np.random.default_rng(4)
        assert model.draw(a, 20.0, 0).size == 0
        b = np.random.default_rng(4)
        assert a.random() == b.random()


@pytest.mark.parametrize(
    "kwargs, match",
    [
        ({"kind": "nb"}, "needs a positive dispersion"),
        ({"kind": "nb", "dispersion": 0.0}, "needs a positive dispersion"),
        ({"kind": "poisson", "dispersion": 2.0}, "no dispersion parameter"),
        ({"kind": "zinb", "dispersion": 2.0, "zero_inflation": 1.0}, "must be in"),
        ({"kind": "zip", "zero_inflation": 0.0}, "disagree about whether"),
        ({"kind": "nb", "dispersion": 2.0, "zero_inflation": 0.3}, "disagree about"),
        ({"kind": "lognormal"}, "unknown kind"),
    ],
)
def test_expression_model_rejects_incoherent_parameters(kwargs, match):
    """THE FAILING INPUTS for the model's own validation.

    The one that matters most is ``poisson`` with a dispersion: silently
    ignoring it would let a sweep report an overdispersion arm that never
    overdispersed, and the table would show invariance because nothing varied.
    """
    with pytest.raises(ExpressionModelError, match=match):
        ExpressionModel("bad", **kwargs)


# ---------------------------------------------------------------------------
# 3. The interval-width curve
# ---------------------------------------------------------------------------


def test_curve_refuses_to_start_below_n_equals_two():
    """THE FAILING INPUT. n=1, where neither interval exists.

    A curve that started at 1 would have a fabricated first point, and it is
    the small-n end that the whole reading is about.
    """
    with pytest.raises(ClosedFormDivergence, match="undefined below n=2"):
        interval_width_curve(n_min=1, n_max=10)


def test_closed_form_guard_fails_on_a_drifted_simulation(committed_simulation):
    """THE FAILING INPUT. The simulated rates shifted by 10 points.

    The closed form is smooth, total and free. That is exactly what makes it
    dangerous: if it stopped describing the estimator the project actually
    runs, the curve would still look like a result. This is the only thing
    standing between the two.
    """
    drifted = committed_simulation.copy()
    drifted.loc[drifted["method"] == "percentile", "false_positive_rate"] += 0.10
    curve = interval_width_curve(n_min=2, n_max=50, simulated=drifted)
    with pytest.raises(ClosedFormDivergence, match="differ by more than"):
        check_closed_form_tracks_simulation(curve)


def test_closed_form_guard_fails_when_there_is_nothing_to_check_against():
    """An unchecked closed form is the object this section is about."""
    curve = interval_width_curve(n_min=2, n_max=50, simulated=None)
    with pytest.raises(ClosedFormDivergence, match="no n in this curve"):
        check_closed_form_tracks_simulation(curve)


def test_closed_form_guard_rejects_a_non_percentile_reference(committed_simulation):
    """Cross-checking against the wrong estimator is not a cross-check.

    ``student_t`` is well calibrated everywhere, so comparing the percentile
    bootstrap's closed form against it would show a large disagreement — or,
    worse, be silently averaged with the percentile rows and dilute one.
    """
    wrong = committed_simulation[committed_simulation["method"] == "student_t"]
    with pytest.raises(ClosedFormDivergence, match="no percentile rows"):
        interval_width_curve(n_min=2, n_max=50, simulated=wrong)


def test_closed_form_guard_passes_on_the_committed_simulation(committed_simulation):
    """POSITIVE CONTROL, and the reported agreement."""
    curve = interval_width_curve(simulated=committed_simulation)
    check_closed_form_tracks_simulation(curve)
    checked = curve[curve["n_simulation_cells"] > 0]
    assert len(checked) == 5, "the committed simulation covers five n values"
    assert checked["abs_diff_vs_median"].max() < CLOSED_FORM_AGREEMENT_TOLERANCE
    # The claim made in the report: agreement within about one point.
    assert checked["abs_diff_vs_median"].max() < 0.01


def test_the_curve_contains_no_data():
    """The point of the whole section, as an executable statement.

    Two calls with no arguments in common except n produce identical curves.
    There is nothing to pass that could change them — no counts matrix, no
    gene, no cohort, no seed — which is why a percentile bootstrap cannot
    report its own miscalibration: the miscalibration is not in the data.
    """
    a = interval_width_curve(n_min=2, n_max=60)
    b = interval_width_curve(n_min=2, n_max=60)
    pd.testing.assert_frame_equal(a, b)
    assert a.loc[a["n_patients"] == 2, "closed_form_false_positive_rate"].iloc[
        0
    ] == pytest.approx(expected_false_positive_rate(2))


def test_the_curve_is_monotone_and_hits_the_reported_landmarks():
    """The published numbers, pinned. Closed form, so exact."""
    curve = interval_width_curve(n_min=2, n_max=200).set_index("n_patients")
    assert curve["width_ratio_vs_t"].is_monotonic_increasing
    assert curve["closed_form_false_positive_rate"].is_monotonic_decreasing
    for n, width, fpr in [
        (2, 0.109, 0.398),
        (4, 0.533, 0.188),
        (10, 0.822, 0.096),
        (20, 0.913, 0.071),
        (44, 0.961, 0.059),
        (100, 0.983, 0.054),
    ]:
        assert curve.loc[n, "width_ratio_vs_t"] == pytest.approx(width, abs=1e-3)
        assert curve.loc[n, "closed_form_false_positive_rate"] == pytest.approx(
            fpr, abs=1e-3
        )
    # Every n below 22 is MISCALIBRATED by the project's own tolerance.
    assert (curve.loc[:21, "verdict"] == "MISCALIBRATED").all()
    assert curve.loc[22, "verdict"] == "CALIBRATED"
