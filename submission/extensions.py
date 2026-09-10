"""Three extensions, each aimed at a limitation the benchmark already admits.

``FINDINGS.md`` closes with six limitations. Three of them are not really
limitations, they are unmeasured quantities, and this module measures them.

1. **The width gate has no calibrated parameter.** Limitation 6 says the width
   criterion fails to bind "partly because the detectable effect the design
   names is large relative to the sampling error". That is a sweepable claim:
   ``VarianceGateKitagawaMethod.DETECTABLE_SHIFT`` is a class attribute and
   nobody calibrated it. :func:`detectable_shift_sweep` sweeps it. The finding
   is that the width gate is a **step function** in that parameter with no
   intermediate regime -- the paper's own thesis, applied to the paper's own
   preferred competitor.

2. **The generator is exactly the model the estimator assumes.** Limitation 1.
   :func:`misspecification_sweep` runs the committed method classes, unchanged,
   against negative-binomial and zero-inflated arms. The finiteness-guard
   result turns out to be invariant, and **for a reason that is provable rather
   than measured** -- see :func:`annihilated_invariance_proof`.

3. **The interval-width curve is evaluated at five points.** ``width_ratio``
   and ``expected_false_positive_rate`` are closed-form functions of ``n``
   alone. The committed job evaluates them only at the n values the real
   cohorts happen to have. :func:`interval_width_curve` evaluates the whole
   thing, which makes visible that the miscalibration of a percentile bootstrap
   is **a pure function of sample size with no data in it** -- so, exactly like
   a recovery curve, it cannot report its own failure.

EVERY GUARD HERE HAS A COMMITTED FAILING INPUT in
``tests/test_bench_extensions.py``. A guard with no such input is untested and
the repository's own rule (``tests/test_checks_can_fail.py``) treats it as
absent.

Everything in this module is synthetic or closed-form. Nothing here is a
result about colorectal cancer.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from src.reference.interval_calibration import (
    NOMINAL_ALPHA,
    calibration_verdict,
    expected_false_positive_rate,
    width_ratio,
)
from submission.bench import BENCH_WORLDS, MEAN_NORMAL, BenchWorld, generate_sample
from submission.competitors import (
    KitagawaNoGateMethod,
    KitagawaPositivityMethod,
    VarianceGateKitagawaMethod,
)
from submission.expression_models import EXPRESSION_MODELS, ExpressionModel


class SweepDoesNotBracketError(ValueError):
    """A sweep that does not contain the transition it claims to locate."""


class InvarianceViolation(ValueError):
    """A quantity asserted to be model-invariant that is not."""


class ClosedFormDivergence(ValueError):
    """A closed form and the simulation it describes have come apart."""


# ---------------------------------------------------------------------------
# 1. The detectable-effect sweep
# ---------------------------------------------------------------------------

#: The value ``competitors.VarianceGateKitagawaMethod`` ships with, taken from
#: the design document's ``s = 0.5``. Read off the class rather than retyped,
#: so the sweep cannot silently stop being anchored to the committed method.
COMMITTED_DETECTABLE_SHIFT = VarianceGateKitagawaMethod.DETECTABLE_SHIFT

#: The grid. Coarse below 0.85 because nothing happens there, then 0.005 steps
#: to 0.99 because the pilot put the transition between 0.90 and 0.95 and a
#: 0.05 grid cannot tell a step from a ramp. Committed as a constant so the
#: table is reproducible from the module rather than from a command line.
DETECTABLE_SHIFT_GRID: tuple[float, ...] = (
    0.50, 0.60, 0.70, 0.80,
) + tuple(round(float(s), 3) for s in np.arange(0.850, 0.9901, 0.005))


def variance_gate_at(detectable_shift: float) -> VarianceGateKitagawaMethod:
    """The committed width-gate method with ``DETECTABLE_SHIFT`` overridden.

    A subclass rather than a monkeypatch. ``fit`` is inherited untouched, so
    every row of the sweep runs the same code path the committed method runs;
    and because the attribute is set on a fresh type, no sweep can leak a value
    into ``VarianceGateKitagawaMethod`` itself and corrupt a later run in the
    same process. That leak is not hypothetical -- it is what setting a class
    attribute in a loop does.
    """
    cls = type(
        f"VarianceGateAt{detectable_shift!r}",
        (VarianceGateKitagawaMethod,),
        {
            "DETECTABLE_SHIFT": float(detectable_shift),
            "name": f"kitagawa+variance-gate@s_detect={detectable_shift:g}",
        },
    )
    return cls()


def critical_detectable_shift(sample) -> float:
    """The exact ``s_detect`` at which the width gate flips, for one replicate.

    The gate refuses when ``half_width > |f_n (s_detect - 1) mean_n|``. For
    ``s_detect < 1`` that rearranges to a threshold on ``s_detect`` alone:

        refuse  <=>  s_detect  >  1 - half_width / (f_n * mean_n)  =:  s*

    So the sweep does not have to find the step by bisection -- there is a
    closed form for it per replicate, and the grid is only there to present it.
    ``-inf`` when no mature tumour cell exists: the finiteness guard refuses at
    every ``s_detect``, including ones that would otherwise never bind.

    The moments come from ``VarianceGateKitagawaMethod._mature_moments``, the
    committed method's own staticmethod, so this is a rearrangement of the
    shipped arithmetic and not a second implementation of it. A second
    implementation would be free to agree with itself.
    """
    n_t = int(sample.n_mature_tumour)
    n_n = int(round(sample.frac_mature_normal * sample.expr_normal.size))
    if n_t < 1 or n_n < 1:
        return float("-inf")
    moments = VarianceGateKitagawaMethod._mature_moments
    _, var_t = moments(sample.expr_tumour, n_t)
    _, var_n = moments(sample.expr_normal, n_n)
    se = sample.frac_mature_normal * float(np.sqrt(var_t / n_t + var_n / n_n))
    half_width = VarianceGateKitagawaMethod.Z * se
    denom = sample.frac_mature_normal * sample.mean_normal
    return 1.0 - half_width / denom


def detectable_shift_sweep(
    *,
    seed: int,
    n_replicates: int = 200,
    grid: tuple[float, ...] = DETECTABLE_SHIFT_GRID,
    worlds: tuple[BenchWorld, ...] = BENCH_WORLDS,
    model: ExpressionModel | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Refusals per (``s_detect``, world), plus the exact per-world step.

    Returns ``(sweep, steps)``. Each replicate's cells are generated ONCE and
    shown to every point on the grid, so a difference between two rows is the
    gate and nothing else. Generating fresh cells per grid point would leave
    the whole sweep confounded with sampling noise at the exact place the step
    is being located.
    """
    from submission.expression_models import POISSON

    model = model or POISSON
    count_gate = KitagawaPositivityMethod()
    gates = {s: variance_gate_at(s) for s in grid}

    sweep_rows, step_rows = [], []
    for world in worlds:
        refusals = {s: 0 for s in grid}
        count_refusals = 0
        criticals = []
        for replicate in range(n_replicates):
            sample = generate_sample(
                world, seed=seed, replicate=replicate, model=model
            )
            criticals.append(critical_detectable_shift(sample))
            count_refusals += int(count_gate.fit(sample).refused)
            for s, gate in gates.items():
                refusals[s] += int(gate.fit(sample).refused)

        crit = np.asarray(criticals, dtype=float)
        finite = crit[np.isfinite(crit)]
        for s in grid:
            sweep_rows.append(
                {
                    "world": world.name,
                    "truth_is_defined": world.truth_is_defined,
                    "s_detect": float(s),
                    "is_committed_value": bool(
                        np.isclose(s, COMMITTED_DETECTABLE_SHIFT)
                    ),
                    "n_replicates": int(n_replicates),
                    "width_gate_refusals": int(refusals[s]),
                    "width_gate_refusal_rate": refusals[s] / n_replicates,
                    "count_gate_refusals": int(count_refusals),
                    "expression_model": model.name,
                }
            )
        ordered = np.sort(finite)
        # The s_detect that reproduces the COUNT gate's refusal count in this
        # world, if one exists. Refusals are a step function of s_detect with
        # jumps at the order statistics of s*, so "the width gate refuses k
        # times" holds exactly on the half-open interval between the k-th and
        # (k+1)-th of them. Reported as a window, because a single midpoint
        # would hide how narrow it is -- which is the finding.
        k = int(count_refusals)
        if 0 < k < ordered.size:
            match_lo, match_hi = float(ordered[k - 1]), float(ordered[k])
        else:
            match_lo = match_hi = float("nan")
        # How steep the step is, in refusals per unit s_detect, across the
        # middle half of it. Defined for every world rather than only where a
        # match exists, because the steepness IS the finding: it says how
        # finely an uncalibrated parameter would have to be tuned to land on
        # any particular refusal count.
        if ordered.size:
            p25, p75 = np.percentile(ordered, [25, 75])
            iqr = float(p75 - p25)
            slope = (0.5 * ordered.size / iqr) if iqr > 0 else float("nan")
        else:
            slope = float("nan")
        step_rows.append(
            {
                "world": world.name,
                "truth_is_defined": world.truth_is_defined,
                "n_replicates": int(n_replicates),
                # -inf criticals are replicates the finiteness guard catches at
                # every s_detect. Counted, never averaged into the step.
                "n_refused_at_every_s_detect": int((~np.isfinite(crit)).sum()),
                "s_star_min": float(ordered.min()) if ordered.size else float("nan"),
                "s_star_median": (
                    float(np.median(ordered)) if ordered.size else float("nan")
                ),
                "s_star_max": float(ordered.max()) if ordered.size else float("nan"),
                "step_width": (
                    float(ordered.max() - ordered.min())
                    if ordered.size
                    else float("nan")
                ),
                "count_gate_refusals": k,
                "s_detect_matching_count_gate_lo": match_lo,
                "s_detect_matching_count_gate_hi": match_hi,
                "matching_window_width": match_hi - match_lo,
                "refusals_per_unit_s_detect": slope,
                "expression_model": model.name,
            }
        )
    steps = pd.DataFrame(step_rows)

    # The number the whole reading turns on, and it is a statement ACROSS
    # worlds rather than within one: the smallest s_detect at which the gate
    # starts refusing in some world where the estimand plainly exists. Carried
    # as a column on every row, repeated, because a reader checking "is there a
    # setting that refuses in depleted_wide but nowhere it should not" needs it
    # next to s_star_max and must not have to recompute it from the raw sweep.
    defines_onset = steps["truth_is_defined"] & (steps["world"] != "depleted_wide")
    onset = (
        float(steps.loc[defines_onset, "s_star_min"].min())
        if bool(defines_onset.any())
        else float("nan")
    )
    steps["collapse_onset_s_detect"] = onset
    # Room between this world's step finishing and the gate starting to refuse
    # somewhere refusing is plainly a false negative. NaN on the four worlds
    # that DEFINE the onset -- there the quantity is self-referential, and a
    # negative number there would read as a finding rather than as arithmetic
    # comparing a set to its own minimum.
    steps["headroom_below_collapse"] = np.where(
        defines_onset, np.nan, onset - steps["s_star_max"]
    )
    return pd.DataFrame(sweep_rows), steps


def check_step_matches_grid(
    sweep: pd.DataFrame, steps: pd.DataFrame, *, seed: int, n_replicates: int,
    grid: tuple[float, ...] = DETECTABLE_SHIFT_GRID,
    worlds: tuple[BenchWorld, ...] = BENCH_WORLDS,
    model: ExpressionModel | None = None,
) -> None:
    """The closed-form step must predict every cell the method actually produced.

    :func:`critical_detectable_shift` is a rearrangement of the shipped gate's
    arithmetic, and the step table is built from it rather than from the grid.
    That is faster and exact, and it is also a second implementation, which is
    the standard way a benchmark ends up agreeing with itself. So: recompute
    every refusal count from the closed form and require it to equal the count
    the grid measured by calling ``fit``, in all cells.

    A disagreement means the rearrangement is wrong -- most likely because the
    gate stopped being monotone in ``s_detect``, or because a branch was added
    ahead of the width comparison -- and the step location, which is the whole
    deliverable, would be fiction.
    """
    from submission.expression_models import POISSON

    model = model or POISSON
    mismatches = {}
    for world in worlds:
        crit = np.array(
            [
                critical_detectable_shift(
                    generate_sample(world, seed=seed, replicate=r, model=model)
                )
                for r in range(n_replicates)
            ]
        )
        block = sweep[sweep["world"] == world.name].set_index("s_detect")
        for s in grid:
            predicted = int((crit < s).sum())
            measured = int(block.loc[float(s), "width_gate_refusals"])
            if predicted != measured:
                mismatches[(world.name, float(s))] = (predicted, measured)
    if mismatches:
        raise SweepDoesNotBracketError(
            "the closed-form critical s_detect does not reproduce the refusals "
            "the gate actually returned, at (world, s_detect) -> (predicted, "
            f"measured): {mismatches}. The step location in the steps table is "
            "derived from that closed form and cannot be trusted while this "
            "disagrees."
        )


def check_sweep_brackets_the_step(sweep: pd.DataFrame) -> None:
    """A sweep must contain both regimes, in a world where the estimand exists.

    Without this the sweep can report "the width gate never binds" when what
    happened is that the grid stopped below the transition. A sweep that only
    ever sees one regime has located nothing, and its table would nevertheless
    look like a result.

    ``annihilated`` is excluded because the finiteness guard refuses there at
    every ``s_detect`` by construction, so it brackets nothing.
    """
    live = sweep[sweep["truth_is_defined"]]
    if live.empty:
        raise SweepDoesNotBracketError(
            "no truth-defined world in the sweep; only the annihilated world "
            "is present, where the finiteness guard refuses at every s_detect "
            "and the grid can therefore not locate any transition."
        )
    bracketed = []
    for world, block in live.groupby("world", observed=True):
        rates = block["width_gate_refusal_rate"]
        if float(rates.min()) <= 0.0 and float(rates.max()) >= 1.0:
            bracketed.append(world)
    if not bracketed:
        spans = {
            str(w): (float(b["width_gate_refusal_rate"].min()),
                     float(b["width_gate_refusal_rate"].max()))
            for w, b in live.groupby("world", observed=True)
        }
        raise SweepDoesNotBracketError(
            "no truth-defined world where the grid spans both never-binding "
            f"(rate 0) and always-binding (rate 1). Observed spans: {spans}. "
            "The grid stopped short of the transition, so any statement about "
            "where the step falls is unsupported by this table."
        )


def check_committed_value_is_in_grid(sweep: pd.DataFrame) -> None:
    """The committed ``DETECTABLE_SHIFT`` must be a row of the sweep.

    A sweep that does not contain the shipped value cannot say what the shipped
    value does, and "the width gate does the job for free" is a claim about the
    shipped value specifically.
    """
    if not bool(sweep["is_committed_value"].any()):
        raise SweepDoesNotBracketError(
            f"the committed DETECTABLE_SHIFT ({COMMITTED_DETECTABLE_SHIFT}) is "
            f"not a point on this grid, so the sweep cannot be read as saying "
            f"anything about the method as shipped. Grid: "
            f"{sorted(sweep['s_detect'].unique())}"
        )


# ---------------------------------------------------------------------------
# 2. Misspecification, and the part of the result that is a theorem
# ---------------------------------------------------------------------------

#: The three gated/ungated variants of identical arithmetic. Deliberately not
#: DEFAULT_METHODS: naive-delta-mean and pseudobulk-de answer a different
#: question, and mixing them into an invariance claim about GATES would make
#: the claim about something else.
GATE_VARIANTS = (
    ("count_gate", KitagawaPositivityMethod()),
    ("width_gate", VarianceGateKitagawaMethod()),
    ("no_gate", KitagawaNoGateMethod()),
)


#: Seeds the misspecification sweep is run at. THREE, not one, and this is the
#: load-bearing decision in this section.
#:
#: Changing the count model changes how many variates the NORMAL arm consumes,
#: which shifts the stream before the tumour arm's maturity draw. So two models
#: do not see the same realised ``n_mature_tumour`` even though the
#: distribution of it -- Binomial(2000, f_t) -- does not depend on the model at
#: all. A single-seed table therefore shows the count gate refusing 88 times
#: under Poisson and 98 under NB(2) and invites exactly one reading, which is
#: the wrong one: the difference is the seed, not the generator.
#:
#: With three seeds the two explanations come apart in the table itself. The
#: count gate's spread across MODELS turns out to sit inside its spread across
#: SEEDS; the width gate's binding under severe overdispersion does not.
MISSPECIFICATION_SEEDS: tuple[int, ...] = (20260829, 20260830, 20260831)


def misspecification_sweep(
    *,
    seed: int | None = None,
    seeds: tuple[int, ...] = MISSPECIFICATION_SEEDS,
    n_replicates: int = 200,
    models: tuple[ExpressionModel, ...] = EXPRESSION_MODELS,
    worlds: tuple[BenchWorld, ...] = BENCH_WORLDS,
) -> pd.DataFrame:
    """The three gates on every world under every count model, at every seed.

    The method objects are the committed classes, constructed once and reused
    across models -- nothing about them is adapted to the generator, which is
    the entire point of the exercise.

    ``seed`` is accepted for symmetry with the other entry points and, when
    given, REPLACES the seed list rather than joining it -- so a caller asking
    for one seed gets one seed and then trips
    :func:`check_seeds_separate_model_from_stream`, rather than quietly getting
    three and believing it measured one.
    """
    seeds = (int(seed),) if seed is not None else tuple(int(s) for s in seeds)
    rows = []
    for run_seed in seeds:
        for model in models:
            for world in worlds:
                refusals = {label: 0 for label, _ in GATE_VARIANTS}
                numbers = {label: 0 for label, _ in GATE_VARIANTS}
                n_mature = []
                for replicate in range(n_replicates):
                    sample = generate_sample(
                        world, seed=run_seed, replicate=replicate, model=model
                    )
                    n_mature.append(sample.n_mature_tumour)
                    for label, method in GATE_VARIANTS:
                        out = method.fit(sample)
                        refusals[label] += int(out.refused)
                        numbers[label] += int(out.intrinsic is not None)
                for label, _ in GATE_VARIANTS:
                    rows.append(
                        {
                            "seed": int(run_seed),
                            "expression_model": model.name,
                            "model_kind": model.kind,
                            "dispersion": (
                                float(model.dispersion)
                                if model.dispersion is not None
                                else float("nan")
                            ),
                            "zero_inflation": float(model.zero_inflation),
                            "theoretical_variance_at_mu": float(
                                model.theoretical_variance(MEAN_NORMAL)
                            ),
                            "world": world.name,
                            "truth_is_defined": bool(world.truth_is_defined),
                            "gate": label,
                            "n_replicates": int(n_replicates),
                            "n_refused": int(refusals[label]),
                            "n_returned_a_number": int(numbers[label]),
                            "refusal_rate": refusals[label] / n_replicates,
                            # Binomial(N_CELLS, f_t) in expectation, and that
                            # distribution does not depend on the count model.
                            # Carried so a reader can see that the count gate's
                            # spread across models is a spread in the REALISED
                            # draw and not in what was being drawn from.
                            "mean_n_mature_tumour": float(np.mean(n_mature)),
                            "max_n_mature_tumour": int(np.max(n_mature)),
                        }
                    )
    return pd.DataFrame(rows)


def model_effect_vs_seed_noise(sweep: pd.DataFrame) -> pd.DataFrame:
    """Is a gate's spread across MODELS bigger than its spread across SEEDS?

    The table this summarises invites one particular misreading -- that the
    count gate is sensitive to the expression model, because its refusal count
    in ``depleted_wide`` moves from 88 to 98 as the model changes. It is not.
    Changing the model shifts the RNG stream before the maturity draw, so each
    model sees a different realisation of the same Binomial. The comparison
    that settles it is against the spread the same gate shows when only the
    seed changes, and that comparison is what this returns.

    The verdict is deliberately three-valued rather than boolean. A bare
    "model spread > seed spread" would have called the count gate
    model-sensitive on 16 against 14, which is not a separation, it is two
    draws from the same thing. The thresholds are 3x to claim the generator and
    1.5x to call it confounded, with an explicit AMBIGUOUS band between them --
    because a comparison forced to come out one way or the other is how a
    benchmark reports a result it does not have.
    """
    rows = []
    for (world, gate), block in sweep.groupby(["world", "gate"], observed=True):
        # Spread across models, worst over seeds -- so a single unlucky seed
        # cannot manufacture a model effect.
        per_seed = block.groupby("seed")["n_refused"]
        model_spread = float((per_seed.max() - per_seed.min()).max())
        # Spread across seeds, worst over models.
        per_model = block.groupby("expression_model")["n_refused"]
        seed_spread = float((per_model.max() - per_model.min()).max())
        rows.append(
            {
                "world": world,
                "gate": gate,
                "n_seeds": int(block["seed"].nunique()),
                "n_models": int(block["expression_model"].nunique()),
                "refusals_min": int(block["n_refused"].min()),
                "refusals_max": int(block["n_refused"].max()),
                "model_spread_worst_seed": model_spread,
                "seed_spread_worst_model": seed_spread,
                "model_to_seed_spread_ratio": _spread_ratio(
                    model_spread, seed_spread
                ),
                "attribution": _attribution(model_spread, seed_spread),
            }
        )
    return pd.DataFrame(rows).sort_values(["world", "gate"]).reset_index(drop=True)


#: Ratio of model spread to seed spread above which the difference is
#: attributed to the generator, and below which it is called confounded.
ATTRIBUTION_RATIO_CLAIM = 3.0
ATTRIBUTION_RATIO_CONFOUNDED = 1.5


def _spread_ratio(model_spread: float, seed_spread: float) -> float:
    """model/seed spread. NaN when neither varies; inf when only seed is 0."""
    if model_spread == 0.0 and seed_spread == 0.0:
        return float("nan")
    if seed_spread == 0.0:
        return float("inf")
    return model_spread / seed_spread


def _attribution(model_spread: float, seed_spread: float) -> str:
    """NO_VARIATION / GENERATOR / AMBIGUOUS / SEED_CONFOUNDED."""
    if model_spread == 0.0 and seed_spread == 0.0:
        return "NO_VARIATION"
    ratio = _spread_ratio(model_spread, seed_spread)
    if ratio >= ATTRIBUTION_RATIO_CLAIM:
        return "GENERATOR"
    if ratio <= ATTRIBUTION_RATIO_CONFOUNDED:
        return "SEED_CONFOUNDED"
    return "AMBIGUOUS"


def check_seeds_separate_model_from_stream(sweep: pd.DataFrame) -> None:
    """A misspecification sweep at one seed may not be published.

    At one seed the expression model and the RNG stream offset are perfectly
    confounded: a different model consumes a different number of variates in
    the normal arm, so it draws a different realisation of the tumour arm's
    maturity mask even though the distribution being drawn from is identical.
    Every difference in the table then has two explanations and the table
    cannot distinguish them -- while looking exactly like a table that can.

    This is the same failure the repository is named after, in the shape a
    benchmark takes it: not a check that cannot fail, but a comparison that
    cannot come out either way.
    """
    n_seeds = int(sweep["seed"].nunique()) if "seed" in sweep.columns else 0
    if n_seeds < 2:
        raise InvarianceViolation(
            f"this sweep has {n_seeds} distinct seed(s). At one seed the count "
            f"model and the RNG stream offset are perfectly confounded -- every "
            f"model draws a different realisation of the same Binomial maturity "
            f"mask -- so no difference in this table can be attributed to the "
            f"generator. Run at least two seeds "
            f"(extensions.MISSPECIFICATION_SEEDS)."
        )


#: Worlds whose refusal behaviour is asserted invariant to the count model.
#: Only the ones with ``frac_mature_tumour`` exactly 0 -- the theorem in
#: ``annihilated_invariance_proof`` is about that parameter and nothing else.
INVARIANT_WORLDS = tuple(
    w.name for w in BENCH_WORLDS if w.frac_mature_tumour == 0.0
)


def annihilated_invariance_proof() -> str:
    """Why the finiteness-guard result cannot depend on the count model.

    Not a summary of the measurement. The measurement is a corollary.

    THE CLAIM. In any world with ``frac_mature_tumour == 0.0`` exactly, the
    refusal decision of every gate in ``GATE_VARIANTS`` is identical for every
    count distribution -- not approximately, not with high probability, but as
    an identity, for any distribution whatsoever including ones nobody has
    written down.

    THE PROOF, in four steps.

    1. ``generate_sample`` labels a tumour cell mature by
       ``rng.random(N_CELLS) < frac_mature_tumour``. ``Generator.random``
       returns values in the half-open interval ``[0, 1)``, so every draw
       satisfies ``x >= 0``. With ``frac_mature_tumour == 0.0`` the comparison
       ``x < 0.0`` is False for every ``x`` in the range -- for every seed,
       every replicate, and every number the generator could ever produce.
       Hence ``n_mature_tumour == 0`` identically. This step uses no property
       of the expression model because the expression model has not been
       consulted yet.

    2. The model is invoked as ``model.draw(rng, mean, int(mature.sum()))``,
       i.e. with ``size = 0``. Whatever the distribution, a draw of zero
       variates is the empty array. The tumour arm is therefore exactly
       ``np.zeros(N_CELLS)``, again for every model.

    3. Each gate's decision on that sample reads only integers and the normal
       arm:

       * ``KitagawaPositivityMethod`` calls ``classify_estimability(0)``, a
         function of the count alone, which returns ``not_estimable``. Refuse.
       * ``VarianceGateKitagawaMethod`` tests ``n_t < 1`` BEFORE computing any
         moment, and returns on the finiteness guard. Refuse. The width
         threshold -- the only place ``DETECTABLE_SHIFT`` appears, and the only
         place a variance appears -- is never reached, which is also why the
         sweep in part 1 leaves this world flat.
       * ``KitagawaNoGateMethod`` substitutes ``0.0`` for the absent mean and
         returns a number. Do not refuse.

    4. None of those three branches reads the tumour expression values or any
       moment of them. Overdispersion, zero inflation and every other
       departure from Poisson are properties of a distribution that was sampled
       zero times.

    WHAT THIS CONVERTS. ``FINDINGS.md`` limitation 1 -- "it tests the gate, not
    robustness to a misspecified model" -- is, for this world, not a limitation
    but a theorem: the headline false-confidence result is robust to
    misspecification because it does not depend on the model at all. The gate
    fires on the ABSENCE of a cell, and an absent cell has no expression
    distribution to misspecify.

    WHAT IT DOES NOT COVER, and this matters more than the theorem. The
    argument is specific to ``frac_mature_tumour == 0``. In ``depleted_wide``
    mature cells exist, the width gate reaches its threshold, and the threshold
    is a function of the sampled variance -- so the width gate's *non-binding*
    is not invariant and is not protected by anything above. The measurement in
    :func:`misspecification_sweep` is what says where the two cases separate.
    """
    return annihilated_invariance_proof.__doc__ or ""


def check_annihilated_refusal_is_model_invariant(sweep: pd.DataFrame) -> None:
    """Refusal counts in a ``frac_mature_tumour == 0`` world must not move.

    This is the mechanised form of :func:`annihilated_invariance_proof`. It is
    a guard rather than a test because the theorem is about the *generator and
    the gates as they are*, and either could be edited: a generator that used
    ``<=`` instead of ``<``, or a width gate that computed its moments before
    checking finiteness, would break the proof silently while every existing
    test still passed.
    """
    block = sweep[sweep["world"].isin(INVARIANT_WORLDS)]
    if block.empty:
        raise InvarianceViolation(
            f"no rows from {list(INVARIANT_WORLDS)} in this sweep -- the "
            f"invariance claim is about those worlds and there is nothing here "
            f"to check it on."
        )
    models = set(block["expression_model"])
    if len(models) < 2:
        raise InvarianceViolation(
            f"only {sorted(models)} present. An invariance claim needs at least "
            f"two count models; with one, the check cannot fail and is worse "
            f"than no check."
        )
    offenders = {}
    for (world, gate), rows in block.groupby(["world", "gate"], observed=True):
        distinct = sorted(set(rows["n_refused"].astype(int)))
        if len(distinct) > 1:
            offenders[f"{world}/{gate}"] = dict(
                zip(rows["expression_model"], rows["n_refused"].astype(int))
            )
    if offenders:
        raise InvarianceViolation(
            "refusal counts moved with the count model in a world where no "
            "mature tumour cell exists. The expression distribution is sampled "
            "zero times there, so this cannot happen unless the generator or a "
            f"gate has changed shape. Offenders: {offenders}"
        )


# ---------------------------------------------------------------------------
# 3. The interval-width curve, over all n rather than five of them
# ---------------------------------------------------------------------------

#: The full curve's domain. 2 because neither interval is defined below it;
#: 200 because it is comfortably past the point where the two agree and the
#: reader should be able to see the curve flatten rather than be told it does.
CURVE_N_MIN, CURVE_N_MAX = 2, 200

#: How far the closed form may sit from the simulated rate, at a given n, on
#: the MEDIAN across the simulation's cells. The closed form is a floor -- it
#: assumes the per-patient values are normal enough that the bootstrap mean is,
#: and a rare transcript's per-patient delta is skewed -- so the simulation is
#: expected slightly above it at small n, and a divergence is readable rather
#: than fatal. 0.015 sits at about 1.5x the largest observed gap.
CLOSED_FORM_AGREEMENT_TOLERANCE = 0.015


def interval_width_curve(
    *,
    n_min: int = CURVE_N_MIN,
    n_max: int = CURVE_N_MAX,
    alpha: float = NOMINAL_ALPHA,
    simulated: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """``width_ratio`` and ``expected_false_positive_rate`` for every n.

    THE POINT OF EVALUATING THE WHOLE THING. Both functions take ``n`` and
    nothing else. No counts matrix, no gene, no cohort, no seed. The committed
    job evaluates them at the five n values the real cohorts happen to have,
    which presents a property of the arithmetic as if it were a property of
    those cohorts. Over the full domain it is visible for what it is: the
    percentile bootstrap's miscalibration is **a pure function of sample size**,
    fixed before any data is collected. At n=2 a nominal 5% test rejects about
    40% of the time, and no amount of care with the data changes that number.

    That is the same shape as the recovery-curve finding one layer over: a
    statistic that cannot report its own failure, because the failure is not in
    the data it is given.

    ``simulated``, when passed, is the committed calibration frame; its
    percentile rows are joined in per n so the closed form can be read against
    the measurement rather than beside it.
    """
    if n_min < 2:
        raise ClosedFormDivergence(
            f"n_min={n_min}: width_ratio is undefined below n=2, where neither "
            f"interval exists. A curve starting at 1 would be a curve with a "
            f"fabricated first point."
        )
    rows = []
    for n in range(int(n_min), int(n_max) + 1):
        fpr = expected_false_positive_rate(n, alpha)
        rows.append(
            {
                "n_patients": int(n),
                "width_ratio_vs_t": float(width_ratio(n, alpha)),
                "nominal_alpha": float(alpha),
                "closed_form_false_positive_rate": float(fpr),
                "excess_over_nominal": float(fpr - alpha),
                "inflation_factor": float(fpr / alpha),
                "verdict": calibration_verdict(fpr, alpha),
            }
        )
    curve = pd.DataFrame(rows)

    for column in (
        "simulated_false_positive_rate_median",
        "simulated_false_positive_rate_min",
        "simulated_false_positive_rate_max",
        "abs_diff_vs_median",
        "abs_diff_worst_cell",
    ):
        curve[column] = float("nan")
    curve["n_simulation_cells"] = 0

    if simulated is not None:
        pct = simulated[simulated["method"] == "percentile"]
        if pct.empty:
            raise ClosedFormDivergence(
                "the simulation frame has no percentile rows. The closed form "
                "describes the percentile bootstrap specifically; cross-checking "
                "it against bca or student_t rows would be comparing it to a "
                "different estimator and calling the disagreement agreement."
            )
        agg = pct.groupby("n_patients")["false_positive_rate"].agg(
            ["median", "min", "max", "count"]
        )
        index = curve.set_index("n_patients")
        for n, row in agg.iterrows():
            if n not in index.index:
                continue
            closed = float(index.loc[n, "closed_form_false_positive_rate"])
            cells = pct.loc[pct["n_patients"] == n, "false_positive_rate"]
            index.loc[n, "simulated_false_positive_rate_median"] = float(row["median"])
            index.loc[n, "simulated_false_positive_rate_min"] = float(row["min"])
            index.loc[n, "simulated_false_positive_rate_max"] = float(row["max"])
            index.loc[n, "abs_diff_vs_median"] = abs(float(row["median"]) - closed)
            index.loc[n, "abs_diff_worst_cell"] = float((cells - closed).abs().max())
            index.loc[n, "n_simulation_cells"] = int(row["count"])
        curve = index.reset_index()
    return curve


def check_closed_form_tracks_simulation(
    curve: pd.DataFrame, tolerance: float = CLOSED_FORM_AGREEMENT_TOLERANCE
) -> None:
    """The arithmetic and the simulation must describe one phenomenon.

    If they diverge, one of two things is true and both are worth stopping
    for: the closed form does not describe the estimator the project actually
    runs, or the simulation has drifted. Either way the curve stops being a
    statement about the committed bootstrap, and it is the curve -- being
    cheap, smooth and total -- that would be believed.

    Checked on the MEDIAN across the simulation's cells. The individual cells
    vary by abundance and by tau, and the closed form is a floor that skew
    pushes the cells above; policing every cell would fire on the skew, which
    is a real effect and is reported in ``abs_diff_worst_cell`` rather than
    treated as an error.
    """
    checked = curve[curve["n_simulation_cells"] > 0]
    if checked.empty:
        raise ClosedFormDivergence(
            "no n in this curve has a simulated rate to check against. An "
            "unchecked closed form is exactly the object this module is about: "
            "smooth, cheap, and unable to report its own failure."
        )
    bad = checked[checked["abs_diff_vs_median"] > tolerance]
    if not bad.empty:
        detail = {
            int(r.n_patients): (
                round(float(r.closed_form_false_positive_rate), 4),
                round(float(r.simulated_false_positive_rate_median), 4),
            )
            for r in bad.itertuples()
        }
        raise ClosedFormDivergence(
            f"closed form and simulated median differ by more than {tolerance} "
            f"at n = {sorted(detail)} (closed, simulated): {detail}. Either the "
            f"closed form is not describing the percentile bootstrap this "
            f"project runs, or the simulation has drifted."
        )
