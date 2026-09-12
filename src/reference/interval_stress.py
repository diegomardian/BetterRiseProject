"""The interval calibration under generators it was not built from.

    from src.reference.interval_stress import REGIMES, stress_calibration_table

WHY THIS EXISTS. ``interval_calibration`` established that the percentile
bootstrap over-rejects at small n and the Student-t interval does not, under a
generator whose per-patient effect is **normal** in the log fold change. The
audit decision record names the limit of that result and asks for this job:
``DECISION_2026-09-11_scope_and_pivot.md`` §4, week 2, *"Add skewed/heavy-tailed
and boundary regimes ... Quantify Monte Carlo uncertainty."* §3 of the same
record already says the committed result *"does not establish that BCa is
generally worse ... or that Student-t is universally calibrated."* This is the
measurement that says how far it does go.

WHAT VARIES. Only the distribution of the per-patient effect, through the same
Poisson thinning ``interval_calibration.simulate_deltas_from_log_fc`` uses:

``gaussian``          effect ~ Normal(0, tau). The committed control.
``skewed``            a centred, standardised log-normal. Real between-patient
                      effects are not symmetric in the log fold change.
``heavy_tailed``      effect ~ standardised Student-t. A few patients carry most
                      of the spread.
``spike_slab``        half the patients at exactly zero effect, half at a larger
                      one, scaled so the overall SD is still tau. The shape a
                      rare gene's patient-to-patient spread actually has.
``boundary_rare``     the gaussian generator at MLH1's abundance, where the
                      detection rate is near 0 and ``cloglog_rate``'s boundary
                      rule is doing work.
``boundary_saturated``the gaussian generator where detection is near 1.

Every regime is a **true null**: mean effect zero, so any rejection is a false
positive and the measured rate is a false-positive rate. The comparison is
between generators, not between hypotheses.

THREE FAMILIES OF CHECK, ON THE SAME SIMULATED DATA. The audit asks which
validation strategy would catch a miscalibrated interval:

``input_validation``  finite values, at least two patients, non-zero spread, a
                      finite interval. Ordinary hygiene.
``original_check``    the diagnostic the project ran before the audit: compare
                      the measured percentile rate to the closed form. It is
                      defined for the percentile bootstrap only.
``audit_check``       recalibrate under the regime and flag when the rate
                      exceeds nominal + tolerance, with the result required to
                      be above nominal on the Wilson interval as well.

The counts are recorded per cell, not asserted: a family that cannot fail shows
up as a detection rate of zero across every regime rather than as a sentence.

MONTE CARLO UNCERTAINTY. A rate over ``n_trials`` has a binomial standard error,
and every cell is repeated across seeds so seed-to-seed spread is visible rather
than assumed away. Both travel with every row.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats

from src.reference.interval_calibration import (
    INTERVAL_METHODS,
    MISCALIBRATION_TOLERANCE,
    NOMINAL_ALPHA,
    CalibrationError,
    calibration_verdict,
    expected_false_positive_rate,
    rejection_rate,
    simulate_deltas_from_log_fc,
    width_ratio,
)

#: Regimes whose effect is drawn by this module. ``kind`` selects the shape;
#: ``cp10k`` and the shape parameters are fixed per regime so a row's generator
#: is fully specified by its name.
_KINDS: frozenset[str] = frozenset(
    {"gaussian", "skewed", "heavy_tailed", "spike_slab"}
)


@dataclass(frozen=True)
class StressRegime:
    """One generator, fully specified by its name.

    ``cp10k`` is the gene's per-cell mean in the normal arm; it sets the
    detection rate. It is carried per regime rather than swept separately so
    that ``boundary_rare`` is not confused with a control at the same abundance.
    """

    name: str
    kind: str
    description: str
    cp10k: float
    tau_scale: float = 1.0
    df: float | None = None
    sigma: float | None = None
    spike_prob: float | None = None

    def __post_init__(self) -> None:
        if self.kind not in _KINDS:
            raise CalibrationError(
                f"regime {self.name!r} has unknown kind {self.kind!r}; "
                f"known: {sorted(_KINDS)}"
            )
        if self.kind == "heavy_tailed" and (self.df is None or self.df <= 2):
            raise CalibrationError(
                f"{self.name!r} is heavy-tailed and needs df > 2; with df <= 2 "
                f"the standardised effect has no finite variance, so 'SD = tau' "
                f"would be a statement about nothing."
            )
        if self.kind == "skewed" and (self.sigma is None or self.sigma <= 0):
            raise CalibrationError(f"{self.name!r} is skewed and needs sigma > 0")
        if self.kind == "spike_slab" and not (
            self.spike_prob is not None and 0 < self.spike_prob < 1
        ):
            raise CalibrationError(
                f"{self.name!r} is a spike-slab and needs 0 < spike_prob < 1"
            )


#: The grids the stress table runs. Ordered so the control reads first and the
#: two boundaries last.
REGIMES: tuple[StressRegime, ...] = (
    StressRegime(
        "gaussian", "gaussian",
        "effect ~ Normal(0, tau); the committed generator, as the control",
        cp10k=3.0,
    ),
    StressRegime(
        "skewed", "skewed",
        "centred standardised log-normal effect; asymmetric patient spread",
        cp10k=3.0, sigma=0.75,
    ),
    StressRegime(
        "heavy_tailed", "heavy_tailed",
        "standardised Student-t(4) effect; a few patients carry the spread",
        cp10k=3.0, df=4.0,
    ),
    StressRegime(
        "spike_slab", "spike_slab",
        "half the patients at zero effect, half larger; overall SD still tau",
        cp10k=3.0, spike_prob=0.5,
    ),
    StressRegime(
        "boundary_rare", "gaussian",
        "gaussian effect at MLH1's abundance; detection near 0 and the "
        "boundary rule is load-bearing",
        cp10k=0.039,
    ),
    StressRegime(
        "boundary_saturated", "gaussian",
        "gaussian effect where detection is near 1; no room to fall",
        cp10k=40.0,
    ),
)

#: The between-patient SD of the effect used across the grid. 0.2 is what
#: Pelka's control genes measure (``interval_heterogeneity``); 0 would make every
#: generator identical and the stress would measure nothing.
STRESS_TAU: float = 0.2

#: Detection family names, in the order a reader should meet them.
CHECK_FAMILIES: tuple[str, ...] = ("input_validation", "original_check", "audit_check")


def standardised_effects(
    regime: StressRegime, n: int, rng: np.random.Generator
) -> np.ndarray:
    """``n`` draws with mean 0 and variance 1, shaped by the regime.

    Standardising to unit variance is what makes ``tau`` mean the same thing
    under every shape: the only thing that differs between regimes is the
    shape, not the size of the patient-to-patient spread.
    """
    if regime.kind == "gaussian":
        return rng.normal(0.0, 1.0, n)
    if regime.kind == "skewed":
        sigma = float(regime.sigma)
        raw = np.exp(sigma * rng.normal(0.0, 1.0, n))
        mean = np.exp(sigma**2 / 2)
        sd = np.sqrt((np.exp(sigma**2) - 1) * np.exp(sigma**2))
        return (raw - mean) / sd
    if regime.kind == "heavy_tailed":
        df = float(regime.df)
        return rng.standard_t(df, n) / np.sqrt(df / (df - 2.0))
    if regime.kind == "spike_slab":
        pi = float(regime.spike_prob)
        active = rng.binomial(1, pi, n)
        return active * rng.normal(0.0, 1.0, n) / np.sqrt(pi)
    raise CalibrationError(f"no draw for kind {regime.kind!r}")


def draw_log_fc(
    regime: StressRegime, tau: float, n: int, rng: np.random.Generator
) -> np.ndarray:
    """Per-patient log fold changes: unit-variance shape scaled by ``tau``.

    ``regime.tau_scale`` lets a regime that is about the *boundary* rather than
    the effect shape still participate without a second definition of tau.
    """
    return tau * regime.tau_scale * standardised_effects(regime, n, rng)


def wilson_interval(
    successes: int, n: int, alpha: float = NOMINAL_ALPHA
) -> tuple[float, float]:
    """Wilson score interval for a binomial rate.

    Used instead of the normal approximation because the rates here sit near 5%
    and the Wald interval puts mass below zero, which would make a lower bound
    read as "worse than impossible".
    """
    if n < 1:
        raise CalibrationError("a Wilson interval needs at least one trial")
    p = successes / n
    z = float(stats.norm.ppf(1 - alpha / 2))
    denom = 1.0 + z**2 / n
    centre = (p + z**2 / (2 * n)) / denom
    half = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / denom
    return max(0.0, centre - half), min(1.0, centre + half)


def input_validation_problems(values: np.ndarray) -> list[str]:
    """What ordinary input validation would reject, if anything.

    This is the family the audit's stop rule excludes by name. It is run on the
    simulated values rather than asserted to pass, so that a generator which
    produced a degenerate sample would be caught here as well as by the audit.
    """
    values = np.asarray(values, dtype=float)
    problems = []
    if values.size < 2:
        problems.append("fewer than two patients")
    if not np.all(np.isfinite(values)):
        problems.append("non-finite value")
    if values.size >= 2 and float(values.std(ddof=1)) == 0.0:
        problems.append("no spread")
    return problems


def original_check_flags(measured: float, closed_form: float) -> bool | None:
    """The pre-audit diagnostic: does the closed form still explain the rate?

    Defined for the percentile bootstrap only -- the closed form is a property
    of the percentile interval, and asking BCa or Student-t to agree with it
    would be a category error. ``None`` records that, rather than a False that
    would read as "the check passed".
    """
    if not np.isfinite(closed_form):
        return None
    return bool(measured - closed_form > MISCALIBRATION_TOLERANCE)


def audit_flags(measured: float, wilson_lo: float) -> bool:
    """The audit's own rule: above tolerance, and above nominal on its interval.

    Requiring the Wilson lower bound to clear nominal is what stops a
    two-percentage-point excess that is really Monte Carlo noise from being
    recorded as a finding. Both conditions are necessary.
    """
    return bool(
        measured > NOMINAL_ALPHA + MISCALIBRATION_TOLERANCE
        and wilson_lo > NOMINAL_ALPHA
    )


def stress_cell(
    *, cohort: str, n_cells: np.ndarray, depth: np.ndarray, regime: StressRegime,
    method: str, seed: int, n_trials: int, tau: float = STRESS_TAU,
) -> dict[str, object]:
    """One (cohort, regime, method, seed) cell of the stress grid."""
    if method not in INTERVAL_METHODS:
        raise CalibrationError(
            f"unknown interval method {method!r}; have {sorted(INTERVAL_METHODS)}"
        )
    n_cells = np.asarray(n_cells, dtype=int)
    depth = np.asarray(depth, dtype=float)
    if n_cells.size < 2:
        raise CalibrationError(
            f"cohort {cohort!r} has {n_cells.size} patients; an interval over "
            f"fewer than two is undefined, so there is no false-positive rate to "
            f"measure."
        )

    problems = 0

    def sampler(r: np.random.Generator) -> np.ndarray:
        nonlocal problems
        log_fc = draw_log_fc(regime, tau, n_cells.size, r)
        deltas = simulate_deltas_from_log_fc(
            n_cells=n_cells, depth=depth, cp10k=regime.cp10k,
            log_fc=log_fc, rng=r,
        )
        if input_validation_problems(deltas):
            problems += 1
        return deltas

    rng = np.random.default_rng(seed)
    measured = rejection_rate(
        sampler, method=method, rng=rng, n_trials=n_trials
    )
    successes = int(round(measured * n_trials))
    wilson_lo, wilson_hi = wilson_interval(successes, n_trials)
    closed_form = (
        expected_false_positive_rate(int(n_cells.size))
        if method == "percentile" else float("nan")
    )
    return {
        "cohort": cohort,
        "n_patients": int(n_cells.size),
        "median_cells_per_arm": float(np.median(n_cells)),
        "regime": regime.name,
        "regime_kind": regime.kind,
        "cp10k_normal": float(regime.cp10k),
        "tau": float(tau),
        "method": method,
        "seed": int(seed),
        "n_trials": int(n_trials),
        "false_positive_rate": float(measured),
        "mc_se": float(np.sqrt(measured * (1 - measured) / n_trials)),
        "wilson_lo": float(wilson_lo),
        "wilson_hi": float(wilson_hi),
        "nominal_alpha": NOMINAL_ALPHA,
        "closed_form_rate": float(closed_form),
        "width_ratio_vs_t": (
            width_ratio(int(n_cells.size)) if method == "percentile"
            else float("nan")
        ),
        "verdict": calibration_verdict(measured),
        "input_validation_failures": int(problems),
        "original_check_flags": original_check_flags(measured, closed_form),
        "audit_flags": audit_flags(measured, wilson_lo),
    }


def stress_calibration_table(
    *, cohorts: dict[str, tuple[np.ndarray, np.ndarray]],
    regimes: tuple[StressRegime, ...] = REGIMES,
    methods: tuple[str, ...] = tuple(INTERVAL_METHODS),
    seeds: tuple[int, ...] = (20260101,), n_trials: int = 1_500,
    tau: float = STRESS_TAU,
) -> pd.DataFrame:
    """The full grid: one row per (cohort, regime, method, seed)."""
    rows = []
    for cohort, (n_cells, depth) in cohorts.items():
        for regime in regimes:
            for method in methods:
                for seed in seeds:
                    rows.append(stress_cell(
                        cohort=cohort, n_cells=n_cells, depth=depth,
                        regime=regime, method=method, seed=seed,
                        n_trials=n_trials, tau=tau,
                    ))
    return pd.DataFrame(rows)


def summarise_stress(frame: pd.DataFrame) -> pd.DataFrame:
    """Per (cohort, regime, method): the rate and what each family detected.

    Seed-to-seed spread is reported beside the binomial uncertainty so a verdict
    that moves with the seed is visible rather than averaged into a point.
    """
    out = (
        frame.groupby(["cohort", "n_patients", "regime", "regime_kind", "method"],
                      observed=True, sort=True)
        .agg(
            n_seeds=("seed", "nunique"),
            n_trials=("n_trials", "first"),
            fpr_median=("false_positive_rate", "median"),
            fpr_min=("false_positive_rate", "min"),
            fpr_max=("false_positive_rate", "max"),
            n_miscalibrated=("verdict", lambda s: int((s == "MISCALIBRATED").sum())),
            input_validation_fired=("input_validation_failures",
                                    lambda s: int((s > 0).sum())),
            # ``original_check_flags`` is object-dtyped because it is None where
            # the check is undefined; ``v is True`` counts only genuine fires and
            # never downcasts None to False through a deprecated fillna.
            original_check_fired=("original_check_flags",
                                  lambda s: int(s.map(lambda v: v is True).sum())),
            audit_fired=("audit_flags", lambda s: int(s.astype(bool).sum())),
        )
        .reset_index()
    )
    for family, column in (
        ("input_validation", "input_validation_fired"),
        ("original_check", "original_check_fired"),
        ("audit_check", "audit_fired"),
    ):
        out[f"detection_rate_{family}"] = out[column] / out["n_seeds"]
    return out.sort_values(
        ["cohort", "regime", "method"], ignore_index=True
    )
