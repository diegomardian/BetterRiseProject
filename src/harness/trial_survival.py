"""The one-line degeneracy check, and what censoring does to it.

``src/harness/trial_recovery`` shows that a recovery curve cannot tell an
estimator apart from its generator when the two are functions of the same
sufficient statistics, and that one column fixes it: compare the estimate
against the *realised* draw rather than against the requested parameter, and
assert the difference is not identically zero.

THIS MODULE IS THE COUNTEREXAMPLE TO THAT FIX. It is not a different failure --
it is the same check, on the endpoint that high-stakes trial simulation actually
uses, where the check quietly stops working.

THE GENERATOR. The strata and the confounded assignment of ``trial_recovery``,
with a survival outcome::

    T ~ Exponential(lambda_g * exp(beta * D)),  lambda = (0.10, 0.16, 0.30)

and administrative censoring at a fixed horizon ``c``: you observe
``(min(T, c), 1{T <= c})`` and nothing else. Sweeping ``c`` moves the censored
fraction from 0% to 78% without touching the effect being estimated.

TWO REALISED TRUTHS, AND THE ENTIRE POINT IS THAT THEY DIVERGE.

``theta_realised_latent``    The standardised log hazard ratio the draw actually
                            took, computed from the UNCENSORED latent event
                            times. This is what a simulator author naturally
                            logs, because the generator knows the latent times
                            and writing them to the truth column costs nothing.

``theta_realised_observed``  The same functional computed from the censored data
                            the analyst is handed -- events and exposure per
                            stratum-arm cell.

With no censoring they are the same number. Under censoring they are not, and
which one you put in the truth column decides whether the check works.

WHAT THE SWEEP SHOWS.

(i)   At 0% censoring the check fires correctly: the exponential MLE reproduces
      the realised truth exactly, to floating point, and is flagged.

(ii)  Censoring SILENCES it against the latent truth. The exponential MLE is
      still, exactly, the generator's own functional of the observed cell
      counts -- its residual against ``theta_realised_observed`` stays
      identically zero at every censoring level -- but its residual against
      ``theta_realised_latent`` grows with the censored fraction, because the
      latent truth is computed from information the estimator was never given.
      The check reads that sampling gap as evidence of independence.

(iii) By roughly half censoring the VERDICT INVERTS. The blind estimator's
      residual against the latent truth exceeds that of the stratified Cox model
      that is genuinely informative. A practitioner ranking estimators by the
      paper's own column would pick the blind one.

THE RULE THIS LICENSES, and it is one line: **in a censored simulator, evaluate
the realised truth on the OBSERVED data, not on the latent times.** The latent
truth is more precise, more available and wrong for this purpose.

THE POSITIVE CONTROL THAT KEEPS (ii) HONEST. Censoring does not silence the
check for everything, and the module includes the case where it does not.
Restricted mean survival to a horizon at or below ``c`` is determined exactly by
the observed data under purely administrative censoring, so ``km-rmst`` is
flagged degenerate against BOTH truths at every censoring level, while the
exponential-model RMST at the same horizon -- which reads the tail through a
fitted rate -- is silenced exactly like the log-hazard-ratio case. The breakdown
is specific to functionals that depend on the unobserved tail, not to censoring
as such.

Estimators are hand-rolled in closed form. ``lifelines`` is deliberately not
added to ``env/w2_harness.yml``: a claim that an estimator is *identically* the
generator's functional needs exact control of the functional, not a library's
choice of tie handling. It is used in ``tests/test_trial_survival.py`` only, as
an independent cross-check of the Cox implementation.
"""

from __future__ import annotations

import argparse
import logging
import math
import sys
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from src.common.io import write_versioned_table

log = logging.getLogger(__name__)

#: Baseline hazard per stratum. Well separated, and not linear in the stratum
#: index, so nothing collapses to an accident of the spacing.
STRATUM_HAZARDS: tuple[float, ...] = (0.10, 0.16, 0.30)

#: P(treated | stratum). Increasing, so assignment is confounded with the
#: covariate -- the situation stratification exists for.
PROPENSITY: tuple[float, ...] = (0.25, 0.50, 0.75)

#: Stratum prevalences.
STRATUM_WEIGHTS: tuple[float, ...] = (0.5, 0.3, 0.2)

#: The requested effect. log(0.6) is a hazard ratio a trial would be powered for.
BETA_REQUESTED: float = math.log(0.6)

#: The horizon standing in for "no administrative censoring". With the slowest
#: cell hazard at 0.06, P(T > 1000) is about 1e-26, so no draw in any sweep this
#: module runs is censored. Finite rather than ``inf`` because the sidecar is
#: JSON and ``Infinity`` is not valid JSON.
NO_CENSORING_TIME: float = 1000.0

#: Censored fractions to hit, marginal over strata and arms.
CENSORING_TARGETS: tuple[float, ...] = (0.0, 0.12, 0.49, 0.78)

#: Restricted-mean horizon. Fixed across every configuration and pre-specified,
#: not chosen per draw: at or below the smallest censoring time in the grid, so
#: RMST is estimable everywhere and the horizon is never a free parameter.
RMST_TAU: float = 2.0

COHORT_SIZES: tuple[int, ...] = (100, 200, 500, 1000, 2000, 3000, 5000)

#: The cohort size the headline table is read at.
HEADLINE_N: int = 3000


# ---------------------------------------------------------------------------
# The censoring grid
# ---------------------------------------------------------------------------


def marginal_censored_fraction(
    censoring_time: float,
    *,
    beta: float = BETA_REQUESTED,
    stratum_hazards: tuple[float, ...] = STRATUM_HAZARDS,
    propensity: tuple[float, ...] = PROPENSITY,
    stratum_weights: tuple[float, ...] = STRATUM_WEIGHTS,
) -> float:
    """P(T > c) in closed form, marginal over strata and arms.

    Closed form rather than simulated, so the grid is a property of the
    generator rather than of whichever seed was used to tune it.
    """
    total = 0.0
    for weight, hazard, p_treated in zip(
        stratum_weights, stratum_hazards, propensity, strict=True
    ):
        treated_rate = hazard * math.exp(beta)
        total += weight * (
            p_treated * math.exp(-treated_rate * censoring_time)
            + (1.0 - p_treated) * math.exp(-hazard * censoring_time)
        )
    return float(total)


def censoring_time_for(target_fraction: float, **kwargs) -> float:
    """The horizon ``c`` whose marginal censored fraction is ``target_fraction``.

    Bisection on a strictly decreasing function. ``target_fraction == 0`` has no
    finite solution and returns :data:`NO_CENSORING_TIME`.
    """
    if not 0.0 <= target_fraction < 1.0:
        raise ValueError(f"censored fraction must be in [0, 1), got {target_fraction}")
    if target_fraction == 0.0:
        return NO_CENSORING_TIME
    low, high = 1e-9, NO_CENSORING_TIME
    for _ in range(200):
        mid = 0.5 * (low + high)
        if marginal_censored_fraction(mid, **kwargs) > target_fraction:
            low = mid
        else:
            high = mid
    return 0.5 * (low + high)


# ---------------------------------------------------------------------------
# The generator
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class LatentCohort:
    """One draw before any horizon is applied. Sorted by latent event time.

    Sorting once here is not cosmetic: ``min(T, c)`` is monotone, so every
    censored view of this cohort is already sorted too, and the Cox risk sets
    below are cumulative sums rather than searches.
    """

    stratum: np.ndarray
    treated: np.ndarray
    latent_time: np.ndarray
    beta_requested: float


def simulate_latent_cohort(
    n_patients: int,
    beta: float = BETA_REQUESTED,
    *,
    rng: np.random.Generator,
    stratum_hazards: tuple[float, ...] = STRATUM_HAZARDS,
    propensity: tuple[float, ...] = PROPENSITY,
    stratum_weights: tuple[float, ...] = STRATUM_WEIGHTS,
) -> LatentCohort:
    """Confounded strata, then exponential event times with a known log-HR."""
    stratum = rng.choice(len(stratum_weights), size=n_patients, p=list(stratum_weights))
    treated = (rng.random(n_patients) < np.asarray(propensity)[stratum]).astype(int)
    rate = np.asarray(stratum_hazards)[stratum] * np.exp(beta * treated)
    latent = rng.exponential(1.0 / rate)

    order = np.argsort(latent, kind="stable")
    return LatentCohort(
        stratum=stratum[order],
        treated=treated[order],
        latent_time=latent[order],
        beta_requested=float(beta),
    )


@dataclass(frozen=True)
class CensoredTrial:
    """One cohort seen through a horizon, plus both realised truths on both scales.

    ``theta_realised_latent`` is computed from ``LatentCohort.latent_time`` --
    the column a simulator author writes because the generator has it.
    ``theta_realised_observed`` is computed from ``records`` alone.
    """

    records: pd.DataFrame
    censoring_time: float
    censored_fraction: float
    beta_requested: float
    theta_realised_latent: dict[str, float]
    theta_realised_observed: dict[str, float]


def censor(
    cohort: LatentCohort, censoring_time: float, *, rmst_tau: float = RMST_TAU
) -> CensoredTrial:
    """Apply an administrative horizon and compute both truths.

    The truths are not separate formulas: each is the natural realised-effect
    functional on its scale, applied once to the latent times and once to the
    observed data. That two applications of one formula disagree is the finding.
    """
    time = np.minimum(cohort.latent_time, censoring_time)
    event = (cohort.latent_time <= censoring_time).astype(int)
    records = pd.DataFrame(
        {
            "stratum": cohort.stratum,
            "treated": cohort.treated,
            "time": time,
            "event": event,
        }
    )
    latent_records = pd.DataFrame(
        {
            "stratum": cohort.stratum,
            "treated": cohort.treated,
            "time": cohort.latent_time,
            "event": np.ones(len(cohort.latent_time), dtype=int),
        }
    )
    return CensoredTrial(
        records=records,
        censoring_time=float(censoring_time),
        censored_fraction=float(1.0 - event.mean()),
        beta_requested=cohort.beta_requested,
        theta_realised_latent={
            "log_hazard_ratio": exponential_mle_standardised(latent_records),
            "rmst_difference": empirical_rmst_standardised(latent_records, tau=rmst_tau),
        },
        theta_realised_observed={
            "log_hazard_ratio": exponential_mle_standardised(records),
            "rmst_difference": km_rmst_standardised(records, tau=rmst_tau),
        },
    )


# ---------------------------------------------------------------------------
# Estimators. Closed form, no library, so the functional is exactly what is
# written here.
# ---------------------------------------------------------------------------


def _standardise(
    stratum: np.ndarray, treated: np.ndarray, cell_value: Callable[[np.ndarray], float]
) -> float:
    """Prevalence-weighted difference of a per-cell summary.

    Returns NaN if any cell has no value -- an arm missing from a stratum is a
    positivity failure, and the estimand has no referent there. It is not zero.
    """
    n = len(stratum)
    if n == 0:
        return float("nan")
    total = 0.0
    for g in np.unique(stratum):
        in_stratum = stratum == g
        weight = in_stratum.sum() / n
        treated_value = cell_value(in_stratum & (treated == 1))
        control_value = cell_value(in_stratum & (treated == 0))
        if not (np.isfinite(treated_value) and np.isfinite(control_value)):
            return float("nan")
        total += weight * (treated_value - control_value)
    return float(total)


def _columns(
    records: pd.DataFrame,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    return (
        records["stratum"].to_numpy(),
        records["treated"].to_numpy(),
        records["time"].to_numpy(dtype=float),
        records["event"].to_numpy(dtype=float),
    )


def exponential_mle_standardised(records: pd.DataFrame) -> float:
    """Per stratum-arm exponential MLE ``d / sum(t)``, standardised over strata.

    DEGENERATE ON THE OBSERVED DATA, BY CONSTRUCTION. For exponential times with
    independent censoring, ``(events, exposure)`` per cell is exactly sufficient,
    and it is exactly what the generator's draw determines -- so this returns
    ``theta_realised_observed["log_hazard_ratio"]`` identically, at every
    censoring level. It is also the estimator anyone writes first when they know
    the generator is exponential.
    """
    stratum, treated, time, event = _columns(records)

    def log_rate(cell: np.ndarray) -> float:
        events = event[cell].sum()
        exposure = time[cell].sum()
        if events <= 0 or exposure <= 0:
            return float("nan")  # no events: the rate is not estimable, not zero
        return float(np.log(events / exposure))

    return _standardise(stratum, treated, log_rate)


def _partial_likelihood(
    beta: float, x: np.ndarray, time: np.ndarray, event: np.ndarray
) -> tuple[float, float, float]:
    """Breslow partial log-likelihood, score and Hessian for one stratum.

    ``time`` must be ascending. Risk sets are reverse cumulative sums, and ties
    at an event time share the Breslow denominator.
    """
    if len(time) == 0:
        return 0.0, 0.0, 0.0
    weight = np.exp(beta * x)
    s0 = np.cumsum(weight[::-1])[::-1]
    s1 = np.cumsum((weight * x)[::-1])[::-1]
    s2 = np.cumsum((weight * x * x)[::-1])[::-1]

    is_event = event > 0
    if not is_event.any():
        return 0.0, 0.0, 0.0
    event_times = time[is_event]
    unique_times, counts = np.unique(event_times, return_counts=True)
    at = np.searchsorted(time, unique_times, side="left")
    d0, d1, d2 = s0[at], s1[at], s2[at]
    sum_x = np.bincount(
        np.searchsorted(unique_times, event_times),
        weights=x[is_event],
        minlength=len(unique_times),
    )
    ratio = d1 / d0
    loglik = float(beta * sum_x.sum() - (counts * np.log(d0)).sum())
    score = float(sum_x.sum() - (counts * ratio).sum())
    hessian = float(-(counts * (d2 / d0 - ratio**2)).sum())
    return loglik, score, hessian


def _cox(
    stratum: np.ndarray, treated: np.ndarray, time: np.ndarray, event: np.ndarray
) -> float:
    """Newton on the stratified Breslow partial likelihood. One covariate.

    Requires ``time`` ascending within each stratum, which every trial in this
    module satisfies because the cohort is sorted at generation and ``min(t, c)``
    preserves the order.

    Returns NaN rather than a number when the partial likelihood carries no
    information about the covariate -- no events, or no within-stratum contrast.
    """
    groups = [
        (treated[stratum == g].astype(float), time[stratum == g], event[stratum == g])
        for g in np.unique(stratum)
    ]
    beta = 0.0
    for _ in range(100):
        loglik = score = hessian = 0.0
        for x, t, e in groups:
            contribution = _partial_likelihood(beta, x, t, e)
            loglik += contribution[0]
            score += contribution[1]
            hessian += contribution[2]
        if hessian >= -1e-12 or not np.isfinite(hessian):
            return float("nan")
        step = float(np.clip(-score / hessian, -5.0, 5.0))
        # Step halving. The partial likelihood is concave so plain Newton
        # converges, but a huge first step off a near-flat Hessian would
        # otherwise be accepted on its way to an overflow.
        for _ in range(40):
            candidate = beta + step
            trial_loglik = sum(
                _partial_likelihood(candidate, x, t, e)[0] for x, t, e in groups
            )
            if np.isfinite(trial_loglik) and trial_loglik >= loglik - 1e-12:
                break
            step *= 0.5
        else:
            return float("nan")
        beta = beta + step
        if abs(step) < 1e-11:
            break
        if abs(beta) > 40.0:
            return float("nan")
    return float(beta)


def stratified_cox_breslow(records: pd.DataFrame) -> float:
    """Cox partial likelihood with the generator's strata, Breslow ties.

    INFORMATIVE. The partial likelihood conditions on the risk sets, so it is
    not a function of the cell exposures the generator's draw fixed; its residual
    against the realised truth is non-zero at every censoring level, which is the
    check working. It is also the estimator the analyst would actually use.
    """
    return _cox(*_columns(records))


def unadjusted_cox(records: pd.DataFrame) -> float:
    """Cox ignoring strata. The confounded control.

    Assignment depends on stratum and the strata differ in baseline hazard, so
    this is biased for ``beta``. Its score test at ``beta = 0`` is the log-rank
    test, which is why no separate log-rank arm is carried: the two are the same
    statistic and would say the same thing about degeneracy here.
    """
    stratum, treated, time, event = _columns(records)
    return _cox(np.zeros_like(stratum), treated, time, event)


def _km_rmst(time: np.ndarray, event: np.ndarray, tau: float) -> float:
    """Restricted mean survival on ``[0, tau]`` from the Kaplan-Meier estimate.

    NaN when ``tau`` exceeds the last observation: the curve is undefined past
    it, and the alternative is an extrapolation the data do not support.
    """
    if len(time) == 0 or tau > time.max():
        return float("nan")
    order = np.argsort(time, kind="stable")
    time, event = time[order], event[order]
    is_event = (event > 0) & (time <= tau)
    if not is_event.any():
        return float(tau)  # no events before tau: S == 1 throughout
    unique_times, counts = np.unique(time[is_event], return_counts=True)
    at_risk = len(time) - np.searchsorted(time, unique_times, side="left")
    survival = np.cumprod(1.0 - counts / at_risk)
    edges = np.concatenate(([0.0], unique_times, [tau]))
    levels = np.concatenate(([1.0], survival))
    return float((levels * np.diff(edges)).sum())


def km_rmst_standardised(records: pd.DataFrame, *, tau: float = RMST_TAU) -> float:
    """Kaplan-Meier RMST difference, standardised over strata.

    DEGENERATE ON THE OBSERVED DATA -- and, under purely administrative
    censoring, degenerate against the LATENT truth too. Censoring at ``c >= tau``
    never bites inside ``[0, tau]``, so this equals the empirical truncated mean
    of the latent times exactly. It is the positive control for the whole
    argument: the check is not silenced by censoring in general, only for
    functionals that depend on the unobserved tail.

    On a different scale from the log-hazard-ratio arms. Its residuals are in
    time units and are not comparable in magnitude with theirs.
    """
    stratum, treated, time, event = _columns(records)
    return _standardise(
        stratum, treated, lambda cell: _km_rmst(time[cell], event[cell], tau)
    )


def empirical_rmst_standardised(
    records: pd.DataFrame, *, tau: float = RMST_TAU
) -> float:
    """RMST difference from the empirical truncated mean of the times.

    The realised RMST truth on uncensored latent times. Identical to
    :func:`km_rmst_standardised` when nothing is censored before ``tau``, which
    is the fact the positive control rests on.
    """
    stratum, treated, time, _event = _columns(records)
    return _standardise(
        stratum,
        treated,
        lambda cell: (
            float(np.minimum(time[cell], tau).mean()) if cell.any() else float("nan")
        ),
    )


def exponential_rmst_standardised(
    records: pd.DataFrame, *, tau: float = RMST_TAU
) -> float:
    """RMST implied by the fitted exponential rate, standardised.

    ``(1 - exp(-rate * tau)) / rate`` per cell. Same scale as
    :func:`km_rmst_standardised` and the same data, but it reaches the horizon
    through a fitted rate, so it depends on the censored tail. It exists so the
    RMST scale has a contrast: the non-parametric arm keeps the check working
    under censoring and this one does not.
    """
    stratum, treated, time, event = _columns(records)

    def rmst(cell: np.ndarray) -> float:
        events = event[cell].sum()
        exposure = time[cell].sum()
        if events <= 0 or exposure <= 0:
            return float("nan")
        rate = events / exposure
        return float((1.0 - np.exp(-rate * tau)) / rate)

    return _standardise(stratum, treated, rmst)


@dataclass(frozen=True)
class Arm:
    """An estimator, the scale its number lives on, and what it is here to show."""

    fn: Callable[[pd.DataFrame], float]
    scale: str
    role: str


ESTIMATORS: dict[str, Arm] = {
    "exponential-mle-standardised": Arm(
        exponential_mle_standardised,
        "log_hazard_ratio",
        "the generator's own functional",
    ),
    "stratified-cox-breslow": Arm(
        stratified_cox_breslow, "log_hazard_ratio", "informative"
    ),
    "unadjusted-cox": Arm(unadjusted_cox, "log_hazard_ratio", "confounded control"),
    "km-rmst-standardised": Arm(
        km_rmst_standardised, "rmst_difference", "the generator's own functional"
    ),
    "exponential-rmst-standardised": Arm(
        exponential_rmst_standardised, "rmst_difference", "informative"
    ),
}

#: Whether the estimate is identically the realised truth computed on the
#: OBSERVED data. Measured against the residual in tests/test_trial_survival.py,
#: not trusted from this dict.
DEGENERATE_ON_OBSERVED: dict[str, bool] = {
    "exponential-mle-standardised": True,
    "stratified-cox-breslow": False,
    "unadjusted-cox": False,
    "km-rmst-standardised": True,
    "exponential-rmst-standardised": False,
}

#: The rule each verdict licenses, carried in the table so a reader of the
#: parquet does not have to come back to this docstring for it.
RULE = (
    "In a censored simulator, evaluate the realised truth on the OBSERVED data. "
    "The latent-time truth is more precise, more available, and wrong for this "
    "purpose: its residual is dominated by what censoring hid, not by whether "
    "the estimator is independent of the generator."
)


# ---------------------------------------------------------------------------
# The sweep
# ---------------------------------------------------------------------------


def run(
    *,
    seed: int,
    n_seeds: int = 6,
    beta: float = BETA_REQUESTED,
    cohort_sizes: tuple[int, ...] = COHORT_SIZES,
    censoring_targets: tuple[float, ...] = CENSORING_TARGETS,
    n_replicates: int = 200,
    rmst_tau: float = RMST_TAU,
) -> pd.DataFrame:
    """One row per (estimator, censoring target, cohort size, seed, replicate).

    The latent cohort is drawn ONCE per (seed, cohort size, replicate) and then
    viewed through each horizon, so the censoring sweep is paired: the four rows
    for one replicate differ only in what was hidden, never in what was drawn.
    """
    horizons = [(target, censoring_time_for(target, beta=beta)) for target in censoring_targets]

    rows: list[dict] = []
    for offset in range(n_seeds):
        run_seed = seed + offset
        for n_patients in cohort_sizes:
            for replicate in range(n_replicates):
                rng = np.random.default_rng([run_seed, n_patients, replicate])
                cohort = simulate_latent_cohort(n_patients, beta, rng=rng)
                for target, horizon in horizons:
                    trial = censor(cohort, horizon, rmst_tau=rmst_tau)
                    for name, arm in ESTIMATORS.items():
                        estimate = arm.fn(trial.records)
                        latent = trial.theta_realised_latent[arm.scale]
                        observed = trial.theta_realised_observed[arm.scale]
                        rows.append(
                            {
                                "estimator": name,
                                "scale": arm.scale,
                                "role": arm.role,
                                "censoring_target": target,
                                "censoring_time": trial.censoring_time,
                                "censored_fraction": trial.censored_fraction,
                                "n_patients": n_patients,
                                "seed": run_seed,
                                "replicate": replicate,
                                "beta_requested": trial.beta_requested,
                                "theta_realised_latent": latent,
                                "theta_realised_observed": observed,
                                "estimate": estimate,
                                "residual_vs_latent": abs(estimate - latent),
                                "residual_vs_observed": abs(estimate - observed),
                            }
                        )
    return pd.DataFrame(rows)


def summarise(runs: pd.DataFrame) -> pd.DataFrame:
    """Per (estimator, censoring target, cohort size): both residual columns.

    ``max`` is reported because that is what a degeneracy claim needs -- one
    non-zero residual anywhere falsifies "identically zero". ``median`` and
    ``q95`` are reported beside it because a maximum over more replicates is
    mechanically larger, so a max is not comparable across replicate counts and
    a median is.
    """
    out = (
        runs.groupby(["estimator", "scale", "role", "censoring_target", "n_patients"])
        .agg(
            n_replicates=("replicate", "count"),
            censored_fraction=("censored_fraction", "mean"),
            max_residual_vs_latent=("residual_vs_latent", "max"),
            median_residual_vs_latent=("residual_vs_latent", "median"),
            q95_residual_vs_latent=("residual_vs_latent", lambda s: s.quantile(0.95)),
            max_residual_vs_observed=("residual_vs_observed", "max"),
            median_residual_vs_observed=("residual_vs_observed", "median"),
        )
        .reset_index()
    )
    out["verdict_vs_latent"] = np.where(
        out["max_residual_vs_latent"] < 1e-9,
        "flagged: shares the generator's functional",
        "passed: looks independent of the generator",
    )
    out["verdict_vs_observed"] = np.where(
        out["max_residual_vs_observed"] < 1e-9,
        "flagged: shares the generator's functional",
        "passed: looks independent of the generator",
    )
    out["rule"] = RULE
    return out.sort_values(
        ["scale", "estimator", "censoring_target", "n_patients"]
    ).reset_index(drop=True)


#: The blind arm and the informative arm whose ordering is the headline.
BLIND = "exponential-mle-standardised"
INFORMATIVE = "stratified-cox-breslow"


def headline(
    runs: pd.DataFrame, *, blind: str = BLIND, informative: str = INFORMATIVE
) -> pd.DataFrame:
    """The ordering inversion, per (censoring target, cohort size).

    Paired within replicate: both estimators see the same draw and the same
    horizon, so ``inversion_rate`` -- the share of replicates on which the blind
    estimator's residual against the LATENT truth is the larger of the two -- is
    a statement about the ordering itself and does not move with the replicate
    count the way a maximum does.
    """
    keys = ["censoring_target", "n_patients", "seed", "replicate"]
    left = runs[runs["estimator"] == blind].set_index(keys)
    right = runs[runs["estimator"] == informative].set_index(keys)
    paired = pd.DataFrame(
        {
            "censored_fraction": left["censored_fraction"],
            "blind_vs_latent": left["residual_vs_latent"],
            "blind_vs_observed": left["residual_vs_observed"],
            "informative_vs_latent": right["residual_vs_latent"],
            "informative_vs_observed": right["residual_vs_observed"],
        }
    ).dropna()
    paired["blind_looks_worse"] = (
        paired["blind_vs_latent"] > paired["informative_vs_latent"]
    )

    out = (
        paired.reset_index()
        .groupby(["censoring_target", "n_patients"])
        .agg(
            n_replicates=("blind_vs_latent", "count"),
            censored_fraction=("censored_fraction", "mean"),
            max_blind_vs_latent=("blind_vs_latent", "max"),
            max_blind_vs_observed=("blind_vs_observed", "max"),
            max_informative_vs_latent=("informative_vs_latent", "max"),
            max_informative_vs_observed=("informative_vs_observed", "max"),
            median_blind_vs_latent=("blind_vs_latent", "median"),
            median_informative_vs_latent=("informative_vs_latent", "median"),
            inversion_rate=("blind_looks_worse", "mean"),
        )
        .reset_index()
    )
    out["ordering_by_latent_truth"] = np.where(
        out["max_blind_vs_latent"] > out["max_informative_vs_latent"],
        "INVERTED: the blind estimator looks worse than the informative one",
        "correct: the blind estimator looks better, as it should",
    )
    out["blind"] = blind
    out["informative"] = informative
    out["rule"] = RULE
    return out.sort_values(["censoring_target", "n_patients"]).reset_index(drop=True)


def _extra_meta(
    *, seeds: list[int], n_replicates: int, runs: pd.DataFrame, head: pd.DataFrame
) -> dict:
    """The sidecar body. Split out so a test can check it for reserved keys."""
    achieved = (
        runs.groupby("censoring_target")["censored_fraction"].mean().round(4).to_dict()
    )
    return {
        "seeds": seeds,
        "n_replicates_per_seed": n_replicates,
        "beta_requested": BETA_REQUESTED,
        "hazard_ratio_requested": float(math.exp(BETA_REQUESTED)),
        "stratum_hazards": list(STRATUM_HAZARDS),
        "propensity": list(PROPENSITY),
        "stratum_weights": list(STRATUM_WEIGHTS),
        "cohort_sizes": list(COHORT_SIZES),
        "censoring_targets": list(CENSORING_TARGETS),
        "censoring_times": {
            str(t): round(censoring_time_for(t), 6) for t in CENSORING_TARGETS
        },
        "censored_fraction_achieved": {str(k): v for k, v in achieved.items()},
        "rmst_tau": RMST_TAU,
        "headline_n_patients": HEADLINE_N,
        "estimators": {k: v.role for k, v in ESTIMATORS.items()},
        "degenerate_on_observed": DEGENERATE_ON_OBSERVED,
        "inversion_rate_at_headline_n": {
            str(row["censoring_target"]): float(row["inversion_rate"])
            for _, row in head[head["n_patients"] == HEADLINE_N].iterrows()
        },
        "rule": RULE,
        "SYNTHETIC": (
            "Simulated trial cohorts with an analytically known hazard ratio and "
            "administrative censoring. Nothing here is a result about any real "
            "trial, patient or endpoint."
        ),
        "what_this_answers": (
            "Whether the residual-against-realised-truth check that catches an "
            "estimator sharing sufficient statistics with its generator still "
            "works once the endpoint is censored. It does when the realised "
            "truth is computed on the observed data, and it does not when the "
            "truth is computed on the latent event times: the residual then "
            "grows with the censored fraction until the blind estimator ranks "
            "worse than an informative one."
        ),
    }


def main(argv: Sequence[str] | None = None) -> int:
    """Write the censoring sweep under ``results/`` with a provenance stamp.

        python -m src.harness.trial_survival
    """
    parser = argparse.ArgumentParser(
        description="What censoring does to the degeneracy check"
    )
    parser.add_argument("--seed", type=int, default=20260909)
    parser.add_argument("--replicates", type=int, default=200)
    parser.add_argument("--seeds", type=int, default=6, help="how many consecutive seeds")
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(message)s")
    runs = run(seed=args.seed, n_seeds=args.seeds, n_replicates=args.replicates)
    table = summarise(runs)
    head = headline(runs)

    at_n = head[head["n_patients"] == HEADLINE_N]
    log.info("n = %d, %d replicates per cell", HEADLINE_N, args.replicates * args.seeds)
    log.info(
        "%-10s %10s %12s %12s %12s %10s",
        "censoring", "achieved", "blind|latent", "blind|obs", "cox|latent", "inversion",
    )
    for _, row in at_n.iterrows():
        log.info(
            "%-10.2f %10.3f %12.4f %12.3g %12.4f %10.3f",
            row["censoring_target"],
            row["censored_fraction"],
            row["max_blind_vs_latent"],
            row["max_blind_vs_observed"],
            row["max_informative_vs_latent"],
            row["inversion_rate"],
        )

    seeds = [args.seed + i for i in range(args.seeds)]
    meta = _extra_meta(
        seeds=seeds, n_replicates=args.replicates, runs=runs, head=head
    )
    path = write_versioned_table(
        table,
        "trial_survival",
        seed=args.seed,
        results_dir=args.results_dir,
        allow_dirty=args.allow_dirty,
        extra_meta=meta,
    )
    log.info("wrote %s (%d rows)", path, len(table))
    head_path = write_versioned_table(
        head,
        "trial_survival_headline",
        seed=args.seed,
        results_dir=args.results_dir,
        allow_dirty=args.allow_dirty,
        extra_meta=meta,
    )
    log.info("wrote %s (%d rows)", head_path, len(head))
    return 0


if __name__ == "__main__":
    sys.exit(main())
