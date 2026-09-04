"""The recovery curve's blind spot, in a trial simulator rather than in cells.

The failure this project hit is not about single-cell data. It is about what a
recovery curve is a function of, and the setup that produces it is the default
one for validating an external control arm. So the demonstration belongs in that
vocabulary, on a generator small enough to reason about in closed form.

THE SETUP. A patient simulator draws covariate strata, assigns treatment with a
stratum-dependent probability -- so the strata confound, which is why anyone
standardises -- and draws outcomes with a specified mean per stratum plus a
treatment effect ``theta``. You then estimate ``theta`` and plot recovered
against requested.

FOUR ESTIMATORS, AND THE POINT IS THAT TWO OF THEM ARE SECRETLY THE SAME.

``gcomp_from_generator``   Standardisation over the strata the generator sampled
                           within, weighting the stratum means the draw
                           realised. This is the *same functional of the same
                           sufficient statistics the generator used*, so it
                           returns the realised effect exactly. Its recovery
                           curve is realised-over-requested: a property of the
                           generator's sampling, containing no estimator.

``ipw_saturated``          Inverse-probability weighting, propensity estimated
                           per stratum from the records. Looks like a different
                           method. It is not: with a saturated propensity model
                           this is algebraically the same functional, so it also
                           returns the realised effect exactly. WE DID NOT
                           EXPECT THIS, and it is the most useful thing here --
                           you cannot tell by reading the code whether your
                           estimator shares statistics with your generator.

``ipw_cross_fitted``       The same estimator with the propensity fitted on one
                           half of the cohort and applied to the other.
                           Consistent, and genuinely a different functional, so
                           it does not reproduce the realised effect. Its curve
                           carries information.

``unadjusted``             The difference in arm means, ignoring strata.
                           Non-degenerate and biased, because assignment depends
                           on stratum. Included so the curve has something it
                           *can* catch.

WHAT THE DEMONSTRATION SHOWS. The recovery curves of the first three are nearly
indistinguishable -- both sit near 1 and both tighten with cohort size. Only one
of them is telling you anything. ``residual_vs_realised`` separates them in one
column: identically zero for the two degenerate estimators and non-zero for the
others, at every cohort size. That is the one-line check, and it is the only
thing here that distinguishes an estimator from its generator.
"""

from __future__ import annotations

import argparse
import logging
import sys
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from src.common.io import write_versioned_table

log = logging.getLogger(__name__)

#: Stratum means of the untreated outcome. Well separated, so ignoring them is a
#: real mistake rather than a rounding error, and deliberately not linear in the
#: stratum index so that nothing collapses to an accident of the spacing.
STRATUM_MEANS: tuple[float, ...] = (10.0, 14.0, 25.0)

#: P(treated | stratum). Increasing, so treatment assignment is confounded with
#: the covariate -- the situation standardisation exists for.
PROPENSITY: tuple[float, ...] = (0.25, 0.50, 0.75)

#: Stratum prevalences.
STRATUM_WEIGHTS: tuple[float, ...] = (0.5, 0.3, 0.2)

OUTCOME_SD: float = 4.0


@dataclass(frozen=True)
class Trial:
    """One simulated cohort, plus both truths."""

    records: pd.DataFrame
    theta_requested: float
    theta_realised: float


def simulate_trial(
    n_patients: int,
    theta: float,
    *,
    rng: np.random.Generator,
    stratum_means: tuple[float, ...] = STRATUM_MEANS,
    propensity: tuple[float, ...] = PROPENSITY,
    stratum_weights: tuple[float, ...] = STRATUM_WEIGHTS,
    outcome_sd: float = OUTCOME_SD,
) -> Trial:
    """Draw a confounded cohort with a known treatment effect.

    ``theta_realised`` is the standardised effect the draw actually took: the
    stratum-weighted difference of realised arm means, using the realised
    stratum weights. It differs from ``theta`` by sampling alone. It is also,
    exactly, what ``gcomp_from_generator`` computes -- which is the whole point.
    """
    stratum = rng.choice(len(stratum_weights), size=n_patients, p=list(stratum_weights))
    treated = rng.random(n_patients) < np.asarray(propensity)[stratum]
    mean = np.asarray(stratum_means)[stratum] + theta * treated
    outcome = rng.normal(mean, outcome_sd)

    records = pd.DataFrame(
        {"stratum": stratum, "treated": treated.astype(int), "outcome": outcome}
    )
    return Trial(
        records=records,
        theta_requested=float(theta),
        theta_realised=_standardised_effect(records),
    )


def _standardised_effect(records: pd.DataFrame) -> float:
    """Stratum-weighted difference of realised arm means.

    Undefined if any stratum is missing an arm, which is a positivity failure --
    the same class of problem the rest of this project is about, and it returns
    NaN rather than a number here too.
    """
    total = len(records)
    effect = 0.0
    for _stratum, group in records.groupby("stratum"):
        treated = group.loc[group["treated"] == 1, "outcome"]
        control = group.loc[group["treated"] == 0, "outcome"]
        if treated.empty or control.empty:
            return float("nan")
        effect += (len(group) / total) * (treated.mean() - control.mean())
    return float(effect)


def gcomp_from_generator(trial: Trial) -> float:
    """Standardisation over the generator's own strata and realised means.

    DEGENERATE BY CONSTRUCTION. This is the same functional of the same
    sufficient statistics the generator used to build the draw, so it returns
    ``trial.theta_realised`` exactly -- not approximately. Everything the
    recovery curve then shows is the generator's sampling behaviour.

    This is not a strawman. It is the shortest path to an arm that runs, and the
    natural first thing to build when the simulator and the estimator are
    written by the same person in the same afternoon.
    """
    return _standardised_effect(trial.records)


def ipw_saturated(trial: Trial) -> float:
    """IPW with a saturated propensity model. DEGENERATE, and not obviously so.

    This looks like a different method from ``gcomp_from_generator``: different
    literature, different code, weights instead of strata means. With the
    propensity estimated per stratum it is algebraically the same functional of
    the same sufficient statistics, so it returns the realised effect exactly
    too, and its recovery curve is just as blind.

    We added this arm expecting it to be the non-degenerate comparison. It is
    the most useful thing in this file: whether an estimator shares sufficient
    statistics with the generator is not visible in the code, only in the
    residual.
    """
    records = trial.records
    weights = np.zeros(len(records))
    for _stratum, group in records.groupby("stratum"):
        p_hat = group["treated"].mean()
        if p_hat in (0.0, 1.0):
            return float("nan")  # positivity fails in this stratum
        idx = group.index.to_numpy()
        weights[idx] = np.where(group["treated"] == 1, 1.0 / p_hat, 1.0 / (1.0 - p_hat))

    treated = records["treated"].to_numpy() == 1
    outcome = records["outcome"].to_numpy()
    n = len(records)
    return float(
        (weights[treated] * outcome[treated]).sum() / n
        - (weights[~treated] * outcome[~treated]).sum() / n
    )


def ipw_cross_fitted(trial: Trial, *, rng: np.random.Generator | None = None) -> float:
    """IPW with the propensity fitted on one fold and applied to the other.

    NON-DEGENERATE. Cross-fitting is what an analyst actually does, and it makes
    the estimator a different functional: the weights applied to a patient come
    from data that patient is not in. Consistent, so its curve still sits near
    1 -- but it no longer reproduces the realised effect, so the curve is a
    statement about the estimator rather than about the generator.
    """
    records = trial.records
    rng = np.random.default_rng(0) if rng is None else rng
    fold = rng.permutation(len(records)) % 2

    total = 0.0
    for held in (0, 1):
        fit = records[fold != held]
        apply_to = records[fold == held]
        if fit.empty or apply_to.empty:
            return float("nan")
        contribution = 0.0
        for stratum, group in apply_to.groupby("stratum"):
            fitted = fit[fit["stratum"] == stratum]
            if fitted.empty:
                return float("nan")
            p_hat = fitted["treated"].mean()
            if p_hat in (0.0, 1.0):
                return float("nan")
            treated = group.loc[group["treated"] == 1, "outcome"].sum() / p_hat
            control = group.loc[group["treated"] == 0, "outcome"].sum() / (1.0 - p_hat)
            contribution += treated - control
        total += contribution
    return float(total / len(records))


def unadjusted(trial: Trial) -> float:
    """Difference in arm means, ignoring strata. Biased, and meant to be.

    Assignment depends on stratum, so this is confounded. It exists to show that
    a recovery curve *does* work when the estimator is not a function of the
    generator's own statistics -- the curve catches this one immediately.
    """
    records = trial.records
    treated = records.loc[records["treated"] == 1, "outcome"]
    control = records.loc[records["treated"] == 0, "outcome"]
    if treated.empty or control.empty:
        return float("nan")
    return float(treated.mean() - control.mean())


def ols_stratum_dummies(trial: Trial) -> float:
    """OLS of outcome on treatment plus additive stratum dummies.

    NON-DEGENERATE, AND NOTHING WAS SPLIT TO MAKE IT SO. This is the estimator
    most people would actually write -- regress the outcome on treatment and the
    covariate -- and for this generator it is correctly specified: the effect is
    additive and homogeneous, so the coefficient on ``treated`` is consistent
    for ``theta``.

    It is still not the generator's functional. Standardisation weights each
    stratum's arm-mean difference by that stratum's prevalence; OLS weights it by
    the within-stratum treatment variance, ``n_g * p_g * (1 - p_g)``. The
    propensities differ across strata, so the two weightings differ, and the
    estimate does not reproduce ``theta_realised``.

    This is here because ``ipw_cross_fitted`` alone leaves an opening: it breaks
    the degeneracy by *sample splitting*, so a reader can conclude the residual
    check merely detects whether folds were used. It does not. This estimator
    sees every record, splits nothing, holds nothing out -- and its residual is
    non-zero at every cohort size. What the check detects is a shared functional,
    not a shared sample.

    The caveat belongs with it: OLS is correctly specified *here* because the
    effect is homogeneous. Under effect heterogeneity its variance-weighted
    estimand is not the standardised one, and its curve would sit off 1 for a
    reason that has nothing to do with this paper.
    """
    records = trial.records
    strata = np.sort(records["stratum"].unique())
    if len(strata) < 2:
        return float("nan")
    columns = [np.ones(len(records)), records["treated"].to_numpy(dtype=float)]
    columns.extend(
        (records["stratum"].to_numpy() == stratum).astype(float)
        for stratum in strata[1:]
    )
    design = np.column_stack(columns)
    if np.linalg.matrix_rank(design) < design.shape[1]:
        return float("nan")
    beta, *_ = np.linalg.lstsq(design, records["outcome"].to_numpy(dtype=float), rcond=None)
    return float(beta[1])


ESTIMATORS: dict[str, Callable[[Trial], float]] = {
    "gcomp-from-generator": gcomp_from_generator,
    "ipw-saturated": ipw_saturated,
    "ipw-cross-fitted": ipw_cross_fitted,
    "ols-stratum-dummies": ols_stratum_dummies,
    "unadjusted": unadjusted,
}

#: Whether the estimator is a function of the generator's own sufficient
#: statistics. The property the recovery curve cannot see.
#: Measured, not asserted: every one of these is checked against the residual in
#: tests/test_trial_recovery.py, which is the only way the claim is worth
#: anything.
IS_DEGENERATE: dict[str, bool] = {
    "gcomp-from-generator": True,
    "ipw-saturated": True,
    "ipw-cross-fitted": False,
    "ols-stratum-dummies": False,
    "unadjusted": False,
}


def run(
    *,
    seed: int,
    theta: float = 3.0,
    cohort_sizes: tuple[int, ...] = (100, 200, 500, 1000, 2000, 5000),
    n_replicates: int = 200,
) -> pd.DataFrame:
    """One row per (estimator, cohort size, replicate).

    Carries both truths and the residual against the realised one, so the
    blindness is a column rather than a claim.
    """
    rows: list[dict] = []
    for n_patients in cohort_sizes:
        for replicate in range(n_replicates):
            rng = np.random.default_rng([seed, n_patients, replicate])
            trial = simulate_trial(n_patients, theta, rng=rng)
            if not np.isfinite(trial.theta_realised):
                continue
            for name, estimator in ESTIMATORS.items():
                estimate = estimator(trial)
                rows.append(
                    {
                        "estimator": name,
                        "shares_sufficient_statistics": IS_DEGENERATE[name],
                        "n_patients": n_patients,
                        "replicate": replicate,
                        "theta_requested": trial.theta_requested,
                        "theta_realised": trial.theta_realised,
                        "estimate": estimate,
                        "residual_vs_realised": abs(estimate - trial.theta_realised),
                        "recovery_ratio": estimate / trial.theta_requested,
                    }
                )
    return pd.DataFrame(rows)


def summarise(runs: pd.DataFrame) -> pd.DataFrame:
    """Per (estimator, cohort size): the curve, and the check that sees through it."""
    out = (
        runs.groupby(["estimator", "shares_sufficient_statistics", "n_patients"])
        .agg(
            n_replicates=("replicate", "count"),
            ratio_median=("recovery_ratio", "median"),
            ratio_q25=("recovery_ratio", lambda s: s.quantile(0.25)),
            ratio_q75=("recovery_ratio", lambda s: s.quantile(0.75)),
            max_residual_vs_realised=("residual_vs_realised", "max"),
        )
        .reset_index()
    )
    return out.sort_values(["estimator", "n_patients"]).reset_index(drop=True)


def main(argv: Sequence[str] | None = None) -> int:
    """Write the demonstration table under ``results/`` with a provenance stamp.

        python -m src.harness.trial_recovery
    """
    parser = argparse.ArgumentParser(description="The recovery curve in a trial simulator")
    parser.add_argument("--seed", type=int, default=20260901)
    parser.add_argument("--replicates", type=int, default=200)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(message)s")
    runs = run(seed=args.seed, n_replicates=args.replicates)
    table = summarise(runs)

    worst = runs.groupby("estimator")["residual_vs_realised"].max()
    for name, value in worst.items():
        log.info("%-22s max |theta_hat - theta_realised| = %.3g  (%s)",
                 name, value, "DEGENERATE" if IS_DEGENERATE[name] else "informative")

    path = write_versioned_table(
        table, "trial_recovery", seed=args.seed,
        results_dir=args.results_dir, allow_dirty=args.allow_dirty,
        extra_meta={
            "n_replicates": args.replicates,
            "theta_requested": 3.0,
            "cohort_sizes": [100, 200, 500, 1000, 2000, 5000],
            "stratum_means": list(STRATUM_MEANS),
            "propensity": list(PROPENSITY),
            "stratum_weights": list(STRATUM_WEIGHTS),
            "outcome_sd": OUTCOME_SD,
            "estimators": {k: ("shares the generator's sufficient statistics"
                               if v else "does not")
                           for k, v in IS_DEGENERATE.items()},
            "max_residual_vs_realised": {k: float(v) for k, v in worst.items()},
            "SYNTHETIC": "Simulated patients with an analytically known treatment "
                         "effect. Nothing here is a result about any real trial.",
            "what_this_answers": (
                "Whether a recovery curve can distinguish an estimator from its "
                "generator. It cannot when the two share sufficient statistics, "
                "and two of these four estimators do — including one that looks "
                "like a different method."
            ),
        },
    )
    log.info("wrote %s (%d rows)", path, len(table))
    return 0


if __name__ == "__main__":
    sys.exit(main())
