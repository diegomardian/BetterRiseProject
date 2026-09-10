"""Degeneracy is a property of the (estimator, truth, design) TRIPLE, not of the
estimator.

``src/harness/trial_recovery.py`` shows that a recovery curve goes blind when the
estimator is a deterministic function of the sufficient statistics the generator
drew: the estimate equals the realised truth identically, the estimator cancels,
and the curve measures the generator. It then ships a dict, ``IS_DEGENERATE``,
keyed by estimator name -- as though blindness were a property an estimator
carries around with it.

It is not, and this module is the counterexample. Four moves, each of which flips
a cell of that dict without touching a line of estimator code.

1. A SECOND REALISED TRUTH. ``theta_realised`` in ``trial_recovery`` is the
   *standardised* effect: stratum differences weighted by realised prevalence.
   That is not what OLS on stratum dummies targets. By Frisch--Waugh, OLS targets
   the variance-weighted contrast, ``w_g = n_g * phat_g * (1 - phat_g)``, and
   against *that* truth OLS is exact to machine precision. The paper's
   "informative" estimator is its most blind one, measured against the estimand
   it actually has. So the residual check needs a truth argument, and the honest
   object is a matrix.

2. THE DESIGN POINT. Under stratified permuted-block randomisation -- what
   ICH E9 tells a real trial to do -- ``phat_g = 0.5`` in every stratum exactly,
   the two weightings become proportional, and the two estimands coincide. OLS
   goes from informative to exactly blind because of how patients were
   *randomised*, with the analysis code untouched. ``propensity_spread`` walks
   continuously between the RCT and the paper's confounded design so the residual
   can be plotted crossing from machine zero to the paper's number. Bernoulli 1:1
   assignment is NOT enough: it gives ``phat_g = 0.5`` only in expectation.

3. TWO MORE ESTIMATORS THAT ARE SECRETLY THE SAME FUNCTIONAL. Saturated AIPW
   without cross-fitting has an augmentation term that is identically zero within
   each stratum, so the doubly-robust estimator collapses onto standardisation.
   Exact stratum matching with M opposite-arm units per unit converges onto the
   same functional as M grows, so blindness arrives *continuously* -- there is no
   line in the code where the estimator becomes degenerate. And IPW whose
   propensity is fitted by ``sklearn.LogisticRegression`` is degenerate in intent
   and never in arithmetic, because the default is L2-penalised: its residual is
   a function of a regularisation constant, and at small ``C`` it scores as more
   informative than cross-fitted IPW. A residual read as a quality score would
   rank the most-degenerate-in-intent estimator top. Turning the penalty off does
   not fix it either -- the residual then stops at the solver's convergence
   tolerance, four orders of magnitude above machine zero, so the number is set
   by two sklearn defaults and no statistical decision at all.

   AND THE ONE WE DID NOT EXPECT. Under block randomisation ``unadjusted`` --
   the raw difference in arm means, which ``trial_recovery`` includes precisely
   as the arm "the curve *can* catch" -- is exactly blind as well, because with
   ``n_g1 = n_g0`` the difference in arm means IS the standardised effect. Under
   the randomisation real trials run, the demonstration's own positive control
   goes blind and the curve has nothing left it can catch.

4. THE INFORMATION RATIO, ``rho = median|theta_hat - theta_realised| /
   median|theta_realised - theta_requested|``. The residual in units of the
   generator's own sampling error. For OLS it is ~0.13 and does not fall with
   cohort size: the OLS recovery curve is mostly generator sampling noise at
   every n. A non-zero residual means the estimator is not *identically* the
   generator; it does not mean the curve is measuring the estimator.

Plus effect heterogeneity, which is cheap and makes the point a fifth time: with
``theta_g = theta + tau * (g - gbar)``, OLS's recovery curve wanders off 1 as tau
grows while its residual against the variance-weighted truth stays at machine
zero. The curve moves for a reason that has nothing to do with estimator quality.

``simulate_trial`` is imported from ``trial_recovery`` unchanged and is still the
reference generator; ``simulate_trial_design`` generalises it and is tested to
reproduce it draw-for-draw at the default design point.
"""

from __future__ import annotations

import argparse
import logging
import sys
from collections.abc import Callable, Sequence
from pathlib import Path

import numpy as np
import pandas as pd

from src.common.io import write_versioned_table
from src.harness.trial_recovery import (
    OUTCOME_SD,
    PROPENSITY,
    STRATUM_MEANS,
    STRATUM_WEIGHTS,
    Trial,
    _standardised_effect,
    gcomp_from_generator,
    ipw_cross_fitted,
    ipw_saturated,
    ols_stratum_dummies,
    simulate_trial,
    unadjusted,
)

log = logging.getLogger(__name__)

__all__ = [
    "DEGENERACY",
    "DEGENERATE_BELOW",
    "DESIGNS",
    "ESTIMATORS",
    "TRUTHS",
    "Design",
    "aipw_saturated",
    "information_ratio",
    "ipw_logistic_l2",
    "matching_exact",
    "run_grid",
    # Re-exported deliberately: this module generalises the reference generator
    # and never replaces it, and every equivalence test imports it from here.
    "simulate_trial",
    "simulate_trial_design",
    "varweighted_effect",
]


# ---------------------------------------------------------------------------
# The second realised truth
# ---------------------------------------------------------------------------


def varweighted_effect(records: pd.DataFrame) -> float:
    """OLS's own finite-sample estimand, in closed form.

    ``sum_g w_g (ybar_g1 - ybar_g0) / sum_g w_g`` with
    ``w_g = n_g * phat_g * (1 - phat_g)``.

    This is exactly the coefficient on ``treated`` from a regression of outcome
    on treatment plus saturated stratum dummies -- Frisch--Waugh, applied to a
    design matrix whose only non-dummy column is the treatment indicator. It is
    an algebraic identity of the realised records, so it holds whatever the
    generator did: any theta, any propensity, homogeneous effect or not.

    Returns NaN under a positivity failure, matching ``_standardised_effect``: a
    stratum with one arm empty has ``w_g = 0`` and an undefined contrast, and
    silently dropping it would change the estimand rather than report a problem.
    """
    numerator = 0.0
    denominator = 0.0
    for _stratum, group in records.groupby("stratum"):
        treated = group.loc[group["treated"] == 1, "outcome"]
        control = group.loc[group["treated"] == 0, "outcome"]
        if treated.empty or control.empty:
            return float("nan")
        n_g = len(group)
        p_hat = len(treated) / n_g
        weight = n_g * p_hat * (1.0 - p_hat)
        numerator += weight * (treated.mean() - control.mean())
        denominator += weight
    if denominator == 0.0:
        return float("nan")
    return float(numerator / denominator)


#: The realised-truth functionals a residual can be taken against. ``theta`` is
#: not among them: the residual against ``theta_requested`` is the recovery
#: curve, which is the thing being audited.
TRUTHS: dict[str, Callable[[Trial], float]] = {
    "standardised": lambda trial: float(trial.theta_realised),
    "varweighted": lambda trial: varweighted_effect(trial.records),
}


# ---------------------------------------------------------------------------
# The generator, with a design knob
# ---------------------------------------------------------------------------


def propensity_at_spread(
    spread: float, propensity: tuple[float, ...] = PROPENSITY
) -> np.ndarray:
    """Interpolate between an RCT and the paper's confounded design.

    ``0.5 + spread * (p_g - 0.5)``. ``spread=0`` is 1:1 randomisation in every
    stratum, ``spread=1`` reproduces the paper's ``(0.25, 0.50, 0.75)`` --
    exactly in binary too, since every step is a halving.
    """
    base = np.asarray(propensity, dtype=float)
    return 0.5 + float(spread) * (base - 0.5)


def _block_sizes(
    n_patients: int, stratum_weights: tuple[float, ...], block_size: int
) -> np.ndarray:
    """Fixed stratum sizes, each a multiple of ``block_size``, summing to n.

    Fixed rather than multinomial: stratified block randomisation allocates
    within a stratum, and leaving the stratum sizes random would put sampling
    noise back into the prevalence weights for no reason. The sizes are the
    design, not a draw.
    """
    if block_size < 2 or block_size % 2 != 0:
        raise ValueError(f"block_size must be even and >= 2, got {block_size}")
    if n_patients % block_size != 0:
        raise ValueError(
            f"n_patients={n_patients} is not a multiple of block_size={block_size}; "
            "exact 1:1 in every stratum is then unattainable, and the design point "
            "would not be the one being claimed"
        )
    weights = np.asarray(stratum_weights, dtype=float)
    n_blocks_total = n_patients // block_size
    raw = weights / weights.sum() * n_blocks_total
    blocks = np.floor(raw).astype(int)
    # Leftover blocks go to the strata with the largest fractional part.
    order = np.argsort(-(raw - blocks))
    for i in range(n_blocks_total - int(blocks.sum())):
        blocks[order[i % len(blocks)]] += 1
    if (blocks <= 0).any():
        raise ValueError(
            f"n_patients={n_patients} with block_size={block_size} leaves a stratum "
            "empty; the design point is not attainable at this cohort size"
        )
    return blocks * block_size


def _block_assignment(
    sizes: np.ndarray, block_size: int, rng: np.random.Generator
) -> tuple[np.ndarray, np.ndarray]:
    """Stratified permuted-block randomisation: exact 1:1 inside every stratum.

    Returns ``(stratum, treated)``. Within each stratum the patients are cut into
    blocks of ``block_size``, each block filled with equal numbers of treated and
    control and then permuted. Exactly half of every stratum is treated, so
    ``phat_g = 0.5`` identically -- not in expectation, which is the whole point
    of blocking and the reason Bernoulli 1:1 is a different design point.
    """
    stratum_parts: list[np.ndarray] = []
    treated_parts: list[np.ndarray] = []
    block = np.repeat([0, 1], block_size // 2)
    for g, size in enumerate(sizes):
        stratum_parts.append(np.full(int(size), g, dtype=int))
        treated_parts.append(
            np.concatenate([rng.permutation(block) for _ in range(int(size) // block_size)])
        )
    return np.concatenate(stratum_parts), np.concatenate(treated_parts)


def _stratum_effects(
    theta: float, tau: float, stratum_weights: tuple[float, ...]
) -> np.ndarray:
    """``theta_g = theta + tau * (g - gbar)``, ``gbar`` the design-weighted mean
    stratum index.

    Centring on the design-weighted mean keeps the population average treatment
    effect equal to ``theta`` for every tau, so a curve that moves under tau is
    moving because the estimand moved and not because the requested effect did.
    """
    weights = np.asarray(stratum_weights, dtype=float)
    index = np.arange(len(weights), dtype=float)
    gbar = float((weights * index).sum() / weights.sum())
    return theta + float(tau) * (index - gbar)


def simulate_trial_design(
    n_patients: int,
    theta: float,
    *,
    rng: np.random.Generator,
    assignment: str = "bernoulli",
    propensity_spread: float = 1.0,
    tau: float = 0.0,
    block_size: int = 4,
    stratum_means: tuple[float, ...] = STRATUM_MEANS,
    propensity: tuple[float, ...] = PROPENSITY,
    stratum_weights: tuple[float, ...] = STRATUM_WEIGHTS,
    outcome_sd: float = OUTCOME_SD,
) -> Trial:
    """``simulate_trial`` with the design exposed.

    At ``assignment="bernoulli", propensity_spread=1.0, tau=0.0`` this consumes
    the same RNG draws in the same order as ``trial_recovery.simulate_trial`` and
    returns an identical cohort. That equivalence is asserted in
    ``tests/test_trial_blindness.py``, because a "generalisation" that quietly
    changed the baseline design point would make every comparison here a
    comparison between two different generators.

    ``theta_realised`` on the returned ``Trial`` is still the standardised
    effect, so every existing estimator and test keeps its meaning. The
    variance-weighted truth is a separate call, ``varweighted_effect``.
    """
    if assignment not in ("bernoulli", "block"):
        raise ValueError(f"assignment must be 'bernoulli' or 'block', got {assignment!r}")

    effects = _stratum_effects(theta, tau, stratum_weights)

    if assignment == "bernoulli":
        p = propensity_at_spread(propensity_spread, propensity)
        stratum = rng.choice(len(stratum_weights), size=n_patients, p=list(stratum_weights))
        treated = rng.random(n_patients) < p[stratum]
    else:
        sizes = _block_sizes(n_patients, stratum_weights, block_size)
        stratum, arm = _block_assignment(sizes, block_size, rng)
        treated = arm.astype(bool)

    mean = np.asarray(stratum_means)[stratum] + effects[stratum] * treated
    outcome = rng.normal(mean, outcome_sd)

    records = pd.DataFrame(
        {"stratum": stratum, "treated": treated.astype(int), "outcome": outcome}
    )
    return Trial(
        records=records,
        theta_requested=float(theta),
        theta_realised=_standardised_effect(records),
    )


class Design:
    """A named point in design space. Immutable, and it names itself."""

    __slots__ = ("name", "assignment", "propensity_spread", "tau", "block_size")

    def __init__(
        self,
        name: str,
        *,
        assignment: str = "bernoulli",
        propensity_spread: float = 1.0,
        tau: float = 0.0,
        block_size: int = 4,
    ) -> None:
        self.name = name
        self.assignment = assignment
        self.propensity_spread = float(propensity_spread)
        self.tau = float(tau)
        self.block_size = int(block_size)

    def draw(self, n_patients: int, theta: float, rng: np.random.Generator) -> Trial:
        return simulate_trial_design(
            n_patients,
            theta,
            rng=rng,
            assignment=self.assignment,
            propensity_spread=self.propensity_spread,
            tau=self.tau,
            block_size=self.block_size,
        )

    def as_meta(self) -> dict:
        return {
            "assignment": self.assignment,
            "propensity_spread": self.propensity_spread,
            "tau": self.tau,
            "block_size": self.block_size,
        }

    def __repr__(self) -> str:  # pragma: no cover - debugging aid
        return (
            f"Design({self.name!r}, assignment={self.assignment!r}, "
            f"propensity_spread={self.propensity_spread}, tau={self.tau})"
        )


#: The three design points the residual matrix is measured on.
#:
#: ``rct-bernoulli`` is not redundant with ``block-randomised``. Both are 1:1 in
#: expectation; only the blocked one is 1:1 in the realised counts, and only the
#: blocked one makes OLS exactly blind. The gap between those two rows is the
#: measurement, and it is why "we randomised 1:1" is not by itself the claim.
DESIGNS: dict[str, Design] = {
    "confounded-bernoulli": Design("confounded-bernoulli", propensity_spread=1.0),
    "rct-bernoulli": Design("rct-bernoulli", propensity_spread=0.0),
    "block-randomised": Design("block-randomised", assignment="block"),
}


# ---------------------------------------------------------------------------
# The new estimators
# ---------------------------------------------------------------------------


def aipw_saturated(trial: Trial) -> float:
    """Augmented IPW / doubly robust, saturated nuisances, no cross-fitting.

    DEGENERATE, and it arrives wearing the most defensive vocabulary in the
    field. The outcome model is ``muhat_gd = ybar_gd`` and the propensity is
    ``ehat_g = phat_g``: both saturated, both fitted on the full sample, which is
    what you get if you write the estimator without cross-fitting and let the
    nuisance models be as flexible as the covariate allows.

    Within stratum g the augmentation over the treated is
    ``(1/phat_g) * sum_{i in g, A=1} (Y_i - ybar_g1) = 0`` identically, and the
    same for the controls. So every augmentation term cancels and the estimator
    collapses onto plain standardisation -- the same functional as
    ``gcomp_from_generator``, reached from the opposite direction. Double
    robustness protects against nuisance misspecification. It does nothing at all
    about sharing sufficient statistics with the generator.
    """
    records = trial.records
    n = len(records)
    total = 0.0
    for _stratum, group in records.groupby("stratum"):
        treated = group.loc[group["treated"] == 1, "outcome"].to_numpy()
        control = group.loc[group["treated"] == 0, "outcome"].to_numpy()
        if treated.size == 0 or control.size == 0:
            return float("nan")
        n_g = len(group)
        e_g = treated.size / n_g
        mu1 = float(treated.mean())
        mu0 = float(control.mean())
        augmentation = (treated - mu1).sum() / e_g - (control - mu0).sum() / (1.0 - e_g)
        total += n_g * (mu1 - mu0) + augmentation
    return float(total / n)


def matching_exact(trial: Trial, *, m: int | None, rng: np.random.Generator) -> float:
    """Exact matching on stratum, ``m`` opposite-arm units per unit.

    Each unit is matched to ``m`` units of the opposite arm drawn with
    replacement from its own stratum, and the estimate is the mean matched
    difference over ALL units -- the ATE version, so it is comparable with the
    standardised truth rather than with the ATT.

    ``m=None`` means every opposite-arm unit in the stratum, at which point each
    treated unit contributes ``y_i - ybar_g0`` and each control ``ybar_g1 - y_j``;
    the sum telescopes to ``n_g (ybar_g1 - ybar_g0)`` and the estimator IS the
    standardised effect. So blindness here is not a property that switches on at
    some line of code: it is approached continuously as ``m`` grows, and there is
    no value of ``m`` at which the method changes.
    """
    if m is not None and m < 1:
        raise ValueError(f"m must be >= 1 or None, got {m}")
    records = trial.records
    n = len(records)
    total = 0.0
    for _stratum, group in records.groupby("stratum"):
        y1 = group.loc[group["treated"] == 1, "outcome"].to_numpy()
        y0 = group.loc[group["treated"] == 0, "outcome"].to_numpy()
        if y1.size == 0 or y0.size == 0:
            return float("nan")
        if m is None:
            total += (y1.size + y0.size) * (float(y1.mean()) - float(y0.mean()))
            continue
        picked0 = y0[rng.integers(0, y0.size, size=(y1.size, m))].mean(axis=1)
        picked1 = y1[rng.integers(0, y1.size, size=(y0.size, m))].mean(axis=1)
        total += float((y1 - picked0).sum() + (picked1 - y0).sum())
    return float(total / n)


def ipw_logistic_l2(trial: Trial, *, C: float, tol: float = 1e-4) -> float:
    """IPW with the propensity from ``sklearn.LogisticRegression`` on saturated
    stratum dummies.

    DEGENERATE IN INTENT, NEVER IN ARITHMETIC. The model is saturated, so the
    unpenalised MLE is exactly ``phat_g`` and this is exactly ``ipw_saturated``.
    But ``LogisticRegression``'s default is L2-penalised with ``C=1.0``, so the
    fitted propensity is shrunk toward the intercept by an amount nobody chose,
    and the residual against the realised truth becomes a smooth function of a
    regularisation constant.

    That makes the residual unusable as a quality score. At small ``C`` this
    estimator's residual exceeds cross-fitted IPW's, so a reader ranking
    estimators by residual would place the one that is algebraically the
    generator above the one that genuinely holds data out. The residual answers
    "is this estimator identically the generator's functional", and only that.

    ``tol`` IS THE SECOND KNOB, AND IT WAS NOT IN THE PLAN. Turning the penalty
    off (``C`` large) does not return the residual to machine zero: it stops at
    about 4e-3, and that floor is ``lbfgs``'s default convergence tolerance
    rather than anything statistical. Tightening ``tol`` to 1e-8 with the same
    ``C`` drops the residual by four more orders of magnitude. So this
    estimator's distance from the generator is set by two numerical constants,
    neither of which is a modelling decision, and both of which a reader would
    have to open the sklearn defaults to find. It is reported at several ``C``
    and two ``tol`` for exactly that reason.
    """
    from sklearn.linear_model import LogisticRegression

    records = trial.records
    strata = np.sort(records["stratum"].unique())
    if len(strata) < 2:
        return float("nan")
    design = np.column_stack(
        [(records["stratum"].to_numpy() == g).astype(float) for g in strata[1:]]
    )
    y = records["treated"].to_numpy(dtype=int)
    if y.min() == y.max():
        return float("nan")
    model = LogisticRegression(
        C=float(C), penalty="l2", solver="lbfgs", max_iter=20000, tol=float(tol)
    )
    model.fit(design, y)
    e = model.predict_proba(design)[:, 1]
    if np.any(e <= 0.0) or np.any(e >= 1.0):
        return float("nan")

    outcome = records["outcome"].to_numpy(dtype=float)
    treated = y == 1
    n = len(records)
    weights = np.where(treated, 1.0 / e, 1.0 / (1.0 - e))
    return float(
        (weights[treated] * outcome[treated]).sum() / n
        - (weights[~treated] * outcome[~treated]).sum() / n
    )


#: Every estimator, on one signature: ``(trial, rng) -> float``.
#:
#: ``ipw_cross_fitted`` is called with its default fold split on purpose, so it
#: is the same estimator as in ``trial_recovery.run`` and the two modules'
#: numbers are comparable rather than two draws of one thing.
ESTIMATORS: dict[str, Callable[[Trial, np.random.Generator], float]] = {
    "gcomp-from-generator": lambda t, rng: gcomp_from_generator(t),
    "ipw-saturated": lambda t, rng: ipw_saturated(t),
    "ipw-cross-fitted": lambda t, rng: ipw_cross_fitted(t),
    "ols-stratum-dummies": lambda t, rng: ols_stratum_dummies(t),
    "unadjusted": lambda t, rng: unadjusted(t),
    "aipw-saturated": lambda t, rng: aipw_saturated(t),
    "matching-1": lambda t, rng: matching_exact(t, m=1, rng=rng),
    "matching-20": lambda t, rng: matching_exact(t, m=20, rng=rng),
    "matching-all": lambda t, rng: matching_exact(t, m=None, rng=rng),
    "ipw-logreg-C0.1": lambda t, rng: ipw_logistic_l2(t, C=0.1),
    "ipw-logreg-C1": lambda t, rng: ipw_logistic_l2(t, C=1.0),
    "ipw-logreg-C100": lambda t, rng: ipw_logistic_l2(t, C=100.0),
    "ipw-logreg-C1e6": lambda t, rng: ipw_logistic_l2(t, C=1e6),
    "ipw-logreg-C1e6-tol1e-8": lambda t, rng: ipw_logistic_l2(t, C=1e6, tol=1e-8),
}

#: The estimators cheap enough to carry through the two continuous sweeps. The
#: sweeps are about the design axis, not the estimator zoo; matching and the four
#: logistic fits cost most of a run and would say nothing extra there.
SWEEP_ESTIMATORS: tuple[str, ...] = (
    "gcomp-from-generator",
    "ipw-saturated",
    "aipw-saturated",
    "ols-stratum-dummies",
    "ipw-cross-fitted",
    "unadjusted",
)

#: Residual below this counts as "identically the generator's functional".
#: Floating-point accumulation over a few thousand records lands at 1e-12 or
#: below; genuinely different functionals here are 1e-2 or above. The measured
#: grid is empty between those two, and that emptiness is asserted in the tests
#: rather than assumed -- a threshold sitting inside a continuum would be a
#: knob, not a criterion.
DEGENERATE_BELOW: float = 1e-9

#: THE CORRECTED TABLE. ``trial_recovery.IS_DEGENERATE`` maps an estimator name
#: to a bool, which presumes blindness is a property of the estimator. It is a
#: property of the (estimator, realised-truth, design-point) triple, and the
#: cells below that disagree with that dict are the demonstration:
#:
#:   ols-stratum-dummies   informative against ``standardised`` under the
#:                         confounded design and EXACTLY BLIND against
#:                         ``varweighted`` -- same code, same data, same run.
#:   ols-stratum-dummies   exactly blind against BOTH truths under
#:                         ``block-randomised``, because blocking makes the two
#:                         truths the same number.
#:   gcomp-from-generator  exactly blind against ``standardised`` everywhere and
#:                         informative against ``varweighted`` under confounding
#:                         -- the paper's blindest estimator has an informative
#:                         cell.
#:   unadjusted            exactly blind under ``block-randomised``. This one was
#:                         NOT anticipated and is the sharpest cell in the grid.
#:                         ``trial_recovery`` includes the raw difference in arm
#:                         means as the arm the curve CAN catch -- "included so
#:                         the curve has something it *can* catch". Under exact
#:                         1:1 allocation in every stratum, ``n_g1 = n_g0 =
#:                         n_g/2``, so the overall treated mean is
#:                         ``sum_g (n_g/n) ybar_g1`` and the difference in arm
#:                         means IS the standardised effect, identically. Under
#:                         the randomisation real trials run, the demonstration's
#:                         own positive control goes blind, and the recovery
#:                         curve has nothing left it can catch.
#:   ipw-logreg-*          exactly blind under ``block-randomised`` at every C.
#:                         With ``phat_g = 0.5`` everywhere the saturated MLE is
#:                         the zero vector, which is also the L2 solution, so the
#:                         regularisation constant that governs this estimator's
#:                         entire residual under confounding stops mattering at
#:                         all. The knob is a property of the design, not of the
#:                         estimator either.
#:
#: Measured, never trusted: ``tests/test_trial_blindness.py`` iterates this whole
#: grid against the residual, and so does the grid test in
#: ``tests/test_trial_recovery.py``.
DEGENERACY: dict[tuple[str, str, str], bool] = {}

#: Estimators that are algebraically the standardised effect whenever the
#: propensity model they use is saturated on the stratum.
_STANDARDISATION = (
    "gcomp-from-generator",
    "ipw-saturated",
    "aipw-saturated",
    "matching-all",
)

#: Estimators that stay informative even under exact blocking, and why:
#: cross-fitting because the weights applied to a patient come from data that
#: patient is not in, and finite-M matching because it resamples the opposite arm
#: instead of averaging all of it.
_INFORMATIVE_EVERYWHERE = ("ipw-cross-fitted", "matching-1", "matching-20")

#: Estimators that are informative under any imbalance and exactly blind under
#: exact 1:1 allocation.
_BLIND_ONLY_UNDER_BLOCKING = (
    "unadjusted",
    "ipw-logreg-C0.1",
    "ipw-logreg-C1",
    "ipw-logreg-C100",
    "ipw-logreg-C1e6",
    "ipw-logreg-C1e6-tol1e-8",
)

for _design_name in DESIGNS:
    _blocked = _design_name == "block-randomised"
    for _truth_name in TRUTHS:
        # Under blocking the two truths are the same number, so an estimator
        # blind to either is blind to both.
        _targets_standardisation = _blocked or _truth_name == "standardised"
        for _name in _STANDARDISATION:
            DEGENERACY[(_name, _truth_name, _design_name)] = _targets_standardisation
        DEGENERACY[("ols-stratum-dummies", _truth_name, _design_name)] = (
            _blocked or _truth_name == "varweighted"
        )
        for _name in _BLIND_ONLY_UNDER_BLOCKING:
            DEGENERACY[(_name, _truth_name, _design_name)] = _blocked
        for _name in _INFORMATIVE_EVERYWHERE:
            DEGENERACY[(_name, _truth_name, _design_name)] = False
del _design_name, _truth_name, _blocked, _targets_standardisation, _name


# ---------------------------------------------------------------------------
# The grid
# ---------------------------------------------------------------------------

#: The RNG stream is a function of these keys in this order, everywhere in this
#: module. Written down because a silently reordered key list is the cheapest way
#: to make "fixed seed" mean nothing.
SEED_KEY_ORDER: tuple[str, ...] = (
    "base_seed",
    "arm",
    "seed_index",
    "design_index",
    "n_patients",
    "replicate",
)

_ARM_CODE: dict[str, int] = {"grid": 0, "spread": 1, "tau": 2, "estimator": 3}


def _rng(base_seed: int, arm: str, *keys: int) -> np.random.Generator:
    return np.random.default_rng([int(base_seed), _ARM_CODE[arm], *(int(k) for k in keys)])


def run_grid(
    *,
    seed: int,
    n_seeds: int = 6,
    theta: float = 3.0,
    cohort_sizes: tuple[int, ...] = (100, 200, 500, 1000, 2000, 5000),
    n_replicates: int = 200,
    designs: Sequence[str] = tuple(DESIGNS),
    estimators: Sequence[str] = tuple(ESTIMATORS),
) -> pd.DataFrame:
    """One row per (design, seed, cohort size, replicate, estimator).

    Both truths are columns, not rows, so the residual matrix is a reshape rather
    than a second pass over the generator: the two residuals for one estimate are
    measured on the *same* draw, which is the only way the comparison in (1) is a
    statement about estimands rather than about sampling.

    ``n_seeds`` independent seed streams, because a machine-zero residual at one
    seed is a coincidence and at six is an identity.
    """
    rows: list[dict] = []
    for design_index, design_name in enumerate(designs):
        design = DESIGNS[design_name]
        for seed_index in range(n_seeds):
            for n_patients in cohort_sizes:
                for replicate in range(n_replicates):
                    draw_rng = _rng(
                        seed, "grid", seed_index, design_index, n_patients, replicate
                    )
                    trial = design.draw(n_patients, theta, draw_rng)
                    truths = {k: fn(trial) for k, fn in TRUTHS.items()}
                    if not all(np.isfinite(v) for v in truths.values()):
                        continue
                    for est_index, name in enumerate(estimators):
                        est_rng = _rng(
                            seed,
                            "estimator",
                            seed_index,
                            design_index,
                            n_patients,
                            replicate * 1000 + est_index,
                        )
                        estimate = ESTIMATORS[name](trial, est_rng)
                        row = {
                            "design": design_name,
                            "assignment": design.assignment,
                            "propensity_spread": design.propensity_spread,
                            "tau": design.tau,
                            "seed_index": seed_index,
                            "n_patients": n_patients,
                            "replicate": replicate,
                            "estimator": name,
                            "theta_requested": trial.theta_requested,
                            "estimate": estimate,
                        }
                        for truth_name, truth_value in truths.items():
                            row[f"truth_{truth_name}"] = truth_value
                            row[f"residual_vs_{truth_name}"] = abs(estimate - truth_value)
                        rows.append(row)
    return pd.DataFrame(rows)


def _long_residuals(runs: pd.DataFrame) -> pd.DataFrame:
    """Melt the two truth columns into a ``truth`` column."""
    parts = []
    for truth_name in TRUTHS:
        part = runs[
            ["design", "seed_index", "n_patients", "replicate", "estimator",
             "theta_requested", "estimate", f"truth_{truth_name}",
             f"residual_vs_{truth_name}"]
        ].copy()
        part = part.rename(
            columns={
                f"truth_{truth_name}": "theta_realised",
                f"residual_vs_{truth_name}": "residual",
            }
        )
        part["truth"] = truth_name
        parts.append(part)
    return pd.concat(parts, ignore_index=True)


def residual_matrix(runs: pd.DataFrame) -> pd.DataFrame:
    """(estimator x truth x design x cohort size) -> residual, and the verdict.

    ``max_residual`` is the column that matters. A median can hide one replicate
    where the identity failed; an identity that fails once is not an identity.
    """
    long = _long_residuals(runs)
    out = (
        long.groupby(["estimator", "truth", "design", "n_patients"])
        .agg(
            n_replicates=("replicate", "count"),
            max_residual=("residual", "max"),
            median_residual=("residual", "median"),
            median_recovery_ratio=("estimate", "median"),
        )
        .reset_index()
    )
    out["median_recovery_ratio"] = out["median_recovery_ratio"] / float(
        long["theta_requested"].iloc[0]
    )
    out["measured_degenerate"] = out["max_residual"] < DEGENERATE_BELOW
    out["declared_degenerate"] = [
        DEGENERACY[(e, t, d)]
        for e, t, d in zip(out["estimator"], out["truth"], out["design"], strict=True)
    ]
    return out.sort_values(["design", "estimator", "truth", "n_patients"]).reset_index(
        drop=True
    )


def information_ratio(runs: pd.DataFrame) -> pd.DataFrame:
    """rho = median|theta_hat - theta_realised| / median|theta_realised - theta|.

    The residual in units of the generator's own sampling error. rho near 1 means
    the estimator's departure from the realised truth is the same size as the
    generator's departure from the requested one, so the recovery curve is
    roughly half estimator and half generator; rho near 0 means the curve is
    almost entirely the generator, whatever the residual's absolute size says.

    ``generator_noise_share = 1 / (1 + rho)`` is a reading aid, not a variance
    decomposition: the two medians are not orthogonal components and do not add.
    It is reported because "rho = 0.13" and "the curve is mostly the generator"
    are the same sentence and only one of them is legible.
    """
    long = _long_residuals(runs)
    long = long.assign(
        generator_error=(long["theta_realised"] - long["theta_requested"]).abs()
    )
    grouped = long.groupby(["estimator", "truth", "design", "n_patients"])
    out = grouped.agg(
        n_replicates=("replicate", "count"),
        median_residual=("residual", "median"),
        median_generator_error=("generator_error", "median"),
    ).reset_index()
    out["rho"] = out["median_residual"] / out["median_generator_error"]
    out["generator_noise_share"] = 1.0 / (1.0 + out["rho"])

    per_seed = (
        long.groupby(["estimator", "truth", "design", "n_patients", "seed_index"])
        .agg(
            median_residual=("residual", "median"),
            median_generator_error=("generator_error", "median"),
        )
        .reset_index()
    )
    per_seed["rho"] = per_seed["median_residual"] / per_seed["median_generator_error"]
    spread = (
        per_seed.groupby(["estimator", "truth", "design", "n_patients"])["rho"]
        .agg(rho_min_over_seeds="min", rho_max_over_seeds="max")
        .reset_index()
    )
    out = out.merge(spread, on=["estimator", "truth", "design", "n_patients"], how="left")
    return out.sort_values(["design", "estimator", "truth", "n_patients"]).reset_index(
        drop=True
    )


def run_sweep(
    *,
    seed: int,
    arm: str,
    values: Sequence[float],
    n_seeds: int = 6,
    theta: float = 3.0,
    n_patients: int = 2000,
    n_replicates: int = 200,
    estimators: Sequence[str] = SWEEP_ESTIMATORS,
) -> pd.DataFrame:
    """One continuous axis: ``arm="spread"`` or ``arm="tau"``.

    ``spread`` walks the propensities from 1:1 to the paper's confounded design;
    ``tau`` walks the effect from homogeneous to strongly stratum-dependent. Both
    return the same columns as ``run_grid`` so the summarisers are shared.
    """
    if arm not in ("spread", "tau"):
        raise ValueError(f"arm must be 'spread' or 'tau', got {arm!r}")
    rows: list[dict] = []
    for value_index, value in enumerate(values):
        design = (
            Design(f"spread={value:g}", propensity_spread=float(value))
            if arm == "spread"
            else Design(f"tau={value:g}", tau=float(value))
        )
        for seed_index in range(n_seeds):
            for replicate in range(n_replicates):
                draw_rng = _rng(seed, arm, seed_index, value_index, n_patients, replicate)
                trial = design.draw(n_patients, theta, draw_rng)
                truths = {k: fn(trial) for k, fn in TRUTHS.items()}
                if not all(np.isfinite(v) for v in truths.values()):
                    continue
                for est_index, name in enumerate(estimators):
                    est_rng = _rng(
                        seed,
                        "estimator",
                        seed_index,
                        value_index,
                        n_patients,
                        replicate * 1000 + est_index,
                    )
                    estimate = ESTIMATORS[name](trial, est_rng)
                    row = {
                        "design": design.name,
                        "sweep_arm": arm,
                        "sweep_value": float(value),
                        "seed_index": seed_index,
                        "n_patients": n_patients,
                        "replicate": replicate,
                        "estimator": name,
                        "theta_requested": trial.theta_requested,
                        "estimate": estimate,
                    }
                    for truth_name, truth_value in truths.items():
                        row[f"truth_{truth_name}"] = truth_value
                        row[f"residual_vs_{truth_name}"] = abs(estimate - truth_value)
                    rows.append(row)
    return pd.DataFrame(rows)


def summarise_sweep(runs: pd.DataFrame) -> pd.DataFrame:
    """Per (sweep value, estimator, truth): the residual and the curve.

    ``median_recovery_ratio`` is carried alongside so the heterogeneity table can
    show the two moving in opposite directions -- the curve wandering off 1 while
    the residual against the variance-weighted truth stays at machine zero.
    """
    long = _long_residuals(runs.assign(design=runs["design"]))
    long = long.merge(
        runs[["design", "sweep_arm", "sweep_value"]].drop_duplicates(),
        on="design",
        how="left",
    )
    out = (
        long.groupby(["sweep_arm", "sweep_value", "estimator", "truth"])
        .agg(
            n_replicates=("replicate", "count"),
            max_residual=("residual", "max"),
            median_residual=("residual", "median"),
            median_estimate=("estimate", "median"),
            median_truth=("theta_realised", "median"),
        )
        .reset_index()
    )
    out["median_recovery_ratio"] = out["median_estimate"] / float(
        runs["theta_requested"].iloc[0]
    )
    out["measured_degenerate"] = out["max_residual"] < DEGENERATE_BELOW
    return out.sort_values(["sweep_arm", "estimator", "truth", "sweep_value"]).reset_index(
        drop=True
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

#: Copied verbatim from ``trial_recovery.main``. Same generator, same disclaimer;
#: paraphrasing it would let the two tables drift apart in what they claim.
SYNTHETIC_NOTE = (
    "Simulated patients with an analytically known treatment effect. Nothing "
    "here is a result about any real trial."
)

SPREAD_VALUES: tuple[float, ...] = (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
#: The drift in OLS's recovery curve is linear in tau and MODEST: the
#: variance weights w_g proportional to n_g p_g (1 - p_g) are not far from
#: the prevalence weights n_g under this design, so tau=4 moves the curve
#: only about 3.6%. The grid runs out to 16 so the effect is larger than the
#: replicate noise and the direction is unambiguous, not because tau=16 is a
#: plausible trial.
TAU_VALUES: tuple[float, ...] = (0.0, 1.0, 2.0, 4.0, 8.0, 16.0)


def main(argv: Sequence[str] | None = None) -> int:
    """Write the four blindness tables under ``results/`` with a provenance stamp.

        python -m src.harness.trial_blindness
    """
    parser = argparse.ArgumentParser(
        description="Degeneracy as a property of (estimator, truth, design)"
    )
    parser.add_argument("--seed", type=int, default=20260909)
    parser.add_argument("--replicates", type=int, default=200)
    parser.add_argument("--seeds", type=int, default=6, help="independent seed streams")
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--allow-dirty", action="store_true")
    parser.add_argument(
        "--quick",
        action="store_true",
        help="a small grid for smoke-testing the wiring; not a result",
    )
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(message)s")

    cohort_sizes = (200, 2000) if args.quick else (100, 200, 500, 1000, 2000, 5000)
    replicates = 5 if args.quick else args.replicates
    n_seeds = 2 if args.quick else args.seeds
    spread_values = (0.0, 1.0) if args.quick else SPREAD_VALUES
    tau_values = (0.0, 16.0) if args.quick else TAU_VALUES

    common = {
        "seed": args.seed,
        "results_dir": args.results_dir,
        "allow_dirty": args.allow_dirty,
    }
    shared_meta = {
        "sweep_seed": args.seed,
        "n_seed_streams": n_seeds,
        "n_replicates_per_seed": replicates,
        "theta_requested": 3.0,
        "stratum_means": list(STRATUM_MEANS),
        "propensity": list(PROPENSITY),
        "stratum_weights": list(STRATUM_WEIGHTS),
        "outcome_sd": OUTCOME_SD,
        "degenerate_below": DEGENERATE_BELOW,
        "seed_key_order": list(SEED_KEY_ORDER),
        "SYNTHETIC": SYNTHETIC_NOTE,
    }

    log.info("grid: %d designs x %d seeds x %d sizes x %d reps x %d estimators",
             len(DESIGNS), n_seeds, len(cohort_sizes), replicates, len(ESTIMATORS))
    runs = run_grid(
        seed=args.seed,
        n_seeds=n_seeds,
        cohort_sizes=cohort_sizes,
        n_replicates=replicates,
    )

    matrix = residual_matrix(runs)
    disagreements = matrix[matrix["measured_degenerate"] != matrix["declared_degenerate"]]
    if len(disagreements):
        log.warning(
            "DECLARED TABLE DISAGREES WITH MEASUREMENT in %d cells:\n%s",
            len(disagreements),
            disagreements[
                ["estimator", "truth", "design", "n_patients", "max_residual",
                 "declared_degenerate"]
            ].to_string(index=False),
        )

    path = write_versioned_table(
        matrix, "trial_blindness_residual_matrix",
        extra_meta=shared_meta | {
            "cohort_sizes": list(cohort_sizes),
            "designs": {k: v.as_meta() for k, v in DESIGNS.items()},
            "truths": {
                "standardised": "sum_g (n_g/n) (ybar_g1 - ybar_g0) -- "
                                "trial_recovery.theta_realised",
                "varweighted": "sum_g w_g (ybar_g1 - ybar_g0) / sum_g w_g, "
                               "w_g = n_g phat_g (1 - phat_g) -- OLS's own "
                               "finite-sample estimand, exact by Frisch-Waugh",
            },
            "estimators": list(ESTIMATORS),
            "declared_vs_measured_disagreements": int(len(disagreements)),
            "what_this_answers": (
                "Whether 'degenerate' is a property of an estimator. It is not: "
                "the same estimator on the same draw is exactly blind against one "
                "realised-truth functional and informative against another, and "
                "OLS flips from informative to exactly blind when the trial is "
                "block-randomised instead of confounded. The honest object is "
                "this matrix, not a dict keyed by estimator name."
            ),
        },
        **common,
    )
    log.info("wrote %s (%d rows)", path, len(matrix))

    rho = information_ratio(runs)
    path = write_versioned_table(
        rho, "trial_blindness_information_ratio",
        extra_meta=shared_meta | {
            "cohort_sizes": list(cohort_sizes),
            "rho_definition": (
                "median|theta_hat - theta_realised| / "
                "median|theta_realised - theta_requested|"
            ),
            "what_this_answers": (
                "How much of a recovery curve is the estimator, in units of the "
                "generator's own sampling error. A non-zero residual is not "
                "evidence that a curve is informative: OLS's rho sits near 0.13 "
                "at every cohort size and does not fall with n, so its curve is "
                "mostly generator noise even where it is not identically zero."
            ),
        },
        **common,
    )
    log.info("wrote %s (%d rows)", path, len(rho))

    spread = summarise_sweep(
        run_sweep(seed=args.seed, arm="spread", values=spread_values,
                  n_seeds=n_seeds, n_replicates=replicates)
    )
    path = write_versioned_table(
        spread, "trial_blindness_propensity_spread",
        extra_meta=shared_meta | {
            "sweep_values": list(spread_values),
            "sweep_n_patients": 2000,
            "sweep_estimators": list(SWEEP_ESTIMATORS),
            "what_this_answers": (
                "Where OLS stops being informative. Its residual against the "
                "standardised truth falls continuously to machine zero as the "
                "propensities collapse to 1:1, so 'informative' is a statement "
                "about the design point and not about the estimator."
            ),
        },
        **common,
    )
    log.info("wrote %s (%d rows)", path, len(spread))

    tau = summarise_sweep(
        run_sweep(seed=args.seed, arm="tau", values=tau_values,
                  n_seeds=n_seeds, n_replicates=replicates)
    )
    path = write_versioned_table(
        tau, "trial_blindness_heterogeneity",
        extra_meta=shared_meta | {
            "sweep_values": list(tau_values),
            "sweep_n_patients": 2000,
            "sweep_estimators": list(SWEEP_ESTIMATORS),
            "theta_g": "theta + tau * (g - gbar), gbar the design-weighted mean "
                       "stratum index, so the population ATE stays theta",
            "what_this_answers": (
                "Whether a recovery curve moving off 1 means the estimator got "
                "worse. Under effect heterogeneity OLS's curve wanders while its "
                "residual against the estimand it actually targets stays at "
                "machine zero: the curve moved because the estimand moved."
            ),
        },
        **common,
    )
    log.info("wrote %s (%d rows)", path, len(tau))
    return 0


if __name__ == "__main__":
    sys.exit(main())
