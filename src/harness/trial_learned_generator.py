"""The residual audit when the generator is LEARNED rather than written down.

``paper/wmhs/sections/full/refdesign.tex`` closes with a prediction the paper
flags as unrun: that replacing the parametric patient simulator with a learned
generative model degrades the equality ``theta_hat(D) == T(D)`` from exact to
APPROXIMATE, so that ``np.allclose`` at default tolerance returns True and sees
nothing while the residual scale ``rho`` still does. This module runs it.

The design is pre-registered in ``docs/prereg_learned_generator.md``, committed
before this file existed, together with seven falsifiers and one amendment.
Read that document first: it fixes what would show the prediction is wrong, and
nothing here is allowed to weaken it.

THE SHAPE. A *training cohort* is drawn from ``trial_recovery.simulate_trial``,
imported unchanged -- covariate strata, confounded assignment, a known requested
effect ``theta``. A generator is fitted to that cohort and asked for a
*synthetic cohort* of the same size. The estimator is re-fit on the synthetic
cohort and never sees the training one. That two-stage wiring -- real draw, fit,
synthetic draw, estimate -- is what a virtual control arm is.

THE OBJECTION THIS MODULE EXISTS TO MEASURE. The paper's stated mechanism for
its prediction is that "a learned generator and a re-fit outcome model share
most but not all of their inputs". That describes the relationship between the
fitted generator and the outcome model -- but the screen does not compare those
two objects. It compares ``T(D~)`` and ``theta_hat(D~)``, BOTH OF WHICH ARE
FUNCTIONALS OF THE SAME SYNTHETIC COHORT ``D~``. Whether they coincide is a
property of the estimator's functional form. The procedure that produced ``D~``
cannot enter it. So the prediction should fail as stated, and the parametric
control sits in the same table at the same cohort sizes so that a reader can
check that rather than take it on assertion.

WHERE APPROXIMATE EQUALITY ACTUALLY COMES FROM, AND IT IS NOT THE GENERATOR.
Two arms here are degenerate in intent:

``gcomp-gmm-outcome-K*``   G-computation whose outcome model is a per-cell
                           Gaussian mixture. At ANY EM fixed point -- indeed
                           after any M-step -- ``sum_k pi_k mu_k = ybar``
                           exactly, because responsibilities sum to one. So a
                           mixture-based outcome model IS the cell mean, to
                           machine precision, and ``tol`` is not a floor. The
                           pre-registration predicted a tol-set floor here and
                           was WRONG; Amendment 1 records that.

``gcomp-ridge-outcome-a*`` G-computation whose outcome model is the same
                           saturated cell-indicator design fitted by penalised
                           least squares. On an orthogonal binary design the
                           solution is ``ybar_gd * n_gd / (n_gd + alpha)``, so
                           the residual is a smooth CHOSEN function of a
                           regularisation constant. These arms exist as the
                           forcing input for the screen's blind band: they can
                           be dialled into the gap between the two things a
                           practitioner might type, which proves the gap is
                           reachable and says what reaches it.

THE TWO THINGS A PRACTITIONER MIGHT TYPE, AND THEY DISAGREE BY 10^3.

    np.allclose(residual, 0.0)        -> |r| <= atol + rtol*|0| = 1e-8
    np.allclose(estimate, reference)  -> |r| <= atol + rtol*|T| ~ 3e-5 at theta=3

The ``rtol`` term multiplies the second argument. Screening the residual
against zero therefore throws the relative tolerance away and is three orders
of magnitude stricter than screening the estimate against the reference. Both
are reported, side by side, in the primary table, because which line gets typed
decides what the screen can see.

NOTHING HERE IS A RESULT ABOUT ANY REAL TRIAL.
"""

from __future__ import annotations

import argparse
import logging
import sys
import warnings
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from src.common.io import write_versioned_table
from src.harness.trial_blindness import DEGENERATE_BELOW
from src.harness.trial_recovery import (
    OUTCOME_SD,
    PROPENSITY,
    STRATUM_MEANS,
    STRATUM_WEIGHTS,
    Trial,
    _standardised_effect,
    gcomp_from_generator,
    ipw_saturated,
    ols_stratum_dummies,
    simulate_trial,
    unadjusted,
)

log = logging.getLogger(__name__)

__all__ = [
    "ALLCLOSE_ATOL",
    "ALLCLOSE_RTOL",
    "ESTIMATORS",
    "GENERATORS",
    "GeneratorFitError",
    "REFERENCES",
    "SATURATED_BELOW",
    "capacity_table",
    "flip_tolerances",
    "primary_table",
    "run_grid",
    "screen_band",
    "trend_table",
]

#: The requested effect. Matched to ``trial_recovery`` and ``trial_blindness``
#: so the parametric control in this table is comparable with the published one.
THETA: float = 3.0

#: numpy's ``allclose`` defaults, written down because the entire (b) question
#: is what they do. ``np.allclose(a, b)`` tests ``|a-b| <= atol + rtol*|b|``,
#: so the second argument is the one ``rtol`` scales.
ALLCLOSE_RTOL: float = 1e-5
ALLCLOSE_ATOL: float = 1e-8

#: A generator whose implied cell means reproduce the TRAINING cohort's cell
#: means below this is the parametric generator in disguise (prereg 2.3). It is
#: reported as a column, never used as a filter.
SATURATED_BELOW: float = 1e-9

#: Draw size for the Monte-Carlo estimate of ``T_model`` where no closed form
#: survives the generator's discretisation step. Its standard error is carried
#: beside it; for the mixture generators it lands near 2e-2, which is the
#: quantitative content of the paper's claim (c).
ORACLE_DRAW: int = 200_000

N_STRATA: int = len(STRATUM_WEIGHTS)
CELLS: tuple[tuple[int, int], ...] = tuple(
    (g, a) for g in range(N_STRATA) for a in (0, 1)
)


class GeneratorFitError(RuntimeError):
    """The training cohort cannot support this generator.

    A stratum-arm cell with no records leaves a plug-in mean undefined. That is
    a positivity failure, and it is counted and reported rather than patched
    with a fallback that would quietly change what the generator is.
    """


# ---------------------------------------------------------------------------
# Seeding
# ---------------------------------------------------------------------------

#: The RNG stream is a function of these keys in this order, everywhere in this
#: module. Written down for the same reason ``trial_blindness`` writes its own
#: down: a silently reordered key list is the cheapest way to make "fixed seed"
#: mean nothing.
SEED_KEY_ORDER: tuple[str, ...] = (
    "seed",
    "arm_code",
    "seed_index",
    "generator_index",
    "n_patients",
    "replicate",
    "estimator_index",
)

_ARM_CODE: dict[str, int] = {
    "train": 0,
    "fit": 1,
    "synth": 2,
    "oracle": 3,
    "estimator": 4,
}


def _rng(base_seed: int, arm: str, *keys: int) -> np.random.Generator:
    return np.random.default_rng(
        [int(base_seed), _ARM_CODE[arm], *(int(k) for k in keys)]
    )


def _seed_int(rng: np.random.Generator) -> int:
    """A reproducible ``random_state`` for sklearn, drawn from our own stream."""
    return int(rng.integers(0, 2**31 - 1))


# ---------------------------------------------------------------------------
# Cohort helpers
# ---------------------------------------------------------------------------


def _records(
    stratum: np.ndarray, treated: np.ndarray, outcome: np.ndarray
) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "stratum": np.asarray(stratum, dtype=int),
            "treated": np.asarray(treated, dtype=int),
            "outcome": np.asarray(outcome, dtype=float),
        }
    )


def cell_means(records: pd.DataFrame) -> dict[tuple[int, int], float]:
    """``ybar_gd`` for every observed cell. Missing cells are simply absent."""
    out: dict[tuple[int, int], float] = {}
    for (g, a), group in records.groupby(["stratum", "treated"]):
        out[(int(g), int(a))] = float(group["outcome"].mean())
    return out


def _as_trial(records: pd.DataFrame, theta: float) -> Trial:
    """Wrap a synthetic cohort so every imported estimator keeps its meaning.

    ``theta_realised`` is ``T_draw``: the standardised effect of the synthetic
    cohort, computed by the SAME function the parametric arm uses. Imported, not
    re-implemented, so the two tables cannot drift apart in what they mean.
    """
    return Trial(
        records=records,
        theta_requested=float(theta),
        theta_realised=_standardised_effect(records),
    )


def _onehot(stratum: np.ndarray) -> np.ndarray:
    return np.column_stack(
        [(np.asarray(stratum) == g).astype(float) for g in range(N_STRATA)]
    )


def _standardise_cell_means(
    records: pd.DataFrame, mu: dict[tuple[int, int], float]
) -> float:
    """``sum_g (n_g/n) (mu_g1 - mu_g0)`` using the cohort's REALISED weights.

    Every G-computation arm in this module differs only in how ``mu`` was
    obtained, so the standardisation step is written once. That is deliberate:
    it means a non-zero residual cannot be an artefact of two arms weighting
    their strata differently.
    """
    n = len(records)
    total = 0.0
    for g, group in records.groupby("stratum"):
        g = int(g)
        if (g, 1) not in mu or (g, 0) not in mu:
            return float("nan")
        total += (len(group) / n) * (mu[(g, 1)] - mu[(g, 0)])
    return float(total)


# ---------------------------------------------------------------------------
# The generators
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class FittedGenerator:
    """A generator that has seen its training cohort and can be sampled.

    ``implied_cell_means`` and ``implied_effect`` are the generator's own
    population quantities, which is what the saturation check (prereg 2.3) and
    the secondary reference ``T_model`` (prereg 3.2) are about. For a learned
    generator neither is free: where the sampling step has no closed form they
    are Monte-Carlo estimates carrying ``implied_effect_mc_se``, and that cost
    is the paper's claim (c) made into a number.
    """

    name: str
    is_learned: bool
    sample: Callable[[int, np.random.Generator], pd.DataFrame]
    implied_cell_means: dict[tuple[int, int], float]
    implied_effect: float
    implied_effect_mc_se: float
    implied_is_monte_carlo: bool


def _gaussian_mixture(n_components: int, seed: int, **kwargs):
    from sklearn.mixture import GaussianMixture

    params = dict(
        n_components=n_components,
        covariance_type="full",
        reg_covar=1e-6,
        tol=1e-3,
        max_iter=200,
        n_init=1,
        init_params="kmeans",
        random_state=seed,
    )
    params.update(kwargs)
    return GaussianMixture(**params)


def _mlp(seed: int, hidden: tuple[int, ...] = (16,)):
    from sklearn.neural_network import MLPRegressor

    return MLPRegressor(
        hidden_layer_sizes=hidden,
        solver="lbfgs",
        max_iter=300,
        random_state=seed,
    )


def fit_parametric_oracle(
    training: pd.DataFrame, *, n_patients: int, theta: float, rng: np.random.Generator
) -> FittedGenerator:
    """THE CONTROL. Not learned, not fitted: it re-draws from the true generator.

    Included so that "exact" and "approximate" are visible in one table rather
    than asserted across two documents. Its population standardised effect is
    ``theta`` in closed form, because the stratum means cancel.
    """
    del training, n_patients, rng

    def sample(n: int, rng_: np.random.Generator) -> pd.DataFrame:
        return simulate_trial(n, theta, rng=rng_).records

    mu = {(g, a): float(STRATUM_MEANS[g]) + float(theta) * a for g, a in CELLS}
    return FittedGenerator(
        name="parametric-oracle",
        is_learned=False,
        sample=sample,
        implied_cell_means=mu,
        implied_effect=float(theta),
        implied_effect_mc_se=0.0,
        implied_is_monte_carlo=False,
    )


def fit_parametric_plugin(
    training: pd.DataFrame, *, n_patients: int, theta: float, rng: np.random.Generator
) -> FittedGenerator:
    """MLE plug-in. SATURATED BY CONSTRUCTION, and included to calibrate that.

    Its implied cell means ARE the training cohort's cell means, so the
    saturation column reads exactly zero on it. A learned generator whose
    column reads the same is this one wearing different vocabulary.
    """
    del theta
    n = len(training)
    counts = training["stratum"].value_counts()
    pi = np.array([counts.get(g, 0) / n for g in range(N_STRATA)], dtype=float)

    mu: dict[tuple[int, int], float] = {}
    sd: dict[tuple[int, int], float] = {}
    p_treat = np.zeros(N_STRATA)
    for g in range(N_STRATA):
        block = training[training["stratum"] == g]
        if block.empty:
            raise GeneratorFitError(f"stratum {g} empty in the training cohort")
        p_treat[g] = float(block["treated"].mean())
        for a in (0, 1):
            y = block.loc[block["treated"] == a, "outcome"].to_numpy()
            if y.size == 0:
                raise GeneratorFitError(f"cell ({g},{a}) empty in the training cohort")
            mu[(g, a)] = float(y.mean())
            sd[(g, a)] = float(y.std(ddof=1)) if y.size > 1 else float(OUTCOME_SD)

    def sample(n_draw: int, rng_: np.random.Generator) -> pd.DataFrame:
        stratum = rng_.choice(N_STRATA, size=n_draw, p=pi / pi.sum())
        treated = (rng_.random(n_draw) < p_treat[stratum]).astype(int)
        means = np.array([mu[(int(g), int(a))] for g, a in zip(stratum, treated)])
        sds = np.array([sd[(int(g), int(a))] for g, a in zip(stratum, treated)])
        return _records(stratum, treated, rng_.normal(means, sds))

    implied = float(sum(pi[g] * (mu[(g, 1)] - mu[(g, 0)]) for g in range(N_STRATA)))
    return FittedGenerator(
        name="parametric-plugin",
        is_learned=False,
        sample=sample,
        implied_cell_means=mu,
        implied_effect=implied,
        implied_effect_mc_se=0.0,
        implied_is_monte_carlo=False,
    )


def _mc_summary(
    records: pd.DataFrame,
) -> tuple[dict[tuple[int, int], float], float, float]:
    """Cell means, standardised effect and its MC standard error, from a draw."""
    mu = cell_means(records)
    effect = _standardised_effect(records)
    n = len(records)
    var = 0.0
    for g, group in records.groupby("stratum"):
        w = len(group) / n
        for a in (0, 1):
            y = group.loc[group["treated"] == a, "outcome"].to_numpy()
            if y.size < 2:
                return mu, float(effect), float("nan")
            var += (w**2) * float(y.var(ddof=1)) / y.size
    return mu, float(effect), float(np.sqrt(var))


def _fit_gmm_joint(
    training: pd.DataFrame,
    *,
    n_components: int,
    rng: np.random.Generator,
) -> FittedGenerator:
    """A LEARNED joint density over ``[onehot(stratum), treated, outcome]``.

    Fitted with ``GaussianMixture`` and sampled with ``.sample()``, the sampled
    rows discretised back: ``stratum = argmax`` over the one-hot columns,
    ``treated = 1`` where that column exceeds 0.5. This is the standard "fit a
    mixture to a table and sample rows" recipe, and its crudeness is part of
    what is being measured.

    Sampled rows are SHUFFLED. ``GaussianMixture.sample`` returns them grouped
    by component, and a cohort whose row order encodes the latent component is
    an artefact of the library rather than a property of the generator.

    Nothing about the discretisation survives into closed form, so this
    generator's implied cell means and implied effect are Monte-Carlo, with the
    standard error reported beside them.
    """
    X = np.column_stack(
        [
            _onehot(training["stratum"].to_numpy()),
            training["treated"].to_numpy(dtype=float),
            training["outcome"].to_numpy(dtype=float),
        ]
    )
    if len(X) < n_components:
        raise GeneratorFitError(
            f"{len(X)} training rows cannot fit {n_components} mixture components"
        )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = _gaussian_mixture(n_components, _seed_int(rng)).fit(X)

    def draw(n_draw: int, rng_: np.random.Generator) -> pd.DataFrame:
        if n_draw < n_components:
            raise GeneratorFitError("cannot sample fewer rows than components")
        model.random_state = _seed_int(rng_)
        z, _ = model.sample(n_draw)
        z = z[rng_.permutation(len(z))]
        stratum = np.argmax(z[:, :N_STRATA], axis=1)
        treated = (z[:, N_STRATA] > 0.5).astype(int)
        return _records(stratum, treated, z[:, N_STRATA + 1])

    oracle = draw(ORACLE_DRAW, _rng(0, "oracle", _seed_int(rng)))
    mu, effect, se = _mc_summary(oracle)
    return FittedGenerator(
        name=f"gmm-joint-K{n_components}",
        is_learned=True,
        sample=draw,
        implied_cell_means=mu,
        implied_effect=effect,
        implied_effect_mc_se=se,
        implied_is_monte_carlo=True,
    )


def fit_mlp_conditional_gaussian(
    training: pd.DataFrame, *, n_patients: int, theta: float, rng: np.random.Generator
) -> FittedGenerator:
    """A LEARNED conditional mean with Gaussian noise.

    ``E[Y|g,a]`` from an ``MLPRegressor`` on a saturated design, residual SD
    from the fit, and ``(stratum, treated)`` resampled from the training
    cohort's empirical joint. The outcome is standardised before fitting and
    un-standardised after: an unscaled target of magnitude 25 is a scaling
    accident rather than a modelling choice, and leaving it in would make this
    arm look bad for a reason that has nothing to do with the paper.

    Its implied quantities ARE closed form -- the fitted network evaluated at
    the six cells -- which is worth noticing: "learned" does not by itself mean
    "no closed form for T_model". That depends on the sampling step, not on
    whether a network was involved.
    """
    del n_patients, theta
    stratum = training["stratum"].to_numpy()
    treated = training["treated"].to_numpy(dtype=float)
    onehot = _onehot(stratum)
    design = np.column_stack([onehot, treated, onehot * treated[:, None]])
    y = training["outcome"].to_numpy(dtype=float)
    centre, scale = float(y.mean()), float(y.std(ddof=1))
    if not np.isfinite(scale) or scale == 0.0:
        raise GeneratorFitError("training outcomes have no spread")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = _mlp(_seed_int(rng)).fit(design, (y - centre) / scale)

    grid = np.array([[g, a] for g, a in CELLS], dtype=float)
    grid_onehot = _onehot(grid[:, 0].astype(int))
    grid_design = np.column_stack(
        [grid_onehot, grid[:, 1], grid_onehot * grid[:, 1][:, None]]
    )
    predicted = model.predict(grid_design) * scale + centre
    mu = {cell: float(v) for cell, v in zip(CELLS, predicted)}

    fitted = model.predict(design) * scale + centre
    resid_sd = float(np.std(y - fitted, ddof=1))
    pairs = np.column_stack([stratum, treated.astype(int)])
    counts = training["stratum"].value_counts()
    pi = np.array([counts.get(g, 0) / len(training) for g in range(N_STRATA)])

    def sample(n_draw: int, rng_: np.random.Generator) -> pd.DataFrame:
        idx = rng_.integers(0, len(pairs), size=n_draw)
        picked = pairs[idx]
        means = np.array([mu[(int(g), int(a))] for g, a in picked])
        return _records(picked[:, 0], picked[:, 1], rng_.normal(means, resid_sd))

    implied = float(sum(pi[g] * (mu[(g, 1)] - mu[(g, 0)]) for g in range(N_STRATA)))
    return FittedGenerator(
        name="mlp-conditional-gaussian",
        is_learned=True,
        sample=sample,
        implied_cell_means=mu,
        implied_effect=implied,
        implied_effect_mc_se=0.0,
        implied_is_monte_carlo=False,
    )


def _gmm_factory(k: int):
    def fit(training, *, n_patients, theta, rng):
        del n_patients, theta
        return _fit_gmm_joint(training, n_components=k, rng=rng)

    return fit


#: Every generator, on one signature. The two parametric arms are the control
#: and sit in the same table at the same cohort sizes and seeds.
GENERATORS: dict[str, Callable[..., FittedGenerator]] = {
    "parametric-oracle": fit_parametric_oracle,
    "parametric-plugin": fit_parametric_plugin,
    "gmm-joint-K1": _gmm_factory(1),
    "gmm-joint-K2": _gmm_factory(2),
    "gmm-joint-K4": _gmm_factory(4),
    "gmm-joint-K8": _gmm_factory(8),
    "mlp-conditional-gaussian": fit_mlp_conditional_gaussian,
}

#: Whether the generator is learned from the training cohort. ``parametric-*``
#: are not; they are the contrast the paper's prediction is stated against.
IS_LEARNED: dict[str, bool] = {
    "parametric-oracle": False,
    "parametric-plugin": False,
    "gmm-joint-K1": True,
    "gmm-joint-K2": True,
    "gmm-joint-K4": True,
    "gmm-joint-K8": True,
    "mlp-conditional-gaussian": True,
}

#: Mixture components, for the capacity axis. ``None`` where the generator has
#: no such knob.
CAPACITY: dict[str, float] = {
    "parametric-oracle": float("nan"),
    "parametric-plugin": float("nan"),
    "gmm-joint-K1": 1.0,
    "gmm-joint-K2": 2.0,
    "gmm-joint-K4": 4.0,
    "gmm-joint-K8": 8.0,
    "mlp-conditional-gaussian": float("nan"),
}


# ---------------------------------------------------------------------------
# The estimators, all re-fit on the synthetic cohort
# ---------------------------------------------------------------------------


def gcomp_saturated(trial: Trial, rng: np.random.Generator) -> float:
    """THE NATURAL WIRING: compute the effect on the cohort you just generated.

    Imported ``gcomp_from_generator`` unchanged, which evaluates
    ``_standardised_effect`` on the synthetic records -- the same function, on
    the same input, that produced ``T_draw``. The residual is therefore bitwise
    zero by construction, for a learned generator exactly as for a parametric
    one, and that is the point: the generator's learnedness never entered.
    """
    del rng
    return gcomp_from_generator(trial)


def gcomp_gmm_outcome(
    trial: Trial, *, n_components: int, rng: np.random.Generator
) -> float:
    """G-computation whose outcome model is a per-cell Gaussian mixture.

    DEGENERATE, and the pre-registration expected it to be only approximately
    so. It was wrong. At any EM fixed point -- after any M-step --
    ``sum_k pi_k mu_k = ybar`` exactly: the M-step sets
    ``mu_k = sum_i r_ik y_i / sum_i r_ik`` and ``pi_k = sum_i r_ik / n``, so the
    mixture mean telescopes to ``(1/n) sum_i y_i sum_k r_ik = ybar`` because
    responsibilities sum to one. Convergence is irrelevant, so ``tol`` sets no
    floor and this arm sits at machine precision.

    Returns NaN when a cell holds fewer records than components, which is a
    fitting impossibility rather than a statistical statement; those draws are
    counted in ``n_nonfinite_estimate``.
    """
    records = trial.records
    mu: dict[tuple[int, int], float] = {}
    for (g, a), group in records.groupby(["stratum", "treated"]):
        y = group["outcome"].to_numpy(dtype=float)
        if y.size < n_components:
            return float("nan")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model = _gaussian_mixture(n_components, _seed_int(rng)).fit(
                y.reshape(-1, 1)
            )
        mu[(int(g), int(a))] = float(model.weights_ @ model.means_.ravel())
    return _standardise_cell_means(records, mu)


def gcomp_ridge_outcome(
    trial: Trial, *, alpha: float, rng: np.random.Generator
) -> float:
    """G-computation whose outcome model is the saturated design, PENALISED.

    THE FORCING INPUT FOR THE SCREEN'S BLIND BAND (Amendment 1). One indicator
    per ``(stratum, treated)`` cell and no intercept, so the design is
    orthogonal and ridge returns ``ybar_gd * n_gd / (n_gd + alpha)`` in closed
    form. The residual against ``T_draw`` is then a smooth function of
    ``alpha``, dialled by the analyst, and ``alpha`` in 1e-5..1e-1 brackets the
    gap between ``allclose(r, 0)`` and ``allclose(estimate, reference)``.

    It is a positive control for the band, NOT evidence for the paper's
    prediction (a): what lands an estimator here is a regularisation constant
    somebody typed, not the generator having been learned.
    """
    from sklearn.linear_model import Ridge

    del rng
    records = trial.records
    observed = sorted(
        {(int(g), int(a)) for g, a in zip(records["stratum"], records["treated"])}
    )
    if len(observed) < 2 * N_STRATA:
        return float("nan")
    index = {cell: j for j, cell in enumerate(observed)}
    design = np.zeros((len(records), len(observed)))
    for i, (g, a) in enumerate(zip(records["stratum"], records["treated"])):
        design[i, index[(int(g), int(a))]] = 1.0
    model = Ridge(alpha=float(alpha), fit_intercept=False, solver="cholesky")
    model.fit(design, records["outcome"].to_numpy(dtype=float))
    mu = {cell: float(model.coef_[j]) for cell, j in index.items()}
    return _standardise_cell_means(records, mu)


def gcomp_mlp_outcome(trial: Trial, rng: np.random.Generator) -> float:
    """G-computation whose outcome model is a neural network, re-fit on D~.

    Saturated in principle -- the design carries a treatment-by-stratum
    interaction, so the network CAN represent the six cell means exactly -- and
    not saturated in practice, because lbfgs stops where it stops. Its residual
    is optimisation error, which is the same species as the ``tol`` floor in
    ``trial_blindness.ipw_logistic_l2`` and is likewise not a statistical
    property of anything.
    """
    records = trial.records
    stratum = records["stratum"].to_numpy()
    treated = records["treated"].to_numpy(dtype=float)
    onehot = _onehot(stratum)
    design = np.column_stack([onehot, treated, onehot * treated[:, None]])
    y = records["outcome"].to_numpy(dtype=float)
    centre, scale = float(y.mean()), float(y.std(ddof=1))
    if not np.isfinite(scale) or scale == 0.0:
        return float("nan")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = _mlp(_seed_int(rng)).fit(design, (y - centre) / scale)

    grid = np.array([[g, a] for g, a in CELLS], dtype=float)
    grid_onehot = _onehot(grid[:, 0].astype(int))
    grid_design = np.column_stack(
        [grid_onehot, grid[:, 1], grid_onehot * grid[:, 1][:, None]]
    )
    predicted = model.predict(grid_design) * scale + centre
    mu = {cell: float(v) for cell, v in zip(CELLS, predicted)}
    return _standardise_cell_means(records, mu)


#: Every estimator, on one signature: ``(trial, rng) -> float``. All are re-fit
#: on the synthetic cohort and none of them ever sees the training cohort.
ESTIMATORS: dict[str, Callable[[Trial, np.random.Generator], float]] = {
    "gcomp-saturated": gcomp_saturated,
    "ipw-saturated": lambda t, rng: ipw_saturated(t),
    "gcomp-gmm-outcome-K1": lambda t, rng: gcomp_gmm_outcome(
        t, n_components=1, rng=rng
    ),
    "gcomp-gmm-outcome-K4": lambda t, rng: gcomp_gmm_outcome(
        t, n_components=4, rng=rng
    ),
    "gcomp-ridge-outcome-a1e-5": lambda t, rng: gcomp_ridge_outcome(
        t, alpha=1e-5, rng=rng
    ),
    "gcomp-ridge-outcome-a1e-3": lambda t, rng: gcomp_ridge_outcome(
        t, alpha=1e-3, rng=rng
    ),
    "gcomp-ridge-outcome-a1e-1": lambda t, rng: gcomp_ridge_outcome(
        t, alpha=1e-1, rng=rng
    ),
    "gcomp-mlp-outcome": gcomp_mlp_outcome,
    "ols-stratum-dummies": lambda t, rng: ols_stratum_dummies(t),
    "unadjusted": lambda t, rng: unadjusted(t),
}

#: Whether the estimator is the reference's own functional, as an ALGEBRAIC
#: claim made before the run. Measured against the residual in
#: ``tests/test_learned_generator_can_fail.py``; a declared value that the
#: measurement contradicts is a finding, not a typo to be edited away.
SHARES_FUNCTIONAL_WITH_T_DRAW: dict[str, bool] = {
    "gcomp-saturated": True,
    "ipw-saturated": True,
    "gcomp-gmm-outcome-K1": True,
    "gcomp-gmm-outcome-K4": True,
    "gcomp-ridge-outcome-a1e-5": False,
    "gcomp-ridge-outcome-a1e-3": False,
    "gcomp-ridge-outcome-a1e-1": False,
    "gcomp-mlp-outcome": False,
    "ols-stratum-dummies": False,
    "unadjusted": False,
}

#: The negative control that carries the experiment's internal validity
#: (prereg 4.1): consistent for theta, and NOT the reference's functional. If
#: its residual is not clearly larger than the near-zero arms, the screen is
#: measuring nothing and no other row may be read.
CONTROL_ESTIMATOR: str = "ols-stratum-dummies"
CONTROL_MUST_EXCEED: float = 1e-2

#: The references a residual may be taken against. ``theta_requested`` is not
#: among them: that residual is the recovery curve, which is the thing being
#: audited.
REFERENCES: tuple[str, ...] = ("T_draw", "T_model")


# ---------------------------------------------------------------------------
# The grid
# ---------------------------------------------------------------------------


def run_grid(
    *,
    seed: int,
    n_seeds: int = 6,
    theta: float = THETA,
    cohort_sizes: tuple[int, ...] = (100, 200, 500, 1000, 2000, 5000),
    n_replicates: int = 30,
    generators: Sequence[str] = tuple(GENERATORS),
    estimators: Sequence[str] = tuple(ESTIMATORS),
) -> pd.DataFrame:
    """One row per (generator, seed, cohort size, replicate, estimator).

    Both references are columns rather than rows, so the two residuals for one
    estimate are measured on the SAME synthetic cohort. Anything else would make
    the comparison between them a comparison between two draws.
    """
    rows: list[dict] = []
    for gen_index, gen_name in enumerate(generators):
        fit = GENERATORS[gen_name]
        for seed_index in range(n_seeds):
            for n_patients in cohort_sizes:
                for replicate in range(n_replicates):
                    keys = (seed_index, gen_index, n_patients, replicate)
                    training = simulate_trial(
                        n_patients, theta, rng=_rng(seed, "train", *keys)
                    ).records
                    try:
                        fitted = fit(
                            training,
                            n_patients=n_patients,
                            theta=theta,
                            rng=_rng(seed, "fit", *keys),
                        )
                        synthetic = fitted.sample(n_patients, _rng(seed, "synth", *keys))
                    except GeneratorFitError as exc:
                        rows.append(
                            {
                                "generator": gen_name,
                                "is_learned": IS_LEARNED[gen_name],
                                "capacity": CAPACITY[gen_name],
                                "seed_index": seed_index,
                                "n_patients": n_patients,
                                "replicate": replicate,
                                "estimator": None,
                                "generator_failed": True,
                                "generator_failure": str(exc),
                            }
                        )
                        continue

                    trial = _as_trial(synthetic, theta)
                    training_means = cell_means(training)
                    deviation = max(
                        (
                            abs(fitted.implied_cell_means[cell] - training_means[cell])
                            for cell in CELLS
                            if cell in fitted.implied_cell_means
                            and cell in training_means
                        ),
                        default=float("nan"),
                    )
                    common = {
                        "generator": gen_name,
                        "is_learned": IS_LEARNED[gen_name],
                        "capacity": CAPACITY[gen_name],
                        "generator_failed": False,
                        "generator_failure": None,
                        "seed_index": seed_index,
                        "n_patients": n_patients,
                        "replicate": replicate,
                        "theta_requested": float(theta),
                        "T_draw": trial.theta_realised,
                        "T_model": fitted.implied_effect,
                        "T_model_mc_se": fitted.implied_effect_mc_se,
                        "T_model_is_monte_carlo": fitted.implied_is_monte_carlo,
                        "max_cellmean_deviation": float(deviation),
                        "generator_is_saturated": bool(deviation < SATURATED_BELOW),
                    }
                    for est_index, est_name in enumerate(estimators):
                        estimate = ESTIMATORS[est_name](
                            trial, _rng(seed, "estimator", *keys, est_index)
                        )
                        row = dict(common)
                        row["estimator"] = est_name
                        row["shares_functional_with_T_draw"] = (
                            SHARES_FUNCTIONAL_WITH_T_DRAW[est_name]
                        )
                        row["estimate"] = estimate
                        row["difference_vs_requested"] = estimate - theta
                        for reference in REFERENCES:
                            row[f"residual_vs_{reference}"] = abs(
                                estimate - row[reference]
                            )
                        rows.append(row)
    return pd.DataFrame(rows)


def _long_residuals(runs: pd.DataFrame) -> pd.DataFrame:
    """Melt the two reference columns into a ``reference`` column."""
    runs = runs[~runs["generator_failed"].astype(bool)]
    parts = []
    keep = [
        "generator",
        "is_learned",
        "capacity",
        "seed_index",
        "n_patients",
        "replicate",
        "estimator",
        "shares_functional_with_T_draw",
        "theta_requested",
        "estimate",
        "difference_vs_requested",
        "max_cellmean_deviation",
        "generator_is_saturated",
        "T_model_mc_se",
        "T_model_is_monte_carlo",
    ]
    for reference in REFERENCES:
        part = runs[keep + [reference, f"residual_vs_{reference}"]].copy()
        part = part.rename(
            columns={reference: "truth", f"residual_vs_{reference}": "residual"}
        )
        part["reference"] = reference
        parts.append(part)
    return pd.concat(parts, ignore_index=True)


# ---------------------------------------------------------------------------
# The screens
# ---------------------------------------------------------------------------


def flip_tolerances(residual: np.ndarray, truth: np.ndarray) -> dict[str, float]:
    """The tolerance at which each of the two ``allclose`` spellings flips.

    ``np.allclose(r, 0.0)`` tests ``|r| <= atol + rtol*|0|``, so the relative
    term vanishes and the screen is pure ``atol``. It flips True->False exactly
    at ``max|r|``.

    ``np.allclose(estimate, reference)`` tests ``|r| <= atol + rtol*|T|``, so
    each pair carries its own allowance and the screen flips at
    ``max(|r| - rtol*|T|)``, clipped below at zero.

    The two differ by ``rtol*|T| ~ 3e-5`` at theta=3 -- three orders of
    magnitude of blindness that depend only on which line was typed.
    """
    residual = np.asarray(residual, dtype=float)
    truth = np.asarray(truth, dtype=float)
    if residual.size == 0:
        return {"flip_atol": float("nan"), "flip_atol_paired": float("nan")}
    paired = np.max(residual - ALLCLOSE_RTOL * np.abs(truth))
    return {
        "flip_atol": float(np.max(residual)),
        "flip_atol_paired": float(max(paired, 0.0)),
    }


def screen_band(allclose_vs_zero: bool, allclose_paired: bool, n_valid: int) -> str:
    """Which of the two spellings, if either, sees this cell.

    ``blind-band`` is the regime prediction (b) is about: the residual is above
    ``atol`` so screening it against zero fires, and below ``atol + rtol*|T|``
    so screening the estimate against the reference does not. A cell landing
    here is invisible to one practitioner and obvious to another, with no
    statistical difference between them.
    """
    if n_valid == 0:
        return "no-valid-pairs"
    if allclose_vs_zero and allclose_paired:
        return "invisible-to-both"
    if not allclose_vs_zero and allclose_paired:
        return "blind-band"
    return "caught-by-both"


def primary_table(runs: pd.DataFrame) -> pd.DataFrame:
    """THE PRIMARY TABLE (prereg 7): one row per
    (generator, estimator, reference, cohort size).

    ``max_residual`` is the column that matters and it is a MAXIMUM, because
    section blind is explicit that "a zero median can conceal nonzero residuals,
    which is why the equality screen uses the maximum".
    """
    long = _long_residuals(runs)
    out_rows: list[dict] = []
    grouped = long.groupby(
        ["generator", "estimator", "reference", "n_patients"], dropna=False
    )
    for (generator, estimator, reference, n_patients), block in grouped:
        residual = block["residual"].to_numpy(dtype=float)
        truth = block["truth"].to_numpy(dtype=float)
        estimate = block["estimate"].to_numpy(dtype=float)
        valid = np.isfinite(residual) & np.isfinite(truth) & np.isfinite(estimate)
        r, t = residual[valid], truth[valid]

        n_valid = int(valid.sum())
        if n_valid == 0:
            # An empty array makes np.allclose return True, which is how a
            # screen reports "clean" because everything was filtered out.
            max_residual = float("nan")
            allclose_zero = False
            allclose_paired = False
            rho = float("nan")
            flips = {"flip_atol": float("nan"), "flip_atol_paired": float("nan")}
            exactly_zero = False
            median_residual = float("nan")
        else:
            max_residual = float(np.max(r))
            allclose_zero = bool(np.allclose(max_residual, 0.0))
            allclose_paired = bool(
                np.allclose(estimate[valid], t, rtol=ALLCLOSE_RTOL, atol=ALLCLOSE_ATOL)
            )
            generator_error = np.abs(
                t - block["theta_requested"].to_numpy(dtype=float)[valid]
            )
            denominator = float(np.median(generator_error))
            rho = float(np.median(r) / denominator) if denominator > 0 else float("nan")
            flips = flip_tolerances(r, t)
            exactly_zero = bool(max_residual == 0.0)
            median_residual = float(np.median(r))

        out_rows.append(
            {
                "generator": generator,
                "is_learned": bool(block["is_learned"].iloc[0]),
                "capacity": float(block["capacity"].iloc[0]),
                "estimator": estimator,
                "shares_functional_with_T_draw": bool(
                    block["shares_functional_with_T_draw"].iloc[0]
                ),
                "reference": reference,
                "n_patients": int(n_patients),
                "n_draws": int(len(block)),
                "n_valid_pairs": n_valid,
                "n_nonfinite_reference": int((~np.isfinite(truth)).sum()),
                "n_nonfinite_estimate": int((~np.isfinite(estimate)).sum()),
                "max_residual": max_residual,
                "median_residual": median_residual,
                "is_exactly_zero": exactly_zero,
                "allclose_residual_vs_zero": allclose_zero,
                "allclose_estimate_vs_reference": allclose_paired,
                "flip_atol": flips["flip_atol"],
                "flip_atol_paired": flips["flip_atol_paired"],
                "screen_band": screen_band(allclose_zero, allclose_paired, n_valid),
                "rho": rho,
                "generator_noise_share": (
                    1.0 / (1.0 + rho) if np.isfinite(rho) else float("nan")
                ),
                "median_estimate": (
                    float(np.median(estimate[valid])) if n_valid else float("nan")
                ),
                "median_recovery_ratio": (
                    float(np.median(estimate[valid])) / THETA
                    if n_valid
                    else float("nan")
                ),
                "median_difference_vs_requested": (
                    float(
                        np.median(
                            block["difference_vs_requested"].to_numpy(dtype=float)[valid]
                        )
                    )
                    if n_valid
                    else float("nan")
                ),
                "generator_is_saturated": bool(block["generator_is_saturated"].all()),
                "max_cellmean_deviation": float(
                    np.nanmax(block["max_cellmean_deviation"].to_numpy(dtype=float))
                ),
                "T_model_mc_se": float(
                    np.nanmedian(block["T_model_mc_se"].to_numpy(dtype=float))
                ),
                "T_model_is_monte_carlo": bool(block["T_model_is_monte_carlo"].iloc[0]),
            }
        )
    out = pd.DataFrame(out_rows)
    return out.sort_values(
        ["reference", "generator", "estimator", "n_patients"]
    ).reset_index(drop=True)


def trend_table(primary: pd.DataFrame) -> pd.DataFrame:
    """Does the residual shrink, grow or stay flat as the cohort grows?

    The question a practitioner scaling up a synthetic cohort actually has. A
    residual pinned at machine precision has no trend to report and says so
    rather than fitting a slope to floating-point noise.
    """
    rows: list[dict] = []
    for (generator, estimator, reference), block in primary.groupby(
        ["generator", "estimator", "reference"]
    ):
        block = block.sort_values("n_patients")
        residual = block["max_residual"].to_numpy(dtype=float)
        n = block["n_patients"].to_numpy(dtype=float)
        usable = np.isfinite(residual) & (residual > 0)
        finite = residual[np.isfinite(residual)]
        largest = float(np.max(finite)) if finite.size else float("nan")
        if usable.sum() < 3 or (np.isfinite(largest) and largest < DEGENERATE_BELOW):
            verdict = (
                "at-machine-precision"
                if np.isfinite(largest) and largest < DEGENERATE_BELOW
                else "too-few-points"
            )
            slope = float("nan")
        else:
            slope = float(
                np.polyfit(np.log10(n[usable]), np.log10(residual[usable]), 1)[0]
            )
            verdict = "shrinks" if slope < -0.2 else "grows" if slope > 0.2 else "flat"
        rows.append(
            {
                "generator": generator,
                "estimator": estimator,
                "reference": reference,
                "log10_slope_vs_log10_n": slope,
                "trend": verdict,
                "max_residual_at_smallest_n": float(residual[0]),
                "max_residual_at_largest_n": float(residual[-1]),
                "smallest_n": int(n[0]),
                "largest_n": int(n[-1]),
            }
        )
    return (
        pd.DataFrame(rows)
        .sort_values(["reference", "generator", "estimator"])
        .reset_index(drop=True)
    )


def capacity_table(primary: pd.DataFrame) -> pd.DataFrame:
    """Residual against mixture components, for the generators that have them.

    Prediction (a) is usually narrated as "more capacity, more approximation".
    Under ``T_draw`` capacity cannot enter, because capacity changes the
    synthetic cohort's distribution and the screen compares two functionals OF
    that cohort. Under ``T_model`` it can. Both are here.
    """
    out = primary[np.isfinite(primary["capacity"].to_numpy(dtype=float))].copy()
    return (
        out[
            [
                "reference",
                "estimator",
                "capacity",
                "generator",
                "n_patients",
                "max_residual",
                "is_exactly_zero",
                "allclose_residual_vs_zero",
                "allclose_estimate_vs_reference",
                "screen_band",
                "rho",
                "T_model_mc_se",
                "max_cellmean_deviation",
                "generator_is_saturated",
            ]
        ]
        .sort_values(["reference", "estimator", "capacity", "n_patients"])
        .reset_index(drop=True)
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

#: Copied verbatim from ``trial_recovery.main`` and ``trial_blindness``. Same
#: generator underneath, same disclaimer; paraphrasing it would let the tables
#: drift apart in what they claim.
SYNTHETIC_NOTE = (
    "Simulated patients with an analytically known treatment effect. Nothing "
    "here is a result about any real trial."
)

#: What ``env/*.yml`` pins, against what the interpreter actually loaded. These
#: differ on the machine this was run on, and the gap is recorded rather than
#: papered over: a residual floor set by a solver's default tolerance is exactly
#: the kind of number that moves between library versions.
PINNED_ENV: dict[str, str] = {
    "scikit-learn": "1.5.1",
    "numpy": "1.26.4",
    "pandas": "2.2.2",
    "scipy": "1.13.1",
}


def _live_env() -> dict[str, str]:
    import scipy
    import sklearn

    return {
        "scikit-learn": sklearn.__version__,
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
    }


def main(argv: Sequence[str] | None = None) -> int:
    """Write the learned-generator tables under ``results/`` with provenance.

        python -m src.harness.trial_learned_generator
    """
    parser = argparse.ArgumentParser(
        description="The residual audit under a learned generative patient model"
    )
    parser.add_argument("--seed", type=int, default=20260913)
    parser.add_argument("--replicates", type=int, default=30)
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
    replicates = 3 if args.quick else args.replicates
    n_seeds = 2 if args.quick else args.seeds

    log.info(
        "grid: %d generators x %d seeds x %d sizes x %d reps x %d estimators",
        len(GENERATORS),
        n_seeds,
        len(cohort_sizes),
        replicates,
        len(ESTIMATORS),
    )
    runs = run_grid(
        seed=args.seed,
        n_seeds=n_seeds,
        cohort_sizes=cohort_sizes,
        n_replicates=replicates,
    )
    primary = primary_table(runs)
    trend = trend_table(primary)
    capacity = capacity_table(primary)

    draw = primary[primary["reference"] == "T_draw"]
    control = draw[draw["estimator"] == CONTROL_ESTIMATOR]
    control_ok = bool((control["max_residual"] > CONTROL_MUST_EXCEED).all())
    if not control_ok:
        log.warning(
            "F6 FIRED: the control estimator %s does not exceed %g at every "
            "cohort size, so the screen has no demonstrated ability to separate "
            "shared-functional from merely-accurate, and NO OTHER ROW MAY BE READ.",
            CONTROL_ESTIMATOR,
            CONTROL_MUST_EXCEED,
        )
    bands = {k: int(v) for k, v in draw["screen_band"].value_counts().items()}
    log.info("screen bands under T_draw: %s", bands)
    for name, block in draw.groupby("estimator"):
        log.info(
            "%-28s max|theta_hat - T_draw| = %.3g  band=%s",
            name,
            block["max_residual"].max(),
            "/".join(sorted(set(block["screen_band"]))),
        )

    common = {
        "seed": args.seed,
        "results_dir": args.results_dir,
        "allow_dirty": args.allow_dirty,
    }
    live = _live_env()
    shared_meta = {
        "sweep_seed": args.seed,
        "n_seed_streams": n_seeds,
        "n_replicates_per_seed": replicates,
        "theta_requested": THETA,
        "cohort_sizes": list(cohort_sizes),
        "is_quick_smoke_test_not_a_result": bool(args.quick),
        "stratum_means": list(STRATUM_MEANS),
        "propensity": list(PROPENSITY),
        "stratum_weights": list(STRATUM_WEIGHTS),
        "outcome_sd": OUTCOME_SD,
        "generators": {
            k: {"learned": IS_LEARNED[k], "capacity": CAPACITY[k]} for k in GENERATORS
        },
        "estimators": list(ESTIMATORS),
        "references": {
            "T_draw": "sum_g (n_g/n)(ybar_g1 - ybar_g0) on the SYNTHETIC cohort, "
            "via trial_recovery._standardised_effect -- the paper's "
            "observed-data reference, PRIMARY",
            "T_model": "the fitted generator's own implied population "
            "standardised effect; closed form for parametric-* and "
            "mlp-conditional-gaussian, Monte-Carlo at n=200000 with a "
            "reported standard error for gmm-joint-*, SECONDARY",
        },
        "allclose_rtol": ALLCLOSE_RTOL,
        "allclose_atol": ALLCLOSE_ATOL,
        "allclose_semantics": (
            "np.allclose(a,b) tests |a-b| <= atol + rtol*|b|. Screening the "
            "residual against 0.0 kills the rtol term and is a pure atol=1e-8 "
            "screen; screening the estimate against the reference allows "
            "atol + rtol*|T| ~ 3e-5 at theta=3. Both are reported."
        ),
        "saturated_below": SATURATED_BELOW,
        "degenerate_below": DEGENERATE_BELOW,
        "oracle_draw": ORACLE_DRAW,
        "seed_key_order": list(SEED_KEY_ORDER),
        "control_estimator": CONTROL_ESTIMATOR,
        "control_must_exceed": CONTROL_MUST_EXCEED,
        "control_estimator_passed": control_ok,
        "prereg": "docs/prereg_learned_generator.md",
        "pinned_env": PINNED_ENV,
        "live_env": live,
        "env_version_gap": {
            k: {"pinned": v, "live": live[k]}
            for k, v in PINNED_ENV.items()
            if live[k] != v
        },
        "SYNTHETIC": SYNTHETIC_NOTE,
    }

    path = write_versioned_table(
        primary,
        "trial_learned_generator_primary",
        extra_meta=shared_meta
        | {
            "screen_bands_under_T_draw": bands,
            "what_this_answers": (
                "Whether a learned generative patient model degrades the "
                "T(D) equality from exact to approximate, as "
                "paper/wmhs/sections/full/refdesign.tex predicts. It does not: "
                "T(D~) and theta_hat(D~) are both functionals of the SAME "
                "synthetic cohort, so whether they coincide is a property of "
                "the estimator's functional form and the generator's "
                "learnedness never enters. Approximate equality does appear, "
                "but from a solver's convergence constant or a regularisation "
                "constant -- numerical decisions, not statistical ones."
            ),
        },
        **common,
    )
    log.info("wrote %s (%d rows)", path, len(primary))

    path = write_versioned_table(
        trend,
        "trial_learned_generator_trend",
        extra_meta=shared_meta
        | {
            "what_this_answers": (
                "Whether the residual shrinks, grows or is flat in cohort size, "
                "which decides whether the problem gets better or worse with "
                "more data."
            ),
        },
        **common,
    )
    log.info("wrote %s (%d rows)", path, len(trend))

    path = write_versioned_table(
        capacity,
        "trial_learned_generator_capacity",
        extra_meta=shared_meta
        | {
            "what_this_answers": (
                "Where -- and whether -- exact becomes approximate as the "
                "learned generator gains capacity. Under T_draw it cannot: "
                "capacity changes the synthetic cohort's distribution, and the "
                "screen compares two functionals OF that cohort."
            ),
        },
        **common,
    )
    log.info("wrote %s (%d rows)", path, len(capacity))
    return 0


if __name__ == "__main__":
    sys.exit(main())
