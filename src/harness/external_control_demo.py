"""Where an observed-data reference does HARM: a plasmode external control arm.

THE OBJECTION THIS ANSWERS. A reviewer of the workshop paper: *"The equality
audit only protects against something when a study uses an observed-data
reference in place of theta as ground truth. The paper never shows that
happening, in its own pipeline or in the seven audited repositories."* Our own
prevalence audit found 0 of 7. So the paper proposed a safeguard against a
mistake nobody was shown making, with no demonstrated consequence.

THE REVIEWER'S ASSUMPTION IS THAT ``theta`` IS AVAILABLE. In plasmode and
resampling designs it is not, and those are the standard settings for synthetic
and external control arms. The cohort is built by resampling real records, so
there is no parametric truth -- only an effect *induced* on real covariates. When
that effect is modified by the covariates, the MARGINAL estimand is a functional
of the realised sample: which records enrolment happened to draw, in what mix,
carrying what latent responsiveness. No scalar in the config equals it. The
analyst must construct a reference from the draw, and there is more than one way
to do that.

THE WORKFLOW. A single-arm trial with an external control arm from a disease
registry (ICH E10 territory). Before touching the real trial the statistician
validates the planned estimator by plasmode simulation on the registry:

    1. a fixed finite registry POOL of records, each with a stratum (line of
       therapy), its own outcome Y0, and a latent responsiveness r fixed to the
       record;
    2. enrol ``n_trial`` records WITHOUT replacement under eligibility weights --
       a new agent enters in refractory disease, so the trial arm is enriched for
       heavily pretreated patients. THIS IS THE CASE-MIX SHIFT, and it is the
       defining nuisance of external control arms;
    3. induce the effect: Y1 = Y0 + DELTA[g] * r. The record keeps its own Y0.
       The benefit is effect-modified -- large in treatment-naive patients, near
       zero in heavily pretreated ones;
    4. draw ``CONTROL_RATIO * n_trial`` external controls uniformly from the rest
       of the pool, keeping their own Y0;
    5. stack them into an analysis file. Arm membership depends on stratum, so
       the strata confound, which is why anybody weights.

The estimand is the ATT: the benefit IN THE ENROLLED POPULATION. That is what the
label would claim and what the confirmatory trial would be powered on.

THE DEFECT. The estimator targets the POOLED analysis population instead.
Saturated IPW with ATE weights ``1/e`` and ``1/(1-e)`` where the estimand needs
ATT weights ``1`` and ``e/(1-e)``. One line. Reachable by a second, wholly
different route -- g-computation standardised over every row in the analysis file
rather than over the trial arm -- which is why this is a way of thinking rather
than a typo. Because the control arm is three times the trial arm, the pooled mix
is dominated by the REGISTRY, so the estimator silently answers "what would this
agent do in the registry?" while the analyst reads "what does it do in the
patients we enrolled?".

WHY THE ANALYST DOES NOT CATCH IT. Positivity is fine (0.083/0.217/0.571, no
stratum empty, no trimming rule fires). The post-weighting covariate balance
table is clean, because ATE weights balance the two ARMS against each other and
the defect lives in balance against the ENROLLED population, which that table
does not look at. The number has the right sign and magnitude and tightens with
n. And the validation confirms it -- which is the part this module is about.

THE FIVE REFERENCES, and the experiment is the gap between them.

``obs-pooled``    THE MISTAKE. sum_g (n_g/n) (Ybar_1g - Ybar_0g) over the whole
                  analysis file. "The average treatment effect in my cohort",
                  computed the way the phrase reads. It is the same functional
                  ``trial_recovery`` ships as ``_standardised_effect`` and the
                  same one the paper analyses. The defective estimator
                  reproduces it exactly, so the validation cannot fail.
``obs-trial``     the same contrast at the TRIAL arm's realised mix. Observed
                  data, right population. The correct estimators reproduce THIS
                  one exactly -- which is the point that equality is a property
                  of the (estimator, reference) pair and not a verdict.
``po-trial``      THE HONEST REFERENCE. mean over the trial arm of (Y1 - Y0):
                  the realised individual causal effects among the enrolled.
                  Available in plasmode because we built Y1 ourselves; not
                  reproducible by any estimator, which sees one arm per record.
                  This is PlasmodeSimulation's criterion 3b, kept.
``po-pooled``     CONTROL. The same potential-outcome contrast over ALL records.
                  Breaks the equality, keeps the wrong population -- and does
                  NOT reveal the defect. Without this arm the demonstration would
                  license "any non-reproduced reference will do", which is false.
``induced-trial`` the parameter-only construction, ignoring r. Secondary anchor;
                  it is NOT the realised truth, which is what makes the point
                  that the reviewer's ``theta`` does not exist here.

Pre-registered with falsifiers before the run: ``docs/prereg_external_control_demo.md``.
"""

from __future__ import annotations

import argparse
import logging
import sys
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from src.common.io import write_versioned_table

log = logging.getLogger(__name__)

# --------------------------------------------------------------------------
# Design constants. Every one of these is fixed in prereg §4, before this file
# existed. Changing one is an amendment to that document, not an edit here.
# --------------------------------------------------------------------------

#: Registry scale. Large enough that drawing both arms without replacement does
#: not measurably deplete any stratum at the largest cohort size.
N_POOL: int = 50_000

#: Line of therapy: treatment-naive / one prior line / heavily pretreated.
STRATUM_NAMES: tuple[str, ...] = ("naive", "one-prior", "heavily-pretreated")

#: The registry's own case mix. Registries are dominated by earlier-line patients.
POOL_MIX: tuple[float, ...] = (0.55, 0.30, 0.15)

#: The enrolled trial's case mix. A new agent enters in refractory disease, so
#: the trial is enriched for the stratum the registry has least of. This is the
#: case-mix shift, and at ``shift = 0`` below it collapses onto POOL_MIX.
TRIAL_MIX: tuple[float, ...] = (0.15, 0.25, 0.60)

#: Untreated PFS in months, falling with line of therapy.
POOL_MEANS: tuple[float, ...] = (12.0, 8.0, 4.0)
POOL_SD: float = 3.0

#: The induced benefit, in months, per stratum. EFFECT-MODIFIED: concentrated
#: where disease biology is least eroded. This is what makes the marginal effect
#: depend on the case mix, and therefore on which records were drawn.
DELTA: tuple[float, ...] = (3.0, 1.5, 0.25)

#: Per-record latent responsiveness, fixed to the record when the pool is built.
#: Mean 1, sd 0.5. Its job is to make the realised marginal truth a property of
#: WHICH RECORDS WERE DRAWN rather than of the parameters plus realised counts.
#: Falsifier F5 checks it actually does that job.
RESPONSIVENESS_SHAPE: float = 4.0
RESPONSIVENESS_SCALE: float = 0.25

#: External control arms are larger than the trial they support.
CONTROL_RATIO: float = 3.0

#: The pre-specified clinically meaningful difference. The decision threshold the
#: analyst would act on.
MCID: float = 1.5

#: The pool is built once, from this seed, and is the same object in every
#: replicate. That is the plasmode mechanic: a fixed finite set of records,
#: resampled, each keeping its own outcome.
POOL_SEED: int = 20260915

#: Equality tolerance, in outcome units (months), chosen explicitly as §blind
#: requires rather than inherited from a library default.
EQUALITY_TOL: float = 1e-10

COHORT_SIZES: tuple[int, ...] = (50, 100, 200, 400, 800, 1600)

#: The conventional post-weighting covariate-balance threshold.
SMD_THRESHOLD: float = 0.1


# --------------------------------------------------------------------------
# The registry, and the plasmode draw
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class Registry:
    """A fixed finite pool of records. Built once; resampled, never re-drawn.

    ``outcome`` is each record's own value and is carried into the analysis file
    unchanged. ``responsiveness`` is latent: it scales the induced effect for
    that record and is never visible to any estimator.
    """

    records: pd.DataFrame

    def __len__(self) -> int:
        return len(self.records)


def build_registry(
    *,
    n_pool: int = N_POOL,
    pool_mix: tuple[float, ...] = POOL_MIX,
    pool_means: tuple[float, ...] = POOL_MEANS,
    pool_sd: float = POOL_SD,
    seed: int = POOL_SEED,
) -> Registry:
    """The registry pool. Deterministic in ``seed``, and that is the point.

    THE POOL IS SYNTHETIC and prereg §10 says so without argument: no real
    registry is available to this repository. What is preserved is the plasmode
    *mechanic* -- fixed finite pool, resampled, each record keeping its own
    outcome, only the effect induced -- not the provenance of the numbers.
    """
    rng = np.random.default_rng(seed)
    stratum = rng.choice(len(pool_mix), size=n_pool, p=list(pool_mix))
    outcome = rng.normal(np.asarray(pool_means)[stratum], pool_sd)
    responsiveness = rng.gamma(RESPONSIVENESS_SHAPE, RESPONSIVENESS_SCALE, n_pool)
    return Registry(
        records=pd.DataFrame(
            {
                "record_id": np.arange(n_pool),
                "stratum": stratum,
                "outcome": outcome,
                "responsiveness": responsiveness,
            }
        )
    )


def shifted_trial_mix(shift: float, *, trial_mix: tuple[float, ...] = TRIAL_MIX,
                      pool_mix: tuple[float, ...] = POOL_MIX) -> tuple[float, ...]:
    """Interpolate the enrolled case mix from the registry's own to the enriched one.

    ``shift = 0`` is the experiment's NULL: the trial enrols the registry's own
    case mix, the ATT and the pooled ATE coincide, and the alleged harm must
    vanish. Falsifier F4 is exactly the requirement that it does.
    """
    return tuple(
        (1.0 - shift) * pool_mix[g] + shift * trial_mix[g] for g in range(len(pool_mix))
    )


@dataclass(frozen=True)
class PlasmodeDraw:
    """One simulated external-control study, plus every reference it admits.

    ``records`` is the analysis file -- what an estimator sees: stratum, arm,
    one outcome per patient. ``trial_truth`` and ``pooled_truth`` are the
    potential-outcome references, which no estimator can see because they read
    both arms for the same record.
    """

    records: pd.DataFrame
    trial_index: np.ndarray
    references: dict[str, float]


def simulate_study(
    n_trial: int,
    *,
    registry: Registry,
    rng: np.random.Generator,
    shift: float = 1.0,
    delta: tuple[float, ...] = DELTA,
    control_ratio: float = CONTROL_RATIO,
) -> PlasmodeDraw:
    """Enrol, induce, and draw the external control arm.

    Enrolment is WITHOUT replacement and under stratum-specific eligibility
    weights, so the trial arm's case mix is the enriched one and the external
    controls that follow carry the registry's own. Both arms are disjoint sets of
    records, which is what "contemporaneous but distinct patients" means.
    """
    pool = registry.records
    n_control = int(round(control_ratio * n_trial))
    if n_trial + n_control > len(pool):
        raise ValueError(
            f"asked for {n_trial + n_control} records from a pool of {len(pool)}"
        )

    target_mix = shifted_trial_mix(shift)
    observed_mix = np.asarray(
        [float((pool["stratum"] == g).mean()) for g in range(len(target_mix))]
    )
    # Eligibility weights: how much more likely a stratum is to enrol than its
    # registry share. This is what an inclusion criterion does.
    eligibility = np.asarray(target_mix) / np.where(observed_mix > 0, observed_mix, 1.0)
    weights = eligibility[pool["stratum"].to_numpy()]
    weights = weights / weights.sum()

    enrolled = rng.choice(len(pool), size=n_trial, replace=False, p=weights)
    remaining = np.setdiff1d(np.arange(len(pool)), enrolled, assume_unique=False)
    controls = rng.choice(remaining, size=n_control, replace=False)

    stratum_all = pool["stratum"].to_numpy()
    y0_all = pool["outcome"].to_numpy()
    effect_all = np.asarray(delta)[stratum_all] * pool["responsiveness"].to_numpy()

    index = np.concatenate([enrolled, controls])
    treated = np.concatenate([np.ones(n_trial, int), np.zeros(n_control, int)])
    y0 = y0_all[index]
    individual_effect = effect_all[index]
    y1 = y0 + individual_effect
    observed = np.where(treated == 1, y1, y0)

    records = pd.DataFrame(
        {
            "stratum": stratum_all[index],
            "treated": treated,
            "outcome": observed,
        }
    )
    trial_rows = np.arange(n_trial)

    references = {
        # THE MISTAKE: the observed arm contrast standardised to the pooled
        # analysis population. This is the functional the defective estimator is.
        "obs-pooled": standardised_contrast(records, weights_from="pooled"),
        # Observed data, right population.
        "obs-trial": standardised_contrast(records, weights_from="trial"),
        # THE HONEST REFERENCE: realised individual causal effects among the
        # enrolled. Reads both arms of the same record, so no estimator can be it.
        "po-trial": float(np.mean(individual_effect[trial_rows])),
        # CONTROL: equality broken, population still wrong.
        "po-pooled": float(np.mean(individual_effect)),
        # Parameters plus realised counts, ignoring responsiveness. NOT the
        # realised truth -- F5 checks that it differs.
        "induced-trial": float(
            np.mean(np.asarray(delta)[stratum_all[enrolled]])
        ),
    }
    return PlasmodeDraw(records=records, trial_index=index[:n_trial], references=references)


# --------------------------------------------------------------------------
# References and estimators
# --------------------------------------------------------------------------


def standardised_contrast(records: pd.DataFrame, *, weights_from: str) -> float:
    """Stratum contrasts of observed arm means, standardised to a chosen population.

    ``weights_from="pooled"`` weights by each stratum's share of the WHOLE
    analysis file. ``weights_from="trial"`` weights by its share of the TREATED
    arm. The two differ by exactly the case-mix shift, and that difference is the
    entire subject of this module.

    Returns NaN when a stratum is missing an arm -- a positivity failure, which
    is not a number and is not reported as one.
    """
    if weights_from == "pooled":
        basis = records
    elif weights_from == "trial":
        basis = records[records["treated"] == 1]
    else:  # pragma: no cover - guarded by the caller's literal
        raise ValueError(f"unknown weights_from {weights_from!r}")
    total = len(basis)
    if total == 0:
        return float("nan")

    effect = 0.0
    for stratum, group in records.groupby("stratum"):
        treated = group.loc[group["treated"] == 1, "outcome"]
        control = group.loc[group["treated"] == 0, "outcome"]
        if treated.empty or control.empty:
            return float("nan")
        share = float((basis["stratum"] == stratum).sum()) / total
        effect += share * (treated.mean() - control.mean())
    return float(effect)


def att_standardisation(draw: PlasmodeDraw) -> float:
    """CORRECT. Standardise the stratum contrasts to the ENROLLED case mix."""
    return standardised_contrast(draw.records, weights_from="trial")


def ate_standardisation(draw: PlasmodeDraw) -> float:
    """THE DEFECT, route one. "Average the stratum contrasts over my analysis file."

    No weight formula was mistyped here. The analyst wrote the sentence they meant
    and standardised over ``len(df)``, which is the registry-dominated pooled
    population rather than the patients they enrolled. This arm exists because
    falsifier F3a says the defect must be reachable without the IPW weight slip,
    or it is a weight-formula curiosity rather than a way of thinking.
    """
    return standardised_contrast(draw.records, weights_from="pooled")


def _saturated_propensity(records: pd.DataFrame) -> dict[int, float] | None:
    """P(treated | stratum) from the records. None if any stratum is degenerate."""
    out: dict[int, float] = {}
    for stratum, group in records.groupby("stratum"):
        p_hat = float(group["treated"].mean())
        if p_hat <= 0.0 or p_hat >= 1.0:
            return None
        out[int(stratum)] = p_hat
    return out


def att_ipw(draw: PlasmodeDraw) -> float:
    """CORRECT. Saturated ATT weights: treated 1, control e/(1-e).

    Algebraically the same functional as ``att_standardisation`` -- which is the
    paper's own point, made again at the right target population. Two estimators
    from different literatures, indistinguishable by reading the code, identical
    in the residual.
    """
    records = draw.records
    e = _saturated_propensity(records)
    if e is None:
        return float("nan")
    n_treated = int((records["treated"] == 1).sum())
    if n_treated == 0:
        return float("nan")

    treated_mask = records["treated"].to_numpy() == 1
    outcome = records["outcome"].to_numpy()
    strata = records["stratum"].to_numpy()
    e_i = np.asarray([e[int(g)] for g in strata])

    w = np.where(treated_mask, 1.0, e_i / (1.0 - e_i))
    return float(
        (w[treated_mask] * outcome[treated_mask]).sum() / n_treated
        - (w[~treated_mask] * outcome[~treated_mask]).sum() / n_treated
    )


def ate_ipw(draw: PlasmodeDraw) -> float:
    """THE DEFECT, route two. Saturated ATE weights: 1/e and 1/(1-e).

    The estimand is the ATT and the weights are the ATE's. One line. It is the
    single most common slip in this design, and it is invisible in every
    diagnostic an external-control protocol actually runs: positivity holds, the
    post-weighting balance table is clean (``balance_table`` measures that), and
    the estimate is tight and plausible.
    """
    records = draw.records
    e = _saturated_propensity(records)
    if e is None:
        return float("nan")
    n = len(records)

    treated_mask = records["treated"].to_numpy() == 1
    outcome = records["outcome"].to_numpy()
    e_i = np.asarray([e[int(g)] for g in records["stratum"].to_numpy()])

    w = np.where(treated_mask, 1.0 / e_i, 1.0 / (1.0 - e_i))
    return float(
        (w[treated_mask] * outcome[treated_mask]).sum() / n
        - (w[~treated_mask] * outcome[~treated_mask]).sum() / n
    )


def ate_ipw_logistic(draw: PlasmodeDraw, *, tol: float = 1e-10,
                     max_iter: int = 200) -> float:
    """THE DEFECT with the propensity FITTED, not counted. Prereg §6.1.

    Saturated-discrete design is what makes the equality exact, so this arm
    attacks our own claim: an analyst fits a propensity model, they do not tally
    stratum proportions by hand. Unpenalised logistic regression on the stratum
    dummies is saturated in intent and therefore recovers the empirical
    proportions -- but only to the solver's convergence tolerance, not to machine
    zero. Falsifier F3c is the question of whether the equality survives that.

    Newton-Raphson rather than sklearn, deliberately: sklearn's default is
    L2-penalised, and ``trial_blindness`` already documents that its residual is
    then a function of a regularisation constant and two library defaults rather
    than of any statistical decision.
    """
    records = draw.records
    strata = np.sort(records["stratum"].unique())
    if len(strata) < 2:
        return float("nan")
    design = np.column_stack(
        [np.ones(len(records))]
        + [(records["stratum"].to_numpy() == g).astype(float) for g in strata[1:]]
    )
    y = records["treated"].to_numpy(dtype=float)

    beta = np.zeros(design.shape[1])
    for _ in range(max_iter):
        eta = design @ beta
        p = 1.0 / (1.0 + np.exp(-eta))
        w = p * (1.0 - p)
        if np.any(w < 1e-12):
            return float("nan")
        hessian = design.T @ (design * w[:, None])
        try:
            step = np.linalg.solve(hessian, design.T @ (y - p))
        except np.linalg.LinAlgError:  # pragma: no cover - singular design
            return float("nan")
        beta = beta + step
        if np.max(np.abs(step)) < tol:
            break

    e_i = 1.0 / (1.0 + np.exp(-(design @ beta)))
    if np.any(e_i <= 0.0) or np.any(e_i >= 1.0):
        return float("nan")
    treated_mask = y == 1
    outcome = records["outcome"].to_numpy()
    n = len(records)
    w = np.where(treated_mask, 1.0 / e_i, 1.0 / (1.0 - e_i))
    return float(
        (w[treated_mask] * outcome[treated_mask]).sum() / n
        - (w[~treated_mask] * outcome[~treated_mask]).sum() / n
    )


def att_matching_all(draw: PlasmodeDraw) -> float:
    """CORRECT. Each enrolled patient against every stratum-matched control.

    The most common external-control method there is. With all controls in the
    stratum used, it is the same functional as ATT standardisation again -- a
    third literature, a third vocabulary, the same arithmetic.
    """
    records = draw.records
    treated = records[records["treated"] == 1]
    if treated.empty:
        return float("nan")
    total = 0.0
    for stratum, group in treated.groupby("stratum"):
        matched = records[(records["stratum"] == stratum) & (records["treated"] == 0)]
        if matched.empty:
            return float("nan")
        total += len(group) * (group["outcome"].mean() - matched["outcome"].mean())
    return float(total / len(treated))


def ols_stratum_dummies(draw: PlasmodeDraw) -> float:
    """A third population again: the variance-weighted one, by Frisch-Waugh.

    What most people would actually type. Its recovery curve sits off 1 here for
    a reason that is neither the defect nor an error -- effect heterogeneity makes
    its estimand a different weighted average -- which is precisely why a
    departure from a reference is not a quality verdict.
    """
    records = draw.records
    strata = np.sort(records["stratum"].unique())
    if len(strata) < 2:
        return float("nan")
    columns = [np.ones(len(records)), records["treated"].to_numpy(dtype=float)]
    columns.extend(
        (records["stratum"].to_numpy() == g).astype(float) for g in strata[1:]
    )
    design = np.column_stack(columns)
    if np.linalg.matrix_rank(design) < design.shape[1]:
        return float("nan")
    beta, *_ = np.linalg.lstsq(
        design, records["outcome"].to_numpy(dtype=float), rcond=None
    )
    return float(beta[1])


def unadjusted(draw: PlasmodeDraw) -> float:
    """Difference in arm means. Confounded, and included so the curve has prey."""
    records = draw.records
    treated = records.loc[records["treated"] == 1, "outcome"]
    control = records.loc[records["treated"] == 0, "outcome"]
    if treated.empty or control.empty:
        return float("nan")
    return float(treated.mean() - control.mean())


ESTIMATORS = {
    "att-standardisation": att_standardisation,
    "att-ipw": att_ipw,
    "att-matching-all": att_matching_all,
    "ate-standardisation": ate_standardisation,
    "ate-ipw": ate_ipw,
    "ate-ipw-logistic": ate_ipw_logistic,
    "ols-stratum-dummies": ols_stratum_dummies,
    "unadjusted": unadjusted,
}

#: Which target population each estimator actually answers for. This is the
#: property the recovery curve cannot see, and the one the defect corrupts.
TARGETS: dict[str, str] = {
    "att-standardisation": "enrolled",
    "att-ipw": "enrolled",
    "att-matching-all": "enrolled",
    "ate-standardisation": "pooled",
    "ate-ipw": "pooled",
    "ate-ipw-logistic": "pooled",
    "ols-stratum-dummies": "variance-weighted",
    "unadjusted": "none",
}

#: The seven estimators of prereg §6. ``ate-ipw-logistic`` is the §6.1 secondary
#: arm and is excluded from the F1(c) "smallest RMSE of all seven" comparison,
#: because it is an implementation probe of an arm already in the set rather than
#: an eighth method a protocol would choose between.
PRIMARY_ESTIMATORS: tuple[str, ...] = (
    "att-standardisation",
    "att-ipw",
    "att-matching-all",
    "ate-standardisation",
    "ate-ipw",
    "ols-stratum-dummies",
    "unadjusted",
)

REFERENCES: tuple[str, ...] = (
    "obs-pooled",
    "obs-trial",
    "po-trial",
    "po-pooled",
    "induced-trial",
)


# --------------------------------------------------------------------------
# The balance table an analyst would actually run
# --------------------------------------------------------------------------


def balance_table(draw: PlasmodeDraw) -> pd.DataFrame:
    """Post-weighting standardised mean differences, per stratum indicator.

    Two columns matter and they answer different questions.

    ``smd_between_arms``  what a protocol reports: does weighting balance the
                          external controls against the trial arm? Under ATE
                          weights both arms are reweighted to the POOLED mix, so
                          this is ~0 and the table passes.
    ``smd_vs_enrolled``   what nobody reports: does the weighted pseudo-population
                          match the patients actually ENROLLED? Under ATE weights
                          it does not, and that gap IS the defect.

    Falsifier F3e is the charge that the routine diagnostic already catches this.
    It fires if ``smd_between_arms`` exceeds SMD_THRESHOLD under ATE weights.
    """
    records = draw.records
    e = _saturated_propensity(records)
    if e is None:
        return pd.DataFrame()
    treated_mask = records["treated"].to_numpy() == 1
    strata = records["stratum"].to_numpy()
    e_i = np.asarray([e[int(g)] for g in strata])

    schemes = {
        "att": np.where(treated_mask, 1.0, e_i / (1.0 - e_i)),
        "ate": np.where(treated_mask, 1.0 / e_i, 1.0 / (1.0 - e_i)),
        "unweighted": np.ones(len(records)),
    }
    enrolled_mix = np.asarray(
        [float((strata[treated_mask] == g).mean()) for g in range(len(POOL_MIX))]
    )

    rows = []
    for scheme, w in schemes.items():
        for g in range(len(POOL_MIX)):
            indicator = (strata == g).astype(float)
            p_t = float(
                np.average(indicator[treated_mask], weights=w[treated_mask])
            )
            p_c = float(
                np.average(indicator[~treated_mask], weights=w[~treated_mask])
            )
            p_w = float(np.average(indicator, weights=w))
            # Standardised difference for a binary covariate.
            pooled_var = 0.5 * (p_t * (1 - p_t) + p_c * (1 - p_c))
            smd_arms = 0.0 if pooled_var <= 0 else (p_t - p_c) / np.sqrt(pooled_var)
            p_e = float(enrolled_mix[g])
            var_e = p_e * (1 - p_e)
            smd_enrolled = 0.0 if var_e <= 0 else (p_w - p_e) / np.sqrt(var_e)
            rows.append(
                {
                    "weighting": scheme,
                    "stratum": g,
                    "stratum_name": STRATUM_NAMES[g],
                    "share_treated_weighted": p_t,
                    "share_control_weighted": p_c,
                    "share_pseudo_population": p_w,
                    "share_enrolled": p_e,
                    "smd_between_arms": float(smd_arms),
                    "smd_vs_enrolled": float(smd_enrolled),
                }
            )
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# The sweep
# --------------------------------------------------------------------------


def run(
    *,
    seed: int,
    registry: Registry | None = None,
    cohort_sizes: tuple[int, ...] = COHORT_SIZES,
    n_replicates: int = 200,
    shift: float = 1.0,
) -> pd.DataFrame:
    """One row per (estimator, reference, cohort size, replicate)."""
    registry = build_registry() if registry is None else registry
    rows: list[dict] = []
    for n_trial in cohort_sizes:
        for replicate in range(n_replicates):
            rng = np.random.default_rng([seed, n_trial, replicate, int(shift * 1000)])
            draw = simulate_study(n_trial, registry=registry, rng=rng, shift=shift)
            estimates = {name: fn(draw) for name, fn in ESTIMATORS.items()}
            for name, estimate in estimates.items():
                for ref_name in REFERENCES:
                    reference = draw.references[ref_name]
                    rows.append(
                        {
                            "estimator": name,
                            "target_population": TARGETS[name],
                            "reference": ref_name,
                            "n_trial": n_trial,
                            "replicate": replicate,
                            "shift": shift,
                            "estimate": estimate,
                            "reference_value": reference,
                            "residual": abs(estimate - reference),
                            "difference": estimate - reference,
                            "ratio": (
                                estimate / reference
                                if np.isfinite(reference) and reference != 0.0
                                else float("nan")
                            ),
                        }
                    )
    return pd.DataFrame(rows)


def summarise(runs: pd.DataFrame) -> pd.DataFrame:
    """Per (estimator, reference, cohort size): the curve, and the equality screen.

    ``n_excluded`` is reported beside every maximum. §blind is explicit that a
    residual maximum without its exclusion count is not a result.
    """
    rows = []
    for (estimator, reference, n_trial), group in runs.groupby(
        ["estimator", "reference", "n_trial"], sort=False
    ):
        valid = group[np.isfinite(group["residual"])]
        n_excluded = len(group) - len(valid)
        if valid.empty:
            rows.append(
                {
                    "estimator": estimator,
                    "target_population": TARGETS[estimator],
                    "reference": reference,
                    "n_trial": int(n_trial),
                    "n_valid": 0,
                    "n_excluded": n_excluded,
                    "max_residual": float("nan"),
                    "equality_flagged": False,
                    "bias": float("nan"),
                    "rmse": float("nan"),
                    "ratio_median": float("nan"),
                    "ratio_q25": float("nan"),
                    "ratio_q75": float("nan"),
                    "estimate_median": float("nan"),
                    "reference_median": float("nan"),
                }
            )
            continue
        max_residual = float(valid["residual"].max())
        rows.append(
            {
                "estimator": estimator,
                "target_population": TARGETS[estimator],
                "reference": reference,
                "n_trial": int(n_trial),
                "n_valid": int(len(valid)),
                "n_excluded": int(n_excluded),
                "max_residual": max_residual,
                "equality_flagged": bool(max_residual <= EQUALITY_TOL),
                "bias": float(valid["difference"].mean()),
                "rmse": float(np.sqrt(np.mean(np.square(valid["difference"])))),
                "ratio_median": float(valid["ratio"].median()),
                "ratio_q25": float(valid["ratio"].quantile(0.25)),
                "ratio_q75": float(valid["ratio"].quantile(0.75)),
                "estimate_median": float(valid["estimate"].median()),
                "reference_median": float(valid["reference_value"].median()),
            }
        )
    out = pd.DataFrame(rows)
    return out.sort_values(["reference", "estimator", "n_trial"]).reset_index(drop=True)


def run_shift_sweep(
    *,
    seed: int,
    registry: Registry | None = None,
    shifts: tuple[float, ...] = (0.0, 0.25, 0.5, 0.75, 1.0),
    n_trial: int = 1600,
    n_replicates: int = 200,
) -> pd.DataFrame:
    """The case-mix sweep. ``shift = 0`` is the experiment's own null.

    Falsifier F4 requires the alleged harm to VANISH at ``shift = 0``. A
    demonstration whose harm survives its own null is not measuring case-mix
    shift and is measuring something else.

    Falsifier F3d reads the smallest shift at which the MCID decision flips: a
    flip that needs the most extreme mix we wrote down and nothing less is a
    designed coincidence rather than a realistic one.
    """
    registry = build_registry() if registry is None else registry
    rows = []
    for shift in shifts:
        runs = run(
            seed=seed,
            registry=registry,
            cohort_sizes=(n_trial,),
            n_replicates=n_replicates,
            shift=shift,
        )
        for estimator in ESTIMATORS:
            sub = runs[runs["estimator"] == estimator]
            obs_pooled = sub[sub["reference"] == "obs-pooled"]
            po_trial = sub[sub["reference"] == "po-trial"]
            valid_obs = obs_pooled[np.isfinite(obs_pooled["residual"])]
            valid_po = po_trial[np.isfinite(po_trial["residual"])]
            if valid_po.empty:
                continue
            estimate_median = float(valid_po["estimate"].median())
            truth_median = float(valid_po["reference_value"].median())
            rows.append(
                {
                    "shift": shift,
                    "estimator": estimator,
                    "target_population": TARGETS[estimator],
                    "n_trial": n_trial,
                    "trial_mix": str(tuple(round(x, 4) for x in shifted_trial_mix(shift))),
                    "estimate_median": estimate_median,
                    "po_trial_median": truth_median,
                    "gap_vs_po_trial": estimate_median - truth_median,
                    "bias_vs_po_trial": float(valid_po["difference"].mean()),
                    "max_residual_vs_obs_pooled": (
                        float(valid_obs["residual"].max()) if not valid_obs.empty
                        else float("nan")
                    ),
                    "equality_flagged_vs_obs_pooled": bool(
                        not valid_obs.empty
                        and float(valid_obs["residual"].max()) <= EQUALITY_TOL
                    ),
                    "mcid": MCID,
                    "estimate_clears_mcid": bool(estimate_median > MCID),
                    "truth_clears_mcid": bool(truth_median > MCID),
                    "decision_flips": bool(
                        (estimate_median > MCID) != (truth_median > MCID)
                    ),
                }
            )
    return pd.DataFrame(rows)


def run_balance_sweep(
    *,
    seed: int,
    registry: Registry | None = None,
    n_trial: int = 1600,
    n_replicates: int = 50,
) -> pd.DataFrame:
    """Mean post-weighting SMDs over replicates. Falsifier F3e reads this."""
    registry = build_registry() if registry is None else registry
    frames = []
    for replicate in range(n_replicates):
        rng = np.random.default_rng([seed, n_trial, replicate, 1000])
        draw = simulate_study(n_trial, registry=registry, rng=rng, shift=1.0)
        table = balance_table(draw)
        if not table.empty:
            frames.append(table)
    if not frames:  # pragma: no cover - only if every draw fails positivity
        return pd.DataFrame()
    combined = pd.concat(frames, ignore_index=True)
    out = (
        combined.groupby(["weighting", "stratum", "stratum_name"])
        .agg(
            n_replicates=("smd_between_arms", "size"),
            mean_smd_between_arms=("smd_between_arms", "mean"),
            max_abs_smd_between_arms=(
                "smd_between_arms", lambda s: float(np.max(np.abs(s)))
            ),
            mean_smd_vs_enrolled=("smd_vs_enrolled", "mean"),
            max_abs_smd_vs_enrolled=(
                "smd_vs_enrolled", lambda s: float(np.max(np.abs(s)))
            ),
            mean_share_pseudo_population=("share_pseudo_population", "mean"),
            mean_share_enrolled=("share_enrolled", "mean"),
        )
        .reset_index()
    )
    out["smd_threshold"] = SMD_THRESHOLD
    out["between_arms_passes"] = out["max_abs_smd_between_arms"] <= SMD_THRESHOLD
    out["vs_enrolled_passes"] = out["max_abs_smd_vs_enrolled"] <= SMD_THRESHOLD
    return out


# --------------------------------------------------------------------------
# The falsifiers, evaluated as code rather than as prose
# --------------------------------------------------------------------------


def evaluate_falsifiers(
    primary: pd.DataFrame, shift: pd.DataFrame, balance: pd.DataFrame
) -> pd.DataFrame:
    """Prereg §8, each falsifier as a row with the number that decided it.

    Returns one row per falsifier with ``fired`` and the evidence. A falsifier
    that cannot be evaluated returns ``fired = True``, because an unevaluable
    falsifier is not a passed one.
    """
    rows: list[dict] = []

    def record(fid: str, statement: str, fired: bool, evidence: str) -> None:
        rows.append(
            {"falsifier": fid, "statement": statement, "fired": bool(fired),
             "evidence": evidence}
        )

    defective = ("ate-ipw", "ate-standardisation")
    obs_pooled = primary[primary["reference"] == "obs-pooled"]
    po_trial = primary[primary["reference"] == "po-trial"]

    # F1 --------------------------------------------------------------------
    d = obs_pooled[obs_pooled["estimator"].isin(defective)]
    a = bool((d["max_residual"] <= EQUALITY_TOL).all()) and not d.empty
    b = bool(d["ratio_median"].between(0.99, 1.01).all())
    ranked = obs_pooled[obs_pooled["estimator"].isin(PRIMARY_ESTIMATORS)]
    smallest = []
    for _n_trial, group in ranked.groupby("n_trial"):
        best = group.loc[group["rmse"].idxmin(), "estimator"]
        smallest.append(best in defective)
    c = bool(smallest) and all(smallest)
    record(
        "F1", "the validation does not pass, so there is nothing to be misled by",
        not (a and b and c),
        f"equality at every n={a}; ratio_median in [0.99,1.01]={b}; "
        f"defective arm has smallest RMSE at every n={c} "
        + (f"(max residual {d['max_residual'].max():.3g})" if not d.empty
           else "(no obs-pooled rows for the defective arm)"),
    )

    # F2 --------------------------------------------------------------------
    p = po_trial[po_trial["estimator"].isin(defective)]
    big = p[p["n_trial"] >= 400] if not p.empty else p
    small_bias = bool((big["bias"].abs() < 0.2).any()) if not big.empty else True
    if p.empty:
        at_large = at_small = float("nan")
        largest_n = float("nan")
    else:
        largest_n = float(p["n_trial"].max())
        at_large = float(p[p["n_trial"] == largest_n]["bias"].abs().mean())
        at_small = (
            float(p[p["n_trial"] == 100]["bias"].abs().mean())
            if (p["n_trial"] == 100).any() else float("nan")
        )
    shrinks = bool(np.isfinite(at_small) and at_large < 0.5 * at_small)
    record(
        "F2", "the honest reference buys nothing",
        p.empty or small_bias or shrinks,
        "no po-trial rows for the defective arm" if p.empty else (
            f"|bias| vs po-trial at n>=400 in "
            f"[{big['bias'].abs().min():.4f},{big['bias'].abs().max():.4f}]; "
            f"|bias| at n=100 {at_small:.4f} -> at n={largest_n:.0f} "
            f"{at_large:.4f}; shrinks={shrinks}"
        ),
    )

    # F3a -------------------------------------------------------------------
    std_row = obs_pooled[obs_pooled["estimator"] == "ate-standardisation"]
    std_eq = bool((std_row["max_residual"] <= EQUALITY_TOL).all()) and not std_row.empty
    std_bias_col = po_trial[po_trial["estimator"] == "ate-standardisation"]["bias"]
    std_bias = float(std_bias_col.abs().min()) if not std_bias_col.empty else float("nan")
    record(
        "F3a", "the defect is only reachable through the IPW weight slip",
        not (std_eq and np.isfinite(std_bias) and std_bias >= 0.2),
        f"ate-standardisation equality at every n={std_eq}; its smallest "
        f"|bias| vs po-trial {std_bias:.4f}",
    )

    # F3b is asserted in the test suite, not here; recorded for completeness.
    record(
        "F3b", "obs-pooled is not the functional the paper itself analyses",
        False,
        "asserted in tests/test_external_control_demo.py against "
        "trial_recovery._standardised_effect",
    )

    # F3c -------------------------------------------------------------------
    logit = obs_pooled[obs_pooled["estimator"] == "ate-ipw-logistic"]
    at_tol = bool((logit["max_residual"] <= EQUALITY_TOL).all()) and not logit.empty
    at_allclose = bool((logit["max_residual"] <= 1e-8).all()) and not logit.empty
    record(
        "F3c", "the equality is an artefact of hand-counted saturated weights",
        not (at_tol or at_allclose),
        "no fitted-propensity rows" if logit.empty else (
            f"fitted-propensity max residual {logit['max_residual'].max():.3g}; "
            f"flags at 1e-10={at_tol}; flags at allclose atol 1e-8={at_allclose}"
        ),
    )

    # F3d -------------------------------------------------------------------
    flips = (
        shift[(shift["estimator"] == "ate-ipw") & shift["decision_flips"].astype(bool)]
        if not shift.empty else shift
    )
    smallest_flip = float(flips["shift"].min()) if not flips.empty else float("nan")
    record(
        "F3d", "the decision only flips at the most extreme case-mix shift",
        not np.isfinite(smallest_flip) or smallest_flip >= 1.0,
        f"smallest shift at which the MCID decision flips: {smallest_flip}",
    )

    # F3e -------------------------------------------------------------------
    ate_balance = balance[balance["weighting"] == "ate"]
    worst_arms = (
        float(ate_balance["max_abs_smd_between_arms"].max())
        if not ate_balance.empty else float("nan")
    )
    worst_enrolled = (
        float(ate_balance["max_abs_smd_vs_enrolled"].max())
        if not ate_balance.empty else float("nan")
    )
    record(
        "F3e", "the routine post-weighting balance table already catches it",
        not np.isfinite(worst_arms) or worst_arms > SMD_THRESHOLD,
        f"worst |SMD| between arms under ATE weights {worst_arms:.3g} "
        f"(threshold {SMD_THRESHOLD}); worst |SMD| vs the ENROLLED population "
        f"{worst_enrolled:.3g}",
    )

    # F4 --------------------------------------------------------------------
    null = (
        shift[(shift["shift"] == 0.0) & (shift["estimator"] == "ate-ipw")]
        if not shift.empty else shift
    )
    null_bias = (
        float(null["bias_vs_po_trial"].abs().max()) if not null.empty else float("nan")
    )
    record(
        "F4", "the demonstration cannot fail: the harm survives its own null",
        not np.isfinite(null_bias) or null_bias > 0.1,
        f"|bias| of ate-ipw vs po-trial at shift=0: {null_bias:.4f} "
        f"(must be <= 0.1)",
    )

    # F5 --------------------------------------------------------------------
    induced = primary[
        (primary["reference"] == "induced-trial")
        & (primary["estimator"] == "att-standardisation")
    ]
    po = primary[
        (primary["reference"] == "po-trial")
        & (primary["estimator"] == "att-standardisation")
    ]
    gap = float(
        np.max(np.abs(induced["reference_median"].to_numpy()
                      - po["reference_median"].to_numpy()))
    ) if len(induced) == len(po) and len(induced) else float("nan")
    record(
        "F5", "the honest reference is secretly the parametric one",
        not np.isfinite(gap) or gap <= EQUALITY_TOL,
        f"max |median(po-trial) - median(induced-trial)| = {gap:.4g}",
    )

    return pd.DataFrame(rows)


def main(argv: Sequence[str] | None = None) -> int:
    """Write the demonstration tables under ``results/`` with a provenance stamp.

        python -m src.harness.external_control_demo
    """
    parser = argparse.ArgumentParser(
        description="Where an observed-data reference does harm: external control arm"
    )
    parser.add_argument("--seed", type=int, default=20260915)
    parser.add_argument("--replicates", type=int, default=200)
    parser.add_argument("--balance-replicates", type=int, default=50)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(message)s")
    registry = build_registry()
    log.info("registry pool: %d records, mix %s", len(registry),
             np.round(registry.records["stratum"].value_counts(normalize=True)
                      .sort_index().to_numpy(), 4))

    runs = run(seed=args.seed, registry=registry, n_replicates=args.replicates)
    primary = summarise(runs)
    shift = run_shift_sweep(
        seed=args.seed, registry=registry, n_replicates=args.replicates
    )
    balance = run_balance_sweep(
        seed=args.seed, registry=registry, n_replicates=args.balance_replicates
    )
    falsifiers = evaluate_falsifiers(primary, shift, balance)

    for _, row in falsifiers.iterrows():
        log.info("%-4s %-8s %s", row["falsifier"],
                 "FIRED" if row["fired"] else "not fired", row["evidence"])

    shared = {
        "SYNTHETIC": (
            "A fixed finite registry pool of simulated records, resampled. The "
            "plasmode mechanic is preserved; the provenance of the numbers is "
            "not real. Nothing here is a result about any real registry, trial "
            "or agent."
        ),
        "prereg": "docs/prereg_external_control_demo.md",
        "n_pool": N_POOL,
        "pool_mix": list(POOL_MIX),
        "trial_mix": list(TRIAL_MIX),
        "pool_means": list(POOL_MEANS),
        "pool_sd": POOL_SD,
        "induced_delta": list(DELTA),
        "responsiveness_gamma": [RESPONSIVENESS_SHAPE, RESPONSIVENESS_SCALE],
        "control_ratio": CONTROL_RATIO,
        "mcid": MCID,
        "pool_seed": POOL_SEED,
        "equality_tolerance": EQUALITY_TOL,
        "n_replicates": args.replicates,
    }

    paths = []
    paths.append(write_versioned_table(
        primary, "external_control_primary", seed=args.seed,
        results_dir=args.results_dir, allow_dirty=args.allow_dirty,
        extra_meta=shared | {
            "cohort_sizes": list(COHORT_SIZES),
            "estimators": dict(TARGETS),
            "references": list(REFERENCES),
            "what_this_answers": (
                "Whether validating an external-control estimator against an "
                "observed-data reference computed on the analysis file can lead "
                "a competent analyst to deploy an estimator that answers for the "
                "wrong target population. The defective arm reproduces "
                "obs-pooled exactly, so the validation cannot fail; against "
                "po-trial, which reads both potential outcomes of the enrolled "
                "records and which no estimator can reproduce, the same recovery "
                "check reveals the defect."
            ),
        },
    ))
    paths.append(write_versioned_table(
        shift, "external_control_shift", seed=args.seed,
        results_dir=args.results_dir, allow_dirty=args.allow_dirty,
        extra_meta=shared | {
            "what_this_answers": (
                "How the harm depends on the case-mix shift between the registry "
                "and the enrolled population, and the smallest shift at which the "
                "MCID decision flips. shift=0 is the experiment's own null: the "
                "harm must vanish there (falsifier F4)."
            ),
        },
    ))
    paths.append(write_versioned_table(
        balance, "external_control_balance", seed=args.seed,
        results_dir=args.results_dir, allow_dirty=args.allow_dirty,
        extra_meta=shared | {
            "balance_replicates": args.balance_replicates,
            "smd_threshold": SMD_THRESHOLD,
            "what_this_answers": (
                "Whether the post-weighting covariate balance table an "
                "external-control protocol actually runs would have caught the "
                "defect. It reports balance BETWEEN ARMS; the defect is in "
                "balance against the ENROLLED population, which it does not look "
                "at."
            ),
        },
    ))
    paths.append(write_versioned_table(
        falsifiers, "external_control_falsifiers", seed=args.seed,
        results_dir=args.results_dir, allow_dirty=args.allow_dirty,
        extra_meta=shared | {
            "what_this_answers": (
                "Prereg §8, evaluated as code. One row per falsifier, with the "
                "number that decided it."
            ),
        },
    ))
    for path in paths:
        log.info("wrote %s", path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
