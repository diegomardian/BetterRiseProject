"""W3.7 — baseline survival models. Clinical covariates only.

**This is a sanity check on the pipeline, not a scientific result.** The brief is
blunt about what it is for:

    "If stage does not come out prognostic, something is wrong upstream — wrong
     endpoint, bad join, censoring bug. Investigate rather than reporting the
     result."

So the module's real output is not a hazard ratio, it is a verdict:
:func:`stage_sanity_check` asks whether stage IV carries a higher hazard than
stage I and whether stage as a whole is significant. Everything else is
supporting evidence for that one question.

No project variables. No expression. Those arrive in Stage 4.

WHAT IS PRE-SPECIFIED ELSEWHERE
-------------------------------
Every modelling choice — which covariates, which endpoints, strata, tie
handling, missing-data policy — comes from ``config/covariate_set.yaml`` through
``src/bulk/covariates.py``. Nothing here hard-codes a formula, and
:func:`fit_cox` calls ``require_locked`` before fitting anything, so a model
cannot be run against a covariate set that is still under discussion.

ON THE PROPORTIONAL-HAZARDS ASSUMPTION
--------------------------------------
Checked per term with Schoenfeld residuals and reported whether or not it holds.
A violated PH assumption is a finding — most often that an effect is
time-varying — and not a reason to quietly switch models until one behaves.
"""

from __future__ import annotations

import warnings
from typing import Any

import numpy as np
import pandas as pd

from src.bulk.covariates import (
    covariate_names,
    require_locked,
    total_df,
)

#: Schoenfeld test p below this flags a term as violating proportional hazards.
PH_ALPHA = 0.05

#: Conventional floor for events per degree of freedom.
EVENTS_PER_DF_FLOOR = 10


class SurvivalError(RuntimeError):
    """A survival model could not be fitted or interpreted."""


def _design_matrix(
    design: pd.DataFrame, names: list[str], spec: dict[str, Any]
) -> pd.DataFrame:
    """One-hot the categoricals at the reference levels the config names.

    Reference levels are pre-specified so a hazard ratio always means the same
    comparison. Letting pandas pick alphabetically would make stage IV's
    coefficient silently relative to whichever level sorted first.
    """
    by_name = {c["name"]: c for c in spec["covariates"]}
    frame = pd.DataFrame(index=design.index)
    for name in names:
        covariate = by_name[name]
        if covariate["type"] == "continuous":
            frame[name] = pd.to_numeric(design[name], errors="coerce")
            continue
        levels = list(covariate["levels"])
        reference = covariate["reference"]
        present = set(design[name].dropna().unique())
        unexpected = present - set(levels)
        if unexpected:
            raise SurvivalError(
                f"{name} contains level(s) {sorted(unexpected)} not in the locked set "
                f"{levels}. Either the clinical table drifted or the config is stale."
            )
        for level in levels:
            if level == reference:
                continue
            frame[f"{name}[{level}]"] = (design[name] == level).astype(float)
    return frame


def fit_cox(
    design: pd.DataFrame,
    spec: dict[str, Any],
    *,
    endpoint: str,
    context: str = "clinical_baseline",
) -> tuple[Any, pd.DataFrame]:
    """Fit a stratified Cox model. Returns (fitter, tidy coefficient table).

    Stratified by whatever ``model.strata`` names — project, so COAD and READ
    get separate baseline hazards. That is the config's answer to the
    pool-or-stratify question and it costs no degrees of freedom.
    """
    from lifelines import CoxPHFitter

    require_locked(spec)
    names = covariate_names(spec, endpoint=endpoint, context=context)
    strata = list(spec["model"].get("strata") or [])

    frame = _design_matrix(design, names, spec)
    frame[f"{endpoint}.time"] = pd.to_numeric(design[f"{endpoint}.time"])
    frame[endpoint] = pd.to_numeric(design[endpoint]).astype(int)
    for stratum in strata:
        frame[stratum] = design[stratum].to_numpy()

    fitter = CoxPHFitter()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fitter.fit(
            frame,
            duration_col=f"{endpoint}.time",
            event_col=endpoint,
            strata=strata or None,
        )

    summary = fitter.summary.reset_index().rename(columns={"covariate": "term"})
    tidy = pd.DataFrame(
        {
            "term": summary["term"],
            "hazard_ratio": summary["exp(coef)"].round(4),
            "ci_low": summary["exp(coef) lower 95%"].round(4),
            "ci_high": summary["exp(coef) upper 95%"].round(4),
            "p": summary["p"].round(5),
        }
    )
    tidy["endpoint"] = endpoint
    tidy["context"] = context
    tidy["n"] = int(len(frame))
    tidy["n_events"] = int(frame[endpoint].sum())
    return fitter, tidy


def proportional_hazards_check(
    fitter: Any, design: pd.DataFrame, spec: dict[str, Any], *, endpoint: str, context: str
) -> pd.DataFrame:
    """Schoenfeld residual test per term. Reported whether or not it passes."""
    names = covariate_names(spec, endpoint=endpoint, context=context)
    strata = list(spec["model"].get("strata") or [])

    frame = _design_matrix(design, names, spec)
    frame[f"{endpoint}.time"] = pd.to_numeric(design[f"{endpoint}.time"])
    frame[endpoint] = pd.to_numeric(design[endpoint]).astype(int)
    for stratum in strata:
        frame[stratum] = design[stratum].to_numpy()

    # NOT fitter.check_assumptions: that returns matplotlib axes, so parsing it
    # for p-values silently yields nothing and every term reads as passing.
    from lifelines.statistics import proportional_hazard_test

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = proportional_hazard_test(fitter, frame, time_transform=["km", "rank"])

    table = result.summary.reset_index()
    term_col, transform_col = table.columns[0], table.columns[1]
    rows = [
        {
            "endpoint": endpoint,
            "term": str(row[term_col]),
            "time_transform": str(row[transform_col]),
            "test_statistic": round(float(row["test_statistic"]), 4),
            "p": round(float(row["p"]), 5),
            "violates_ph": bool(float(row["p"]) < PH_ALPHA),
        }
        for _, row in table.iterrows()
    ]
    return pd.DataFrame(rows)


def stage_sanity_check(tidy: pd.DataFrame, *, endpoint: str) -> dict[str, Any]:
    """Did the pipeline reproduce the single best-known fact in colorectal cancer?

    Stage IV must carry a substantially higher hazard than stage I. If it does
    not, the brief says to investigate upstream rather than report — so this
    returns a verdict, not a p-value to interpret generously.
    """
    stage_terms = tidy.loc[tidy["term"].str.startswith("stage[")]
    if stage_terms.empty:
        raise SurvivalError("no stage terms in the model — stage is the sanity check")

    four = stage_terms.loc[stage_terms["term"] == "stage[IV]"]
    if four.empty:
        raise SurvivalError("stage[IV] missing; cannot run the sanity check")
    hr = float(four["hazard_ratio"].iloc[0])
    p = float(four["p"].iloc[0])
    monotone = _is_monotone(stage_terms)

    passed = hr > 1.0 and p < 0.05
    return {
        "endpoint": endpoint,
        "stage_IV_hazard_ratio": round(hr, 3),
        "stage_IV_ci_low": float(four["ci_low"].iloc[0]),
        "stage_IV_ci_high": float(four["ci_high"].iloc[0]),
        "stage_IV_p": p,
        "hazard_increases_with_stage": monotone,
        "verdict": "stage is prognostic" if passed else "STAGE NOT PROGNOSTIC — investigate",
        "passed": passed,
    }


def _is_monotone(stage_terms: pd.DataFrame) -> bool:
    """Do hazard ratios rise across stage II -> III -> IV? Informative, not required.

    Stage II versus III often overlaps in TCGA, so non-monotonicity is common
    and is not on its own a pipeline fault. Reported so a reader can see it.
    """
    order = ["stage[II]", "stage[III]", "stage[IV]"]
    values = [
        float(stage_terms.loc[stage_terms["term"] == term, "hazard_ratio"].iloc[0])
        for term in order
        if (stage_terms["term"] == term).any()
    ]
    return all(b >= a for a, b in zip(values, values[1:], strict=False))


def kaplan_meier_by_stage(
    design: pd.DataFrame, *, endpoint: str
) -> tuple[pd.DataFrame, float]:
    """Per-stage survival curves plus the multivariate log-rank p."""
    from lifelines import KaplanMeierFitter
    from lifelines.statistics import multivariate_logrank_test

    rows = []
    for stage, group in design.groupby("stage", observed=True):
        fitter = KaplanMeierFitter()
        fitter.fit(group[f"{endpoint}.time"], group[endpoint], label=str(stage))
        median = fitter.median_survival_time_
        rows.append(
            {
                "endpoint": endpoint,
                "stage": str(stage),
                "n": int(len(group)),
                "n_events": int((group[endpoint] == 1).sum()),
                "median_survival_days": None if not np.isfinite(median) else float(median),
                "survival_at_1y": round(float(fitter.predict(365.0)), 4),
                "survival_at_3y": round(float(fitter.predict(1095.0)), 4),
            }
        )
    test = multivariate_logrank_test(
        design[f"{endpoint}.time"], design["stage"], design[endpoint]
    )
    return pd.DataFrame(rows).sort_values("stage").reset_index(drop=True), float(test.p_value)


def power_note(design: pd.DataFrame, spec: dict[str, Any], endpoint: str, context: str) -> str:
    """One line on whether this model had enough events to be believed."""
    names = covariate_names(spec, endpoint=endpoint, context=context)
    df = total_df(spec, names)
    events = int((design[endpoint] == 1).sum())
    ratio = events / df if df else float("nan")
    verdict = "adequate" if ratio >= EVENTS_PER_DF_FLOOR else "UNDERPOWERED"
    return f"{events} events / {df} df = {ratio:.1f} per df — {verdict}"
