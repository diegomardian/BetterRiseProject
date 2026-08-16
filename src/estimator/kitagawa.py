"""Kitagawa decomposition. W4.

Kitagawa (1955) demographic standardisation. NOT regression Oaxaca-Blinder,
which decomposes by regression coefficients and answers a different question.

    compositional = Δ(mature fraction) × normal per-cell mean
    intrinsic     = tumour mature fraction × Δ(per-cell mean)
    interaction   = Δ(mature fraction) × Δ(per-cell mean)

The split is not unique: normal-weighted and tumour-weighted give different
answers and the difference lives in the interaction term. Report both plus
doubly-robust. CLAUDE.md invariant 7: the interaction term is reported
separately and never folded into either arm.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
import pandas as pd

from src.harness.positivity import CUTPOINTS, Cutpoints, classify_estimability
from src.schema import WEIGHTINGS

Weighting = Literal["normal", "tumour", "doubly_robust"]


@dataclass(frozen=True)
class Decomposition:
    """One (patient, gene, rung, axis, weighting) result, pre-schema.

    ``intrinsic`` is None when the estimate is not identifiable. It is never
    0.0 in that case — CLAUDE.md invariant 1. ``src.schema.validate_results``
    will reject the frame if this is violated, but the honest value starts here.
    """

    compositional: float | None
    intrinsic: float | None
    interaction: float | None
    weighting: Weighting
    n_cells_mature: int


def decompose(
    frac_mature_normal: float,
    frac_mature_tumour: float,
    mean_normal: float,
    mean_tumour: float,
    *,
    n_cells_mature: int,
    weighting: Weighting = "normal",
) -> Decomposition:
    """Two-term Kitagawa split with the interaction reported separately.

    The identity holds exactly for the normal weighting::

        total = compositional + intrinsic + interaction

    Estimability is NOT decided here — call
    ``src.harness.classify_estimability(n_cells_mature)`` and null out the
    intrinsic term yourself, so the decision is visible at the call site rather
    than buried in the arithmetic.

    W4 — the scalar identity below is correct and unit-tested; the real work is
    the per-patient version over the AnnData, the doubly-robust reweighting
    (weeks 3-4) and the patient-level bootstrap (weeks 4-5).
    """
    d_frac = frac_mature_tumour - frac_mature_normal
    d_mean = mean_tumour - mean_normal

    if weighting == "normal":
        compositional = d_frac * mean_normal
        intrinsic = frac_mature_normal * d_mean
        # Sign convention: normal-weighted puts the cross term positive,
        # tumour-weighted has already absorbed it, so it is reported negated.
        interaction = d_frac * d_mean
    elif weighting == "tumour":
        compositional = d_frac * mean_tumour
        intrinsic = frac_mature_tumour * d_mean
        interaction = -d_frac * d_mean
    elif weighting == "doubly_robust":
        compositional, intrinsic, interaction = _doubly_robust_split(
            frac_mature_normal, frac_mature_tumour, mean_normal, mean_tumour
        )
    else:
        raise ValueError(f"unknown weighting {weighting!r}")

    return Decomposition(
        compositional=compositional,
        intrinsic=intrinsic,
        interaction=interaction,
        weighting=weighting,
        n_cells_mature=n_cells_mature,
    )


def _doubly_robust_split(
    frac_mature_normal: float,
    frac_mature_tumour: float,
    mean_normal: float,
    mean_tumour: float,
) -> tuple[float, float, float]:
    """Pooled-reference split. First-pass answer for weighting="doubly_robust".

    Kline (2011), "Oaxaca-Blinder as a Reweighting Estimator," AER P&P 101(3):
    532-537, shows that using the POOLED (average) reference weights instead of
    either group's own values is numerically equivalent to a reweighting
    estimator — which is where the "doubly robust" framing for this weighting
    comes from. It also has a property neither plain weighting has: summing
    compositional and intrinsic reproduces total exactly, with the interaction
    term at zero by construction (algebra below), so it does not depend on
    which side is picked as the reference — the asymmetry invariant 7 exists to
    police is removed rather than chosen away.

    Proof compositional + intrinsic == total, for pooled p=(p_N+p_T)/2 and
    m=(m_N+m_T)/2:
        (p_T-p_N)*m + p*(m_T-m_N)
      = [(p_T-p_N)(m_N+m_T) + (p_N+p_T)(m_T-m_N)] / 2
      = (2 p_T m_T - 2 p_N m_N) / 2 = p_T*m_T - p_N*m_N = total.

    This is a defensible, citable first cut, not the last word. A genuine
    covariate-level AIPW estimator (propensity model over cell-level
    covariates, paired with an outcome regression) needs cell-level Lee data
    to fit and validate — replace this once that data is in hand, and keep the
    "agreement with the plain version, quantified" comparison from the week
    3-4 deliverable in src/estimator/README.md.
    """
    pooled_frac = (frac_mature_normal + frac_mature_tumour) / 2.0
    pooled_mean = (mean_normal + mean_tumour) / 2.0
    d_frac = frac_mature_tumour - frac_mature_normal
    d_mean = mean_tumour - mean_normal
    compositional = d_frac * pooled_mean
    intrinsic = pooled_frac * d_mean
    interaction = 0.0
    return compositional, intrinsic, interaction


#: Columns a per-cohort summary frame must carry — one row per
#: (patient, study, gene, rung, axis), the key the schema itself uses minus
#: ``weighting`` (decompose_cohort fans each row out across all three).
_SUMMARY_COLUMNS = (
    "patient_id",
    "study_id",
    "gene",
    "granularity_rung",
    "labeling_axis",
    "frac_mature_normal",
    "frac_mature_tumour",
    "mean_normal",
    "mean_tumour",
    "n_cells_mature",
)


def decompose_cohort(summary: pd.DataFrame, *, cutpoints: Cutpoints = CUTPOINTS) -> pd.DataFrame:
    """Kitagawa split for a whole cohort, all three weightings, schema-shaped.

    ``summary`` has one row per (patient_id, study_id, gene, granularity_rung,
    labeling_axis) carrying the four scalar summary statistics ``decompose()``
    needs — computed upstream from labelled AnnData, one axis/rung/gene at a
    time. This function does the fan-out across weightings and the
    estimability call; it does not touch cells.

    Estimability is decided once per key, not once per weighting — the
    mature-cell count does not depend on which weighting is used to report the
    split, and this is the one and only place ``intrinsic`` gets nulled out
    (CLAUDE.md invariant 1; do not reimplement ``classify_estimability`` — see
    src/estimator/README.md). ``cutpoints`` defaults to W2's current
    ``CUTPOINTS`` and is threaded through explicitly, not looked up freshly per
    row, so a whole cohort run is reported against one calibration.

    Returns a frame with the schema's columns except ``ci_low``/``ci_high``,
    which are left as ``None`` — ``bootstrap_over_patients`` fills those in.
    Coerce with ``src.schema.coerce_results`` and merge in the bootstrap CIs
    before calling ``write_results``.
    """
    missing = [c for c in _SUMMARY_COLUMNS if c not in summary.columns]
    if missing:
        raise ValueError(f"summary is missing column(s): {missing}")

    rows: list[dict] = []
    for record in summary.to_dict("records"):
        n = int(record["n_cells_mature"])
        estimability = classify_estimability(n, cutpoints)
        for weighting in WEIGHTINGS:
            d = decompose(
                record["frac_mature_normal"],
                record["frac_mature_tumour"],
                record["mean_normal"],
                record["mean_tumour"],
                n_cells_mature=n,
                weighting=weighting,
            )
            intrinsic = None if estimability == "not_estimable" else d.intrinsic
            rows.append(
                {
                    "patient_id": record["patient_id"],
                    "study_id": record["study_id"],
                    "gene": record["gene"],
                    "granularity_rung": record["granularity_rung"],
                    "labeling_axis": record["labeling_axis"],
                    "n_cells_mature": n,
                    "compositional": d.compositional,
                    "intrinsic": intrinsic,
                    "interaction": d.interaction,
                    "weighting": weighting,
                    "estimability": estimability,
                    "ci_low": None,
                    "ci_high": None,
                }
            )
    return pd.DataFrame(rows)


_BOOTSTRAP_GROUP_COLS = ("study_id", "gene", "granularity_rung", "labeling_axis")
_BOOTSTRAP_TERMS = ("compositional", "intrinsic", "interaction")


def bootstrap_over_patients(
    summary: pd.DataFrame,
    *,
    n_boot: int = 2000,
    seed: int,
    cutpoints: Cutpoints = CUTPOINTS,
    alpha: float = 0.05,
) -> pd.DataFrame:
    """Patient-level bootstrap CIs. CLAUDE.md invariant 5 — patients, not cells.

    Resamples whole patients with replacement, independently within each
    (study_id, gene, granularity_rung, labeling_axis) group, and re-runs
    ``decompose_cohort`` on each resample, averaging across the resampled
    patients within the group on every draw. Resampling cells instead inflates
    n by roughly the number of cells per patient and produces intervals wrong
    by an order of magnitude, in the flattering direction — that mistake is
    why this function takes a patient-indexed ``summary`` frame and never
    touches a cell count directly.

    ``seed`` is required, not defaulted, and drives one
    ``numpy.random.default_rng`` used for every resample — global RNG state
    (``set_global_seeds``) is deliberately not used here; the docstring on
    that function says why.

    Returns one row per (study_id, gene, granularity_rung, labeling_axis,
    weighting, term) with a **cohort-level** ``ci_low``/``ci_high`` at the
    percentile bounds implied by ``alpha`` (default: a 95% interval). This is
    deliberately a separate frame from ``decompose_cohort``'s per-patient
    point estimates, not merged into it here — the schema has one
    ``ci_low``/``ci_high`` slot per patient-row, and which quantity belongs
    there (this cohort-level bootstrap band on ``intrinsic`, a per-patient
    interval from within-patient cell-count uncertainty, or a value read off
    the hierarchical model instead) is an open call for whoever owns the
    week 4-5 hierarchical-model deliverable to make deliberately. See
    ``attach_intrinsic_ci`` below for the one choice this module commits to.
    """
    missing = [c for c in _SUMMARY_COLUMNS if c not in summary.columns]
    if missing:
        raise ValueError(f"summary is missing column(s): {missing}")
    if n_boot < 1:
        raise ValueError(f"n_boot={n_boot} must be positive")

    rng = np.random.default_rng(seed)
    lo_q, hi_q = 100 * (alpha / 2), 100 * (1 - alpha / 2)

    records: list[dict] = []
    for key, group in summary.groupby(list(_BOOTSTRAP_GROUP_COLS)):
        patients = group["patient_id"].unique()
        n_patients = len(patients)
        draws = {f"{w}::{t}": [] for w in WEIGHTINGS for t in _BOOTSTRAP_TERMS}

        for _ in range(n_boot):
            resampled_ids = rng.choice(patients, size=n_patients, replace=True)
            resampled = pd.concat(
                [group[group["patient_id"] == pid] for pid in resampled_ids], ignore_index=True
            )
            result = decompose_cohort(resampled, cutpoints=cutpoints)
            for weighting in WEIGHTINGS:
                sub = result[result["weighting"] == weighting]
                for term in _BOOTSTRAP_TERMS:
                    values = sub[term].dropna()
                    draws[f"{weighting}::{term}"].append(
                        float(values.mean()) if len(values) else np.nan
                    )

        study_id, gene, rung, axis = key
        for weighting in WEIGHTINGS:
            for term in _BOOTSTRAP_TERMS:
                values = np.asarray(draws[f"{weighting}::{term}"], dtype=float)
                values = values[~np.isnan(values)]
                ci_low, ci_high = (
                    (None, None)
                    if len(values) == 0
                    else (float(np.percentile(values, lo_q)), float(np.percentile(values, hi_q)))
                )
                records.append(
                    {
                        "study_id": study_id,
                        "gene": gene,
                        "granularity_rung": rung,
                        "labeling_axis": axis,
                        "weighting": weighting,
                        "term": term,
                        "ci_low": ci_low,
                        "ci_high": ci_high,
                        "n_boot": n_boot,
                    }
                )
    return pd.DataFrame(records)


def attach_intrinsic_ci(point: pd.DataFrame, ci: pd.DataFrame) -> pd.DataFrame:
    """Fill the schema's single ci_low/ci_high slot from the intrinsic band.

    ``point`` is ``decompose_cohort``'s output (per patient); ``ci`` is
    ``bootstrap_over_patients``'s output (per cohort/gene/rung/axis/weighting,
    one row per term). The schema carries one CI pair per row but three point
    estimates (compositional, intrinsic, interaction) — this function commits
    to reporting the ``intrinsic`` term's cohort-level bootstrap band, since
    that is the estimand under identifiability scrutiny (estimability is
    defined for intrinsic, not for compositional or interaction). Every
    patient row within a (study, gene, rung, axis, weighting) group carries
    the same cohort-level band; this is a population interval broadcast onto
    each patient row, not a per-patient interval — say so if you present it.
    """
    intrinsic_cols = [*_BOOTSTRAP_GROUP_COLS, "weighting", "ci_low", "ci_high"]
    intrinsic_ci = ci.loc[ci["term"] == "intrinsic", intrinsic_cols]
    merged = point.drop(columns=["ci_low", "ci_high"]).merge(
        intrinsic_ci, on=[*_BOOTSTRAP_GROUP_COLS, "weighting"], how="left"
    )
    return merged
