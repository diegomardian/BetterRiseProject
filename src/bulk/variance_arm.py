"""Stage 4's variance question, exactly as the locked pre-specification states it.

    How much variance in bulk GUCA2A and CDX2 is explained by mature-colonocyte
    fraction alone?

Never as a raw R-squared, and never as two raw R-squared values compared between
genes. Issue #54 showed why: R-squared is a share of variance, so a gene near the
assay floor loses it to measurement noise whatever its biology, and under PURE
composition the 1.4%-of-normal gene returns 0.891 -> 0.000 as the floor rises
while the 94.7% gene holds at ~0.86. Both arms of the original prediction
satisfied with nothing intrinsic in the data. The statistic is therefore a
PERCENTILE within an abundance-matched null, and the null is drawn by the rule
committed in the spec.

THE ONE THING THAT MAKES THE PERCENTILE MEAN ANYTHING is that the target gene
and every null gene go through the *same* model. A target adjusted for purity
and plate, compared against a null that was not, measures the covariates. So
there is one function, :func:`gene_r_squared`, and both sides call it -- not two
code paths that agree by review.

WHERE THE LOCKED SPEC IS SILENT, AND WHAT IS DONE ABOUT IT. `variance_reported`
asks for both the partial R-squared (adjusted, "the honest quantity") and the
marginal one ("the one §6.2's sentence is about"), and the prediction arms say
"R-squared" without saying which. That is a real gap in the lock. Both are
computed, both get a percentile, and both verdicts are reported. If they
disagree, that is reported as a disagreement -- picking one after seeing them is
the exact move the pre-specification exists to prevent.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)


class VarianceArmError(RuntimeError):
    """The variance arm cannot be run on these inputs."""


#: Plate is `strata_or_random_effect` in the covariate lock. For OLS it enters as
#: fixed dummies: the lock's "29 levels over 624 samples is too many for fixed
#: effects" is an events-per-degree-of-freedom argument about Cox models, and a
#: linear model on 624 samples has ~590 residual df left after them. Recorded
#: here and in the sidecar rather than decided quietly.
PLATE_AS = "fixed_dummies"


@dataclass
class Attrition:
    """Who was dropped and why. Invariant 1's reporting half."""

    rows: list[dict] = field(default_factory=list)

    def record(self, stage: str, before: int, after: int, reason: str) -> None:
        self.rows.append({
            "stage": stage, "n_before": before, "n_after": after,
            "n_dropped": before - after, "reason": reason,
        })

    def to_frame(self) -> pd.DataFrame:
        return pd.DataFrame(self.rows)


def build_design(
    fractions: pd.DataFrame,
    covariates: pd.DataFrame | None,
    *,
    covariate_names: list[str],
    attrition: Attrition,
) -> pd.DataFrame:
    """Predictor plus covariates, one row per sample, complete cases only.

    A patient whose mature fraction is not estimable contributes `None` and is
    DROPPED, never entered as a zero fraction. The locked spec says so and the
    reason is that entering them as zero manufactures exactly the compositional
    signal the analysis is testing for.
    """
    before = len(fractions)
    design = fractions[["sample_id", "mature_colonocyte_fraction"]].copy()

    estimable = design["mature_colonocyte_fraction"].notna()
    design = design[estimable]
    attrition.record(
        "estimability", before, len(design),
        "mature_colonocyte_fraction is None (not estimable). Invariant 1: not "
        "entered as 0.0, which would manufacture the compositional signal.",
    )

    if covariates is not None and covariate_names:
        missing = [c for c in covariate_names if c not in covariates.columns]
        if missing:
            raise VarianceArmError(
                f"the covariate frame is missing {missing}. The locked spec's "
                f"`expression_models` context requires them; running without is "
                f"a different model than the one pre-specified."
            )
        before_join = len(design)
        design = design.merge(
            covariates[["sample_id", *covariate_names]], on="sample_id", how="inner"
        )
        attrition.record(
            "covariate_join", before_join, len(design),
            "no covariate record for this sample",
        )
        before_complete = len(design)
        design = design.dropna()
        attrition.record(
            "complete_case", before_complete, len(design),
            "missing a covariate value; OLS is complete-case",
        )
    if len(design) < 20:
        raise VarianceArmError(
            f"only {len(design)} samples survive the design. An R-squared on "
            f"this many is not interpretable and its matched-null percentile "
            f"less so."
        )
    return design


def _model_matrix(design: pd.DataFrame, covariate_names: list[str]) -> np.ndarray:
    """Covariates as a numeric matrix with an intercept. Categoricals -> dummies."""
    if not covariate_names:
        return np.ones((len(design), 1))
    frame = design[covariate_names]
    encoded = pd.get_dummies(frame, drop_first=True, dummy_na=False)
    matrix = np.column_stack([np.ones(len(design)), encoded.to_numpy(dtype=float)])
    return matrix


def _r_squared(y: np.ndarray, x: np.ndarray) -> float:
    """OLS R-squared of y on x (x already carries an intercept column)."""
    coefs, *_ = np.linalg.lstsq(x, y, rcond=None)
    resid = y - x @ coefs
    ss_res = float(resid @ resid)
    centred = y - y.mean()
    ss_tot = float(centred @ centred)
    if ss_tot <= 0:
        return float("nan")
    return 1.0 - ss_res / ss_tot


@dataclass(frozen=True)
class GeneFit:
    """One gene's two R-squared values, on one design."""

    gene: str
    n: int
    marginal_r2: float
    partial_r2: float
    covariate_r2: float

    def value(self, kind: str) -> float:
        if kind == "marginal":
            return self.marginal_r2
        if kind == "partial":
            return self.partial_r2
        raise VarianceArmError(f"unknown R-squared kind {kind!r}")


def gene_r_squared(
    expression: pd.Series,
    design: pd.DataFrame,
    covariate_names: list[str],
) -> GeneFit:
    """Marginal and partial R-squared for one gene. Target and null both use this.

    ``marginal`` is the gene on fraction alone -- the quantity §6.2's sentence
    asks about. ``partial`` is the fraction's share over and above the
    covariates, ``(R2_full - R2_cov) / (1 - R2_cov)``, which is the adjusted one.

    Both sides of the comparison call this same function, which is the only
    reason the percentile is a comparison at all.
    """
    y = expression.to_numpy(dtype=float)
    fraction = design["mature_colonocyte_fraction"].to_numpy(dtype=float)
    intercept = np.ones(len(design))

    marginal = _r_squared(y, np.column_stack([intercept, fraction]))
    covariate_matrix = _model_matrix(design, covariate_names)
    r2_cov = _r_squared(y, covariate_matrix)
    r2_full = _r_squared(y, np.column_stack([covariate_matrix, fraction]))
    partial = (
        float("nan") if not np.isfinite(r2_cov) or r2_cov >= 1.0
        else (r2_full - r2_cov) / (1.0 - r2_cov)
    )
    return GeneFit(
        gene=str(expression.name), n=len(design), marginal_r2=marginal,
        partial_r2=partial, covariate_r2=r2_cov,
    )


@dataclass(frozen=True)
class NullComparison:
    """A gene's R-squared against the null drawn at its own abundance."""

    gene: str
    kind: str                  # "marginal" | "partial"
    r2: float
    n_null: int
    null_median: float
    null_p05: float
    null_p95: float
    percentile: float
    excess: float
    exceeds_null: bool         # percentile > 0.95, the spec's threshold


def compare_to_null(fit: GeneFit, null_fits: list[GeneFit], kind: str) -> NullComparison:
    """Where this gene's R-squared sits in its abundance-matched null."""
    values = np.array([f.value(kind) for f in null_fits], dtype=float)
    values = values[np.isfinite(values)]
    if values.size < 20:
        raise VarianceArmError(
            f"only {values.size} null genes returned a finite {kind} R-squared "
            f"for {fit.gene}. A percentile within a null this small is not one."
        )
    r2 = fit.value(kind)
    percentile = float((values <= r2).mean())
    return NullComparison(
        gene=fit.gene, kind=kind, r2=r2, n_null=int(values.size),
        null_median=float(np.median(values)),
        null_p05=float(np.percentile(values, 5)),
        null_p95=float(np.percentile(values, 95)),
        percentile=percentile, excess=r2 - float(np.median(values)),
        exceeds_null=bool(percentile > 0.95),
    )


# ---------------------------------------------------------------------------
# The pre-registered verdicts


def primary_verdict(guca2a: NullComparison, cdx2: NullComparison) -> tuple[str, str]:
    """The locked primary arm, and its own `disconfirmed_if`, applied verbatim.

        confirmed     GUCA2A percentile <= 0.95 AND CDX2 percentile > 0.95
        disconfirmed  GUCA2A > 0.95 AND CDX2 <= 0.95
        indeterminate both, or neither -- "reported in those words rather than
                      resolved toward either arm"

    The indeterminate branch is why total instrument failure does not read as
    confirmation: a constant predictor gives both genes an R-squared of zero,
    neither exceeds, and this returns indeterminate. It still needs
    `check_predictor` to say the number came from a broken instrument rather
    than from the biology -- indeterminate for the right reason and
    indeterminate for the wrong one are the same string here.
    """
    if guca2a.kind != cdx2.kind:
        raise VarianceArmError("cannot compare a marginal to a partial R-squared")
    if not guca2a.exceeds_null and cdx2.exceeds_null:
        return "confirmed", (
            "GUCA2A falls within its abundance-matched null while CDX2 rises "
            "above its own, in the direction PR #49 predicts."
        )
    if guca2a.exceeds_null and not cdx2.exceeds_null:
        return "disconfirmed", (
            "GUCA2A exceeds its matched null while CDX2 does not. That is the "
            "compositional reading and it contradicts #49."
        )
    both = "both" if guca2a.exceeds_null else "neither"
    return "indeterminate", (
        f"{both} gene exceeds its abundance-matched null. The locked spec says "
        f"this is reported in those words rather than resolved toward either arm."
    )


def secondary_verdict(guca2a: NullComparison, cdx2: NullComparison) -> tuple[str, str]:
    """CDX2's excess over its null median greater than GUCA2A's, or not."""
    if guca2a.excess >= cdx2.excess:
        return "disconfirmed", (
            f"GUCA2A's excess ({guca2a.excess:+.4f}) is at least CDX2's "
            f"({cdx2.excess:+.4f}), which the locked spec names as disconfirming."
        )
    return "confirmed", (
        f"CDX2's excess ({cdx2.excess:+.4f}) exceeds GUCA2A's "
        f"({guca2a.excess:+.4f})."
    )


def negative_control_verdict(
    housekeeping: list[GeneFit], null_median: float, kind: str,
) -> tuple[str, str]:
    """Both negative controls from the locked spec, reported together.

    High-abundance: R-squared(ACTB, GAPDH ~ fraction) < 0.10, or every
    R-squared is an upper bound. Low-abundance: the matched-null median at
    GUCA2A's level < 0.10, or fraction predicts floor expression generally and
    the primary is indeterminate rather than confirmed.
    """
    breached = [f for f in housekeeping if np.isfinite(f.value(kind)) and f.value(kind) >= 0.10]
    floor_breached = np.isfinite(null_median) and null_median >= 0.10
    if not breached and not floor_breached:
        return "clean", "both negative controls are below 0.10"
    parts = []
    if breached:
        parts.append(
            "housekeeping tracks the fraction ("
            + ", ".join(f"{f.gene} R2={f.value(kind):.3f}" for f in breached)
            + "): the model is picking up library size, purity or composition "
              "structure, and every R-squared here is an UPPER BOUND"
        )
    if floor_breached:
        parts.append(
            f"the matched-null median is {null_median:.3f} >= 0.10: fraction "
            f"predicts low-abundance expression in general, most likely a "
            f"detection-rate artefact, so the primary is INDETERMINATE rather "
            f"than confirmed"
        )
    return "breached", "; ".join(parts)


def benjamini_hochberg(p_values: dict[str, float]) -> dict[str, float]:
    """BH-adjusted values, as the locked spec's `multiplicity` asks.

    The percentile is a permutation p-value in disguise: a gene at percentile
    0.97 has 3% of its own null above it. So `p = 1 - percentile`, adjusted
    across the outcome genes within a cohort and rung.
    """
    if not p_values:
        return {}
    names = list(p_values)
    values = np.array([p_values[n] for n in names], dtype=float)
    order = np.argsort(values)
    ranked = values[order]
    m = len(ranked)
    adjusted = ranked * m / np.arange(1, m + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1].clip(0.0, 1.0)
    out = dict(zip([names[i] for i in order], adjusted, strict=True))
    return {name: float(out[name]) for name in names}


def resolve_r_squared_kinds(spec: dict[str, Any]) -> tuple[str, ...]:
    """Which R-squared the verdicts are stated on. Both, because the lock is silent.

    `variance_reported` asks for both and the prediction arms say "R-squared"
    without saying which. Choosing one after seeing them is precisely what the
    pre-specification prevents, so both are carried through to the verdict and a
    disagreement is reported as a disagreement.
    """
    reported = str(spec["model"]["variance_reported"]).lower()
    kinds = tuple(k for k in ("partial", "marginal") if k in reported)
    if not kinds:
        raise VarianceArmError(
            "the locked spec's `variance_reported` names neither a partial nor "
            "a marginal R-squared, so there is nothing to compute."
        )
    return kinds
