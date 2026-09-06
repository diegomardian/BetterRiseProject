"""Detection rates are not comparable across genes. This is the scale that is.

WHY THIS EXISTS. ``coexpression_silencing`` scores a per-cell DETECTION rate --
the share of cells carrying >= 1 UMI -- and that choice is well argued where it
is made: detection survives a depth difference once the arms are matched, and a
mean over mostly-zero counts is dominated by the few cells that fired. All of
that is about comparing ONE gene between two arms, and for that it is right.

``specificity()`` then did something the statistic does not support: it compared
the detection DELTAS OF DIFFERENT GENES to each other. Those genes sit at
baseline detection rates from 0.36 to 0.98 on the adenoma cohort. A proportion
at 0.98 has fifteen times less variance -- and correspondingly less room to move
-- than one at 0.44, so a delta of -0.02 and a delta of -0.17 can be the same
underlying fold change. Ranking genes by raw delta ranks them substantially by
abundance.

This is the repository's own signature defect, one layer up. ``SATURATION_CEILING``
already exists in ``coexpression_silencing`` because a control pinned at 1.00
"cannot fall", and ``premise_holds`` switches such a control to log2 expression.
The specificity table never switched, and its guard could not have noticed: the
test fixture synthesised ``delta_detect`` directly, with no baseline rate
anywhere in it, so a gene at 0.98 and a gene at 0.44 were the same object to it.

THE FIX. Model detection as Poisson thinning -- a cell is detected when at least
one of its ``mu`` expected UMIs lands, so ``p = 1 - exp(-mu)``. Then

    cloglog(p) = log(-log(1 - p)) = log(mu)

and the between-arm difference is ``log(mu_T / mu_N)``: a LOG FOLD CHANGE, on the
same scale for every gene whatever its abundance. Under a uniform down-regulation
of the whole programme this difference is identical across genes, which is
exactly the null the raw deltas cannot express.

It is a monotone transform of the statistic the project already chose, so it
changes the scale without changing the measurement. ``log2_cp10k_ratio`` -- which
is already on every row -- is the independent corroborating statistic: it reaches
the same quantity by a different route, and where the two agree the conclusion
does not rest on this model being right.

MEASURED AGREEMENT, on the 377 adenoma rows carrying both: Pearson r = 0.851
against ``ln(CP10K fold change)``. The slope is 0.62 rather than the 1.0 the
Poisson model predicts, because real counts are overdispersed and detection
therefore under-responds. That attenuation is why this is reported BESIDE the
log2 statistic and not instead of it.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

#: Detection is undefined for cloglog at exactly 0 and exactly 1: ``log(0)`` and
#: ``log(-log(0))``. THE BOUNDARY RULE, stated once and applied everywhere.
#:
#: A rate is turned back into its cell counts, given the Jeffreys-style
#: correction ``(k + 1/2) / (n + 1)``, and transformed. So the correction is a
#: function of how many cells the estimate rests on, not a fixed epsilon -- a
#: patient with 15 mature cells and a patient with 500 are pulled toward the
#: middle by very different amounts, which is the behaviour wanted. A fixed
#: epsilon would let a 15-cell zero and a 500-cell zero claim the same evidence.
#:
#: The correction is what makes the transform total. ``MIN_CELLS_PER_ARM`` upstream
#: keeps ``n`` at 15 or more, so the induced bias is bounded and small; it is
#: pinned by ``test_the_boundary_rule_is_a_function_of_cell_count``.
BOUNDARY_PSEUDOCOUNT = 0.5
BOUNDARY_DENOMINATOR_PAD = 1.0

#: How much of the cross-gene spread in detection deltas a single uniform
#: thinning may explain before the ordering stops being readable as biology.
#:
#: This is a ceiling on a NULL's explanatory power, so a high value is the bad
#: outcome: if one factor with no gene-specificity in it reproduces the pattern,
#: the pattern is not evidence of gene-specific anything. Set at 0.50 -- half
#: the cross-gene variance -- because below that the null is leaving the
#: majority of the structure unexplained.
UNIFORM_THINNING_R2_CEILING = 0.50

#: Baseline rates must actually differ before the null can be fitted. If every
#: gene sits at the same abundance the null has nothing to discriminate on and
#: any R² it reports is an artefact of six numbers and one free parameter.
#: Below this spread the diagnostic ABSTAINS rather than reporting clean, which
#: is the failure mode this repository is named after.
MIN_BASELINE_SPREAD = 0.05


class ScaleError(ValueError):
    """A detection rate cannot be placed on the log fold-change scale."""


def cloglog_rate(p, n):
    """``log(mu)`` for a detection rate ``p`` resting on ``n`` cells.

    Inverts the Poisson detection model ``p = 1 - exp(-mu)``. The returned value
    is a log expected-counts-per-cell, so a DIFFERENCE of two of them is a log
    fold change and is comparable across genes of any abundance.

    Applies the boundary rule above, which is what makes this total on [0, 1]
    rather than only on the open interval.
    """
    p = np.asarray(p, dtype=float)
    n = np.asarray(n, dtype=float)
    if np.any(n < 1):
        raise ScaleError("a detection rate resting on fewer than one cell")
    if np.any((p < 0) | (p > 1)):
        raise ScaleError("detection rates must lie in [0, 1]")
    k = np.round(p * n)
    adjusted = (k + BOUNDARY_PSEUDOCOUNT) / (n + BOUNDARY_DENOMINATOR_PAD)
    return np.log(-np.log1p(-adjusted))


def delta_cloglog(frame: pd.DataFrame) -> pd.Series:
    """Per-row log fold change, from the two arms' detection rates and counts.

    Needs ``detect_normal``/``n_normal`` and ``detect_tumour``/``n_tumour``. The
    cell counts are not optional: they are what the boundary rule is a function
    of, and taking them as read is how a 15-cell arm gets treated like a
    500-cell one.
    """
    needed = {"detect_normal", "n_normal", "detect_tumour", "n_tumour"}
    missing = needed - set(frame.columns)
    if missing:
        raise ScaleError(
            f"cannot place detection on the log fold-change scale without "
            f"{sorted(missing)}. The cell counts are what the boundary rule "
            f"scales with; substituting a fixed epsilon would let a 15-cell "
            f"arm claim the evidence of a 500-cell one."
        )
    return pd.Series(
        cloglog_rate(frame["detect_tumour"], frame["n_tumour"])
        - cloglog_rate(frame["detect_normal"], frame["n_normal"]),
        index=frame.index,
    )


def _predicted_delta(baseline: np.ndarray, factor: float) -> float:
    """Mean detection delta a common thinning by ``factor`` would produce."""
    return float(np.mean((1.0 - (1.0 - baseline) ** factor) - baseline))


def uniform_thinning_null(
    deltas: pd.DataFrame, *, grid: np.ndarray | None = None
) -> tuple[pd.DataFrame, dict]:
    """Could ONE common fold change, with no gene-specificity, produce this?

    THE DIAGNOSTIC THE SPECIFICITY TABLE WAS MISSING. A gradient of detection
    deltas across genes looks like a graded biological effect and is also what a
    single uniform down-regulation produces, because the detection curve's
    sensitivity ``dp/dlog(mu) = (1-p)·log(1/(1-p))`` peaks near p = 0.63 and
    falls away at both ends. So the two have to be told apart rather than
    assumed apart.

    Fits one factor ``c`` to all genes at once -- the most gene-unspecific model
    there is -- and reports how much of the cross-gene spread it accounts for.
    A HIGH R² is the bad outcome: it means the ordering carries no information
    the abundances did not already carry.

    Returns ``(per-gene table, verdict dict)``.
    """
    needed = {"gene", "detect_normal", "delta_detect"}
    missing = needed - set(deltas.columns)
    if missing:
        raise ScaleError(f"the thinning null needs {sorted(missing)}")

    baselines, observed, genes = {}, {}, []
    for gene, block in deltas.groupby("gene", observed=True):
        rows = block[["detect_normal", "delta_detect"]].dropna()
        if rows.empty:
            continue
        genes.append(str(gene))
        baselines[str(gene)] = rows["detect_normal"].to_numpy(dtype=float)
        observed[str(gene)] = float(rows["delta_detect"].mean())

    if len(genes) < 3:
        return pd.DataFrame(), {
            "verdict": "UNDEFINED",
            "detail": (
                f"{len(genes)} gene(s) with a baseline rate. The null fits one "
                f"parameter across genes and needs at least 3 to say anything."
            ),
        }

    mean_baselines = np.array([baselines[g].mean() for g in genes])
    spread = float(mean_baselines.max() - mean_baselines.min())
    if spread < MIN_BASELINE_SPREAD:
        return pd.DataFrame(), {
            "verdict": "UNDEFINED",
            "baseline_spread": spread,
            "detail": (
                f"baseline detection spans only {spread:.3f}, below "
                f"{MIN_BASELINE_SPREAD}. With every gene at the same abundance "
                f"the null has nothing to discriminate on, and an R² from one "
                f"free parameter over {len(genes)} near-identical points would "
                f"be an artefact. ABSTAINING rather than reporting clean."
            ),
        }

    if grid is None:
        grid = np.arange(0.05, 3.0005, 0.001)
    target = np.array([observed[g] for g in genes])
    errors = np.array([
        float(np.sum((np.array([_predicted_delta(baselines[g], c) for g in genes])
                      - target) ** 2))
        for c in grid
    ])
    factor = float(grid[int(np.argmin(errors))])

    predicted = np.array([_predicted_delta(baselines[g], factor) for g in genes])
    residual = target - predicted
    ss_res = float(np.sum(residual ** 2))
    ss_tot = float(np.sum((target - target.mean()) ** 2))
    r_squared = float(1.0 - ss_res / ss_tot) if ss_tot > 0 else float("nan")

    table = pd.DataFrame({
        "gene": genes,
        "baseline_detection": mean_baselines,
        "observed_delta": target,
        "null_predicts": predicted,
        "residual": residual,
    }).sort_values("baseline_detection", ignore_index=True)

    order_rho = _rank_correlation(target, predicted)
    explains = np.isfinite(r_squared) and r_squared >= UNIFORM_THINNING_R2_CEILING
    verdict = {
        "verdict": "GRADIENT IS ABUNDANCE" if explains else "STRUCTURE SURVIVES",
        "common_factor": factor,
        "variance_explained": r_squared,
        "rank_correlation_with_null": order_rho,
        "baseline_spread": spread,
        "ceiling": UNIFORM_THINNING_R2_CEILING,
        "detail": (
            f"one common fold change of {factor:.3f}, with no gene-specificity "
            f"in it, explains {100 * r_squared:.1f}% of the cross-gene spread "
            f"in detection deltas (rank correlation {order_rho:+.2f}). "
            + (
                "At or above the ceiling: the ordering of the raw deltas is "
                "not readable as biology, and any tier structure must be shown "
                "on the log fold-change scale instead."
                if explains else
                "Below the ceiling: abundance alone does not reproduce the "
                "pattern, so structure beyond a uniform thinning survives. This "
                "licenses the ORDERING as more than abundance; it does not "
                "license reading the raw delta MAGNITUDES across genes."
            )
        ),
    }
    return table, verdict


def _rank_correlation(a: np.ndarray, b: np.ndarray) -> float:
    """Spearman, computed here so the diagnostic carries no new dependency."""
    if len(a) < 2:
        return float("nan")
    ra, rb = _ranks(a), _ranks(b)
    if ra.std() == 0 or rb.std() == 0:
        return float("nan")
    return float(np.corrcoef(ra, rb)[0, 1])


def _ranks(values: np.ndarray) -> np.ndarray:
    order = np.argsort(np.argsort(values))
    return order.astype(float)
