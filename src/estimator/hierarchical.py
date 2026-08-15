"""Hierarchical (mixed-effects) uncertainty for the Kitagawa building blocks. W4.

Complementary to ``kitagawa.bootstrap_over_patients``, not a replacement — the
week 4-5 deliverable in src/estimator/README.md lists both together:
"Patient-level bootstrap; hierarchical model with patient as grouping factor —
CIs that reflect the real unit of inference." The bootstrap resamples patients
nonparametrically; this module fits two mixed-effects models with a random
intercept per patient (statsmodels MixedLM) to get patient-clustered standard
errors on the two Kitagawa inputs — Δ(mature fraction) and Δ(per-cell mean) —
and propagates them into compositional/intrinsic through the same arithmetic
as ``kitagawa.decompose()``.

This is NOT regression Oaxaca-Blinder (README, CLAUDE.md invariant on this):
the mixed models here estimate the uncertainty of the two summary statistics
Kitagawa needs. They do not define the split — that stays exactly
``kitagawa.decompose()``'s arithmetic, unchanged.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from src.estimator.kitagawa import Weighting, decompose

_REQUIRED_CELL_COLUMNS = {"patient_id", "tissue", "mature", "expression"}


def _fit_tissue_effect(cells: pd.DataFrame, outcome: str) -> tuple[float, float]:
    """Random-intercept-per-patient model of ``outcome`` ~ tissue.

    Returns ``(coefficient, standard_error)`` for the tumour-vs-normal
    contrast. Cells within a patient are correlated; a plain per-cell average
    or OLS treats each cell as an independent observation and understates the
    uncertainty by roughly the number of cells per patient — the same failure
    mode CLAUDE.md invariant 5 calls out for a naive cell-level bootstrap.

    For ``outcome="mature"`` (a 0/1 indicator) this is a linear probability
    model with a random intercept, not a mixed logistic regression —
    statsmodels has no mature binomial-mixed-model fitter outside
    ``BinomialBayesMixedGLM``'s variational fit, which is a different fitting
    paradigm entirely and a bigger change to make without real data to
    validate it against. Expect convergence warnings from statsmodels here,
    especially near extreme fractions; they are visible, not silent, and are
    the first thing to revisit once real Lee data replaces synthetic
    fixtures. Imports statsmodels lazily so ``kitagawa.py``'s core arithmetic
    never needs it — only this module does.
    """
    import statsmodels.formula.api as smf

    data = cells.copy()
    data["tissue_tumour"] = (data["tissue"] == "tumour").astype(float)
    fit = smf.mixedlm(f"{outcome} ~ tissue_tumour", data, groups=data["patient_id"]).fit(reml=True)
    return float(fit.params["tissue_tumour"]), float(fit.bse["tissue_tumour"])


def hierarchical_intrinsic_ci(
    cells: pd.DataFrame,
    *,
    weighting: Weighting = "normal",
    seed: int,
    n_draws: int = 2000,
    alpha: float = 0.05,
) -> dict:
    """Model-based, patient-clustered CI for one (gene, rung, axis, weighting).

    ``cells`` is long-format, one row per cell, for a single gene/rung/axis:
    columns ``patient_id``, ``tissue`` ('normal'/'tumour'), ``mature`` (0/1),
    ``expression`` (float). Two random-intercept models are fit:

      1. ``mature ~ tissue``                          -> Δ(mature fraction)
      2. ``expression ~ tissue``, mature cells only    -> Δ(per-cell mean)

    Reference points (the normal- or tumour-group mature fraction and mean,
    depending on ``weighting``) are taken as the observed group means — only
    the *deltas* carry model-based uncertainty here. That is a deliberate
    simplification, and it is conservative in one specific sense: it will
    generally understate the true interval relative to also propagating
    uncertainty in the reference point, and that direction of error was
    chosen over the alternative. If it needs tightening once real Lee data is
    in hand, tighten it there, not by widening this function silently.

    Draws ``n_draws`` samples from Normal(coefficient, standard_error) for
    each delta independently, recombines each draw through
    ``kitagawa.decompose()``, and reports the ``alpha``-level percentile
    interval on ``intrinsic`` — the same shape of output as
    ``bootstrap_over_patients``, so the two are directly comparable as a
    cross-check on each other, per the week 4-5 "done when."
    """
    missing = _REQUIRED_CELL_COLUMNS - set(cells.columns)
    if missing:
        raise ValueError(f"cells is missing column(s): {sorted(missing)}")
    bad_tissue = set(cells["tissue"].unique()) - {"normal", "tumour"}
    if bad_tissue:
        raise ValueError(f"tissue must be 'normal' or 'tumour', found {sorted(bad_tissue)}")
    if n_draws < 1:
        raise ValueError(f"n_draws={n_draws} must be positive")

    frac_normal = float(cells.loc[cells["tissue"] == "normal", "mature"].mean())
    frac_tumour = float(cells.loc[cells["tissue"] == "tumour", "mature"].mean())
    mature_cells = cells[cells["mature"] == 1]
    mean_normal = float(mature_cells.loc[mature_cells["tissue"] == "normal", "expression"].mean())
    mean_tumour = float(mature_cells.loc[mature_cells["tissue"] == "tumour", "expression"].mean())
    n_cells_mature = int((mature_cells["tissue"] == "tumour").sum())

    _frac_coef, frac_se = _fit_tissue_effect(cells, "mature")
    _mean_coef, mean_se = _fit_tissue_effect(mature_cells, "expression")

    rng = np.random.default_rng(seed)
    d_frac_draws = rng.normal(frac_tumour - frac_normal, frac_se, size=n_draws)
    d_mean_draws = rng.normal(mean_tumour - mean_normal, mean_se, size=n_draws)

    intrinsic_draws = np.empty(n_draws)
    for i in range(n_draws):
        d = decompose(
            frac_normal,
            frac_normal + d_frac_draws[i],
            mean_normal,
            mean_normal + d_mean_draws[i],
            n_cells_mature=n_cells_mature,
            weighting=weighting,
        )
        intrinsic_draws[i] = d.intrinsic

    point = decompose(
        frac_normal,
        frac_tumour,
        mean_normal,
        mean_tumour,
        n_cells_mature=n_cells_mature,
        weighting=weighting,
    )
    lo_q, hi_q = 100 * (alpha / 2), 100 * (1 - alpha / 2)
    return {
        "compositional": point.compositional,
        "intrinsic": point.intrinsic,
        "interaction": point.interaction,
        "n_cells_mature": n_cells_mature,
        "ci_low": float(np.percentile(intrinsic_draws, lo_q)),
        "ci_high": float(np.percentile(intrinsic_draws, hi_q)),
        "delta_frac_se": frac_se,
        "delta_mean_se": mean_se,
    }
