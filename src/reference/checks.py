"""G1: is the residual signal ambient RNA rather than biology? W1, week 4.

Ambient counts are enriched for whatever is abundant. So **if a gene's
behaviour tracks its abundance, that behaviour is a property of the soup and
not of the tumour.** G1 is the check that catches it.

**The thresholds in this module were committed before it was written.**
Decision #17, ``docs/open_decisions.md``. Nothing here was tuned against a G1
result, because none existed when the numbers were chosen — and after this runs
once, that ordering can never be recovered. Do not edit the constants to make a
run pass; the whole value of a pre-registered threshold is that it is allowed to
fail.

**Two statistics, both reported.** ``execution_plan.md`` §4 names G1 as
*post-correction retention vs total abundance*. Decision #17 as first written
substituted *apparent loss vs abundance* — a different measurement — without
saying so. The correction is on the record and the resolution is to run both:

- :func:`retention_correlation` — the plan's, and the **named gate criterion**.
  Asks whether *the correction* is abundance-driven.
- :func:`loss_correlation` — asks whether *the project's signal* is
  abundance-driven, which is what G1's stated consequence is about.

Where the two disagree, the disagreement is the finding, not a tie to break.
:func:`compare_statistics` reports it.
"""

from __future__ import annotations

from typing import Any, Final

import numpy as np
import pandas as pd

from src.common.panel import tier_of

__all__ = [
    "G1Error",
    "G1_TIERS",
    "MIN_GENES_PER_TIER",
    "MIN_TIER_SEPARATION",
    "TIER_D_MAX_RHO",
    "compare_statistics",
    "g1_verdict",
    "loss_correlation",
    "retention_correlation",
    "tier_correlations",
]


class G1Error(Exception):
    """Raised when a G1 input cannot be read as asked."""


#: The tiers G1 compares. A is the compositional target set, B the intrinsic
#: set, D the negative control — genes chosen to carry no differentiation story.
#: Tier C is deliberately absent: it is not part of the falsification design.
G1_TIERS: Final[tuple[str, ...]] = ("A", "B", "D")

#: **Pre-committed, decision #17.** Above this, tier D's abundance relationship
#: has no biological reading left — tier D genes have no differentiation story,
#: so a strong correlation there is the soup, measured.
#:
#: 0.5 is the same rank-correlation line #16 already uses for method agreement,
#: so the project carries one meaning of "these two things track each other"
#: rather than a different one per test.
TIER_D_MAX_RHO: Final[float] = 0.5

#: **Pre-committed, decision #17.** If A, B and D all fall within this of one
#: another, the panel is measuring abundance and the tier structure — the whole
#: falsification design — carries no information.
#:
#: 0.2 is the smallest separation that survives n≈8 genes per tier. Below it,
#: tier differences are not distinguishable from noise at this panel size, and
#: pretending otherwise would manufacture a pass.
MIN_TIER_SEPARATION: Final[float] = 0.2

#: Below this many genes a within-tier Spearman is not worth reading, and the
#: tier is reported ``not_estimable`` rather than given a number.
#:
#: **Not a pass.** CLAUDE.md invariant 1 in its general form: an unestimable
#: term is `None`, never a convenient default. A tier that cannot be measured
#: must not clear a gate by silence.
MIN_GENES_PER_TIER: Final[int] = 4


def _spearman(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """Spearman rho and its p-value.

    **Spearman, not Pearson**, and for a stated reason: abundance spans orders
    of magnitude, so a Pearson coefficient here is decided by a handful of very
    high genes. The same reasoning as the retention comparison in #16.
    """
    from scipy.stats import spearmanr

    result = spearmanr(x, y)
    return float(result.statistic), float(result.pvalue)


def tier_correlations(
    frame: pd.DataFrame,
    *,
    value_column: str,
    abundance_column: str = "abundance",
    gene_column: str = "gene",
) -> pd.DataFrame:
    """Spearman of abundance against `value_column`, **within each tier**.

    One row per tier in :data:`G1_TIERS`, always all three, whether or not each
    was estimable — a tier missing from the output would read as "fine" to
    anything downstream that only looks at the rows present.

    Genes outside the panel are dropped, and genes whose tier is not one of
    A/B/D are dropped with them. Tier assignment comes from the frozen panel
    (``config/panel.yaml``), never from the caller, so a tier cannot be
    relabelled to move a result.
    """
    for column in (gene_column, abundance_column, value_column):
        if column not in frame.columns:
            raise G1Error(
                f"input is missing column {column!r}; has {sorted(frame.columns)}"
            )
    if frame.empty:
        raise G1Error("input frame is empty")

    work = frame[[gene_column, abundance_column, value_column]].copy()
    work["tier"] = work[gene_column].map(tier_of)
    work = work[work["tier"].isin(G1_TIERS)]
    work = work.replace([np.inf, -np.inf], np.nan).dropna(
        subset=[abundance_column, value_column]
    )

    rows: list[dict[str, Any]] = []
    for tier in G1_TIERS:
        genes = work[work["tier"] == tier]
        n = int(len(genes))
        if n < MIN_GENES_PER_TIER:
            rows.append({
                "tier": tier,
                "n_genes": n,
                "rho": None,
                "p_value": None,
                "estimability": "not_estimable",
                "reason": (
                    f"{n} usable gene(s), below the {MIN_GENES_PER_TIER} a "
                    f"within-tier Spearman needs"
                ),
            })
            continue
        # Constant input makes rho undefined rather than zero. scipy returns
        # nan and warns; say so explicitly instead of emitting nan downstream.
        x = genes[abundance_column].to_numpy(dtype=float)
        y = genes[value_column].to_numpy(dtype=float)
        if np.ptp(x) == 0 or np.ptp(y) == 0:
            rows.append({
                "tier": tier, "n_genes": n, "rho": None, "p_value": None,
                "estimability": "not_estimable",
                "reason": "abundance or value is constant within the tier",
            })
            continue
        rho, p = _spearman(x, y)
        rows.append({
            "tier": tier, "n_genes": n, "rho": rho, "p_value": p,
            "estimability": "estimated", "reason": "",
        })

    out = pd.DataFrame(rows)
    out["statistic"] = value_column
    return out


def loss_correlation(
    frame: pd.DataFrame, *, loss_column: str = "loss", **kwargs: Any
) -> pd.DataFrame:
    """Abundance against **apparent loss** (Δ per-cell mean, tumour − normal).

    Decision #17's substituted statistic. Asks whether *the project's signal*
    is abundance-driven. Reported, but **not** the named gate criterion — see
    the module docstring.
    """
    return tier_correlations(frame, value_column=loss_column, **kwargs)


def retention_correlation(
    frame: pd.DataFrame, *, retention_column: str = "retention", **kwargs: Any
) -> pd.DataFrame:
    """Abundance against **post-correction retention** (`after / before`).

    ``execution_plan.md`` §4's statistic, and **the named gate criterion**.
    Asks whether *the correction* is abundance-driven.

    Expect this one to look bad on this cohort, near-tautologically: soup is
    enriched for abundant genes, so any method that removes soup strips abundant
    genes hardest. That is precisely why it is reported next to
    :func:`loss_correlation` rather than alone.
    """
    return tier_correlations(frame, value_column=retention_column, **kwargs)


def g1_verdict(correlations: pd.DataFrame) -> dict[str, Any]:
    """Apply decision #17's pre-committed thresholds. **Do not retune.**

    Returns ``verdict`` in ``{"pass", "fail", "indeterminate",
    "not_estimable"}`` with the reasons that produced it.

    ``fail`` if either pre-committed condition holds:

    1. ``|rho| > TIER_D_MAX_RHO`` in tier D.
    2. All three tier correlations fall within :data:`MIN_TIER_SEPARATION` of
       one another.

    ``pass`` if tier D is flat **and** A and B each separate from D by more than
    :data:`MIN_TIER_SEPARATION`.

    **``indeterminate`` is a real outcome here, not a defect in this code.**
    The rule as written in #17 leaves a gap: tier D flat, the overall range
    wider than 0.2, but only *one* of A and B separated from D satisfies
    neither the fail conditions nor the pass condition. Resolving that gap by
    picking a side would be choosing a threshold after seeing the data, which is
    the one thing #17 exists to prevent — so it is surfaced and left to the
    team. #17 already says the rule should be ratified before this runs.
    """
    required = {"tier", "rho", "estimability"}
    missing = required - set(correlations.columns)
    if missing:
        raise G1Error(f"correlations is missing column(s): {sorted(missing)}")

    by_tier = correlations.set_index("tier")
    unestimable = [
        t for t in G1_TIERS
        if t not in by_tier.index
        or by_tier.loc[t, "estimability"] != "estimated"
    ]
    if unestimable:
        return {
            "verdict": "not_estimable",
            "reasons": [
                f"tier(s) {', '.join(unestimable)} could not be estimated; the "
                f"gate is not evaluable rather than passed"
            ],
            "rho": {t: None for t in G1_TIERS},
        }

    rho = {t: float(by_tier.loc[t, "rho"]) for t in G1_TIERS}
    reasons: list[str] = []

    d_flat = abs(rho["D"]) <= TIER_D_MAX_RHO
    if not d_flat:
        reasons.append(
            f"tier D |rho| = {abs(rho['D']):.3f} > {TIER_D_MAX_RHO}: tier D "
            f"genes carry no differentiation story, so an abundance "
            f"relationship there has no biological reading left"
        )

    spread = max(rho.values()) - min(rho.values())
    tiers_alike = spread < MIN_TIER_SEPARATION
    if tiers_alike:
        reasons.append(
            f"all three tiers fall within {spread:.3f} < "
            f"{MIN_TIER_SEPARATION}: the panel is measuring abundance and the "
            f"tier structure carries no information"
        )

    if reasons:
        return {"verdict": "fail", "reasons": reasons, "rho": rho}

    separated = {
        t: abs(rho[t] - rho["D"]) > MIN_TIER_SEPARATION for t in ("A", "B")
    }
    if d_flat and all(separated.values()):
        return {
            "verdict": "pass",
            "reasons": [
                f"tier D flat (|rho| = {abs(rho['D']):.3f}) and tiers A and B "
                f"both separate from it by more than {MIN_TIER_SEPARATION}"
            ],
            "rho": rho,
        }

    unseparated = [t for t, ok in separated.items() if not ok]
    return {
        "verdict": "indeterminate",
        "reasons": [
            f"tier D is flat and the tiers are not all alike, but tier(s) "
            f"{', '.join(unseparated)} sit within {MIN_TIER_SEPARATION} of D. "
            f"#17's rule does not cover this case; it needs the team, not a "
            f"choice made here after seeing the number"
        ],
        "rho": rho,
    }


def compare_statistics(
    named: pd.DataFrame, secondary: pd.DataFrame
) -> dict[str, Any]:
    """Run both verdicts and report whether they agree.

    `named` is the gate criterion (retention vs abundance); `secondary` is the
    loss version. **Disagreement is the finding**, not something to resolve by
    preferring one — the gate is decided by `named`, and the disagreement is
    reported beside it.
    """
    named_verdict = g1_verdict(named)
    secondary_verdict = g1_verdict(secondary)
    agree = named_verdict["verdict"] == secondary_verdict["verdict"]
    return {
        "gate_verdict": named_verdict["verdict"],
        "named": named_verdict,
        "secondary": secondary_verdict,
        "agree": agree,
        "note": (
            "Both statistics reach the same verdict."
            if agree else
            "The two statistics disagree. The gate is decided by the "
            "pre-registered one (retention vs abundance); the disagreement is "
            "reported, not resolved — see decision #17."
        ),
    }
