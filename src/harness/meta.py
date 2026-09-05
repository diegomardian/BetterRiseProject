"""Random-effects meta-analysis. The estimator invariant 4 has always required.

CLAUDE.md invariant 4, execution_plan.md §6.1 and README.md all say the same
thing -- "estimate per study, then meta-analyse, never pool" -- and until now
nothing in ``src/`` implemented the second half. The WMHS appendix says so
plainly: *"Our design named meta-analytic intervals as an output and the data
cannot supply them. A meta-analysis needs two quotable (study, rung) cells and
this project has one."* The ICBI atlas supplies fourteen, which is the first
time the rule has had anything to operate on.

DerSimonian-Laird, because it is the two-stage estimator the design document
names and because its assumptions are legible: each study estimates its own
quantity, those quantities are drawn from a distribution with variance tau^2,
and the pooled estimate weights by within-study variance plus tau^2.

WHY HETEROGENEITY GATES THE VERDICT RATHER THAN ANNOTATING IT.
The premise this is applied to -- "are the two arms of this study comparable?"
-- is a property OF A STUDY. Pool three studies whose controls sit at -0.7,
+0.1 and +0.6 and the random-effects estimate lands near zero with an interval
that may well clear the tolerance, reporting "the premise holds" when what is
actually true is that no two studies agree on what the premise is. The pooled
number is arithmetically correct and substantively meaningless.

So :func:`premise_verdict` refuses to read a pooled estimate when I^2 exceeds
``MAX_I_SQUARED``. That is not a caveat attached to an answer; it is a third
answer, and it has the same shape as every other three-state verdict in this
project -- satisfied, refused, and undecided.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np

log = logging.getLogger(__name__)


class MetaError(RuntimeError):
    """The inputs cannot support a meta-analysis."""


#: Fewer studies than this and tau^2 is not estimable in any useful sense --
#: DerSimonian-Laird's Q has k-1 degrees of freedom, so k=2 gives one.
MIN_STUDIES = 3

#: Above this, the pooled estimate is not read. Cochrane's conventional
#: "considerable heterogeneity" boundary, committed here rather than chosen
#: after seeing a value.
MAX_I_SQUARED = 0.75

#: Normal quantile for a 95% interval. The per-study inputs are bootstrap
#: intervals over patients, so this is a second-stage normal approximation on
#: already-resampled quantities.
Z_95 = 1.959963984540054


def se_from_interval(ci_low: float, ci_high: float, *, z: float = Z_95) -> float:
    """Standard error implied by a symmetric 95% interval.

    The per-study inputs are percentile bootstrap intervals, which are not
    exactly symmetric. Taking the half-width over ``z`` is the standard
    two-stage move and it is an approximation -- recorded here so it is visible
    rather than buried in a caller.
    """
    if not (np.isfinite(ci_low) and np.isfinite(ci_high)):
        return float("nan")
    if ci_high < ci_low:
        raise MetaError(f"interval [{ci_low}, {ci_high}] is inverted")
    return (ci_high - ci_low) / (2.0 * z)


@dataclass(frozen=True)
class MetaResult:
    """A pooled estimate and everything needed to distrust it."""

    k: int
    pooled: float
    se: float
    ci_low: float
    ci_high: float
    tau_squared: float
    i_squared: float
    q: float
    df: int
    #: Fixed-effect pooled estimate, for comparison. When tau^2 is 0 the two
    #: agree exactly, which is a cheap check that the arithmetic is right.
    pooled_fixed: float

    @property
    def homogeneous(self) -> bool:
        return bool(np.isfinite(self.i_squared) and self.i_squared <= MAX_I_SQUARED)

    def as_row(self) -> dict:
        return {
            "k_studies": self.k, "pooled": self.pooled, "se": self.se,
            "ci_low": self.ci_low, "ci_high": self.ci_high,
            "tau_squared": self.tau_squared, "i_squared": self.i_squared,
            "cochran_q": self.q, "df": self.df,
            "pooled_fixed_effect": self.pooled_fixed,
            "homogeneous": self.homogeneous,
        }


def meta_analyse(
    estimates: np.ndarray | list[float],
    standard_errors: np.ndarray | list[float],
) -> MetaResult:
    """DerSimonian-Laird random-effects pooling.

    ``estimates`` are per-study effects on a common scale; ``standard_errors``
    their within-study SEs. Studies carrying a non-finite estimate or SE are
    DROPPED with a warning rather than silently contributing -- a NaN in a
    weighted mean propagates to the whole pooled value, and this project has
    already published one number that turned out to be a coerced missing.
    """
    y = np.asarray(estimates, dtype=float)
    s = np.asarray(standard_errors, dtype=float)
    if y.shape != s.shape:
        raise MetaError(f"{y.size} estimates against {s.size} standard errors")

    usable = np.isfinite(y) & np.isfinite(s) & (s > 0)
    if usable.sum() < y.size:
        log.warning("meta: dropping %d of %d studies with a non-finite "
                    "estimate or a non-positive SE", int((~usable).sum()), y.size)
    y, s = y[usable], s[usable]
    k = y.size
    if k < MIN_STUDIES:
        raise MetaError(
            f"{k} usable study/studies, below the {MIN_STUDIES} minimum. "
            f"Cochran's Q has k-1 degrees of freedom, so tau^2 from fewer is "
            f"not an estimate of anything. Report the studies side by side."
        )

    v = s ** 2
    w = 1.0 / v
    sum_w = w.sum()
    pooled_fixed = float((w * y).sum() / sum_w)

    q = float((w * (y - pooled_fixed) ** 2).sum())
    df = k - 1
    c = float(sum_w - (w ** 2).sum() / sum_w)
    # tau^2 is a variance and cannot be negative; Q below its expectation means
    # less between-study variation than sampling alone predicts, which is
    # tau^2 = 0 rather than a negative number.
    tau_squared = max(0.0, (q - df) / c) if c > 0 else 0.0
    i_squared = max(0.0, (q - df) / q) if q > 0 else 0.0

    w_star = 1.0 / (v + tau_squared)
    sum_w_star = w_star.sum()
    pooled = float((w_star * y).sum() / sum_w_star)
    se = float(np.sqrt(1.0 / sum_w_star))

    return MetaResult(
        k=k, pooled=pooled, se=se,
        ci_low=pooled - Z_95 * se, ci_high=pooled + Z_95 * se,
        tau_squared=tau_squared, i_squared=i_squared, q=q, df=df,
        pooled_fixed=pooled_fixed,
    )


def premise_verdict(result: MetaResult, tolerance: float) -> tuple[str, str]:
    """The premise at the meta level: satisfied, refused, or undecided.

    Applied to a CONTROL statistic, where the question is whether a
    housekeeping gene moved between arms by more than ``tolerance``. The three
    states mirror the per-study check in
    ``src/reference/jobs/coexpression_silencing.py`` exactly, and gain a fourth
    route into "undecided" that only exists at this level: the studies may not
    be describing a common quantity at all.
    """
    if not result.homogeneous:
        return "UNRESOLVED", (
            f"I^2 = {result.i_squared:.1%} exceeds the pre-committed "
            f"{MAX_I_SQUARED:.0%} ceiling (tau^2 = {result.tau_squared:.4f}, "
            f"Q = {result.q:.2f} on {result.df} df). The studies are not "
            f"estimating a common quantity, so the pooled value "
            f"({result.pooled:+.3f}) is arithmetically correct and "
            f"substantively meaningless. Report the per-study estimates; do not "
            f"read this one."
        )

    reach = max(abs(result.ci_low), abs(result.ci_high))
    floor = 0.0 if result.ci_low <= 0.0 <= result.ci_high else min(
        abs(result.ci_low), abs(result.ci_high)
    )
    interval = f"{result.pooled:+.3f} [{result.ci_low:+.3f}, {result.ci_high:+.3f}]"

    if reach <= tolerance:
        return "HOLDS", (
            f"pooled control {interval} over k = {result.k}, entirely within "
            f"the {tolerance} tolerance (I^2 = {result.i_squared:.1%})."
        )
    if floor > tolerance:
        return "REFUSED", (
            f"pooled control {interval} over k = {result.k}, entirely beyond "
            f"the {tolerance} tolerance (I^2 = {result.i_squared:.1%}). The arms "
            f"are not comparable and the detection reading is not licensed."
        )
    return "UNRESOLVED", (
        f"pooled control {interval} over k = {result.k} straddles the "
        f"{tolerance} tolerance (I^2 = {result.i_squared:.1%}). Not refused and "
        f"not satisfied -- undecided, and reported as such."
    )
