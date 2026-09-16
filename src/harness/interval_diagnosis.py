"""Why the within-patient interval over-rejects. Causes tested by construction.

The WMHS submission reports null rejection of 12-30% at 50 mature cells and up
to 9.5% at 800, against a nominal 5%, and explains the five-cell case by
all-zero tumour draws. A reviewer's objection is that this explanation does not
reach 50 or 800 cells, and names three candidate causes: **zero inflation**,
**skew**, and the **fixed-fraction resampling** design.

This module tests all three by *substituting the thing* rather than by arguing
about it. Each cause is removed from the resampling population -- or from the
design -- and the null rejection is re-measured with everything else held
identical. What a cause is worth is then the difference between two measured
rates, in percentage points, and not an opinion.

THE POOL FAMILIES
-----------------
Under ``shift = 1.0`` the generator's ``_apply_shift`` returns early, so the
mature cells of both arms are i.i.d. draws with replacement from one finite pool
``F``: the held-out patients' mature-cell counts of the target gene. Every
family below replaces ``F`` and changes nothing else.

``empirical``
    ``F`` itself. The baseline.
``zeros_removed_mean_matched``
    ``F`` restricted to strictly positive values, then multiplied by a constant
    so its mean returns to ``mean(F)``. Removes the zero atom; keeps the mean
    and the shape of the positive part. The difference from ``empirical`` is
    what zero inflation is worth.
``skew_matched_no_zeros``
    A continuous, strictly positive distribution matched to ``mean(F)``,
    ``var(F)`` and ``skew(F)``: a three-parameter (shifted) gamma where the
    shift is non-negative, and a two-parameter gamma matched on mean and
    variance where it is not. Which one was used, and the realised moments, are
    recorded on every row -- a moment match that silently did not hold would
    make the skew attribution meaningless.
``gaussian_matched``
    ``Normal(mean(F), var(F))``. No zeros, no discreteness, no skew. This is the
    floor, and what survives here is the ``z``-vs-``t`` and plug-in-sd term whose
    closed form ``reference.interval_calibration.expected_false_positive_rate``
    gives. Comparing the two is falsifier F5 of the pre-registration: if the
    Gaussian floor does not match its own arithmetic, nothing else in this
    module can be trusted.

THE DESIGN FACTOR
-----------------
The committed sweep varies a mature *fraction* against a fixed 2,000 cells, so
the reference arm sits at 800 mature cells at every grid point while the tumour
arm is starved. ``balanced`` sets ``n_n = n_t`` instead. Note in advance that
the two designs coincide at 800 cells by construction -- the design can only
matter where the arms are unbalanced -- so a non-zero design effect at 800 would
be a bug, and is checked as one.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Final

import numpy as np
from scipy import stats

from src.harness.interval_repair import (
    ArmResamples,
    excludes_zero,
    mcse,
    percentile_interval,
)

#: Substituted for the empirical pool, one at a time.
POOL_FAMILIES: Final[tuple[str, ...]] = (
    "empirical",
    "zeros_removed_mean_matched",
    "skew_matched_no_zeros",
    "gaussian_matched",
)

#: ``fixed_fraction`` is what the committed sweep does.
DESIGNS: Final[tuple[str, ...]] = ("fixed_fraction", "balanced")

#: Reference-arm mature cells under ``fixed_fraction``: ``0.40 * 2000``.
N_REFERENCE_FIXED: Final[int] = 800

#: The scale factor on the intrinsic term. Constant across the whole sweep.
FRAC_MATURE_NORMAL: Final[float] = 0.40


class MomentMatch:
    """The realised moments of a substituted pool, carried with it.

    A family called ``skew_matched_no_zeros`` that did not in fact match the
    skew would make the attribution in the results table a fiction. So the
    realised moments travel with the sampler and land in the table, and
    ``skew_matched`` is a column a reader can check rather than a promise in a
    docstring.
    """

    __slots__ = ("family", "mean", "var", "skew", "form", "skew_matched", "min_support")

    def __init__(self, family, mean, var, skew, form, skew_matched, min_support):
        self.family = family
        self.mean = float(mean)
        self.var = float(var)
        self.skew = float(skew)
        self.form = form
        self.skew_matched = bool(skew_matched)
        self.min_support = float(min_support)

    def as_dict(self) -> dict:
        return {
            "family": self.family,
            "pool_mean": self.mean,
            "pool_var": self.var,
            "pool_skew": self.skew,
            "pool_form": self.form,
            "skew_matched": self.skew_matched,
            "pool_min_support": self.min_support,
        }


def _moments(values: np.ndarray) -> tuple[float, float, float]:
    """Population mean, variance and skewness of a finite pool."""
    values = np.asarray(values, dtype=float)
    mean = float(values.mean())
    var = float(values.var(ddof=0))
    if var <= 0:
        return mean, var, 0.0
    skew = float(np.mean((values - mean) ** 3) / var**1.5)
    return mean, var, skew


class PoolSampler:
    """Draws from one substituted population. ``None`` when it cannot exist."""

    __slots__ = ("_kind", "_values", "_args", "moments")

    def __init__(self, kind: str, moments: MomentMatch, values=None, args=None):
        self._kind = kind
        self._values = values
        self._args = args
        self.moments = moments

    def draw(self, n: int, rng: np.random.Generator) -> np.ndarray:
        if self._kind == "discrete":
            return rng.choice(self._values, size=n, replace=True)
        if self._kind == "gamma":
            shift, shape, scale = self._args
            return shift + rng.gamma(shape, scale, size=n)
        if self._kind == "normal":
            mean, sd = self._args
            return rng.normal(mean, sd, size=n)
        raise ValueError(f"unknown sampler kind {self._kind!r}")


def build_pool(values: np.ndarray, family: str) -> PoolSampler | None:
    """The substituted population for one family, or ``None`` if undefined.

    ``None`` rather than a fallback: a family that cannot be built for a given
    pool (no positive values, zero variance) has no rate, and inventing one by
    quietly substituting a different distribution is how a diagnosis becomes
    decorative. The caller counts these.
    """
    if family not in POOL_FAMILIES:
        raise ValueError(f"unknown family {family!r}; known: {list(POOL_FAMILIES)}")
    values = np.asarray(values, dtype=float)
    if values.size < 2:
        return None
    mean, var, skew = _moments(values)

    if family == "empirical":
        return PoolSampler(
            "discrete",
            MomentMatch(family, mean, var, skew, "empirical", True, float(values.min())),
            values=values,
        )

    if var <= 0:
        return None

    if family == "zeros_removed_mean_matched":
        positive = values[values > 0]
        if positive.size < 2 or positive.mean() <= 0:
            return None
        scaled = positive * (mean / positive.mean())
        m, v, s = _moments(scaled)
        return PoolSampler(
            "discrete",
            MomentMatch(family, m, v, s, "empirical_positive_rescaled", False,
                        float(scaled.min())),
            values=scaled,
        )

    if family == "skew_matched_no_zeros":
        if skew <= 0:
            return None
        # Three-parameter gamma: mean = shift + k*theta, var = k*theta^2,
        # skew = 2/sqrt(k). Exact on all three moments when the shift lands at
        # or above zero, which is what keeps the support strictly positive.
        shape = 4.0 / skew**2
        scale = float(np.sqrt(var / shape))
        shift = mean - shape * scale
        if shift >= 0.0:
            return PoolSampler(
                "gamma",
                MomentMatch(family, mean, var, skew, "shifted_gamma", True, shift),
                args=(shift, shape, scale),
            )
        # The shift would be negative, so matching all three moments would put
        # mass below zero and the family would no longer be what it is named.
        # Match mean and variance instead, and SAY SO on the row rather than
        # reporting a skew match that did not happen.
        shape = mean**2 / var
        scale = var / mean
        if not (shape > 0 and scale > 0):
            return None
        realised_skew = 2.0 / float(np.sqrt(shape))
        return PoolSampler(
            "gamma",
            MomentMatch(family, mean, var, realised_skew, "gamma_mean_var_only",
                        False, 0.0),
            args=(0.0, shape, scale),
        )

    # gaussian_matched
    return PoolSampler(
        "normal",
        MomentMatch(family, mean, var, 0.0, "normal", True, float("-inf")),
        args=(mean, float(np.sqrt(var))),
    )


def arm_sizes(count: int, design: str) -> tuple[int, int]:
    """``(n_normal, n_tumour)`` under one design."""
    if design == "fixed_fraction":
        return N_REFERENCE_FIXED, count
    if design == "balanced":
        return count, count
    raise ValueError(f"unknown design {design!r}; known: {list(DESIGNS)}")


def null_rejection(
    pool_values: Sequence[np.ndarray],
    *,
    count: int,
    family: str,
    design: str,
    n_boot: int,
    seeds: Sequence[int],
    frac_mature_normal: float = FRAC_MATURE_NORMAL,
) -> dict:
    """Measure null rejection for one (family, design, count).

    ``pool_values`` is one empirical pool per replicate -- the held-out
    patients' mature-cell counts for that replicate, so the pool varies from
    replicate to replicate exactly as it does in the real sweep. Substituting a
    single pooled population instead would remove the between-holdout variation
    that is part of what is being measured.

    The truth is exactly zero: under ``shift = 1.0`` both arms are i.i.d. draws
    from the same population, whatever that population is. So a rejection is a
    false positive by construction and no ground-truth bookkeeping is needed.
    """
    if len(pool_values) != len(seeds):
        raise ValueError(
            f"{len(pool_values)} pools for {len(seeds)} seeds; one pool per "
            f"replicate is what makes the holdout variation real"
        )
    n_n, n_t = arm_sizes(count, design)
    rejects = 0
    attempted = 0
    undefined_pool = 0
    widths: list[float] = []
    moment_rows: list[dict] = []

    for values, seed in zip(pool_values, seeds, strict=True):
        sampler = build_pool(values, family)
        if sampler is None:
            undefined_pool += 1
            continue
        moment_rows.append(sampler.moments.as_dict())
        rng = np.random.default_rng(seed)
        normal = sampler.draw(n_n, rng)
        tumour = sampler.draw(n_t, rng)
        boot_n = ArmResamples(normal, n_boot=n_boot, rng=rng)
        boot_t = ArmResamples(tumour, n_boot=n_boot, rng=rng)
        interval = percentile_interval(
            boot_n, boot_t, frac_mature_normal=frac_mature_normal
        )
        verdict = excludes_zero(interval)
        if verdict is None:
            continue
        attempted += 1
        rejects += int(verdict)
        widths.append(interval[1] - interval[0])

    rate = rejects / attempted if attempted else float("nan")
    row = {
        "family": family,
        "design": design,
        "n_cells_mature": count,
        "n_normal": n_n,
        "n_tumour": n_t,
        "n_replicates": len(seeds),
        "n_scored": attempted,
        "n_pool_undefined": undefined_pool,
        "null_rejection": rate,
        "null_rejection_mcse": mcse(rate, attempted) if attempted else float("nan"),
        "median_ci_width": float(np.median(widths)) if widths else float("nan"),
    }
    if moment_rows:
        row |= {
            "pool_mean": float(np.mean([m["pool_mean"] for m in moment_rows])),
            "pool_var": float(np.mean([m["pool_var"] for m in moment_rows])),
            "pool_skew": float(np.mean([m["pool_skew"] for m in moment_rows])),
            "pool_form": moment_rows[0]["pool_form"],
            "skew_matched": all(m["skew_matched"] for m in moment_rows),
        }
    return row


# --------------------------------------------------------------------------
# guards. Each has a committed failing input in tests/test_checks_can_fail.py
# --------------------------------------------------------------------------


class DiagnosisError(ValueError):
    """A diagnosis that does not describe the thing it claims to diagnose."""


def gaussian_floor(n_normal: int, n_tumour: int, alpha: float = 0.05) -> float:
    """Closed-form false-positive rate of the TWO-SAMPLE percentile bootstrap.

    ``interval_calibration.expected_false_positive_rate`` is the one-sample
    statement: ``P(|T_{n-1}| > z*sqrt((n-1)/n))``. The interval under test is a
    two-sample one with unequal arms, so the one-sample form has no single ``n``
    to be evaluated at, and using it anyway -- at ``n_t``, or at some effective
    n -- would be comparing a measurement to an arithmetic statement about a
    different quantity.

    The two-sample form, for arms of equal variance (which is exactly true in
    the ``gaussian_matched`` family, since both arms are drawn from one normal
    distribution):

    * the bootstrap resamples each arm at its own n with the plug-in,
      divide-by-n variance, so its SD is
      ``sigma * sqrt((n_n-1)/n_n^2 + (n_t-1)/n_t^2)``,
    * the interval with 95% coverage is ``t(v) * sigma * sqrt(1/n_n + 1/n_t)``
      at the Welch degrees of freedom ``v``,

    so the percentile interval is ``shrink * z / t(v)`` times the width it
    claims, and its false-positive rate is ``P(|T_v| > z * shrink)``.

    Reduces to the one-sample form's ``sqrt((n-1)/n)`` shrink at ``n_n = n_t``,
    where ``v`` becomes ``2(n-1)``; that identity is asserted in the tests.
    """
    n_n, n_t = int(n_normal), int(n_tumour)
    if n_n < 2 or n_t < 2:
        return float("nan")
    boot_var = (n_n - 1) / n_n**2 + (n_t - 1) / n_t**2
    true_var = 1.0 / n_n + 1.0 / n_t
    shrink = float(np.sqrt(boot_var / true_var))
    df = true_var**2 / (
        1.0 / (n_n**2 * (n_n - 1)) + 1.0 / (n_t**2 * (n_t - 1))
    )
    z = float(stats.norm.ppf(1 - alpha / 2))
    return float(2 * stats.t.sf(z * shrink, df))


def check_gaussian_floor_matches_its_arithmetic(
    rows, *, tolerance_mcse: float = 3.0
) -> None:
    """Falsifier F5. The Gaussian floor must obey the closed form.

    ``gaussian_matched`` removes zeros, discreteness and skew, so what is left
    is the plug-in-variance and ``z``-vs-``t`` term, and :func:`gaussian_floor`
    says exactly what that costs at a given pair of arm sizes. If the measured
    floor and the arithmetic disagree by more than three Monte-Carlo standard
    errors, the apparatus is not measuring what this module says it measures and
    every attribution downstream is void.

    Applies to BOTH designs, because the two-sample closed form does not need
    the arms to be equal.
    """
    bad = []
    for row in rows:
        if row["family"] != "gaussian_matched":
            continue
        if not np.isfinite(row["null_rejection"]):
            continue
        expected = gaussian_floor(row["n_normal"], row["n_tumour"])
        se = row["null_rejection_mcse"]
        if not np.isfinite(se) or se <= 0 or not np.isfinite(expected):
            continue
        if abs(row["null_rejection"] - expected) > tolerance_mcse * se:
            bad.append(
                f"{row['design']} n=({row['n_normal']},{row['n_tumour']}): "
                f"measured {row['null_rejection']:.4f} vs closed form "
                f"{expected:.4f} "
                f"({abs(row['null_rejection'] - expected) / se:.1f} MCSE)"
            )
    if bad:
        raise DiagnosisError(
            "the Gaussian floor does not match the closed form it must obey "
            "(pre-registration falsifier F5), so the apparatus is not measuring "
            "what this module claims and no attribution below it is meaningful: "
            + "; ".join(bad)
        )


def check_design_effect_vanishes_where_the_arms_coincide(rows) -> None:
    """``fixed_fraction`` and ``balanced`` ARE the same design at 800 cells.

    ``arm_sizes(800, "fixed_fraction") == arm_sizes(800, "balanced") == (800, 800)``,
    so any difference between them at that count is Monte-Carlo noise at most.
    A systematic one would mean the two designs are not being run on the same
    thing, which would silently inflate the design's share of the attribution --
    the one number the reviewer specifically asks about.
    """
    by_key: dict[tuple, dict[str, dict]] = {}
    for row in rows:
        if row["n_normal"] != row["n_tumour"]:
            continue
        by_key.setdefault((row["family"], row["n_cells_mature"]), {})[row["design"]] = row

    bad = []
    for (family, count), designs in by_key.items():
        if set(designs) != set(DESIGNS):
            continue
        a, b = designs["fixed_fraction"], designs["balanced"]
        if not (np.isfinite(a["null_rejection"]) and np.isfinite(b["null_rejection"])):
            continue
        se = float(np.hypot(a["null_rejection_mcse"], b["null_rejection_mcse"]))
        if se <= 0:
            continue
        if abs(a["null_rejection"] - b["null_rejection"]) > 4.0 * se:
            bad.append(
                f"{family} at n={count}: {a['null_rejection']:.4f} vs "
                f"{b['null_rejection']:.4f} ({abs(a['null_rejection'] - b['null_rejection']) / se:.1f} MCSE)"
            )
    if bad:
        raise DiagnosisError(
            "the two designs give different null rejection where they are the "
            "same design (both arms at the same n), so they are not being run "
            "on the same thing and the design's share of the attribution is not "
            "trustworthy: " + "; ".join(bad)
        )


def attribute(rows, *, count: int) -> dict:
    """Percentage points attributable to each cause, at one count.

    Reports BOTH a sequential attribution and the single-factor deltas, plus the
    residual left over. The two disagree exactly when the causes are not
    additive, which is a finding about the mechanism and is reported as one
    rather than resolved by choosing an order.

    WHAT IS DELIBERATELY NOT CREDITED. The reviewer named three causes, and only
    those three are credited. The step from ``zeros_removed_mean_matched`` to
    ``skew_matched_no_zeros`` replaces a discrete empirical pool with a
    continuous one matched on its first three moments, so it removes
    DISCRETENESS and every moment above the third. That is a fourth cause, it is
    nobody's hypothesis, and crediting it to "skew" would make the skew number
    an overstatement. It is reported as ``pp_discreteness_and_higher_moments``
    and it falls into the residual, which is where an unnamed cause belongs.
    """
    lookup = {
        (r["family"], r["design"]): r["null_rejection"]
        for r in rows
        if r["n_cells_mature"] == count
    }
    missing = [
        (f, "fixed_fraction") for f in POOL_FAMILIES
        if (f, "fixed_fraction") not in lookup
    ]
    if missing:
        raise DiagnosisError(f"no rate for {missing} at n={count}")

    p0 = lookup[("empirical", "fixed_fraction")]
    p1 = lookup[("zeros_removed_mean_matched", "fixed_fraction")]
    p2 = lookup[("skew_matched_no_zeros", "fixed_fraction")]
    p3 = lookup[("gaussian_matched", "fixed_fraction")]
    balanced = lookup.get(("empirical", "balanced"), float("nan"))

    zero_inflation = p0 - p1
    skew_given_no_zeros = p2 - p3
    # Not credited to any named cause; see the docstring. Lands in the residual.
    discreteness = p1 - p2
    shape_total = p1 - p3
    design = p0 - balanced
    floor = p3 - 0.05
    explained = zero_inflation + skew_given_no_zeros + design + floor
    excess = p0 - 0.05

    return {
        "n_cells_mature": count,
        "observed_null_rejection": p0,
        "excess_over_nominal": excess,
        "pp_zero_inflation": zero_inflation,
        "pp_skew_given_no_zeros": skew_given_no_zeros,
        "pp_discreteness_and_higher_moments": discreteness,
        "pp_non_gaussian_shape_total": shape_total,
        "pp_fixed_fraction_design": design,
        "pp_floor_z_vs_t": floor,
        "pp_explained_sequential": explained,
        "pp_residual_unexplained": excess - explained,
        "share_unexplained": (excess - explained) / excess if excess else float("nan"),
        "rate_empirical": p0,
        "rate_zeros_removed": p1,
        "rate_skew_matched": p2,
        "rate_gaussian": p3,
        "rate_empirical_balanced": balanced,
    }


def closure_verdict(
    attribution: dict,
    *,
    max_unexplained_share: float = 0.40,
    negligible_excess: float = 0.005,
) -> str:
    """Falsifier F3, as a named verdict rather than a bare boolean.

    ``over_explained`` is not a pedantic extra case. A sequential attribution
    whose parts sum to MORE than the excess describes the mechanism no better
    than one that leaves most of it open -- the causes are interacting, or one
    of the substituted pools is not the thing it is named -- and a bare
    "closes / does not close" would have reported the first as success.

    ``no_excess_to_explain`` is the case where the observed rate is already at
    nominal. There is then nothing to attribute, and calling that a closed
    diagnosis would be claiming a result from an empty question.
    """
    excess = attribution["excess_over_nominal"]
    if not np.isfinite(excess) or abs(excess) <= negligible_excess:
        return "no_excess_to_explain"
    share = attribution["share_unexplained"]
    if not np.isfinite(share):
        return "leaves_most_unexplained"
    if share > max_unexplained_share:
        return "leaves_most_unexplained"
    if share < -max_unexplained_share:
        return "over_explained"
    return "closes"


def diagnosis_closes(attribution: dict, *, max_unexplained_share: float = 0.40) -> bool:
    """``True`` only when the named causes account for the excess.

    Not an exception. An unexplained residual is a legitimate and publishable
    result -- the pre-registration says so in those words -- so this returns a
    verdict to be written into the table, not an error that stops the run.
    """
    return closure_verdict(
        attribution, max_unexplained_share=max_unexplained_share
    ) == "closes"


def normal_quantile_shortfall(n: int, alpha: float = 0.05) -> float:
    """``z*sqrt((n-1)/n)/t(n-1)`` -- the width the percentile interval is missing.

    Re-exported here so a reader of the diagnosis does not have to go to
    ``reference.interval_calibration`` to find out how small this term is. It is
    0.966 at n=50 and 0.9993 at n=800: a floor, not an explanation.
    """
    z = float(stats.norm.ppf(1 - alpha / 2))
    t = float(stats.t.ppf(1 - alpha / 2, max(n - 1, 1)))
    return z * float(np.sqrt((n - 1) / n)) / t
