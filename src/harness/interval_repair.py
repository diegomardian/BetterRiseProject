"""Candidate repairs for the within-patient intrinsic interval. W2.

WHAT THIS IS FOR
----------------
``src/harness/interval.py`` builds a **percentile bootstrap over cells** for one
patient's intrinsic term, and the WMHS submission measures its null rejection at
37-67% (5 cells), 12-30% (50 cells) and up to 9.5% (800 cells) against a nominal
5%. A reviewer's strongest objection to that paper is that it diagnoses the
defect and never tries to repair it. This module is the repair attempt: five
candidate intervals plus the status-quo baseline, every one of them applied to
the *same* draws so that a difference between them is a difference between
methods and nothing else.

The candidate list, the pass criterion and the falsifiers are pre-registered in
``docs/prereg_repaired_interval.md`` and were committed before any result table.
Nothing may be added here without carrying ``post_hoc`` into the results.

WHICH BOOTSTRAP THIS IS
-----------------------
The one over **cells**, within a patient. ``src/reference/interval_calibration``
analyses the one over **patients**, within a cohort, and proves a closed form:
the percentile interval is ``z*sqrt((n-1)/n)/t(n-1)`` times the width it claims.
That arithmetic is about a mean of ``n`` exchangeable values and so applies here
too, at the cell count — but it is worth 0.8 percentage points at n=50 and 0.07
at n=800. It is a floor, not the explanation, and this module does not offer it
as one. See the pre-registration §0.

THE ESTIMAND, WHICH IS AFFINE AND THAT MATTERS
----------------------------------------------
Under ``weighting="normal"`` -- what the sweep and the committed diagnostic run
-- ``kitagawa.decompose`` gives ``intrinsic = f_n * (mean_t - mean_n)``, and
``within_patient_intrinsic_ci`` holds ``f_n`` fixed across resamples. ``f_t``
does not enter this term at all. So every interval below is an affine image of a
two-sample interval on a difference of means, with the constant ``f_n``.

This module therefore computes ``f_n * (mean_t - mean_n)`` directly rather than
calling ``decompose`` once per bootstrap draw, because at B=2000 resamples times
2000 replicates times 32 settings the Python-level call is the entire runtime.
That is a performance shortcut around the real estimator, which is exactly the
kind of shortcut that silently goes stale, so it is pinned:
``tests/test_interval_repair.py::test_affine_shortcut_equals_the_real_estimator``
asserts elementwise identity against ``decompose()``. If W4 changes the
estimator, that test fails rather than these results quietly ceasing to describe
it.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Final

import numpy as np
from scipy import stats

#: Nominal level every interval here claims.
NOMINAL_ALPHA: Final[float] = 0.05

#: Resamples per interval. Held identical across candidates: B is not a free
#: parameter to be tuned per method. The sweep's own 200 is too few for a BCa
#: tail -- its 2.5% adjusted quantile would rest on about five draws.
DEFAULT_N_BOOT: Final[int] = 2000

#: The frozen candidate list. ``percentile`` is the status quo and a baseline,
#: not a candidate. Order is the order they are reported in.
CANDIDATES: Final[tuple[str, ...]] = (
    "percentile",
    "studentised",
    "bca",
    "welch_t",
    "patient_cluster",
    "pseudobulk_patient_t",
)

#: Candidates that need cluster labels (patient of origin per drawn cell).
NEEDS_PATIENTS: Final[frozenset[str]] = frozenset(
    {"patient_cluster", "pseudobulk_patient_t"}
)

#: Returned in place of an interval. ``None``, never ``0.0`` -- CLAUDE.md
#: invariant 1. An abstention is counted in its own column upstream and is
#: never folded into a rate.
ABSTAIN: Final[tuple[None, None]] = (None, None)

Interval = tuple[float | None, float | None]


# --------------------------------------------------------------------------
# the point estimate and its standard error
# --------------------------------------------------------------------------


def intrinsic_point(
    mean_normal: float, mean_tumour: float, *, frac_mature_normal: float
) -> float:
    """``f_n * (mean_t - mean_n)`` -- the normal-weighted intrinsic term.

    The affine shortcut described in the module docstring. Pinned against
    ``kitagawa.decompose`` by test, not by comment.
    """
    return float(frac_mature_normal) * (float(mean_tumour) - float(mean_normal))


def _unbiased_var(values: np.ndarray) -> float:
    """Sample variance with ddof=1, or 0.0 for a single observation."""
    return float(values.var(ddof=1)) if values.size > 1 else 0.0


def delta_method_se(
    normal: np.ndarray, tumour: np.ndarray, *, frac_mature_normal: float
) -> float:
    """Welch standard error of the scaled difference of means.

    ``f_n * sqrt(s_t^2/n_t + s_n^2/n_n)``. Zero when both arms are constant,
    which callers must treat as *undefined* rather than as certainty.
    """
    n_n, n_t = normal.size, tumour.size
    if n_n < 1 or n_t < 1:
        return float("nan")
    var = _unbiased_var(tumour) / n_t + _unbiased_var(normal) / n_n
    return float(frac_mature_normal) * float(np.sqrt(var))


# --------------------------------------------------------------------------
# resampling, done once and shared by every bootstrap candidate
# --------------------------------------------------------------------------


class ArmResamples:
    """Bootstrap means and variances for one arm. Computed once, reused.

    Three of the six candidates resample cells, and resampling the 800-cell
    reference arm is most of the cost of the whole study. Doing it once per
    replicate and sharing it is not only faster: it is what makes the
    candidates comparable, because they then see *identical* draws.
    """

    __slots__ = ("means", "variances", "n")

    def __init__(self, values: np.ndarray, *, n_boot: int, rng: np.random.Generator):
        values = np.asarray(values, dtype=float)
        n = values.size
        self.n = n
        if n == 0:
            self.means = np.empty(0)
            self.variances = np.empty(0)
            return
        idx = rng.integers(0, n, size=(n_boot, n))
        drawn = values[idx]
        self.means = drawn.mean(axis=1)
        # ddof=1 to match ``delta_method_se``; the studentised pivot divides by
        # an estimate of the same quantity on each side, and mixing ddof
        # between the two would put a systematic factor in the pivot.
        self.variances = drawn.var(axis=1, ddof=1) if n > 1 else np.zeros(n_boot)


# --------------------------------------------------------------------------
# the candidates
# --------------------------------------------------------------------------


def percentile_interval(
    boot_normal: ArmResamples,
    boot_tumour: ArmResamples,
    *,
    frac_mature_normal: float,
    alpha: float = NOMINAL_ALPHA,
) -> Interval:
    """R0. The status quo: percentiles of the bootstrap draws.

    Reproduces ``interval.within_patient_intrinsic_ci`` exactly given the same
    resamples; that equivalence is asserted in the tests.
    """
    if boot_normal.n == 0 or boot_tumour.n == 0:
        return ABSTAIN
    draws = frac_mature_normal * (boot_tumour.means - boot_normal.means)
    lo, hi = np.percentile(draws, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return float(lo), float(hi)


def studentised_interval(
    normal: np.ndarray,
    tumour: np.ndarray,
    boot_normal: ArmResamples,
    boot_tumour: ArmResamples,
    *,
    frac_mature_normal: float,
    alpha: float = NOMINAL_ALPHA,
) -> Interval:
    """R1. Bootstrap-t: pivot on the studentised statistic, not on the estimate.

    ``t* = (theta* - theta_hat) / SE*`` with ``SE*`` recomputed on each
    resample, then ``theta_hat -/+ q(t*) * SE``. Second-order accurate where the
    percentile method is only first-order, which is the textbook reason to
    expect it to help with skew.

    Abstains when the observed ``SE`` is zero -- both arms constant gives an
    undefined pivot, and a zero-width interval around a non-zero estimate would
    reject the null with probability one for no reason. Resamples whose own
    ``SE*`` is zero are dropped from the quantile and counted against a majority
    rule: if more than half the pivots are degenerate there is no distribution
    left to take a quantile of.
    """
    if boot_normal.n == 0 or boot_tumour.n == 0:
        return ABSTAIN
    se = delta_method_se(normal, tumour, frac_mature_normal=frac_mature_normal)
    if not np.isfinite(se) or se <= 0.0:
        return ABSTAIN

    theta = intrinsic_point(
        float(normal.mean()), float(tumour.mean()), frac_mature_normal=frac_mature_normal
    )
    draws = frac_mature_normal * (boot_tumour.means - boot_normal.means)
    se_star = frac_mature_normal * np.sqrt(
        boot_tumour.variances / boot_tumour.n + boot_normal.variances / boot_normal.n
    )
    usable = se_star > 0
    if usable.sum() <= usable.size // 2:
        return ABSTAIN
    pivots = (draws[usable] - theta) / se_star[usable]
    q_lo, q_hi = np.percentile(pivots, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    # Note the crossed subtraction: the upper pivot quantile sets the LOWER
    # bound. Writing this the intuitive way round is the classic bootstrap-t
    # sign error and it produces an interval that looks plausible and is
    # reflected about the estimate.
    return float(theta - q_hi * se), float(theta - q_lo * se)


def _jackknife_centred(values: np.ndarray) -> np.ndarray:
    """``mean(jackknife) - jackknife`` for the mean, in closed form.

    ``jack_i = (S - x_i)/(n-1)``, so this is ``(x_i - xbar)/(n-1)``. The naive
    ``np.delete`` loop in ``interval_calibration.bca_interval`` is O(n^2) and
    fine at ten patients; here n reaches 800 cells and it is not. The two agree
    to floating point, which is asserted in the tests.
    """
    n = values.size
    if n < 2:
        return np.zeros(n)
    return (values - values.mean()) / (n - 1)


def bca_interval(
    normal: np.ndarray,
    tumour: np.ndarray,
    boot_normal: ArmResamples,
    boot_tumour: ArmResamples,
    *,
    frac_mature_normal: float,
    alpha: float = NOMINAL_ALPHA,
) -> Interval:
    """R2. Bias-corrected and accelerated, two-sample.

    ``z0`` from the share of bootstrap draws below the estimate; acceleration
    from a jackknife pooled over the cells of BOTH arms, which is the two-sample
    form (Efron & Tibshirani §14.3) rather than the one-sample formula applied
    to one side.

    ``interval_calibration`` found BCa *slightly worse* than the plain
    percentile over ten patients, because both corrections are estimated from
    the same ten numbers. Over 50 or 800 cells they are estimated from 50 or 800
    numbers, so that objection does not transfer and this is a real candidate
    rather than a formality.
    """
    if boot_normal.n == 0 or boot_tumour.n == 0:
        return ABSTAIN
    draws = frac_mature_normal * (boot_tumour.means - boot_normal.means)
    theta = intrinsic_point(
        float(normal.mean()), float(tumour.mean()), frac_mature_normal=frac_mature_normal
    )

    n_boot = draws.size
    share_below = float(np.mean(draws < theta))
    # Clamp off the open ends: ppf(0) and ppf(1) are infinite, and an infinite
    # z0 silently returns an extreme order statistic as a bound.
    share_below = min(max(share_below, 1.0 / n_boot), 1.0 - 1.0 / n_boot)
    z0 = float(stats.norm.ppf(share_below))

    # theta is decreasing in the normal arm's cells and increasing in the
    # tumour arm's, so the normal arm's influence values carry a sign flip.
    centred = np.concatenate(
        [-_jackknife_centred(normal), _jackknife_centred(tumour)]
    )
    denom = 6.0 * float(np.sum(centred**2)) ** 1.5
    accel = float(np.sum(centred**3) / denom) if denom > 0 else 0.0

    out = []
    for z in (stats.norm.ppf(alpha / 2), stats.norm.ppf(1 - alpha / 2)):
        shifted = z0 + z
        adjusted = z0 + shifted / (1 - accel * shifted)
        if not np.isfinite(adjusted):
            return ABSTAIN
        out.append(float(np.percentile(draws, 100 * stats.norm.cdf(adjusted))))
    lo, hi = out
    if hi < lo:
        return ABSTAIN
    return lo, hi


def welch_t_interval(
    normal: np.ndarray,
    tumour: np.ndarray,
    *,
    frac_mature_normal: float,
    alpha: float = NOMINAL_ALPHA,
) -> Interval:
    """R3. The non-bootstrap reference: Welch's t on the cells.

    ``interval_calibration`` found the Student-t interval calibrated at every
    patient count it measured, where both bootstraps were not. This is the
    cell-level, two-sample analogue of that finding, and it is in the candidate
    list because "the sophisticated method is the wrong one here" has already
    been true once in this repository.

    Welch-Satterthwaite degrees of freedom, so the arms are not assumed to have
    equal variance -- they emphatically do not, since the sweep starves one of
    them.
    """
    n_n, n_t = normal.size, tumour.size
    if n_n < 2 or n_t < 2:
        return ABSTAIN
    v_n, v_t = _unbiased_var(normal) / n_n, _unbiased_var(tumour) / n_t
    se = frac_mature_normal * float(np.sqrt(v_n + v_t))
    if se <= 0.0:
        return ABSTAIN
    denom = v_n**2 / (n_n - 1) + v_t**2 / (n_t - 1)
    if denom <= 0.0:
        return ABSTAIN
    df = (v_n + v_t) ** 2 / denom
    crit = float(stats.t.ppf(1 - alpha / 2, df))
    theta = intrinsic_point(
        float(normal.mean()), float(tumour.mean()), frac_mature_normal=frac_mature_normal
    )
    return theta - crit * se, theta + crit * se


def _cluster_means(
    values: np.ndarray, patients: np.ndarray
) -> tuple[list[np.ndarray], np.ndarray]:
    """Per-patient cell blocks and the distinct patient labels, in a fixed order."""
    labels = np.unique(patients)
    return [values[patients == p] for p in labels], labels


def patient_cluster_interval(
    normal: np.ndarray,
    tumour: np.ndarray,
    normal_patients: np.ndarray,
    tumour_patients: np.ndarray,
    *,
    frac_mature_normal: float,
    n_boot: int = DEFAULT_N_BOOT,
    rng: np.random.Generator,
    alpha: float = NOMINAL_ALPHA,
) -> Interval:
    """R4. Cluster bootstrap: resample PATIENTS, not cells.

    Within each arm the distinct patients are resampled with replacement and
    all of a resampled patient's cells come with it. This is the candidate the
    reviewer names, and it is the one CLAUDE.md invariant 5 would reach for --
    but invariant 5 is about a claim over a *population* of patients, and the
    number of patients available inside one synthetic sample is the number held
    out by the generator, which the committed sweep sets to **two**.

    A bootstrap over two clusters has three reachable resamples per arm. That is
    not a small sample, it is a categorically different object, and if this
    candidate fails that is the reason rather than anything about cluster
    bootstraps. The driver runs it at ``n_held_out`` of 2 and 5 so the two
    explanations are separated by measurement.

    The arms are resampled independently, because in this generator they *are*
    independent draws from the held-out pool.

    Abstains with fewer than two distinct patients in either arm.
    """
    blocks_n, labels_n = _cluster_means(normal, normal_patients)
    blocks_t, labels_t = _cluster_means(tumour, tumour_patients)
    if labels_n.size < 2 or labels_t.size < 2:
        return ABSTAIN

    draws = np.empty(n_boot)
    k_n, k_t = labels_n.size, labels_t.size
    pick_n = rng.integers(0, k_n, size=(n_boot, k_n))
    pick_t = rng.integers(0, k_t, size=(n_boot, k_t))
    for b in range(n_boot):
        m_n = np.concatenate([blocks_n[i] for i in pick_n[b]]).mean()
        m_t = np.concatenate([blocks_t[i] for i in pick_t[b]]).mean()
        draws[b] = frac_mature_normal * (m_t - m_n)
    lo, hi = np.percentile(draws, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return float(lo), float(hi)


def pseudobulk_patient_t_interval(
    normal: np.ndarray,
    tumour: np.ndarray,
    normal_patients: np.ndarray,
    tumour_patients: np.ndarray,
    *,
    frac_mature_normal: float,
    alpha: float = NOMINAL_ALPHA,
) -> Interval:
    """R5. Pseudobulk: one number per patient per arm, then Student-t.

    Each patient contributes the mean of their own mature cells in each arm;
    the interval is a **paired** Student-t on the per-patient differences,
    scaled by ``f_n``. Paired is the right form here and not a choice: both arms
    are drawn from the same held-out patients, so a patient appears on both
    sides and pairing removes the between-patient level that is common to them.

    This is the most defensible pseudobulk reading of the estimand -- it makes
    the patient the unit, which is what CLAUDE.md invariant 5 asks for -- and it
    is also the one most exposed to the cohort size, because the number of
    patients inside one sample is the generator's ``n_held_out``. At two, the
    critical value is ``t(1) = 12.7``.

    Abstains with fewer than two patients present in both arms.
    """
    shared = np.intersect1d(np.unique(normal_patients), np.unique(tumour_patients))
    if shared.size < 2:
        return ABSTAIN
    diffs = np.array(
        [
            tumour[tumour_patients == p].mean() - normal[normal_patients == p].mean()
            for p in shared
        ]
    )
    scaled = frac_mature_normal * diffs
    n = scaled.size
    se = float(scaled.std(ddof=1)) / np.sqrt(n)
    if not np.isfinite(se) or se <= 0.0:
        return ABSTAIN
    crit = float(stats.t.ppf(1 - alpha / 2, n - 1))
    mean = float(scaled.mean())
    return mean - crit * se, mean + crit * se


# --------------------------------------------------------------------------
# all six on one sample
# --------------------------------------------------------------------------


def all_intervals(
    normal: np.ndarray,
    tumour: np.ndarray,
    *,
    frac_mature_normal: float,
    seed: int,
    n_boot: int = DEFAULT_N_BOOT,
    alpha: float = NOMINAL_ALPHA,
    normal_patients: np.ndarray | None = None,
    tumour_patients: np.ndarray | None = None,
    candidates: Sequence[str] = CANDIDATES,
) -> dict[str, Interval]:
    """Every requested candidate on ONE sample, sharing one set of resamples.

    Sharing is the point. Two candidates measured on different draws differ by
    Monte-Carlo noise as well as by method, and at the rates in question --
    5% against 9.5% -- that noise is the same size as the effect.

    Cluster candidates abstain rather than raise when patient labels are not
    supplied, so a caller may request the whole list against a sample that
    cannot support all of it, and the abstention is then visible in the results
    instead of being an exception nobody sees.
    """
    unknown = [c for c in candidates if c not in CANDIDATES]
    if unknown:
        raise ValueError(f"unknown candidate(s) {unknown}; known: {list(CANDIDATES)}")

    normal = np.asarray(normal, dtype=float)
    tumour = np.asarray(tumour, dtype=float)
    rng = np.random.default_rng(seed)

    out: dict[str, Interval] = {}
    if normal.size == 0 or tumour.size == 0:
        return {c: ABSTAIN for c in candidates}

    needs_cells = {"percentile", "studentised", "bca"} & set(candidates)
    boot_n = boot_t = None
    if needs_cells:
        boot_n = ArmResamples(normal, n_boot=n_boot, rng=rng)
        boot_t = ArmResamples(tumour, n_boot=n_boot, rng=rng)

    for name in candidates:
        if name == "percentile":
            out[name] = percentile_interval(
                boot_n, boot_t, frac_mature_normal=frac_mature_normal, alpha=alpha
            )
        elif name == "studentised":
            out[name] = studentised_interval(
                normal, tumour, boot_n, boot_t,
                frac_mature_normal=frac_mature_normal, alpha=alpha,
            )
        elif name == "bca":
            out[name] = bca_interval(
                normal, tumour, boot_n, boot_t,
                frac_mature_normal=frac_mature_normal, alpha=alpha,
            )
        elif name == "welch_t":
            out[name] = welch_t_interval(
                normal, tumour, frac_mature_normal=frac_mature_normal, alpha=alpha
            )
        elif name == "patient_cluster":
            out[name] = (
                ABSTAIN
                if normal_patients is None or tumour_patients is None
                else patient_cluster_interval(
                    normal, tumour,
                    np.asarray(normal_patients), np.asarray(tumour_patients),
                    frac_mature_normal=frac_mature_normal,
                    n_boot=n_boot, rng=rng, alpha=alpha,
                )
            )
        elif name == "pseudobulk_patient_t":
            out[name] = (
                ABSTAIN
                if normal_patients is None or tumour_patients is None
                else pseudobulk_patient_t_interval(
                    normal, tumour,
                    np.asarray(normal_patients), np.asarray(tumour_patients),
                    frac_mature_normal=frac_mature_normal, alpha=alpha,
                )
            )
    return out


def excludes_zero(interval: Interval) -> bool | None:
    """Whether the interval excludes zero. ``None`` for an abstention.

    ``None``, not ``False``. An abstention did not fail to reject; it declined
    to answer, and ``calibration.coverage_and_discrimination`` documents at
    length what happens when those two are folded together.
    """
    lo, hi = interval
    if lo is None or hi is None:
        return None
    return bool(lo > 0 or hi < 0)


def covers(interval: Interval, truth: float) -> bool | None:
    """Whether the interval covers ``truth``. ``None`` for an abstention."""
    lo, hi = interval
    if lo is None or hi is None:
        return None
    return bool(lo <= truth <= hi)


def width(interval: Interval) -> float:
    """Interval width, or ``nan`` for an abstention."""
    lo, hi = interval
    if lo is None or hi is None:
        return float("nan")
    return float(hi - lo)


def mcse(rate: float, n: int) -> float:
    """Binomial Monte-Carlo standard error of a measured rate.

    Reported beside every rate in every table this module feeds. A rate without
    one is not a measurement, and the difference this study exists to resolve --
    5% against 9.5% -- is under two standard errors at the 200 replicates the
    committed diagnostic used.
    """
    if n <= 0:
        return float("nan")
    return float(np.sqrt(max(rate, 0.0) * max(1.0 - rate, 0.0) / n))
