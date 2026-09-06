"""An interval that overstates its own confidence, and how to tell.

WHY THIS EXISTS. Every reading in this project ends in a percentile bootstrap
over patients: ``premise_holds``, ``summarise``, ``specificity`` and
``control_log2_interval`` all resample the per-patient values ``N_BOOTSTRAP``
times and take the 2.5th and 97.5th percentiles of the mean. At n=44 patients
that is fine. The MLH1 positive control would run at **n=10**, and at n=10 it is
not fine: under a true null -- both arms drawn from the same distribution, no
effect of any kind -- the interval excludes zero about **9-11%** of the time
against a nominal 5%.

This is the same shape as the defect the repository is named after, one layer
over. A guard that cannot fail reports success; an interval that is narrower
than it claims reports significance. Both turn an absence of evidence into a
green light, and both do it silently, because nothing raises. The difference is
that this one is measurable in advance, which is what this module does.

WHAT WAS TRIED AND DID NOT WORK. The obvious repairs are the sophisticated
ones, and they are worse:

* **BCa** -- bias-corrected and accelerated -- is *slightly worse* than the
  plain percentile at n=10. Its bias correction ``z0`` and acceleration ``a``
  are themselves estimated from the same ten numbers, and a jackknife
  acceleration over ten points is not a stable quantity.
* **More bootstrap replicates do not help.** ``N_BOOTSTRAP`` is 10,000 already.
  The error is not Monte Carlo error in the resampling; it is that the
  bootstrap distribution of a mean over ten skewed values is not the sampling
  distribution of that mean. Raising B estimates the wrong thing more precisely.

The plain **Student-t interval is well calibrated everywhere measured** -- 4-5%
at n=10, 20 and 44, at both a rare gene's abundance and a common one's. The
sophisticated method is the wrong one here, and that is not a general statement
about bootstraps: it is a statement about this n.

THE SECOND DEFECT, WHICH IS SUBTLER AND IS WHY ``power_curve`` IS SHAPED AS IT
IS. A power figure and a false-positive rate are properties of an *interval
method*, not of a design. It is possible -- and it happened, in the analysis
this module was written to support -- to quote power from the percentile
bootstrap while planning to report the Student-t interval. That overstates
power, because the interval that fires 11% of the time under the null also
fires more often under the alternative. The two numbers must come from the same
method or neither means anything.

So ``power_curve`` does not let them be separated: every row it emits carries
the false-positive rate of **its own method, measured in the same call, on the
same generator**. There is no code path that produces a power number without
its calibration attached, and ``check_power_carries_its_own_calibration``
refuses a frame where they have drifted apart. The failing input is committed
in ``tests/test_checks_can_fail.py``.

THE GENERATIVE MODEL, stated once. A cell is detected when at least one of its
``mu`` expected UMIs lands, so ``p = 1 - exp(-mu)`` -- the same Poisson thinning
``detection_scale`` inverts, used forwards. Two things are modelled that the
first version of this simulation left out:

1. **Per-patient cell counts and depths are taken from the real cohort**, not
   idealised. Pelka's methylated arm runs 79 to 460 mature cells; a simulation
   at the median is a simulation of a cohort that does not exist, and it is
   optimistic because the noisiest patients are the ones dropped.
2. **Between-patient heterogeneity, ``tau``.** The per-patient true log fold
   change is drawn from ``Normal(log(fc), tau^2)`` rather than being identical
   for everyone. A power calculation whose only noise is binomial cannot come
   out underpowered from patient-to-patient variation, because there is none in
   it -- which is a check that cannot fail, again. ``heterogeneity_tau``
   measures ``tau`` from committed per-patient tables instead of assuming it.

OWNERSHIP. CONTRIBUTING §2: this is W1's, and it lives here rather than in
``src/harness/`` because what it calibrates is W1's statistic -- the per-patient
bootstrap in ``coexpression_silencing`` and ``specificity``. It has no
dependency on the harness and the harness has none on it.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
import pandas as pd
from scipy import stats

from src.reference.detection_scale import cloglog_rate

#: The confidence level every interval in this project claims.
NOMINAL_ALPHA = 0.05

#: How far a measured false-positive rate may sit above nominal before the
#: interval is called MISCALIBRATED.
#:
#: 0.02 -- so 7% against a nominal 5%. At the replicate counts used here the
#: Monte Carlo standard error on a 5% rate is about 0.5 percentage points, so
#: this bar sits roughly 4 SE above nominal and is not tripped by simulation
#: noise. It is deliberately loose: the point is to catch an interval that is
#: wrong by a factor, not to police the third decimal.
MISCALIBRATION_TOLERANCE = 0.02

#: Bootstrap replicates inside one interval. Smaller than the 10,000 the
#: analysis jobs use, because this measures a rate over many trials and the
#: per-interval Monte Carlo error averages out across them. Raising it does not
#: move the calibration verdict -- the error is not Monte Carlo error. See the
#: module docstring.
DEFAULT_N_BOOT = 2_000

#: Trials per calibration cell. 1,500 gives ~0.5pp standard error at 5%.
DEFAULT_N_TRIALS = 1_500


class CalibrationError(ValueError):
    """A calibration or power measurement that cannot be read as stated."""


# ---------------------------------------------------------------------------
# Interval methods. One signature, so a caller cannot quietly use two.
# ---------------------------------------------------------------------------


def percentile_interval(
    values: np.ndarray, *, rng: np.random.Generator,
    alpha: float = NOMINAL_ALPHA, n_boot: int = DEFAULT_N_BOOT,
) -> tuple[float, float]:
    """The interval this project uses everywhere. Resample, take percentiles."""
    values = np.asarray(values, dtype=float)
    draws = rng.choice(values, size=(n_boot, values.size), replace=True).mean(axis=1)
    lo, hi = np.percentile(draws, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return float(lo), float(hi)


def bca_interval(
    values: np.ndarray, *, rng: np.random.Generator,
    alpha: float = NOMINAL_ALPHA, n_boot: int = DEFAULT_N_BOOT,
) -> tuple[float, float]:
    """Bias-corrected and accelerated. The textbook repair, measured not assumed.

    Both corrections are estimated from the same sample the interval is about.
    At n=10 the jackknife acceleration is a function of ten leave-one-out means
    and is not stable, which is why this comes out no better than the plain
    percentile. Kept because "we tried the sophisticated one" has to be a
    measurement rather than a claim.
    """
    values = np.asarray(values, dtype=float)
    n = values.size
    theta = float(values.mean())
    draws = rng.choice(values, size=(n_boot, n), replace=True).mean(axis=1)

    share_below = float(np.mean(draws < theta))
    # Clamp off the open ends: ppf(0) and ppf(1) are infinite, and an infinite
    # z0 would silently return the extreme order statistic as a bound.
    share_below = min(max(share_below, 1.0 / n_boot), 1.0 - 1.0 / n_boot)
    z0 = float(stats.norm.ppf(share_below))

    jack = np.array([np.delete(values, i).mean() for i in range(n)])
    centred = jack.mean() - jack
    denom = 6.0 * float(np.sum(centred ** 2)) ** 1.5
    accel = float(np.sum(centred ** 3) / denom) if denom > 0 else 0.0

    out = []
    for z in (stats.norm.ppf(alpha / 2), stats.norm.ppf(1 - alpha / 2)):
        adjusted = z0 + (z0 + z) / (1 - accel * (z0 + z))
        out.append(float(np.percentile(draws, 100 * stats.norm.cdf(adjusted))))
    return out[0], out[1]


def student_t_interval(
    values: np.ndarray, *, rng: np.random.Generator | None = None,
    alpha: float = NOMINAL_ALPHA, n_boot: int = DEFAULT_N_BOOT,
) -> tuple[float, float]:
    """Mean ± t·SE. No resampling; ``rng`` is accepted and ignored.

    The signature matches the others so that the method can be selected by name
    from a table. A measurement that switches methods must not also switch call
    shapes -- that is how a comparison ends up not being one.
    """
    values = np.asarray(values, dtype=float)
    n = values.size
    if n < 2:
        raise CalibrationError(
            f"a t-interval needs at least 2 observations, got {n}. "
            f"With one patient the standard error is undefined and any interval "
            f"reported would be a statement about nothing."
        )
    se = float(values.std(ddof=1)) / np.sqrt(n)
    crit = float(stats.t.ppf(1 - alpha / 2, n - 1))
    mean = float(values.mean())
    return mean - crit * se, mean + crit * se


INTERVAL_METHODS: dict[str, Callable[..., tuple[float, float]]] = {
    "percentile": percentile_interval,
    "bca": bca_interval,
    "student_t": student_t_interval,
}

#: What the MLH1 reading reports. Chosen by the measurement in this module, not
#: by convention: it is the only one of the three calibrated at n=10.
CALIBRATED_METHOD = "student_t"


def excludes_zero(interval: tuple[float, float]) -> bool:
    lo, hi = interval
    return bool(lo > 0 or hi < 0)


# ---------------------------------------------------------------------------
# The generative model
# ---------------------------------------------------------------------------


def simulate_deltas(
    *, n_cells: np.ndarray, depth: np.ndarray, cp10k: float,
    fold_change: float, tau: float, rng: np.random.Generator,
) -> np.ndarray:
    """Per-patient cloglog detection deltas under a stated truth.

    Parameters
    ----------
    n_cells, depth:
        Per-patient mature cells per arm and median UMIs per cell, **taken from
        the real cohort**. Both arms get the same values because the analysis
        depth-matches them, which is what makes that legitimate here.
    cp10k:
        The gene's per-cell mean in the normal arm, counts per 10k. This is what
        sets the detection rate and therefore the whole power problem.
    fold_change:
        Tumour mean over normal mean. 1.0 is the null.
    tau:
        Between-patient SD of the true log fold change. 0.0 is the homogeneous
        model, which is optimistic; see ``heterogeneity_tau``.

    Returns the per-patient deltas on the log fold-change scale, through the
    same ``cloglog_rate`` the analysis uses -- boundary rule included, so the
    attenuation it induces at low counts is in the measurement rather than
    argued around.
    """
    n_cells = np.asarray(n_cells, dtype=int)
    depth = np.asarray(depth, dtype=float)
    if n_cells.shape != depth.shape:
        raise CalibrationError(
            f"n_cells has {n_cells.shape} and depth has {depth.shape}. These are "
            f"per-patient vectors of the same cohort; a mismatch means one of "
            f"them is not the cohort you think it is."
        )
    if tau < 0:
        raise CalibrationError(f"tau={tau} is a standard deviation")

    mu_normal = cp10k / 1e4 * depth
    log_fc = np.full(n_cells.size, np.log(fold_change))
    if tau > 0:
        log_fc = log_fc + rng.normal(0.0, tau, n_cells.size)

    p_normal = 1.0 - np.exp(-mu_normal)
    p_tumour = 1.0 - np.exp(-mu_normal * np.exp(log_fc))
    k_normal = rng.binomial(n_cells, p_normal)
    k_tumour = rng.binomial(n_cells, p_tumour)
    return np.asarray(
        cloglog_rate(k_tumour / n_cells, n_cells)
        - cloglog_rate(k_normal / n_cells, n_cells)
    )


def rejection_rate(
    sampler: Callable[[np.random.Generator], np.ndarray],
    *, method: str, rng: np.random.Generator,
    n_trials: int = DEFAULT_N_TRIALS, alpha: float = NOMINAL_ALPHA,
) -> float:
    """Share of trials whose interval excludes zero.

    Under a null sampler this is the false-positive rate; under an alternative
    it is power. Deliberately one function: they are the same measurement and
    splitting them is how they end up computed with different methods.
    """
    if method not in INTERVAL_METHODS:
        raise CalibrationError(
            f"unknown interval method {method!r}; have {sorted(INTERVAL_METHODS)}"
        )
    fn = INTERVAL_METHODS[method]
    hits = 0
    for _ in range(n_trials):
        values = sampler(rng)
        hits += excludes_zero(fn(values, rng=rng, alpha=alpha))
    return hits / n_trials


# ---------------------------------------------------------------------------
# Calibration
# ---------------------------------------------------------------------------


def width_ratio(n_patients: int, alpha: float = NOMINAL_ALPHA) -> float:
    """How much narrower the percentile bootstrap is than the correct interval.

    THE MECHANISM, WHICH IS ARITHMETIC AND CONTAINS NO DATA. The bootstrap
    distribution of a mean is asymptotically normal with standard deviation
    ``s * sqrt((n-1)/n)`` -- the plug-in, divide-by-n standard deviation, not
    the unbiased one -- so the percentile interval is approximately

        mean  ±  z * s * sqrt((n-1)/n) / sqrt(n)

    while the interval that actually has 95% coverage for a mean estimated with
    an unknown variance is

        mean  ±  t(n-1) * s / sqrt(n)

    The ratio is ``z * sqrt((n-1)/n) / t(n-1)``: a function of ``n`` alone. It
    is 0.53 at n=4, 0.82 at n=10, 0.91 at n=20, 0.96 at n=44 and 0.98 at n=100.

    So the percentile bootstrap does not fail because a transcript is rare, or
    because ``cloglog`` is nonlinear, or because of anything about this project.
    It fails because it uses a normal quantile where a t quantile is needed and
    a biased standard deviation where an unbiased one is needed, and both errors
    run the same way. Everything else in this module is a check that the
    arithmetic describes the real thing.
    """
    if n_patients < 2:
        raise CalibrationError(
            f"n={n_patients}: neither interval is defined on fewer than two "
            f"observations, so there is no ratio between them."
        )
    z = float(stats.norm.ppf(1 - alpha / 2))
    t = float(stats.t.ppf(1 - alpha / 2, n_patients - 1))
    return z * np.sqrt((n_patients - 1) / n_patients) / t


def expected_false_positive_rate(n_patients: int,
                                 alpha: float = NOMINAL_ALPHA) -> float:
    """Closed-form false-positive rate of the percentile bootstrap of a mean.

    ``P(|T| > z * sqrt((n-1)/n))`` for ``T ~ t(n-1)``: 18.8% at n=4, 9.6% at
    n=10, 7.1% at n=20, 5.9% at n=44.

    **This is a floor, not a prediction.** It assumes the per-patient values are
    normal enough that the bootstrap mean is; a rare transcript's per-patient
    delta is also skewed, and the measured rates sit slightly above the closed
    form at small n for that reason. The two agreeing to within a point or so is
    what says the simulation and the arithmetic are describing one phenomenon --
    and where they diverge, the divergence is the skew and is worth reading.
    """
    z = float(stats.norm.ppf(1 - alpha / 2))
    return float(2 * stats.t.sf(z * np.sqrt((n_patients - 1) / n_patients),
                                n_patients - 1))


def calibration_verdict(false_positive_rate: float,
                        alpha: float = NOMINAL_ALPHA) -> str:
    """CALIBRATED / MISCALIBRATED / CONSERVATIVE, on the measured rate."""
    if false_positive_rate > alpha + MISCALIBRATION_TOLERANCE:
        return "MISCALIBRATED"
    if false_positive_rate < alpha - MISCALIBRATION_TOLERANCE:
        return "CONSERVATIVE"
    return "CALIBRATED"


def calibration_table(
    *, cohorts: dict[str, tuple[np.ndarray, np.ndarray]],
    abundances: dict[str, float], taus: tuple[float, ...] = (0.0, 0.2),
    methods: tuple[str, ...] = tuple(INTERVAL_METHODS),
    seed: int, n_trials: int = DEFAULT_N_TRIALS,
) -> pd.DataFrame:
    """False-positive rate of each interval method under a true null.

    ``cohorts`` maps a label to that cohort's ``(n_cells, depth)`` vectors, so
    the sweep runs at the n values the project actually reports at rather than
    at round numbers.
    """
    rows = []
    for cohort, (n_cells, depth) in cohorts.items():
        for gene_label, cp10k in abundances.items():
            for tau in taus:
                for method in methods:
                    rng = np.random.default_rng(seed)

                    def null(r, *, n=n_cells, d=depth, c=cp10k, t=tau):
                        return simulate_deltas(
                            n_cells=n, depth=d, cp10k=c,
                            fold_change=1.0, tau=t, rng=r,
                        )

                    rate = rejection_rate(
                        null, method=method, rng=rng, n_trials=n_trials
                    )
                    n = int(len(n_cells))
                    rows.append({
                        "cohort": cohort,
                        "n_patients": n,
                        "median_cells_per_arm": float(np.median(n_cells)),
                        "abundance": gene_label,
                        "cp10k_normal": float(cp10k),
                        "tau": float(tau),
                        "method": method,
                        "nominal_alpha": NOMINAL_ALPHA,
                        "false_positive_rate": rate,
                        "verdict": calibration_verdict(rate),
                        # The closed form applies to the percentile bootstrap
                        # only. Carried on every row so a reader can see that
                        # the simulated rate is not a lone number: where the two
                        # agree the mechanism is understood, and where they
                        # part the gap is the skew a normal approximation drops.
                        "closed_form_rate": (
                            expected_false_positive_rate(n)
                            if method == "percentile" else float("nan")
                        ),
                        "width_ratio_vs_t": (
                            width_ratio(n) if method == "percentile"
                            else float("nan")
                        ),
                        "n_trials": int(n_trials),
                    })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Power, which may not be separated from its calibration
# ---------------------------------------------------------------------------


def power_curve(
    *, n_cells: np.ndarray, depth: np.ndarray, cp10k: float,
    fold_changes: tuple[float, ...], taus: tuple[float, ...],
    method: str = CALIBRATED_METHOD, seed: int,
    n_trials: int = DEFAULT_N_TRIALS, cohort: str = "",
) -> pd.DataFrame:
    """Power at each fold change, **carrying its own false-positive rate**.

    THE STRUCTURAL POINT. Power and calibration are properties of an interval
    method. Quoting power from one method while reporting another's interval
    overstates the design, because the interval that over-rejects under the null
    also over-rejects under the alternative. That is not hypothetical: the first
    version of the MLH1 power statement read 86% from the percentile bootstrap
    while planning to report the Student-t interval, whose power at the same
    fold change is 74%.

    So there is no way to call this and get a power number alone. Every row
    carries ``false_positive_rate`` measured on the same generator, the same
    method and the same cohort, and ``check_power_carries_its_own_calibration``
    refuses a frame where the two have come apart.
    """
    rows = []
    for tau in taus:
        rng = np.random.default_rng(seed)

        def null(r, *, t=tau):
            return simulate_deltas(
                n_cells=n_cells, depth=depth, cp10k=cp10k,
                fold_change=1.0, tau=t, rng=r,
            )

        fpr = rejection_rate(null, method=method, rng=rng, n_trials=n_trials)

        for fold_change in fold_changes:
            rng = np.random.default_rng(seed)

            def alternative(r, *, fc=fold_change, t=tau):
                return simulate_deltas(
                    n_cells=n_cells, depth=depth, cp10k=cp10k,
                    fold_change=fc, tau=t, rng=r,
                )

            rows.append({
                "cohort": cohort,
                "n_patients": int(len(n_cells)),
                "median_cells_per_arm": float(np.median(n_cells)),
                "cp10k_normal": float(cp10k),
                "tau": float(tau),
                "method": method,
                "fold_change": float(fold_change),
                "silencing_pct": float(100 * (1 - fold_change)),
                "power": rejection_rate(
                    alternative, method=method, rng=rng, n_trials=n_trials
                ),
                "false_positive_rate": fpr,
                "calibration_verdict": calibration_verdict(fpr),
                "n_trials": int(n_trials),
            })
    return pd.DataFrame(rows)


def check_power_carries_its_own_calibration(frame: pd.DataFrame) -> None:
    """Refuse a power table whose power and calibration are from different runs.

    THE GUARD FOR THE DEFECT THIS MODULE EXISTS TO PREVENT. A power figure is
    only interpretable beside the false-positive rate of the same interval, and
    the two are easy to assemble from different sources -- one from a
    simulation someone ran last week, one from a table. This asserts they were
    not: within a (cohort, tau, cp10k) cell every row must name one method and
    one false-positive rate.

    A frame that passes this is not thereby correct; it is merely not
    self-inconsistent in the one way that has already happened.
    """
    required = {"cohort", "tau", "cp10k_normal", "method",
                "power", "false_positive_rate"}
    missing = required - set(frame.columns)
    if missing:
        raise CalibrationError(
            f"a power table must carry {sorted(required)}; missing "
            f"{sorted(missing)}. A power number without the false-positive rate "
            f"of the interval that produced it is not a claim about a design."
        )
    if frame.empty:
        return
    grouped = frame.groupby(["cohort", "tau", "cp10k_normal"], dropna=False)
    for key, block in grouped:
        if block["method"].nunique() > 1:
            raise CalibrationError(
                f"{key} mixes interval methods {sorted(block['method'].unique())} "
                f"in one power curve. Power and its false-positive rate must "
                f"come from the same interval or neither is readable."
            )
        if block["false_positive_rate"].nunique() > 1:
            raise CalibrationError(
                f"{key} carries {block['false_positive_rate'].nunique()} "
                f"different false-positive rates for one method. The rate is a "
                f"property of the method and the cohort, so more than one means "
                f"they were measured on different generators and the power "
                f"figures beside them are not comparable."
            )


# ---------------------------------------------------------------------------
# Between-patient heterogeneity, measured rather than assumed
# ---------------------------------------------------------------------------


def heterogeneity_tau(
    deltas: pd.DataFrame, *, seed: int, n_boot: int = 4_000,
) -> pd.DataFrame:
    """Per-gene between-patient SD of the delta, net of sampling noise.

    WHY A POWER CALCULATION NEEDS THIS. The observed spread of per-patient
    deltas is sampling noise **plus** real patient-to-patient variation. A
    simulation containing only the first is a simulation of a cohort more
    homogeneous than any that exists, and it reports power that cannot be
    achieved.

    The sampling part is obtained by parametric bootstrap **through
    ``cloglog_rate`` itself**, not through a delta-method formula: the boundary
    rule is a nonlinear function of cell count, and a closed form that ignores
    it understates the variance exactly where the counts are small, which is
    where this matters.

    ``tau`` is the excess, floored at zero. A negative excess means the observed
    spread is smaller than binomial sampling alone predicts -- reported as
    ``tau_squared_raw`` rather than hidden, because a systematically negative
    excess would mean the sampling model is wrong, not that patients agree.
    """
    needed = {"gene", "n_normal", "n_tumour", "detect_normal", "detect_tumour"}
    missing = needed - set(deltas.columns)
    if missing:
        raise CalibrationError(f"heterogeneity needs {sorted(missing)}")

    rng = np.random.default_rng(seed)
    rows = []
    for gene, block in deltas.groupby("gene", observed=True):
        block = block.dropna(subset=["n_normal", "n_tumour",
                                     "detect_normal", "detect_tumour"])
        if len(block) < 3:
            continue
        observed = np.asarray(
            cloglog_rate(block["detect_tumour"], block["n_tumour"])
            - cloglog_rate(block["detect_normal"], block["n_normal"])
        )
        observed_var = float(np.var(observed, ddof=1))

        per_patient = []
        for _, row in block.iterrows():
            n_n, n_t = int(row["n_normal"]), int(row["n_tumour"])
            k_n = rng.binomial(n_n, float(row["detect_normal"]), n_boot)
            k_t = rng.binomial(n_t, float(row["detect_tumour"]), n_boot)
            drawn = (cloglog_rate(k_t / n_t, np.full(n_boot, n_t))
                     - cloglog_rate(k_n / n_n, np.full(n_boot, n_n)))
            per_patient.append(float(np.var(drawn, ddof=1)))
        sampling_var = float(np.mean(per_patient))

        excess = observed_var - sampling_var
        rows.append({
            "gene": str(gene),
            "n_patients": int(len(block)),
            "mean_delta_cloglog": float(observed.mean()),
            "baseline_detection": float(block["detect_normal"].mean()),
            "observed_sd": float(np.sqrt(observed_var)),
            "sampling_sd": float(np.sqrt(sampling_var)),
            "tau_squared_raw": excess,
            "tau": float(np.sqrt(excess)) if excess > 0 else 0.0,
        })
    return pd.DataFrame(rows).sort_values("baseline_detection",
                                          ignore_index=True)
