"""Per-patient intervals for the intrinsic term. W2.

WHY THIS IS NOT A VIOLATION OF INVARIANT 5
------------------------------------------
CLAUDE.md invariant 5 says *bootstrap over patients, not cells*, and this module
resamples cells. Read the estimand before reading that as a contradiction —
they are two different questions:

==================================  =========================  =================
Question                            Sampling unit              Lives in
==================================  =========================  =================
"What is the cohort's intrinsic     **patients**               W4's
loss, and how sure are we?"         (invariant 5)              ``bootstrap_over_patients``
"Does *this* patient have enough    **cells**, by              here
mature cells for their own          construction
estimate to mean anything?"
==================================  =========================  =================

Invariant 5 exists because resampling cells to make a *population* claim inflates
n by the number of cells per patient and produces intervals wrong by an order of
magnitude, in the flattering direction. That failure mode requires the claim to
be about patients. The positivity cutpoint is not: it asks whether one patient's
own per-cell mean is pinned down by the cells that patient has, and there the
cells *are* the sample.

This matters concretely. W4's ``attach_intrinsic_ci`` broadcasts a cohort-level
band onto every patient row and says so explicitly — a patient with 800 mature
cells and one with 21 get the **same** interval. Coverage and discrimination
against ``n_cells_mature`` would then be flat, with no crossing point, and
``calibration.calibrate_cutpoints`` could never fire. A cutpoint on cell count
needs an interval that responds to cell count. This is that interval; its width
scales as roughly 1/sqrt(n_mature).

Both intervals are real and they answer different questions. Neither replaces
the other, and a results table must say which one it is carrying.
"""

from __future__ import annotations

import numpy as np

from src.estimator.kitagawa import decompose

#: Enough for a stable 95% percentile bound without dominating a sweep.
DEFAULT_N_BOOT = 400


def within_patient_intrinsic_ci(
    mature_normal: np.ndarray,
    mature_tumour: np.ndarray,
    *,
    frac_mature_normal: float,
    frac_mature_tumour: float,
    n_boot: int = DEFAULT_N_BOOT,
    seed: int,
    weighting: str = "normal",
    alpha: float = 0.05,
) -> tuple[float | None, float | None]:
    """Percentile interval on one patient's intrinsic term.

    Parameters
    ----------
    mature_normal, mature_tumour:
        Target-gene expression in the mature cells of each arm, one value per
        cell. These are the only cells that carry information about the
        per-cell mean, which is why the interval narrows with their number and
        why the cutpoint is defined on that number.
    frac_mature_normal, frac_mature_tumour:
        Held fixed across resamples. The compositional term's uncertainty is a
        different question and the cutpoint is not defined on it.

    Returns
    -------
    ``(ci_low, ci_high)``, or ``(None, None)`` when either arm has no mature
    cells. **Undefined, not zero** — CLAUDE.md invariant 1. A caller that turns
    this into ``(0.0, 0.0)`` has produced the exact error the third segment
    exists to prevent.
    """
    mature_normal = np.asarray(mature_normal, dtype=float)
    mature_tumour = np.asarray(mature_tumour, dtype=float)
    if n_boot < 1:
        raise ValueError(f"n_boot={n_boot} must be positive")
    if not 0.0 < alpha < 1.0:
        raise ValueError(f"alpha={alpha} must be in (0, 1)")
    if mature_normal.size == 0 or mature_tumour.size == 0:
        return None, None

    rng = np.random.default_rng(seed)
    n_n, n_t = mature_normal.size, mature_tumour.size

    # Resample each arm at its own n — the tumour arm is the one the sweep
    # starves, and its n is what classify_estimability() reads.
    means_n = mature_normal[rng.integers(0, n_n, size=(n_boot, n_n))].mean(axis=1)
    means_t = mature_tumour[rng.integers(0, n_t, size=(n_boot, n_t))].mean(axis=1)

    draws = np.empty(n_boot)
    for i in range(n_boot):
        # Through the real estimator, not a re-derivation of its formula. If
        # W4 changes decompose(), this interval changes with it, which is the
        # behaviour we want from something that calibrates W4's cutpoints.
        draws[i] = decompose(
            frac_mature_normal,
            frac_mature_tumour,
            float(means_n[i]),
            float(means_t[i]),
            n_cells_mature=n_t,
            weighting=weighting,
        ).intrinsic

    lo_q, hi_q = 100 * (alpha / 2), 100 * (1 - alpha / 2)
    return float(np.percentile(draws, lo_q)), float(np.percentile(draws, hi_q))


def ci_width(ci_low: float | None, ci_high: float | None) -> float:
    """Interval width, or ``nan`` when the interval is undefined."""
    if ci_low is None or ci_high is None:
        return float("nan")
    return float(ci_high - ci_low)
