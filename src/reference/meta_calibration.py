"""A weight that has no average, and the heterogeneity verdict built on it.

WHY THIS EXISTS. ``src/harness/meta.py`` pools per-study estimates by
inverse-variance weighting and then reads a verdict off ``I^2`` against a
pre-committed 75% ceiling. The weights come from
``coexpression_meta.per_study_stats``, which is ``sd(ddof=1)/sqrt(n)`` over
patients -- and in the ICBI meta ``n`` runs down to **three**.

DerSimonian-Laird treats those standard errors as known constants. They are
estimated, and what is estimated is a variance, so the weight is an INVERSE
chi-square::

    w_hat / w_true  =  1 / (chi2(n-1) / (n-1))

which gives, exactly::

    E[w_hat] / w_true  =  (n - 1) / (n - 3)      finite only for n >= 4
    Var[w_hat]                                    finite only for n >= 6

**At n = 3 the inverse-variance weight has no finite mean.** Not "is imprecise"
-- ``E[1/chi2(2)]`` diverges, and a simulation of it does not converge: the
running mean over 10^4, 10^5, 10^6, 10^7 draws wanders 9.7, 22.2, 13.0, 16.5.
``weight_inflation`` returns ``inf`` there rather than a number, because there
is no number.

THE FLOOR WAS SET BY THE WRONG CRITERION. ``MIN_PREMISE_PATIENTS = 3`` in
``coexpression_silencing.py`` is documented as "3 needed to tell a shift from
noise" -- a statement about whether the ESTIMATE exists. Nothing asked whether
the WEIGHT built from it is integrable. The numerator's estimability was used to
admit a denominator. That is the defect in one line, and it is the reason this
module's floors are named after moments of the weight rather than after a count
of patients.

WHAT IT DOES TO THE VERDICT. Under perfect homogeneity -- every study
estimating one common effect, nothing between them but sampling noise -- at the
eleven patient counts the ICBI meta actually has:

* Cochran's Q rejects at **32.5%** against a nominal 5%.
* ``I^2`` has median **0.269**, 90th percentile 0.701, 95th percentile 0.785.
* ``P(I^2 > 0.75)`` is **6.9%**.

Give all eleven studies Pelka's n = 29 and the median null ``I^2`` is **0.000**
and ``P(I^2 > 0.75)`` is 0.02%. **The baseline heterogeneity is manufactured by
the small studies, through the weight rather than through the biology.**

THE SAME SHAPE AS ``interval_calibration``, ONE LAYER UP. There, an interval was
``z*sqrt((n-1)/n)/t(n-1)`` times the width it claimed -- a function of ``n``
alone. Here, a weight is ``(n-1)/(n-3)`` times what it claims -- also a function
of ``n`` alone. Both are arithmetic. Both were invisible because nothing raised.
The difference is direction: that one made a check too eager to find an effect,
this one makes a gate too eager to refuse one.

AND THE GUARD, WHICH IS THE POINT OF THE MODULE. ``MAX_I_SQUARED = 0.75`` is
Higgins' rule of thumb for a literature of trials with hundreds of subjects
each. Applied at n = 3 it is a threshold whose null distribution nobody
computed, which is this repository's signature defect wearing a citation. So a
heterogeneity verdict does not exist here without the null it was read against:
``null_i_squared`` and ``calibrated_p`` are emitted in the same row as the
verdict, and ``check_heterogeneity_carries_its_own_null`` refuses a frame where
they have drifted apart. The failing input is committed in
``tests/test_checks_can_fail.py``.

WHAT THIS DOES NOT DO. It does not drop a study, and it must not be used to.
Khaliq_2022 is the most influential point in the KRT8 meta and removing it by
hand is post-hoc outlier removal. The floors below are fixed by the chi-square
degrees-of-freedom condition and were written into
``docs/prereg_meta_weight_calibration.md`` before any floor curve was run.

OWNERSHIP. CONTRIBUTING §2: ``src/harness/meta.py`` is W2's and this module does
not modify it -- it imports ``meta_analyse`` and reads its output. The floor is
applied at the job layer, which is W1's. Changing ``MIN_STUDIES`` or
``MAX_I_SQUARED`` themselves is a W2 PR with two approvals and is out of scope.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats

log = logging.getLogger(__name__)

__all__ = [
    "MetaCalibrationError",
    "NullHeterogeneity",
    "PATIENT_FLOOR",
    "WEAK_PATIENT_FLOOR",
    "STATUS_QUO_PATIENT_FLOOR",
    "FLOORS_REPORTED",
    "calibrated_p",
    "check_heterogeneity_carries_its_own_null",
    "floor_label",
    "min_patients_for_finite_weight_mean",
    "min_patients_for_finite_weight_variance",
    "null_heterogeneity_table",
    "null_i_squared",
    "weight_inflation",
]


class MetaCalibrationError(ValueError):
    """A heterogeneity number was asked for without the null it means anything against."""


#: Draws per floor for the null simulation. The quantity is a tail probability
#: around 0.02-0.07, so the Monte Carlo SE at 200k is ~5e-4 -- two orders below
#: the differences between floors that the verdict turns on.
DEFAULT_N_TRIALS = 200_000

#: Fixed seed, invariant 10.
DEFAULT_SEED = 20260907

#: **Primary floor.** The smallest n at which ``Var[w_hat]`` is finite. A weight
#: with a finite mean but infinite variance is still a weight that can take any
#: value, and DerSimonian-Laird treats it as a constant -- Q, tau^2 and I^2 all
#: assume it is one. Two moments is the least that assumption needs.
PATIENT_FLOOR = 6

#: The weaker sensitivity: the smallest n at which ``E[w_hat]`` is finite at
#: all. Reported alongside, never switched to.
WEAK_PATIENT_FLOOR = 4

#: What ``coexpression_silencing.MIN_PREMISE_PATIENTS`` currently is. Carried so
#: the change is visible rather than asserted.
STATUS_QUO_PATIENT_FLOOR = 3

#: Reported together, always. ``docs/prereg_meta_weight_calibration.md`` §5:
#: the floor is a modelling choice and a modelling choice is reported as a
#: curve, exactly as ``config/labeling_axes.yaml`` requires of the rung axis.
FLOORS_REPORTED = (STATUS_QUO_PATIENT_FLOOR, WEAK_PATIENT_FLOOR, PATIENT_FLOOR)


def min_patients_for_finite_weight_mean() -> int:
    """Smallest ``n`` with ``E[1/se_hat^2]`` finite.

    ``se_hat^2 = s^2/n`` with ``s^2 ~ sigma^2 chi2(n-1)/(n-1)``, so the weight is
    an inverse chi-square on ``n-1`` df. ``E[1/chi2(v)]`` is finite iff ``v > 2``.
    """
    return 4


def min_patients_for_finite_weight_variance() -> int:
    """Smallest ``n`` with ``Var[1/se_hat^2]`` finite. ``E[1/chi2(v)^2]`` needs ``v > 4``."""
    return 6


def weight_inflation(n_patients: int) -> float:
    """``E[w_hat] / w_true = (n-1)/(n-3)``, or ``inf`` where there is no mean.

    A closed form, a function of ``n`` alone, and the direct analogue of
    ``interval_calibration.width_ratio``. Small studies are over-weighted in
    expectation by this factor -- 3.00x at n=4, 2.00x at n=5, 1.29x at n=10,
    1.08x at n=29 -- and at n<=3 the expectation does not exist.
    """
    if n_patients < 2:
        raise MetaCalibrationError(
            f"n = {n_patients}: an SE over patients needs at least 2."
        )
    if n_patients < min_patients_for_finite_weight_mean():
        return float("inf")
    return float((n_patients - 1) / (n_patients - 3))


@dataclass(frozen=True)
class NullHeterogeneity:
    """The distribution of ``I^2`` when every study estimates the same thing.

    Carried next to any observed ``I^2`` that is about to be compared to a
    threshold. Depends only on the patient counts -- never on the estimates --
    which is what makes it safe to compute before looking at an answer.
    """

    n_patients: tuple[int, ...]
    k: int
    median: float
    q90: float
    q95: float
    p_exceeds_ceiling: float
    q_rejection_rate: float
    n_trials: int
    seed: int

    def as_row(self) -> dict[str, float | int | str]:
        return {
            "k": self.k,
            "null_i_squared_median": self.median,
            "null_i_squared_q90": self.q90,
            "null_i_squared_q95": self.q95,
            "null_p_exceeds_ceiling": self.p_exceeds_ceiling,
            "null_q_rejection_rate": self.q_rejection_rate,
            "null_n_trials": self.n_trials,
            "null_seed": self.seed,
        }


def _simulate_null(
    n_patients: np.ndarray, n_trials: int, seed: int
) -> tuple[np.ndarray, np.ndarray]:
    """``I^2`` and ``Q`` under exact homogeneity, at the given patient counts.

    The generative model is the one DL assumes and nothing more: every study
    estimates a common effect (taken as 0 without loss -- Q and I^2 are
    location-invariant), each patient mean is normal with the study's own
    sampling variance, and each study's ``s^2`` is drawn from its own
    chi-square. Between-study variance is EXACTLY ZERO by construction, so any
    ``I^2`` this returns is manufactured entirely by estimating the weights.

    Unit variance is likewise without loss: Q and I^2 are scale-invariant, so a
    common sigma cancels. Studies differing in sigma is a different question and
    is not what the ceiling is being calibrated against.
    """
    rng = np.random.default_rng(seed)
    ns = np.asarray(n_patients, dtype=int)
    k = ns.size
    means = rng.normal(0.0, 1.0 / np.sqrt(ns), size=(n_trials, k))
    s_squared = rng.chisquare(ns - 1, size=(n_trials, k)) / (ns - 1)
    weights = ns / s_squared                      # = 1 / (s^2 / n)
    fixed = (weights * means).sum(1) / weights.sum(1)
    q = (weights * (means - fixed[:, None]) ** 2).sum(1)
    df = k - 1
    i_squared = np.maximum(0.0, (q - df) / np.where(q > 0, q, np.nan))
    return np.nan_to_num(i_squared), q


def null_i_squared(
    n_patients: "np.ndarray | list[int]",
    *,
    ceiling: float | None = None,
    n_trials: int = DEFAULT_N_TRIALS,
    seed: int = DEFAULT_SEED,
) -> NullHeterogeneity:
    """What ``I^2`` does at these patient counts when there is no heterogeneity.

    ``ceiling`` defaults to ``meta.MAX_I_SQUARED`` -- imported rather than
    restated, so the two cannot drift.
    """
    from src.harness.meta import MAX_I_SQUARED, MIN_STUDIES

    ns = np.asarray(list(n_patients), dtype=int)
    if ns.size < MIN_STUDIES:
        raise MetaCalibrationError(
            f"{ns.size} studies, below meta.MIN_STUDIES = {MIN_STUDIES}. "
            f"There is no verdict to calibrate."
        )
    if (ns < 2).any():
        raise MetaCalibrationError(
            f"patient counts {sorted(ns.tolist())} include one below 2; an SE "
            f"over patients does not exist there."
        )
    bar = MAX_I_SQUARED if ceiling is None else ceiling
    i_squared, q = _simulate_null(ns, n_trials, seed)
    df = ns.size - 1
    return NullHeterogeneity(
        n_patients=tuple(int(v) for v in ns),
        k=int(ns.size),
        median=float(np.median(i_squared)),
        q90=float(np.quantile(i_squared, 0.90)),
        q95=float(np.quantile(i_squared, 0.95)),
        p_exceeds_ceiling=float((i_squared > bar).mean()),
        q_rejection_rate=float((q > stats.chi2.ppf(0.95, df)).mean()),
        n_trials=int(n_trials),
        seed=int(seed),
    )


def calibrated_p(
    observed_i_squared: float,
    n_patients: "np.ndarray | list[int]",
    *,
    n_trials: int = DEFAULT_N_TRIALS,
    seed: int = DEFAULT_SEED,
) -> float:
    """``P(I^2 >= observed)`` under homogeneity at these patient counts.

    The number ``MAX_I_SQUARED`` was standing in for and never supplied. Uses
    the ``(r+1)/(B+1)`` convention so a p of exactly zero is never reported from
    a finite simulation.
    """
    ns = np.asarray(list(n_patients), dtype=int)
    i_squared, _ = _simulate_null(ns, n_trials, seed)
    exceed = int((i_squared >= observed_i_squared).sum())
    return float((exceed + 1) / (n_trials + 1))


def floor_label(floor: int) -> str:
    """Why a floor is where it is, in words, for the report."""
    if floor >= min_patients_for_finite_weight_variance():
        return f"n >= {floor}: Var[w] finite"
    if floor >= min_patients_for_finite_weight_mean():
        return f"n >= {floor}: E[w] finite, Var[w] not"
    return f"n >= {floor}: E[w] DIVERGES"


def null_heterogeneity_table(
    n_patients: "np.ndarray | list[int]",
    floors: "tuple[int, ...]" = FLOORS_REPORTED,
    *,
    n_trials: int = DEFAULT_N_TRIALS,
    seed: int = DEFAULT_SEED,
) -> pd.DataFrame:
    """One row per floor: which studies survive it and what the null looks like.

    Depends on patient counts only. Safe to compute -- and in this project's
    order of operations, computed -- before any estimate is read.
    """
    ns = np.asarray(list(n_patients), dtype=int)
    rows = []
    for floor in floors:
        kept = ns[ns >= floor]
        row: dict[str, object] = {
            "patient_floor": int(floor),
            "floor_rationale": floor_label(floor),
            "weight_mean_finite": bool(floor >= min_patients_for_finite_weight_mean()),
            "weight_variance_finite": bool(
                floor >= min_patients_for_finite_weight_variance()
            ),
            "n_studies_dropped": int((ns < floor).sum()),
        }
        try:
            row.update(null_i_squared(kept, n_trials=n_trials, seed=seed).as_row())
        except MetaCalibrationError as exc:
            row.update({"k": int(kept.size), "null_refused": str(exc)})
        rows.append(row)
    return pd.DataFrame(rows)


#: Columns that make a heterogeneity claim, and the ones that calibrate it.
#: Neither group may appear without the other.
_VERDICT_COLUMNS = ("i_squared", "verdict")
_NULL_COLUMNS = ("null_i_squared_median", "null_p_of_observed")


def check_heterogeneity_carries_its_own_null(frame: pd.DataFrame) -> None:
    """Refuse a frame that states a heterogeneity verdict without its null.

    The exact analogue of
    ``interval_calibration.check_power_carries_its_own_calibration``, and it
    exists for the same reason: ``I^2 = 87.6%`` and "exceeds the 75% ceiling"
    are not interpretable side by side unless the reader is also told that at
    these patient counts homogeneity itself produces a median ``I^2`` of 0.269.
    A frame carrying the verdict and not the null lets a threshold that was
    never calibrated read as one that was.

    Raises ``MetaCalibrationError``; the failing input is committed in
    ``tests/test_checks_can_fail.py``.
    """
    if frame.empty:
        return
    has_verdict = [c for c in _VERDICT_COLUMNS if c in frame.columns]
    has_null = [c for c in _NULL_COLUMNS if c in frame.columns]
    if has_verdict and not has_null:
        raise MetaCalibrationError(
            f"frame states {has_verdict} without any of {list(_NULL_COLUMNS)}. "
            f"A heterogeneity verdict is a comparison against a threshold whose "
            f"null distribution depends on the patient counts; at the ICBI n's "
            f"homogeneity alone gives a median I^2 of 0.269 and rejects "
            f"Cochran's Q 32.5% of the time. Emit the null in the same row."
        )
    if not has_verdict:
        return
    missing = [c for c in _NULL_COLUMNS if c not in frame.columns]
    if missing:
        raise MetaCalibrationError(
            f"frame carries a verdict and only part of its calibration; "
            f"missing {missing}."
        )
    for column in _NULL_COLUMNS:
        if frame[column].isna().any():
            bad = frame.loc[frame[column].isna()]
            raise MetaCalibrationError(
                f"{int(len(bad))} row(s) state a verdict with {column} missing. "
                f"Invariant 1: a null that was not computed is not a null of 0."
            )
