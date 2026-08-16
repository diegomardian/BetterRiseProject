"""Deriving the positivity cutpoints from the sweep. W2, week 5.

The provisional cutpoints (n>=50 / 20<=n<50 / n<20) are round numbers. These
are derived from where the estimator stops being able to answer, against three
values pre-registered in ``docs/harness_design_spec.md`` §4 **before** any sweep
was run. That ordering is the whole reason the result is a calibration and not a
choice, so :func:`calibrate_cutpoints` takes the criteria as an argument and
records them in its report rather than reading them off the data.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
import pandas as pd

from src.harness.positivity import PROVISIONAL, Cutpoints


@dataclass(frozen=True)
class CalibrationCriteria:
    """Pre-registered. Fixed before the sweep, not tuned to it."""

    #: The effect we require the estimator to be able to see: a halving of
    #: per-cell output. Mid-band of the published attenuation range.
    detectable_shift: float = 0.5
    #: A nominal-95% interval that covers the truth less often than this is not
    #: an interval.
    coverage_target: float = 0.90
    #: Probability the interval excludes zero at the detectable effect.
    discrimination_target: float = 0.80
    #: Which arm the cutpoints are calibrated on. "oracle" — the cutpoints
    #: govern single-cell decomposition, and the bulk arm is measured for §2.2
    #: rather than used to decide what is estimable.
    arm: str = "oracle"


#: The pre-registered criteria, as a named singleton. Referring to this rather
#: than constructing one per call makes it a single auditable object: if the
#: numbers ever move, they move here, in one place, in a reviewable diff.
PREREGISTERED = CalibrationCriteria()


@dataclass(frozen=True)
class CalibrationReport:
    """The cutpoints plus the evidence, so a reader can check the derivation."""

    cutpoints: Cutpoints
    criteria: CalibrationCriteria
    table: pd.DataFrame
    provisional: Cutpoints = PROVISIONAL

    def comparison(self) -> pd.DataFrame:
        """Provisional vs calibrated, side by side, for the gate memo."""
        return pd.DataFrame(
            [
                {"cutpoint": "ok", "provisional": self.provisional.ok,
                 "calibrated": self.cutpoints.ok},
                {"cutpoint": "wide", "provisional": self.provisional.wide,
                 "calibrated": self.cutpoints.wide},
            ]
        )  # fmt: skip


def _bin_edges(n_cells_mature: pd.Series, n_bins: int) -> np.ndarray:
    """Log-spaced bins over the observed mature-cell counts."""
    positive = n_cells_mature[n_cells_mature > 0]
    if positive.empty:
        return np.array([0, 1])
    lo, hi = max(int(positive.min()), 1), int(positive.max())
    if lo >= hi:
        return np.array([0, lo, hi + 1])
    return np.unique(
        np.concatenate([[0], np.geomspace(lo, hi, num=n_bins).round().astype(int)])
    )


def coverage_and_discrimination(
    sweep: pd.DataFrame,
    criteria: CalibrationCriteria = PREREGISTERED,
    *,
    n_bins: int = 12,
) -> pd.DataFrame:
    """Per mature-count bin: how often the CI covers truth and excludes zero.

    Requires ``ci_low``/``ci_high`` on the sweep rows. ``run_sweep`` leaves them
    null — CIs come from ``kitagawa.bootstrap_over_patients`` or
    ``hierarchical.hierarchical_intrinsic_ci``, both of which are W4's, and
    attaching them is the step that turns a sweep into a calibration.
    """
    rows = sweep[
        (sweep["arm"] == criteria.arm) & (sweep["shift"] == criteria.detectable_shift)
    ].copy()
    if rows.empty:
        raise ValueError(
            f"sweep has no arm={criteria.arm!r} rows at the pre-registered "
            f"detectable shift {criteria.detectable_shift}. Calibrating on a "
            f"different effect than the one registered would make the cutpoints "
            f"chosen rather than derived."
        )
    if rows[["ci_low", "ci_high"]].isna().all().all():
        raise ValueError(
            "sweep carries no confidence intervals, so coverage and "
            "discrimination are undefined. Attach CIs first — see the docstring."
        )

    edges = _bin_edges(rows["n_cells_mature"], n_bins)
    rows["bin"] = pd.cut(rows["n_cells_mature"], bins=edges, include_lowest=True)

    truth = rows["intrinsic_true_realised"]
    covered = (rows["ci_low"] <= truth) & (truth <= rows["ci_high"])
    excludes_zero = (rows["ci_low"] > 0) | (rows["ci_high"] < 0)
    rows = rows.assign(
        _covered=covered, _excludes_zero=excludes_zero,
        _width=rows["ci_high"] - rows["ci_low"],
    )

    out = (
        rows.groupby("bin", observed=True)
        .agg(
            n_replicates=("replicate", "count"),
            n_cells_mature=("n_cells_mature", "median"),
            coverage=("_covered", "mean"),
            discrimination=("_excludes_zero", "mean"),
            median_ci_width=("_width", "median"),
        )
        .reset_index(drop=True)
        .sort_values("n_cells_mature")
    )
    out["shift"] = criteria.detectable_shift
    out["verdict"] = [
        _verdict(c, d, criteria)
        for c, d in zip(out["coverage"], out["discrimination"], strict=True)
    ]
    return out.reset_index(drop=True)


def _verdict(coverage: float, discrimination: float, criteria: CalibrationCriteria) -> str:
    if coverage >= criteria.coverage_target and discrimination >= criteria.discrimination_target:
        return "ok"
    if coverage >= criteria.coverage_target:
        return "wide_interval"
    return "not_estimable"


def calibrate_cutpoints(
    sweep: pd.DataFrame,
    criteria: CalibrationCriteria = PREREGISTERED,
    *,
    n_bins: int = 12,
    source: str | None = None,
) -> CalibrationReport:
    """Smallest mature-cell counts meeting each criterion. Week-5 deliverable.

    ``ok``   — smallest count with coverage >= target AND discrimination >= target
    ``wide`` — smallest count with coverage >= target but discrimination below it

    The returned ``Cutpoints.source`` names the run, so the number that ends up
    in ``positivity.CUTPOINTS`` traces back to the sweep that produced it rather
    than appearing as a constant nobody can account for.

    Raises if no bin qualifies as ``ok``. That is not a failure of this
    function — it is gate criterion **G4** firing, and it means
    non-identifiability is the headline result rather than a caveat.
    """
    table = coverage_and_discrimination(sweep, criteria, n_bins=n_bins)

    ok_bins = table[table["verdict"] == "ok"]
    if ok_bins.empty:
        raise ValueError(
            "no mature-cell count in this sweep meets both the coverage and "
            "discrimination targets. The estimator cannot resolve the "
            "pre-registered effect at any n tested — that is G4, and the "
            "non-identifiability finding becomes the headline result rather "
            "than a caveat (execution_plan.md §5). Do not widen the criteria "
            "to make this go away."
        )
    ok = int(ok_bins["n_cells_mature"].min())

    wide_bins = table[table["verdict"].isin(["ok", "wide_interval"])]
    wide = int(wide_bins["n_cells_mature"].min()) if not wide_bins.empty else ok

    return CalibrationReport(
        cutpoints=Cutpoints(
            ok=ok,
            wide=min(wide, ok),
            source=source or f"calibrated on {len(sweep)} sweep rows, {asdict(criteria)}",
        ),
        criteria=criteria,
        table=table,
    )
