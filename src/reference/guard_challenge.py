"""Constructed inputs run against checks, with the expected outcome sealed first.

Week 1, deliverable 2 of ``docs/DECISION_2026-09-11_scope_and_pivot.md`` §4.

The question this answers is narrow and it is not "how many bugs are there". It
is: **does the audit separate a check that missed a scientific failure from a
check whose 'failure' is an algebraic identity behaving as identities do?** If
it cannot, the audit flags everything, and a method that flags everything ranks
nothing.

So every case declares, before it runs:

``expected``
    ``fires`` or ``passes``. Sealed in the case definition; comparing it to the
    observed outcome is the whole measurement.
``consequence``
    ``decision_relevant`` — a miss would change a scientific conclusion.
    ``input_validation`` — a miss produces a crash or a malformed table, not a
    wrong conclusion. Ordinary software hygiene.
    ``algebraic_identity`` — the check reports something mathematically
    necessary. **Excluded from the stop rule by construction**, because counting
    identities as findings is how an audit inflates its own hit rate.

Only ``decision_relevant`` cases where the guard did **not** behave as it should
count toward the week-1 stop rule.

BLINDING, STATED PLAINLY
------------------------
Every case here has ``saw_ledger=True``. They were written by someone who had
read ``docs/HANDOFF.md`` §3 and the existing guard suite, so they are
**development material and not the held-out evaluation** — the same objection
the decision record makes against counting the 53 existing guard tests as
evidence. ``docs/week1_challenge_protocol.md`` is how an unblinded challenger
produces cases that do count. This file is the harness and the worked examples;
it is not the result the stop rule is read from.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from src.harness.depth_confound import MATURITY_DEPTH_RHO_TOLERANCE, max_attainable_rho

CONSEQUENCES: frozenset[str] = frozenset(
    {"decision_relevant", "input_validation", "algebraic_identity"}
)
OUTCOMES: frozenset[str] = frozenset({"fires", "passes"})


class ChallengeError(ValueError):
    """A case is not well formed, or its expectation was not sealed."""


@dataclass(frozen=True)
class ChallengeCase:
    """One constructed input, one check, one sealed expectation."""

    case_id: str
    target_check: str
    description: str
    expected: str
    consequence: str
    run: Callable[[], bool] = field(repr=False)
    constructed_by: str = "rater_a"
    saw_ledger: bool = True
    ledger_entry: str | None = None
    #: Set when a case's input was changed AFTER its sealed expectation was
    #: falsified. Changing an input to make an expectation come true is the move
    #: that must never be invisible, so it is a column rather than a comment.
    revised_after_falsification: bool = False
    revision_note: str | None = None

    def __post_init__(self) -> None:
        if self.expected not in OUTCOMES:
            raise ChallengeError(f"{self.case_id}: expected must be one of {sorted(OUTCOMES)}")
        if self.consequence not in CONSEQUENCES:
            raise ChallengeError(
                f"{self.case_id}: consequence must be one of {sorted(CONSEQUENCES)}"
            )
        if self.revised_after_falsification and not self.revision_note:
            raise ChallengeError(
                f"{self.case_id}: the input was revised after its expectation was "
                "falsified and no revision_note says so. An unrecorded revision is "
                "how a challenge set quietly becomes a set of cases that pass."
            )


# ---------------------------------------------------------------------------
# The checks under challenge. Each returns True for "this check fired".
# ---------------------------------------------------------------------------


def check_interval_is_calibrated(false_positive_rate: float, *, alpha: float = 0.05) -> bool:
    """Fires when an interval's false-positive rate misses its nominal level.

    The tolerance is the Monte Carlo standard error at the replicate count used
    below, doubled. Fixing it by measurement rather than by taste is the point:
    a tolerance chosen after seeing the rate is the defect this audit is about.
    """
    return abs(false_positive_rate - alpha) > 2 * np.sqrt(alpha * (1 - alpha) / 4000)


def check_recovery_curve_is_informative(residuals: np.ndarray) -> bool:
    """Fires when a recovery curve's residuals are identically zero.

    A curve whose estimator reproduces its comparator exactly measures the
    generator, not the estimator. Whether that is a defect or an identity is
    what ``consequence`` records — this function only reports the arithmetic.
    """
    return bool(np.max(np.abs(residuals)) < 1e-12)


def check_statistic_can_reach_tolerance(ceiling: float, tolerance: float) -> bool:
    """Fires when a statistic's attainable maximum sits below its threshold."""
    return bool(ceiling < tolerance)


def check_detection_can_fall(baseline: float, *, floor: float = 0.05) -> bool:
    """Fires when a proportion is close enough to 1.0 to have nowhere to fall."""
    return bool(1.0 - baseline < floor)


def check_counts_three_outcomes(values: np.ndarray, threshold: float) -> bool:
    """Fires when abstentions are silently folded into failures.

    ``nan > x`` is False, so a naive ``(values > threshold).sum()`` scores an
    abstention exactly as it scores a failure. The check compares the naive
    count with one that separates the three outcomes.
    """
    naive = int((values > threshold).sum())
    separated = int((values[~np.isnan(values)] > threshold).sum())
    n_abstained = int(np.isnan(values).sum())
    return naive == separated and n_abstained > 0


def check_branch_is_directional(excludes_zero: bool, effect: float, predicted_sign: int) -> bool:
    """Fires when a directional prediction is decided by a direction-free test.

    The shape of ledger entry L22: ``excludes_zero`` is true on both sides of
    zero, so a branch reading only that clause confirms a prediction its own
    data contradicts.
    """
    direction_free_says = excludes_zero
    directional_says = excludes_zero and np.sign(effect) == predicted_sign
    return direction_free_says != directional_says


# ---------------------------------------------------------------------------
# Simulations backing the cases. Fixed seeds; no external data.
# ---------------------------------------------------------------------------

N_REPLICATES = 4000
N_PATIENTS = 10


def _null_false_positive_rate(method: str, *, seed: int = 20260101) -> float:
    """How often a 95% interval for a mean excludes zero when the mean is zero."""
    rng = np.random.default_rng(seed)
    excluded = 0
    for _ in range(N_REPLICATES):
        sample = rng.normal(0.0, 1.0, N_PATIENTS)
        if method == "student_t":
            from scipy import stats

            half = stats.t.ppf(0.975, N_PATIENTS - 1) * sample.std(ddof=1) / np.sqrt(N_PATIENTS)
            low, high = sample.mean() - half, sample.mean() + half
        elif method == "percentile":
            idx = rng.integers(0, N_PATIENTS, size=(400, N_PATIENTS))
            boots = sample[idx].mean(axis=1)
            low, high = np.percentile(boots, [2.5, 97.5])
        else:  # pragma: no cover - guarded by the case definitions
            raise ChallengeError(f"unknown interval method {method!r}")
        excluded += int(low > 0 or high < 0)
    return excluded / N_REPLICATES


def _oracle_residuals(*, against: str, seed: int = 20260101) -> np.ndarray:
    """Residuals of a sample mean against realised or parametric truth."""
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(200):
        cells = rng.normal(1.5, 1.0, 300)
        estimate = cells.mean()
        truth = cells.mean() if against == "realised" else 1.5
        out.append(estimate - truth)
    return np.asarray(out)


# ---------------------------------------------------------------------------
# The cases
# ---------------------------------------------------------------------------

CASES: tuple[ChallengeCase, ...] = (
    ChallengeCase(
        case_id="C01",
        target_check="check_interval_is_calibrated",
        description=(
            "CLEAN CONTROL. A Student-t interval for a sample mean at n=10, under a "
            "true null. A correctly calibrated check must NOT fire on a correctly "
            "calibrated interval; if it does, the audit flags everything."
        ),
        expected="passes",
        consequence="decision_relevant",
        ledger_entry=None,
        run=lambda: check_interval_is_calibrated(_null_false_positive_rate("student_t")),
    ),
    ChallengeCase(
        case_id="C02",
        target_check="check_interval_is_calibrated",
        description=(
            "The percentile bootstrap of a mean at the same n, same null, same "
            "replicate count. Differs from C01 only in the interval."
        ),
        expected="fires",
        consequence="decision_relevant",
        ledger_entry="L11",
        run=lambda: check_interval_is_calibrated(_null_false_positive_rate("percentile")),
    ),
    ChallengeCase(
        case_id="C03",
        target_check="check_recovery_curve_is_informative",
        description=(
            "ALGEBRAIC IDENTITY. A sample mean compared with the realised sample "
            "mean. The residual is exactly zero because the two are the same "
            "arithmetic on the same cells. The check fires and it is RIGHT to "
            "fire, and this is still not a missed scientific failure — which is "
            "the distinction ledger entry L01 is disputed over."
        ),
        expected="fires",
        consequence="algebraic_identity",
        ledger_entry="L01",
        run=lambda: check_recovery_curve_is_informative(_oracle_residuals(against="realised")),
    ),
    ChallengeCase(
        case_id="C04",
        target_check="check_recovery_curve_is_informative",
        description=(
            "CLEAN CONTROL, and C03's pair. The same estimator against the "
            "PARAMETRIC truth, which the sample does not reproduce exactly. "
            "Residuals are non-zero and the curve is informative. src/harness/"
            "calibration.py already makes this distinction in a comment; the "
            "WMHS prose is where it goes missing."
        ),
        expected="passes",
        consequence="decision_relevant",
        ledger_entry="L01",
        run=lambda: check_recovery_curve_is_informative(_oracle_residuals(against="parametric")),
    ),
    ChallengeCase(
        case_id="C05",
        target_check="check_statistic_can_reach_tolerance",
        description=(
            "A rare label at prevalence 0.00144 -- the Lee/SMC best4 tumour arm's "
            "8 mature cells in 5,564 -- where the attainable maximum of the "
            "depth-maturity correlation is 0.066 against a 0.20 tolerance. No "
            "dataset can fail this check at this prevalence."
        ),
        expected="fires",
        consequence="decision_relevant",
        ledger_entry="L03",
        revised_after_falsification=True,
        revision_note=(
            "First written at prevalence 0.02 and the sealed expectation was "
            "FALSIFIED on the first run: sqrt(3*0.02*0.98) = 0.2425, which clears "
            "the tolerance. The error was the constructor's, not the guard's -- a "
            "prevalence chosen without computing the ceiling, when "
            "max_attainable_rho's own docstring states the crossing at p = 1.37%. "
            "The input was changed to a prevalence the repository actually "
            "reports. Recorded rather than silently corrected, because a sealed "
            "expectation that is falsified and then edited is indistinguishable "
            "from one that was never sealed."
        ),
        run=lambda: check_statistic_can_reach_tolerance(
            max_attainable_rho(0.00144), MATURITY_DEPTH_RHO_TOLERANCE
        ),
    ),
    ChallengeCase(
        case_id="C06",
        target_check="check_statistic_can_reach_tolerance",
        description=(
            "CLEAN CONTROL. The same check at prevalence 0.5, where the ceiling "
            "clears the tolerance and the comparison is real."
        ),
        expected="passes",
        consequence="decision_relevant",
        ledger_entry="L03",
        run=lambda: check_statistic_can_reach_tolerance(
            max_attainable_rho(0.5), MATURITY_DEPTH_RHO_TOLERANCE
        ),
    ),
    ChallengeCase(
        case_id="C07",
        target_check="check_detection_can_fall",
        description=(
            "A premise control whose baseline detection is 0.995. A real fall in "
            "per-cell output cannot register as a fall in detection."
        ),
        expected="fires",
        consequence="decision_relevant",
        ledger_entry="L06",
        run=lambda: check_detection_can_fall(0.995),
    ),
    ChallengeCase(
        case_id="C08",
        target_check="check_detection_can_fall",
        description="CLEAN CONTROL. The same control at baseline 0.51, where a fall is visible.",
        expected="passes",
        consequence="decision_relevant",
        ledger_entry="L06",
        run=lambda: check_detection_can_fall(0.51),
    ),
    ChallengeCase(
        case_id="C09",
        target_check="check_counts_three_outcomes",
        description=(
            "A coverage vector carrying abstentions as NaN. `nan > x` is False, so "
            "a naive count scores an abstention exactly as it scores a failure."
        ),
        expected="fires",
        consequence="decision_relevant",
        ledger_entry="L04",
        run=lambda: check_counts_three_outcomes(
            np.array([0.9, 0.8, np.nan, np.nan, 0.2]), 0.5
        ),
    ),
    ChallengeCase(
        case_id="C10",
        target_check="check_counts_three_outcomes",
        description=(
            "INPUT VALIDATION. The same vector with no abstentions at all. Nothing "
            "is miscounted, and a check that fired here would be complaining about "
            "well-formed input."
        ),
        expected="passes",
        consequence="input_validation",
        ledger_entry="L04",
        run=lambda: check_counts_three_outcomes(np.array([0.9, 0.8, 0.4, 0.2]), 0.5),
    ),
    ChallengeCase(
        case_id="C11",
        target_check="check_branch_is_directional",
        description=(
            "An interval excluding zero on the side OPPOSITE the prediction. A "
            "direction-free branch reports the prediction replicating; a "
            "directional one reports contradiction. L22's exact shape, in the "
            "abstract."
        ),
        expected="fires",
        consequence="decision_relevant",
        ledger_entry="L22",
        run=lambda: check_branch_is_directional(True, effect=-0.42, predicted_sign=+1),
    ),
    ChallengeCase(
        case_id="C12",
        target_check="check_branch_is_directional",
        description=(
            "CLEAN CONTROL. The same branch on an interval excluding zero in the "
            "PREDICTED direction, where direction-free and directional agree. This "
            "is the case that actually ran in the Chen subtype contrast, which is "
            "why the defect changed no result."
        ),
        expected="passes",
        consequence="decision_relevant",
        ledger_entry="L22",
        run=lambda: check_branch_is_directional(True, effect=+0.523, predicted_sign=+1),
    ),
)


def run_challenges(cases: tuple[ChallengeCase, ...] = CASES) -> pd.DataFrame:
    """Run every case and compare the observed outcome with the sealed one."""
    ids = [c.case_id for c in cases]
    if len(set(ids)) != len(ids):
        raise ChallengeError(f"duplicate case ids: {sorted(ids)}")

    rows = []
    for case in cases:
        fired = bool(case.run())
        observed = "fires" if fired else "passes"
        rows.append(
            {
                "case_id": case.case_id,
                "target_check": case.target_check,
                "ledger_entry": case.ledger_entry,
                "description": case.description,
                "constructed_by": case.constructed_by,
                "saw_ledger": case.saw_ledger,
                "revised_after_falsification": case.revised_after_falsification,
                "revision_note": case.revision_note,
                "expected": case.expected,
                "observed": observed,
                "as_expected": observed == case.expected,
                "consequence": case.consequence,
                "is_clean_control": case.expected == "passes",
            }
        )
    frame = pd.DataFrame(rows)
    frame["counts_toward_stop_rule"] = (
        (~frame["as_expected"]) & (frame["consequence"] == "decision_relevant")
    )
    return frame


def stop_rule_summary(frame: pd.DataFrame) -> dict:
    """The week-1 arithmetic, computed rather than asserted."""
    return {
        "n_cases": int(len(frame)),
        "n_clean_controls": int(frame["is_clean_control"].sum()),
        "n_as_expected": int(frame["as_expected"].sum()),
        "n_unexpected": int((~frame["as_expected"]).sum()),
        "n_algebraic_identity": int((frame["consequence"] == "algebraic_identity").sum()),
        "n_input_validation": int((frame["consequence"] == "input_validation").sum()),
        "n_counting_toward_stop_rule": int(frame["counts_toward_stop_rule"].sum()),
        "n_revised_after_falsification": int(frame["revised_after_falsification"].sum()),
        "any_case_blind_to_the_ledger": bool((~frame["saw_ledger"]).any()),
    }
