"""What the gate costs at the cohort sizes that actually exist. W2.

Everything downstream of the plan was designed for roughly 60 patients. The
cohorts are smaller than that, and nobody had re-costed the gate against the real
numbers — W2 promised W4 this twice.

    GSE178341   36 matched of 62   (~26 have no normal arm)
    SMC         10 paired, not 23
    the plan    ~60

Two different quantities are affected and they are not interchangeable:

1. **G4's verdict.** G4 is a proportion tested against a pre-committed 0.50 line,
   so its precision is binomial in the number of *patients*. This is where the
   shortfall bites hardest, and it is the half nobody had looked at.
2. **The cohort band on the decomposition terms.** ``bootstrap_over_patients``
   resamples patients (CLAUDE.md invariant 5), so its width scales as roughly
   1/sqrt(n_patients). That one is arithmetic, but it is measured here rather
   than asserted.

WHAT THIS MODULE DOES NOT DO
----------------------------
It does not change G4's rule. The 0.50 line is pre-committed (execution_plan.md
§5, decision #19) and moving it after seeing what n does to it is the move this
project refuses everywhere else. What it does is report the interval alongside
the verdict, so that "PASS" and "PASS, and this cohort could not have said
otherwise" stop being the same sentence.
"""

from __future__ import annotations

from math import comb
from typing import Final

import numpy as np
import pandas as pd

from src.harness.positivity import (
    NON_IDENTIFIABILITY_HEADLINE_FRACTION,
    gate_g4_verdict,
    wilson_interval,
)

#: The real cohort sizes, with what each one is. Kept here rather than passed in
#: so that a re-costing cannot quietly be run against the plan's numbers again.
COHORT_SIZES: Final = {
    "smc_paired": 10,
    "gse178341_matched": 36,
    "plan_assumed": 60,
}


# ---------------------------------------------------------------------------
# G4 — a proportion against a pre-committed line
# ---------------------------------------------------------------------------


def g4_pass_probability(n_patients: int, true_fraction_below: float) -> float:
    """P(G4 returns PASS) when the true below-threshold fraction is known.

    Exact binomial rather than a normal approximation: at n = 10 the
    approximation is wrong in the third decimal and this is a decision rule, not
    a plot. G4 passes when the *observed* fraction is strictly below 0.50, so the
    largest passing count is ``ceil(n/2) - 1``.
    """
    if n_patients < 1:
        raise ValueError(f"n_patients={n_patients} must be at least 1")
    if not 0 <= true_fraction_below <= 1:
        raise ValueError(f"true_fraction_below={true_fraction_below} outside [0, 1]")
    largest_passing = -(-n_patients // 2) - 1  # ceil(n/2) - 1
    p = true_fraction_below
    return float(
        sum(
            comb(n_patients, k) * p**k * (1 - p) ** (n_patients - k)
            for k in range(0, largest_passing + 1)
        )
    )


def largest_clean_pass(n_patients: int, *, alpha: float = 0.05) -> int:
    """Most below-threshold patients that still give a PASS the cohort can defend.

    "Defend" meaning the interval on the observed fraction excludes the 0.50
    line, so the verdict is not a point estimate presented as a fact. Returns -1
    if no count qualifies, which is itself the answer for a small enough cohort.
    """
    best = -1
    for k in range(n_patients + 1):
        low, high = wilson_interval(k, n_patients, alpha=alpha)
        passes = k / n_patients < NON_IDENTIFIABILITY_HEADLINE_FRACTION
        resolvable = not (low < NON_IDENTIFIABILITY_HEADLINE_FRACTION < high)
        if passes and resolvable:
            best = k
    return best


def effective_decision_line(n_patients: int, *, alpha: float = 0.05) -> float:
    """The 0.50 line as the cohort can actually enforce it.

    G4's rule says 50%. Read with an interval, a cohort of 36 needs 33.3% and a
    cohort of 10 needs 10%. The pre-committed line does not move; what moves is
    how far below it the data has to land before the call means anything, and
    that gap is the cost of the shortfall.
    """
    best = largest_clean_pass(n_patients, alpha=alpha)
    return float("nan") if best < 0 else best / n_patients


def g4_operating_characteristic(
    *,
    sizes: dict[str, int] | None = None,
    true_fractions: tuple[float, ...] = (0.20, 0.35, 0.40, 0.45, 0.50, 0.60, 0.75),
) -> pd.DataFrame:
    """P(PASS) by cohort size and true fraction. One row per (cohort, fraction).

    Read down a column to see how a cohort behaves as the truth moves; read
    across to see what n buys. Near the line it buys less than people expect,
    which is the honest headline: at a true 0.45 the gate says PASS 67% of the
    time at n=36 and 74% at n=60. Neither is a decision procedure you would
    choose; the pre-committed consequence is what makes it one.
    """
    sizes = dict(sizes or COHORT_SIZES)
    rows = []
    for label, n in sizes.items():
        for fraction in true_fractions:
            rows.append(
                {
                    "cohort": label,
                    "n_patients": n,
                    "true_fraction_below": fraction,
                    "p_gate_says_pass": g4_pass_probability(n, fraction),
                    "largest_clean_pass": largest_clean_pass(n),
                    "effective_decision_line": effective_decision_line(n),
                }
            )
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# The cohort band on the terms themselves
# ---------------------------------------------------------------------------


def synthetic_cohort(
    n_patients: int,
    *,
    rng: np.random.Generator,
    study_id: str = "RECOST",
    gene: str = "GUCA2A",
) -> pd.DataFrame:
    """A patient-level summary frame of the shape ``decompose_cohort`` takes.

    Spread chosen so that the mature fraction and the per-cell mean both vary
    between patients, because a cohort band with no between-patient variance is
    zero-width whatever n is, and would make this measurement say nothing.
    """
    return pd.DataFrame(
        {
            "patient_id": [f"P{i:03d}" for i in range(n_patients)],
            "study_id": study_id,
            "gene": gene,
            "granularity_rung": "lineage",
            "labeling_axis": "stem_pole",
            "frac_mature_normal": rng.uniform(0.30, 0.50, n_patients),
            "frac_mature_tumour": rng.uniform(0.05, 0.30, n_patients),
            "mean_normal": rng.normal(10.0, 1.5, n_patients),
            "mean_tumour": rng.normal(7.0, 1.5, n_patients),
            "n_cells_mature": rng.integers(60, 900, n_patients),
        }
    )


def cohort_ci_width_by_n(
    *,
    seed: int,
    sizes: dict[str, int] | None = None,
    n_boot: int = 200,
    n_replicates: int = 5,
    weighting: str = "normal",
) -> pd.DataFrame:
    """Measured width of the patient-level bootstrap band, by cohort size.

    Runs W4's ``bootstrap_over_patients`` — the real one, not a stand-in — over
    synthetic cohorts of each size. Replicated, because a single cohort draw
    gives one width and the question is about the size, not about that draw.

    The expected scaling is 1/sqrt(n): 60 -> 36 widens by 1.29x, 60 -> 10 by
    2.45x. Measuring it rather than quoting it is the point, since the estimator
    nulls the intrinsic term below the positivity cutpoint and a formula does not
    know that.
    """
    from src.estimator.kitagawa import bootstrap_over_patients

    sizes = dict(sizes or COHORT_SIZES)
    rng = np.random.default_rng(seed)
    rows = []
    for label, n in sizes.items():
        for replicate in range(n_replicates):
            summary = synthetic_cohort(n, rng=rng)
            band = bootstrap_over_patients(summary, n_boot=n_boot, seed=seed + replicate)
            band = band[band["weighting"] == weighting]
            for _, row in band.iterrows():
                if row["ci_low"] is None or row["ci_high"] is None:
                    continue
                rows.append(
                    {
                        "cohort": label,
                        "n_patients": n,
                        "replicate": replicate,
                        "term": row["term"],
                        "ci_width": float(row["ci_high"]) - float(row["ci_low"]),
                        "n_boot": n_boot,
                        "seed": seed + replicate,
                    }
                )
    return pd.DataFrame(rows)


def widening_vs_plan(widths: pd.DataFrame, *, reference: str = "plan_assumed") -> pd.DataFrame:
    """How much wider each cohort's band is than the plan's, per term.

    The number to quote at the gate. ``ratio`` above 1 means the real cohort buys
    a wider interval than everything downstream was designed against.
    """
    median = (
        widths.groupby(["cohort", "n_patients", "term"], observed=True)["ci_width"]
        .median()
        .reset_index()
    )
    base = median[median["cohort"] == reference].set_index("term")["ci_width"]
    median["ratio_vs_plan"] = [
        row.ci_width / base[row.term] if row.term in base.index else float("nan")
        for row in median.itertuples()
    ]
    return median.sort_values(["term", "n_patients"]).reset_index(drop=True)


# ---------------------------------------------------------------------------
# G4 against a real counts frame
# ---------------------------------------------------------------------------


#: A patient with both arms. Decision #19: G4's population is matched patients
#: only, and ``n_unmatched_patients`` is required rather than defaulted, because
#: a default is how the wrong population gets used without anyone choosing it.
def matched_and_unmatched(counts: pd.DataFrame) -> tuple[set, set]:
    """Split patient ids into those with both arms and those with one.

    "Matched" is defined here as *observed in both tissues*, not as a column
    somebody set upstream. The compositional term is a within-patient contrast,
    so a patient with one arm contributes to neither side of it.
    """
    for column in ("patient_id", "tissue"):
        if column not in counts.columns:
            raise ValueError(f"counts is missing {column!r}")
    arms = counts.groupby("patient_id", observed=True)["tissue"].nunique()
    return set(arms[arms >= 2].index), set(arms[arms < 2].index)


#: The two tumour-arm definitions of prereg amendment 1 / decision #15. W1's
#: ``mature_cell_counts`` emits BOTH, and its docstring says "Report both
#: definitions; do not pick one."
TUMOUR_ARM_COLUMN: Final = "tumour_arm"


def g4_over_rungs(counts: pd.DataFrame, *, arm: str = "tumour") -> pd.DataFrame:
    """G4's verdict for every (axis, rung, tumour_arm) in a counts frame.

    One row per group, because **G4 does not have a single answer**. The mature
    population is defined by the rung, so the fraction below the positivity
    threshold is too — and it is defined again by which tumour-arm definition is
    in force, because the filtered arm drops tumour epithelium not called
    malignant.

    IT GROUPS BY ``tumour_arm``, AND THE FIRST VERSION DID NOT
    ----------------------------------------------------------
    That was a real defect and it reached a published table. W1's frame carries
    both definitions stacked — 30 patients x 2 tissues x 2 axes x 4 rungs x
    **2 arms**. Without grouping, the per-patient dedup below collapsed two *arm
    definitions* rather than two duplicate rows, and ``.first()`` silently took
    whichever came first in file order: ``filtered``, in **224 of 224**
    combinations.

    The verdicts differ enormously. At ``stem_pole``/``lineage``: 13 of 28
    patients below the cutpoint under ``filtered``, **1 of 28** under
    ``unfiltered``. Five patients have *zero* epithelium in the filtered arm
    (C119, C137, C152, C165, C170 — no aneuploid population separated), so they
    enter ``fraction_below`` as guaranteed-below on a **malignancy-calling**
    fact. That is exactly the failure :func:`gate_g4_verdict`'s own docstring
    says it exists to prevent — a cohort-design fact reported as a positivity
    finding — reproduced one level up, in the function that consumes it.

    Prereg amendment 1 pre-commits the answer when the two disagree: report both
    and treat disagreement as **not identifiable** rather than choosing.
    :func:`g4_arms_agree` does that reduction; this function does not do it for
    you, because collapsing the arms is the mistake.
    """
    required = {"patient_id", "tissue", "labeling_axis", "granularity_rung", "n_cells_mature"}
    missing = sorted(required - set(counts.columns))
    if missing:
        raise ValueError(f"counts is missing column(s) {missing}")

    matched, unmatched = matched_and_unmatched(counts)
    if not matched:
        raise ValueError("no patient has both arms; G4's population is empty")

    group_cols = ["labeling_axis", "granularity_rung"]
    if TUMOUR_ARM_COLUMN in counts.columns:
        group_cols.append(TUMOUR_ARM_COLUMN)

    rows = []
    for key, group in counts.groupby(group_cols, observed=True):
        axis, rung = key[0], key[1]
        tumour_arm = key[2] if len(key) > 2 else "unspecified"
        sub = group[(group["tissue"] == arm) & (group["patient_id"].isin(matched))]
        if sub.empty:
            continue
        # One row per patient. A patient counted twice is invariant 5's mistake
        # one level up — it was the fourth instance of the bug family this repo
        # keeps finding (issue #36).
        per_patient = sub.groupby("patient_id", observed=True)["n_cells_mature"].first()
        verdict = gate_g4_verdict(
            [int(n) for n in per_patient], n_unmatched_patients=len(unmatched)
        )
        rows.append(
            {
                "labeling_axis": axis,
                "granularity_rung": rung,
                "tumour_arm": tumour_arm,
                "arm": arm,
                "n_patients": verdict["n_patients"],
                "n_unmatched_excluded": verdict["n_unmatched_excluded"],
                "n_below_threshold": verdict["n_below_threshold"],
                "fraction_below": verdict["fraction_below"],
                "fraction_below_ci_low": verdict["fraction_below_ci_low"],
                "fraction_below_ci_high": verdict["fraction_below_ci_high"],
                "threshold": verdict["threshold"],
                "passes": verdict["passes"],
                "resolvable": verdict["resolvable"],
            }
        )
    return (
        pd.DataFrame(rows)
        .sort_values(["labeling_axis", "granularity_rung", "tumour_arm"])
        .reset_index(drop=True)
    )


def g4_arms_agree(verdicts: pd.DataFrame) -> pd.DataFrame:
    """Reduce per-arm verdicts to one, per prereg amendment 1.

    The amendment pre-commits the answer: report both tumour-arm definitions and
    treat disagreement as **not identifiable** rather than choosing between them.
    A rung where ``filtered`` and ``unfiltered`` return different verdicts is not
    a FAIL and not a PASS — it is a rung this cohort cannot answer, and saying so
    is the pre-registered response rather than a hedge.

    A rung where the arms agree but some arm cannot resolve its own interval is
    also ``not_identifiable``: agreeing on an undecidable number is not agreement.
    """
    if TUMOUR_ARM_COLUMN not in verdicts.columns:
        raise ValueError(
            "verdicts has no tumour_arm column, so agreement between the two "
            "definitions cannot be checked. Pass g4_over_rungs output from a "
            "frame that carries both arms."
        )
    present = sorted(set(verdicts[TUMOUR_ARM_COLUMN]))
    if len(present) < 2:
        raise ValueError(
            f"only one tumour-arm definition is present ({present}), so there is "
            f"nothing to check agreement between. Prereg amendment 1 requires "
            f"BOTH — reducing a single arm to a verdict is the choice the "
            f"amendment exists to prevent."
        )
    rows = []
    for (axis, rung), group in verdicts.groupby(
        ["labeling_axis", "granularity_rung"], observed=True
    ):
        passes = set(group["passes"])
        resolvable = bool(group["resolvable"].all())
        if len(passes) > 1 or not resolvable:
            verdict = "not_identifiable"
        else:
            verdict = "PASS" if passes == {True} else "FAIL"
        rows.append(
            {
                "labeling_axis": axis,
                "granularity_rung": rung,
                "verdict": verdict,
                "arms_agree": len(passes) == 1,
                "all_arms_resolvable": resolvable,
                "n_arms": int(len(group)),
            }
        )
    return (
        pd.DataFrame(rows)
        .sort_values(["labeling_axis", "granularity_rung"])
        .reset_index(drop=True)
    )
