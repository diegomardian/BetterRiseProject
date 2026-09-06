"""The decomposition job, against inputs whose answer is fixed independently.

Two kinds of test here. The first kind checks arithmetic the frozen spec already
determines — `epithelial` MUST return a compositional term of exactly zero,
because every epithelial cell is mature at that rung and Δ(mature fraction) is
therefore identically zero. That is not a property of the data; it is what the
rung means, and it is the lower bound the granularity curve is reported against.

The second kind is the repository's usual discipline: the input that forces each
guard to fail. A decomposition that silently defaulted a missing mature fraction
would produce a full-looking table of numbers computed from nothing.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.estimator.kitagawa import decompose_cohort
from src.reference.jobs.adenoma_decomposition import (
    DENOMINATORS,
    PRIMARY_DENOMINATOR,
    DecompositionError,
    compositional_estimability,
    denominator_disagreements,
    student_t_companion,
    summary_frame,
)

GENES = ("ACTB", "CDX2", "GUCA2A")


def _deltas(n_patients: int = 12, *, rung: str = "lineage",
            frac_normal: float = 0.5, frac_tumour: float = 0.3,
            resolved: int = 200, mature: int = 100) -> pd.DataFrame:
    rng = np.random.default_rng(4)
    rows = []
    for i in range(n_patients):
        for gene in GENES:
            loss = 0.2 if gene == "GUCA2A" else 0.95
            base = 20.0 if gene != "CDX2" else 1.5
            rows.append({
                "patient_id": f"P{i}", "study_id": "Chen_2021_Cell",
                "gene": gene, "granularity_rung": rung,
                "cp10k_normal": base * float(rng.uniform(0.9, 1.1)),
                "cp10k_tumour": base * loss * float(rng.uniform(0.9, 1.1)),
                "n_tumour": mature, "n_normal": mature,
                "frac_mature_normal": frac_normal,
                "frac_mature_tumour": frac_tumour,
                "frac_mature_normal_all_epithelial": frac_normal * 0.75,
                "frac_mature_tumour_all_epithelial": frac_tumour * 0.75,
                "n_cells_resolved_normal": resolved,
                "n_cells_resolved_tumour": resolved,
                "unresolved_share_normal": 0.25,
                "unresolved_share_tumour": 0.25,
            })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# What the frozen rung spec determines, independently of any data
# ---------------------------------------------------------------------------


def test_the_epithelial_rung_returns_an_exactly_zero_compositional_term():
    """The lower bound of the granularity curve, and it is arithmetic.

    At `epithelial` every epithelial cell is mature, so both arms' mature
    fraction is 1.0, Δf is exactly 0, and the whole change must land in the
    intrinsic term. `RUNG_SPECS` says this rung "is supposed to look degenerate
    — that is what it demonstrates," and a curve whose lower bound is not
    degenerate is not measuring what it claims.
    """
    deltas = _deltas(rung="epithelial", frac_normal=1.0, frac_tumour=1.0)
    split = decompose_cohort(summary_frame(deltas, denominator=PRIMARY_DENOMINATOR))
    assert (split["compositional"] == 0.0).all()
    assert (split["interaction"] == 0.0).all()
    intrinsic = pd.to_numeric(split["intrinsic"], errors="coerce")
    assert intrinsic.abs().max() > 0, "the change has to land somewhere"


def test_a_starved_rung_writes_none_not_zero_into_intrinsic():
    """Invariant 1, through the real estimator rather than a re-derivation."""
    deltas = _deltas(mature=3)
    split = decompose_cohort(summary_frame(deltas, denominator=PRIMARY_DENOMINATOR))
    starved = split[split["estimability"] == "not_estimable"]
    assert not starved.empty
    assert starved["intrinsic"].isna().all()
    assert starved["compositional"].notna().all(), (
        "an unestimable INTRINSIC term does not make the compositional term "
        "unestimable; the row must not be dropped"
    )


# ---------------------------------------------------------------------------
# The inputs that force each guard
# ---------------------------------------------------------------------------


def test_a_table_without_mature_fractions_is_refused_not_defaulted():
    """THE ONE THAT MATTERS. A defaulted fraction is a fabricated compositional
    term, and it would look exactly like a real one."""
    deltas = _deltas().drop(columns=["frac_mature_normal", "frac_mature_tumour"])
    with pytest.raises(DecompositionError, match="--collect-fractions"):
        summary_frame(deltas, denominator=PRIMARY_DENOMINATOR)


def test_the_sensitivity_denominator_is_refused_when_absent_too():
    deltas = _deltas().drop(columns=["frac_mature_normal_all_epithelial",
                                     "frac_mature_tumour_all_epithelial"])
    with pytest.raises(DecompositionError):
        summary_frame(deltas, denominator="all_epithelial")


def test_compositional_estimability_is_not_inherited_from_the_intrinsic_call():
    """It gates on a different count, under different (non-provisional) cutpoints."""
    deltas = _deltas().drop(columns=["n_cells_resolved_normal",
                                     "n_cells_resolved_tumour"])
    with pytest.raises(DecompositionError, match="not the same question"):
        compositional_estimability(deltas)


def test_compositional_estimability_can_differ_from_the_intrinsic_one():
    """A row may have an estimable fraction and an unestimable intrinsic term.

    That is the ordinary case at a starved rung, and folding the two into one
    verdict would report the compositional term as unavailable when it is not.
    """
    deltas = _deltas(mature=3, resolved=500)
    comp = compositional_estimability(deltas)
    split = decompose_cohort(summary_frame(deltas, denominator=PRIMARY_DENOMINATOR))
    assert (comp["compositional_estimability"] == "ok").all()
    assert (split["estimability"] == "not_estimable").all()


def test_the_denominator_comparison_fires_when_the_answer_flips():
    """The sensitivity analysis that replaced the threshold gate.

    Built so the two denominators give opposite-signed compositional terms:
    the primary has the tumour arm LESS mature, the sensitivity arm MORE.
    """
    deltas = _deltas(frac_normal=0.5, frac_tumour=0.3)
    deltas["frac_mature_normal_all_epithelial"] = 0.3
    deltas["frac_mature_tumour_all_epithelial"] = 0.5
    companions = {}
    for name in DENOMINATORS:
        summary = summary_frame(deltas, denominator=name)
        companions[name] = student_t_companion(
            summary, decompose_cohort(summary), seed=1
        )
    disputed = denominator_disagreements(companions)
    assert not disputed.empty, (
        "a sign flip between denominators must be reported, not averaged"
    )
    assert (disputed["term"] == "compositional").any()


def test_the_denominator_comparison_stays_quiet_when_they_agree():
    """The other half: it must not call everything denominator-dependent."""
    deltas = _deltas()
    deltas["frac_mature_normal_all_epithelial"] = deltas["frac_mature_normal"]
    deltas["frac_mature_tumour_all_epithelial"] = deltas["frac_mature_tumour"]
    companions = {}
    for name in DENOMINATORS:
        summary = summary_frame(deltas, denominator=name)
        companions[name] = student_t_companion(
            summary, decompose_cohort(summary), seed=1
        )
    assert denominator_disagreements(companions).empty


def test_the_companion_carries_the_schema_bands_real_error_rate():
    """So a reader comparing the two intervals never has to leave the table."""
    deltas = _deltas(n_patients=20)
    summary = summary_frame(deltas, denominator=PRIMARY_DENOMINATOR)
    companion = student_t_companion(summary, decompose_cohort(summary), seed=1)
    assert companion["percentile_band_false_positive_rate"].notna().all()
    assert companion["percentile_band_false_positive_rate"].max() == pytest.approx(
        0.071, abs=0.002
    ), "n=20 -> 7.1%, the number the prereg quotes"
