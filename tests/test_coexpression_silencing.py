"""The silencing reading, against the inputs that must break it.

This job exists because the decomposition cannot answer the mechanism question
on this panel. That makes it the kind of analysis that gets believed too easily,
so every claim it makes is tested here against data constructed to defeat it.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.reference.jobs.coexpression_silencing import (
    CONTROL_LOG2_TOLERANCE,
    CONTROL_TOLERANCE,
    GENE_ROLES,
    MIN_PREMISE_PATIENTS,
    premise_holds,
    summarise,
)


def _deltas(n_patients: int = 5, *, detection=0.5, log2=0.0, **genes) -> pd.DataFrame:
    """Per-patient rows, which is what premise_holds takes.

    ``detection`` sets the arm rates, so a test can put a control at ceiling
    (0.99) or leave it room to move (0.5). ``log2`` is the expression ratio the
    premise falls back to when detection is saturated.
    """
    return pd.DataFrame([
        {"gene": g, "patient_id": f"p{i}", "delta_detect": v,
         "detect_normal": detection, "detect_tumour": detection + v,
         "log2_cp10k_ratio": log2}
        for g, v in genes.items() for i in range(n_patients)
    ])


# ---------------------------------------------------------------------------
# The premise control
# ---------------------------------------------------------------------------


def test_the_premise_is_refused_when_a_control_gene_moves():
    """The violating input: the population itself changed.

    If ACTB detection has shifted, the cells scored in the two arms are not the
    same kind of cell, and a marker falling inside them is not silencing. This
    is the failure mode the whole reading turns on, so it must refuse rather
    than report.
    """
    holds, reading = premise_holds(_deltas(ACTB=-0.4, KRT8=-0.01, GUCA2A=-0.5))
    assert not holds
    assert "REFUSED" in reading and "ACTB" in reading


def test_the_premise_holds_when_the_controls_are_still():
    holds, reading = premise_holds(_deltas(ACTB=-0.01, KRT8=0.02, GUCA2A=-0.5))
    assert holds and "holds" in reading


def test_a_control_just_past_the_tolerance_refuses():
    """The boundary is a decision, so it is tested rather than assumed."""
    assert premise_holds(_deltas(ACTB=CONTROL_TOLERANCE - 0.001, KRT8=0.0))[0]
    assert not premise_holds(_deltas(ACTB=CONTROL_TOLERANCE + 0.001, KRT8=0.0))[0]


def test_scoring_no_control_gene_refuses_rather_than_passing():
    """A premise nothing tested is not a premise that held. This is the
    project's rule about checks that cannot fire, applied to itself."""
    holds, reading = premise_holds(_deltas(GUCA2A=-0.5, CDX2=-0.1))
    assert not holds
    assert "UNDEFINED" in reading


def test_every_declared_control_is_a_gene_the_job_scores():
    """A tolerance on a gene that is never measured protects nothing."""
    controls = [g for g, r in GENE_ROLES.items() if r == "control"]
    assert len(controls) >= 2, "one control cannot distinguish a shift from noise"


# ---------------------------------------------------------------------------
# The summary must not pool studies
# ---------------------------------------------------------------------------


def test_studies_are_summarised_separately_and_never_pooled():
    """CLAUDE.md invariant 4. Pooling here is not a stylistic preference: the
    two cohorts disagree about whether GUCA2A falls further than the identity
    markers, and a pooled mean reports a cleaner effect than either study
    supports. The violating input is exactly that disagreement.
    """
    rows = []
    for i in range(6):
        rows.append({"study_id": "A", "gene": "GUCA2A", "patient_id": f"a{i}",
                     "delta_detect": -0.60, "delta_given_conditioner": -0.60})
        rows.append({"study_id": "B", "gene": "GUCA2A", "patient_id": f"b{i}",
                     "delta_detect": -0.05, "delta_given_conditioner": -0.05})
    out = summarise(pd.DataFrame(rows))
    detection = out[out["statistic"] == "detection"].set_index("study_id")

    assert set(detection.index) == {"A", "B"}, "one row per study, not one pooled row"
    assert detection.loc["A", "mean_delta"] == pytest.approx(-0.60, abs=1e-9)
    assert detection.loc["B", "mean_delta"] == pytest.approx(-0.05, abs=1e-9)

    pooled = np.mean([-0.60] * 6 + [-0.05] * 6)
    assert not np.isclose(detection.loc["A", "mean_delta"], pooled)
    assert not np.isclose(detection.loc["B", "mean_delta"], pooled)


def test_a_study_with_too_few_patients_gets_no_interval():
    """Two patients cannot carry a bootstrap over patients. The mean is still
    reported; the interval is not invented."""
    rows = [{"study_id": "A", "gene": "GUCA2A", "patient_id": p,
             "delta_detect": -0.5, "delta_given_conditioner": -0.5} for p in ("p1", "p2")]
    out = summarise(pd.DataFrame(rows))
    row = out[out["statistic"] == "detection"].iloc[0]
    assert row["n_patients"] == 2
    assert not np.isfinite(row["ci_low"]) and not np.isfinite(row["ci_high"])
    assert not row["excludes_zero"]


def test_the_bootstrap_resamples_patients_not_rows():
    """Invariant 5. A gene measured on three patients has three units of
    inference however many cells sit behind them, so the interval must be wide."""
    rows = [{"study_id": "A", "gene": "GUCA2A", "patient_id": p, "delta_detect": d,
             "delta_given_conditioner": d}
            for p, d in zip(("p1", "p2", "p3"), (-0.9, -0.1, -0.5), strict=True)]
    row = summarise(pd.DataFrame(rows)).iloc[0]
    assert row["n_patients"] == 3
    assert row["ci_high"] - row["ci_low"] > 0.3, (
        "three disagreeing patients must not produce a tight interval"
    )


def test_the_premise_is_undefined_on_too_few_patients():
    """The violating input: one patient.

    A control's mean over one patient IS that patient. Comparing a single noisy
    number against a tolerance reports clean whatever the data did, which is the
    exact shape of failure this repository is about. The smoke run on GSE178341
    returned "holds" off one patient with KRT8 at -0.050, three times either Lee
    cohort, and nothing in the guard noticed.
    """
    holds, reading = premise_holds(_deltas(n_patients=1, ACTB=-0.01, KRT8=-0.05))
    assert not holds
    assert "UNDEFINED" in reading and "below the" in reading

    enough = premise_holds(_deltas(n_patients=MIN_PREMISE_PATIENTS, ACTB=-0.01, KRT8=-0.02))
    assert enough[0], "the floor must not refuse a genuinely quiet control"


def test_a_control_at_ceiling_is_assessed_on_expression_not_detection():
    """The defect this caught, and the reason it was invisible.

    ACTB and KRT8 are detected in 0.99-1.00 of cells in BOTH arms of both Lee
    cohorts, so their detection rate cannot fall and reports "no change"
    whatever happens. The first version of this guard passed the premise on
    exactly that, and the underlying expression was moving by 1.5x
    simultaneously. A control at ceiling must be asked a question it can answer.
    """
    # Saturated on detection, and its expression moved a long way: must refuse.
    holds, reading = premise_holds(
        _deltas(detection=0.99, log2=-1.2, ACTB=-0.005, KRT8=-0.004)
    )
    assert not holds
    assert "log2 expression" in reading and "REFUSED" in reading

    # Saturated on detection, expression quiet: the premise may hold, and the
    # reading has to say the controls were judged on expression.
    holds, reading = premise_holds(
        _deltas(detection=0.99, log2=0.05, ACTB=-0.005, KRT8=-0.004)
    )
    assert holds and "saturated on detection" in reading


def test_the_expression_tolerance_is_a_boundary_that_is_tested():
    assert premise_holds(
        _deltas(detection=0.99, log2=CONTROL_LOG2_TOLERANCE - 0.01, ACTB=0.0, KRT8=0.0))[0]
    assert not premise_holds(
        _deltas(detection=0.99, log2=CONTROL_LOG2_TOLERANCE + 0.01, ACTB=0.0, KRT8=0.0))[0]


def test_the_premise_is_undefined_without_the_arm_rates():
    """Whether a control had room to fall cannot be assessed from the delta
    alone, and a premise that cannot be assessed has not been satisfied."""
    bare = pd.DataFrame([{"gene": "ACTB", "patient_id": f"p{i}", "delta_detect": 0.0}
                         for i in range(5)])
    holds, reading = premise_holds(bare)
    assert not holds and "UNDEFINED" in reading


def test_a_control_astride_the_tolerance_is_unresolved_not_refused():
    """The third state, and why it exists.

    GSE144735's ACTB read +0.486 against a 0.5 tolerance under one seed and
    +0.540 under another: the same data, a different draw, the premise flipping
    from held to refused. A verdict that moves with the seed has not been
    reached. Refused and satisfied are both claims; undecided is the honest
    third answer and the project already draws this distinction for cutpoints.
    """
    astride = pd.DataFrame([
        {"gene": "ACTB", "patient_id": f"p{i}", "delta_detect": -0.002,
         "detect_normal": 0.99, "detect_tumour": 0.988, "log2_cp10k_ratio": v}
        for i, v in enumerate([0.1, 0.3, 0.5, 0.7, 0.9, 1.1])
    ])
    holds, reading = premise_holds(astride)
    assert not holds
    assert "UNRESOLVED" in reading and "straddles" in reading

    # Far enough past it that no resample lands inside: a real refusal.
    clear = astride.assign(log2_cp10k_ratio=[2.0, 2.1, 2.2, 2.3, 2.4, 2.5])
    holds, reading = premise_holds(clear)
    assert not holds and "REFUSED" in reading

    # Tightly inside it: genuinely held.
    quiet = astride.assign(log2_cp10k_ratio=[0.01, 0.02, 0.0, -0.01, 0.02, 0.0])
    assert premise_holds(quiet)[0]


def test_the_premise_verdict_does_not_move_with_the_seed_when_it_is_settled():
    """A settled verdict is one the draw cannot change. An unsettled one is
    exactly what UNRESOLVED is for."""
    quiet = pd.DataFrame([
        {"gene": "ACTB", "patient_id": f"p{i}", "delta_detect": -0.002,
         "detect_normal": 0.99, "detect_tumour": 0.988, "log2_cp10k_ratio": 0.01}
        for i in range(6)
    ])
    assert {premise_holds(quiet, seed=s)[0] for s in (1, 2, 3, 20260101, 20260904)} == {True}
