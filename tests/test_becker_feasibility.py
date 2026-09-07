"""B1's gate, against the outcomes it has to be able to return.

The gate exists to STOP the replication cheaply, so the test that matters is
that it can. A gate that returns FULL DESIGN on everything would let the whole
avenue proceed on a cohort where the target gene is undetectable — which is the
specific, physical risk B1 carries, since GUCA2A and MS4A12 are cytoplasmic and
snRNA-seq samples nuclei.

One test per branch of `docs/prereg_becker_replication.md` §3, plus the two
things the job refuses to assume.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.reference.jobs.becker_feasibility import (
    CRITICAL_GENE,
    MIN_DETECTION,
    MIN_PATIENT_SHARE_NONZERO,
    FeasibilityError,
    detection_table,
    gate,
    verdict,
)

PANEL = ("ACTB", "KRT8", "EPCAM", "CDX2", "MS4A12", "GUCA2A")


def _detection(**overrides) -> pd.DataFrame:
    base = dict.fromkeys(PANEL, 0.60)
    base.update(overrides)
    return pd.DataFrame({
        "gene": list(base),
        "detection": list(base.values()),
        "share_patients_nonzero": [1.0] * len(base),
    })


def test_the_gate_stops_the_replication_when_the_critical_gene_fails():
    """THE BRANCH THE GATE EXISTS FOR.

    GUCA2A under the floor ends it: the primary claim is its four cross-block
    contrasts, so there is nothing left to test. And the verdict must say this
    is a measurement about snRNA-seq rather than a negative about biology —
    that distinction is the whole reason the gate runs before the analysis.
    """
    out = verdict(gate(_detection(GUCA2A=0.04)))
    assert out["verdict"] == "CANNOT RUN"
    assert "not a negative result about the biology" in out["detail"]


def test_a_missing_critical_gene_blames_the_identifier_space_first():
    frame = _detection()
    out = verdict(gate(frame[frame.gene != CRITICAL_GENE]))
    assert out["verdict"] == "CANNOT RUN"
    assert "identifier space" in out["detail"]


def test_ms4a12_alone_failing_costs_the_secondary_claim_only():
    out = verdict(gate(_detection(MS4A12=0.05)))
    assert out["verdict"] == "PRIMARY ONLY"
    assert "not gene-specific" in out["detail"]
    assert "single-cohort" in out["detail"]


def test_a_failing_control_shrinks_the_comparator_set():
    out = verdict(gate(_detection(CDX2=0.02)))
    assert out["verdict"] == "REDUCED PANEL"
    assert "fewer than 8" in out["detail"]


def test_a_clean_panel_passes():
    """The other half: the gate must not refuse everything."""
    assert verdict(gate(_detection()))["verdict"] == "FULL DESIGN"


def test_the_patient_share_gate_bites_independently_of_detection():
    """A gene detected overall but concentrated in a few patients fails.

    Pooled detection can clear the floor while most patients carry none of it,
    and a per-patient design cannot use a gene like that.
    """
    frame = _detection()
    frame.loc[frame.gene == CRITICAL_GENE, "share_patients_nonzero"] = 0.4
    out = verdict(gate(frame))
    assert out["verdict"] == "CANNOT RUN"
    assert f"{MIN_PATIENT_SHARE_NONZERO:.0%}" in out["detail"]


def test_the_gate_reports_the_fold_change_against_chen():
    """A gene that clears the floor after dropping 4x is still news."""
    gated = gate(_detection(GUCA2A=0.11))
    row = gated.loc[gated.gene == "GUCA2A"].iloc[0]
    assert row["passes"]
    assert row["fold_vs_chen"] == pytest.approx(0.11 / 0.437, rel=1e-6)


def test_a_detection_table_without_its_columns_is_refused():
    with pytest.raises(FeasibilityError, match="missing"):
        gate(pd.DataFrame({"gene": ["GUCA2A"], "detection": [0.5]}))


# ---------------------------------------------------------------------------
# The counting itself
# ---------------------------------------------------------------------------


def test_detection_and_patient_share_are_counted_correctly():
    counts = np.array([
        [5, 0],   # P0 detects gene 0
        [0, 0],   # P0
        [0, 3],   # P1 detects gene 1 only
        [0, 0],   # P1
    ], dtype=float)
    table = detection_table(counts, {"GUCA2A": 0, "ACTB": 1},
                            ["P0", "P0", "P1", "P1"]).set_index("gene")
    assert table.loc["GUCA2A", "detection"] == pytest.approx(0.25)
    assert table.loc["GUCA2A", "share_patients_nonzero"] == pytest.approx(0.5)
    assert table.loc["ACTB", "share_patients_nonzero"] == pytest.approx(0.5)
    assert table.loc["GUCA2A", "n_patients"] == 2


def test_the_thresholds_are_the_pre_registered_ones():
    """Pinned, because a gate whose numbers drift is not a pre-registration."""
    assert MIN_DETECTION == 0.10
    assert MIN_PATIENT_SHARE_NONZERO == 0.75
    assert CRITICAL_GENE == "GUCA2A"
