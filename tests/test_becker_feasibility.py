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
    out = verdict(gate(_detection(GUCA2A=0.04)), mature_labelled=True)
    assert out["verdict"] == "CANNOT RUN"
    assert "not a negative result about the biology" in out["detail"]


def test_a_missing_critical_gene_blames_the_identifier_space_first():
    frame = _detection()
    out = verdict(gate(frame[frame.gene != CRITICAL_GENE]), mature_labelled=True)
    assert out["verdict"] == "CANNOT RUN"
    assert "identifier space" in out["detail"]


def test_ms4a12_alone_failing_costs_the_secondary_claim_only():
    out = verdict(gate(_detection(MS4A12=0.05)), mature_labelled=True)
    assert out["verdict"] == "PRIMARY ONLY"
    assert "not gene-specific" in out["detail"]
    assert "single-cohort" in out["detail"]


def test_a_failing_control_shrinks_the_comparator_set():
    out = verdict(gate(_detection(CDX2=0.02)), mature_labelled=True)
    assert out["verdict"] == "REDUCED PANEL"
    assert "fewer than 8" in out["detail"]


def test_a_clean_panel_passes():
    """The other half: the gate must not refuse everything."""
    assert verdict(gate(_detection()), mature_labelled=True)["verdict"] == "FULL DESIGN"


def test_the_patient_share_gate_bites_independently_of_detection():
    """A gene detected overall but concentrated in a few patients fails.

    Pooled detection can clear the floor while most patients carry none of it,
    and a per-patient design cannot use a gene like that.
    """
    frame = _detection()
    frame.loc[frame.gene == CRITICAL_GENE, "share_patients_nonzero"] = 0.4
    out = verdict(gate(frame), mature_labelled=True)
    assert out["verdict"] == "CANNOT RUN"
    assert f"{MIN_PATIENT_SHARE_NONZERO:.0%}" in out["detail"]


def test_the_gate_reports_the_fold_change_against_chen():
    """A gene that clears the floor after dropping 4x is still news."""
    gated = gate(_detection(GUCA2A=0.11))
    row = gated.loc[gated.gene == "GUCA2A"].iloc[0]
    assert row["passes"]
    # On the detection scale: cloglog(p) = log(mu), so the reported fold is a
    # ratio of expected UMIs, not of probabilities.
    def mu(p):
        return -np.log1p(-p)

    assert row["fold_mu_vs_chen"] == pytest.approx(mu(0.11) / mu(0.437), rel=1e-6)
    assert row["log_fc_vs_chen"] == pytest.approx(
        np.log(mu(0.11) / mu(0.437)), rel=1e-6)
    # The old number is kept, under a name that says what it is.
    assert row["naive_ratio_of_p"] == pytest.approx(0.11 / 0.437, rel=1e-6)


def test_the_ratio_of_probabilities_reverses_the_gene_ordering():
    """The defect, as a test — and it is not hypothetical, it is the run.

    Detection is bounded at 1 and ACTB enters saturated at 0.984, so dividing
    Becker's small number by it compresses the loss. On Becker's own measured
    values the two scales disagree about which gene the nuclear protocol treated
    worse, and they disagree about the pair the pre-registration is ABOUT:
    the ratio scale says the housekeeping control was better retained than the
    target, the detection scale says the opposite.
    """
    measured = _detection(ACTB=0.199808, GUCA2A=0.033845)
    gated = gate(measured).set_index("gene")
    actb, guca = gated.loc["ACTB"], gated.loc["GUCA2A"]

    assert actb["naive_ratio_of_p"] > guca["naive_ratio_of_p"]   # ACTB "better"
    assert guca["fold_mu_vs_chen"] > actb["fold_mu_vs_chen"]     # GUCA2A better
    assert guca["log_fc_vs_chen"] > actb["log_fc_vs_chen"]


def test_a_failing_lower_bound_does_not_get_to_say_cannot_run():
    """The defect this job shipped twice, as a test.

    §3's gate is the mature cells of the reference arm. Nothing in this job
    labels cells, so what it measures is strictly SMALLER than what §3 is
    about. A failure there is a failure of a lower bound, and emitting §3's
    pre-registered "CANNOT RUN" from it claims a verdict that was never taken.
    """
    unlabelled = verdict(gate(_detection(GUCA2A=0.04)), mature_labelled=False)
    assert unlabelled["verdict"] == "GATE NOT RUN — LOWER BOUND FAILS"
    assert "UNDETERMINED" in unlabelled["detail"]
    # and it must say how far short, so the reader can judge the enrichment
    assert "2.50x" in unlabelled["detail"]

    labelled = verdict(gate(_detection(GUCA2A=0.04)), mature_labelled=True)
    assert labelled["verdict"] == "CANNOT RUN"


def test_a_passing_lower_bound_is_allowed_to_pass():
    """The asymmetry §3 already relies on: enrichment only helps.

    A gene clearing the floor WITHOUT the mature-cell restriction clears it
    with one, so a pass needs no hedge. Only failure is one-sided.
    """
    out = verdict(gate(_detection()), mature_labelled=False)
    assert out["verdict"] == "FULL DESIGN"


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


# ---------------------------------------------------------------------------
# The mature label, and the way it fails
# ---------------------------------------------------------------------------


def test_no_marker_is_on_the_panel():
    """Invariant 2. A target in its own label returns the threshold it was given.

    The whole panel is excluded, not just GUCA2A, so ACTB/KRT8/EPCAM stay
    usable as the depth-audit controls INSIDE the label.
    """
    from src.reference.jobs.becker_feasibility import (
        DEPTH_AUDIT_CONTROLS,
        MATURE_MARKERS,
    )

    assert not set(MATURE_MARKERS) & set(PANEL)
    assert CRITICAL_GENE not in MATURE_MARKERS
    assert set(DEPTH_AUDIT_CONTROLS) <= set(PANEL)


def test_a_nucleus_needs_two_markers_not_one():
    from src.reference.jobs.becker_feasibility import label_mature

    counts = np.array([
        [0, 0, 0],   # nothing
        [3, 0, 0],   # one marker  -> not mature
        [3, 2, 0],   # two markers -> mature
        [1, 1, 1],   # three       -> mature
    ], dtype=float)
    mask = label_mature(counts, {"CA1": 0, "CA2": 1, "AQP8": 2}, min_umi=1)
    assert mask.tolist() == [False, False, True, True]


def test_a_label_with_no_located_marker_raises_rather_than_being_empty():
    """An all-False mask reads as 'this tissue has no mature cells'.

    That is a biological claim produced by a failed lookup, which is the exact
    error this repository has made four times in the other direction.
    """
    from src.reference.jobs.becker_feasibility import label_mature

    with pytest.raises(FeasibilityError, match="identifier space"):
        label_mature(np.zeros((4, 3)), {})


def test_the_audit_catches_a_label_that_only_selected_deeper_nuclei():
    """THE FAILURE THIS AUDIT EXISTS FOR.

    "Detects >= 2 markers" correlates with library size, and a deeper nucleus
    detects everything more. If the mature restriction lifts GUCA2A by the same
    factor it lifts the housekeeping controls, nothing was enriched — the gate
    would be cleared by depth. The target must land BEYOND the control band.
    """
    from src.reference.jobs.becker_feasibility import enrichment_audit

    whole = _detection(ACTB=0.30, KRT8=0.15, EPCAM=0.27, GUCA2A=0.08,
                       CDX2=0.06, MS4A12=0.12)
    # every gene lifted by the same multiple: pure depth, no enrichment
    depth_only = whole.copy()
    depth_only["detection"] = 1 - (1 - depth_only["detection"]) ** 2.0
    audit = enrichment_audit(whole, depth_only).set_index("gene")
    assert not bool(audit.loc[CRITICAL_GENE, "beyond_control_band"])

    # now GUCA2A alone lifted further: real enrichment
    enriched = depth_only.copy()
    enriched.loc[enriched.gene == CRITICAL_GENE, "detection"] = 0.42
    audit2 = enrichment_audit(whole, enriched).set_index("gene")
    assert bool(audit2.loc[CRITICAL_GENE, "beyond_control_band"])


def test_the_audit_reports_the_control_band_it_judged_against():
    """A verdict whose threshold is not in the row cannot be re-checked."""
    from src.reference.jobs.becker_feasibility import enrichment_audit

    whole = _detection()
    mature = whole.copy()
    mature["detection"] = 1 - (1 - mature["detection"]) ** 1.5
    audit = enrichment_audit(whole, mature)
    assert audit["control_band_low"].notna().all()
    assert audit["control_band_high"].notna().all()
    assert (audit["control_band_low"] <= audit["control_band_high"]).all()


def test_a_pass_bought_by_depth_does_not_get_reported_as_full_design():
    """THE DEFECT THE 2026-09-07 RUN SHIPPED, on that run's own numbers.

    The labelled gate returned FULL DESIGN while the audit beside it said
    GUCA2A's enrichment sat inside the control band. The audit warned and did
    not bind, so a floor cleared by library size was reported as a design that
    runs. A check that warns without binding is a check that reports success.
    """
    from src.reference.jobs.becker_feasibility import enrichment_audit

    whole = _detection(ACTB=0.304561, KRT8=0.154186, EPCAM=0.269934,
                       GUCA2A=0.078697, MS4A12=0.125358, CDX2=0.061325)
    mature = _detection(ACTB=0.456126, KRT8=0.295804, EPCAM=0.473717,
                        GUCA2A=0.147823, MS4A12=0.334034, CDX2=0.115213)
    audit = enrichment_audit(whole, mature)
    gated = gate(mature)

    # every gene clears the floor, so the unaudited verdict is the optimistic one
    assert gated["passes"].all()
    assert verdict(gated, mature_labelled=True)["verdict"] == "FULL DESIGN"

    # with the audit in hand it is refused, and it says why
    out = verdict(gated, mature_labelled=True, audit=audit)
    assert out["verdict"] == "CLEARED BY DEPTH — NOT LICENSED"
    assert "library size" in out["detail"]
    assert "MS4A12" in out["detail"], "the genuinely enriched gene is the contrast"


def test_a_target_that_outruns_the_controls_still_passes():
    """The other half, so the guard above is not simply a refusal machine."""
    from src.reference.jobs.becker_feasibility import enrichment_audit

    whole = _detection(ACTB=0.30, KRT8=0.15, EPCAM=0.27, GUCA2A=0.08,
                       MS4A12=0.12, CDX2=0.06)
    mature = whole.copy()
    mature["detection"] = 1 - (1 - mature["detection"]) ** 2.0
    mature.loc[mature.gene == CRITICAL_GENE, "detection"] = 0.55   # beyond depth
    audit = enrichment_audit(whole, mature)
    out = verdict(gate(mature), mature_labelled=True, audit=audit)
    assert out["verdict"] == "FULL DESIGN"
