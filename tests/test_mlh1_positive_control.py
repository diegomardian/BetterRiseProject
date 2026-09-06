"""The positive control's verdict function, against each branch it can take.

``instrument_verdict`` is where the pre-registration's §5 stops being prose and
becomes a decision. That makes it exactly the kind of thing this repository
keeps finding broken: a function with five documented outcomes and, until
somebody constructs the input for each, no evidence that more than one of them
is reachable. A verdict function that can only ever say one thing is a check
that cannot fail wearing a different hat.

So there is one test per branch, each forcing its own outcome, plus one that
asserts the branches are mutually exclusive on the same input.

The job itself needs the 30 GB atlas and cannot run here. What can run is every
decision it makes once the cells have been scored, which is where the
interpretation lives.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.reference.jobs.mlh1_positive_control import (
    MIN_PATIENTS_WITH_SIGNAL,
    PRIMARY_STRATUM,
    SECONDARY_STRATUM,
    TARGET,
    UNDERPOWERED_STRATUM,
    arm_of,
    arm_reading,
    instrument_verdict,
)

HELD = (True, "holds over 10 patients on 2 control(s)")
FAILED = (False, "UNRESOLVED: control ACTB +0.31 [-0.02, +0.60] on detection")


def _primary(**overrides) -> dict:
    base = {
        "gene": TARGET, "arm": PRIMARY_STRATUM, "n_patients": 10,
        "patients_with_signal": 9, "estimability": "estimable",
        "mean_delta_cloglog": -0.8, "ci_low": -1.2, "ci_high": -0.4,
        "excludes_zero": True, "detail": "",
    }
    return base | overrides


# ---------------------------------------------------------------------------
# One test per pre-registered branch
# ---------------------------------------------------------------------------


def test_a_fall_that_excludes_zero_says_the_instrument_sees_silencing():
    out = instrument_verdict(_primary(), HELD)
    assert out["verdict"] == "INSTRUMENT SEES KNOWN SILENCING"
    assert "without" in out["detail"].lower(), (
        "the consequence must carry its own limit: there is no mechanistic "
        "control arm behind this result"
    )


def test_a_rise_that_excludes_zero_is_the_falsifier_and_not_a_success():
    """§5: an instrument firing the wrong way on a known event is not calibrated.

    THE BRANCH MOST AT RISK OF BEING WRITTEN AS A PASS. The interval excludes
    zero and the gene moved, which is what "it worked" looks like from a
    distance. The direction is the whole test.
    """
    out = instrument_verdict(
        _primary(mean_delta_cloglog=+0.9, ci_low=+0.4, ci_high=+1.4), HELD
    )
    assert out["verdict"] == "WRONG DIRECTION"
    assert "not rehabilitated" in out["detail"]


def test_an_interval_containing_zero_says_the_nulls_stay_uninformative():
    out = instrument_verdict(
        _primary(mean_delta_cloglog=-0.1, ci_low=-0.5, ci_high=+0.3,
                 excludes_zero=False),
        HELD,
    )
    assert out["verdict"] == "INSTRUMENT DOES NOT SEE IT"
    assert "moderate silencing" in out["detail"], (
        "a null here is informative against STRONG silencing only, and the "
        "power qualification must travel with the verdict"
    )


def test_a_failed_premise_is_uninterpretable_rather_than_negative():
    """The distinction the whole project turns on: undecided is not negative."""
    out = instrument_verdict(_primary(), FAILED)
    assert out["verdict"] == "UNINTERPRETABLE"
    assert "not a negative result" in out["detail"]


def test_too_few_patients_with_signal_is_not_estimable_and_not_zero():
    """Invariant 1, at the one place it would be easiest to violate.

    With MLH1 undetected in most patients the delta is a difference between two
    boundary corrections -- arithmetic about the pseudocount, not a measurement
    of a gene. That is `not_estimable`, and writing it as "no silencing" is the
    single most likely route to a wrong conclusion in this project.
    """
    out = instrument_verdict(
        _primary(patients_with_signal=MIN_PATIENTS_WITH_SIGNAL - 1), HELD
    )
    assert out["verdict"] == "NOT ESTIMABLE"
    assert "invariant 1" in out["detail"].lower()
    assert "not be written as zero" in out["detail"]


def test_an_unestimable_arm_is_not_estimable_whatever_its_numbers_look_like():
    out = instrument_verdict(
        _primary(estimability="not_estimable", detail="1 patient(s)"), HELD
    )
    assert out["verdict"] == "NOT ESTIMABLE"


def test_the_premise_is_checked_before_anything_else():
    """A failed premise wins over every other branch, including a clean fall."""
    strong_fall = _primary(mean_delta_cloglog=-2.0, ci_low=-2.5, ci_high=-1.5)
    assert instrument_verdict(strong_fall, FAILED)["verdict"] == "UNINTERPRETABLE"
    assert instrument_verdict(strong_fall, HELD)["verdict"] != "UNINTERPRETABLE"


def test_every_branch_is_reachable_and_they_are_distinct():
    """Five documented outcomes, five reachable, no two the same."""
    cases = [
        (_primary(), HELD),
        (_primary(mean_delta_cloglog=+0.9, ci_low=+0.4, ci_high=+1.4), HELD),
        (_primary(excludes_zero=False, ci_low=-0.5, ci_high=+0.3), HELD),
        (_primary(), FAILED),
        (_primary(patients_with_signal=0), HELD),
    ]
    verdicts = {instrument_verdict(p, q)["verdict"] for p, q in cases}
    assert len(verdicts) == 5, f"branches collapsed onto {verdicts}"


# ---------------------------------------------------------------------------
# Arm assignment and the reading itself
# ---------------------------------------------------------------------------


def test_the_arms_are_the_pre_registered_ones():
    assert arm_of("mlh1_methylated") == PRIMARY_STRATUM
    assert arm_of("mlh1_intact_mmrd") == UNDERPOWERED_STRATUM
    assert arm_of("mmr_proficient") == SECONDARY_STRATUM
    assert arm_of("mlh1_deficient_unmethylated") == SECONDARY_STRATUM


def _scored(n_patients: int, detect_normal: float, detect_tumour: float,
            cells: int = 262) -> pd.DataFrame:
    return pd.DataFrame({
        "patient_id": [f"P{i}" for i in range(n_patients)],
        "gene": TARGET,
        "n_normal": cells, "n_tumour": cells,
        "detect_normal": detect_normal, "detect_tumour": detect_tumour,
    })


def test_a_reading_reports_the_calibrated_interval_not_the_bootstrap():
    row = arm_reading(_scored(10, 0.032, 0.012), gene=TARGET,
                      arm=PRIMARY_STRATUM, seed=1)
    assert row["interval_method"] == "student_t"
    assert row["mean_delta_cloglog"] < 0
    assert row["ci_low"] < row["ci_high"]


def test_a_reading_carries_the_rate_of_the_interval_it_is_NOT_using():
    """The n=4 arm's uselessness has to be on its own row, not in a footnote.

    ``percentile_false_positive_rate`` is what the project's usual interval
    would do at this n. It is carried so that a reader who quotes the n=4 arm
    sees, in the same table, why they should not.
    """
    row = arm_reading(_scored(4, 0.023, 0.020), gene=TARGET,
                      arm=UNDERPOWERED_STRATUM, seed=1)
    assert row["percentile_false_positive_rate"] == pytest.approx(0.188, abs=0.002)


def test_an_arm_of_one_patient_is_not_estimable_rather_than_a_point():
    row = arm_reading(_scored(1, 0.03, 0.01), gene=TARGET,
                      arm=UNDERPOWERED_STRATUM, seed=1)
    assert row["estimability"] == "not_estimable"
    assert np.isnan(row["mean_delta_cloglog"])
    assert row["excludes_zero"] is False


def test_patients_with_signal_counts_the_normal_arm_not_the_tumour_arm():
    """Silencing empties the TUMOUR arm; the question is whether the reference
    arm ever had the gene to lose."""
    frame = _scored(10, 0.032, 0.0)
    frame.loc[frame.index[:6], "detect_normal"] = 0.0
    row = arm_reading(frame, gene=TARGET, arm=PRIMARY_STRATUM, seed=1)
    assert row["patients_with_signal"] == 4
    assert instrument_verdict(
        row | {"arm": PRIMARY_STRATUM}, HELD
    )["verdict"] == "NOT ESTIMABLE"


def test_the_strata_merge_does_not_suffix_away_the_patient_id():
    """The bug this would have been: a KeyError after reading 30 GB.

    ``strata_for`` returns a frame carrying its own ``patient_id`` in the SHORT
    form, while the scored deltas carry it in the long ``Study.Cxxx`` form. A
    plain merge on ``short_id`` suffixes both into ``patient_id_x`` /
    ``patient_id_y`` and the bare column ceases to exist -- so every downstream
    ``drop_duplicates("patient_id")`` raises, on the cluster, after the atlas
    has been loaded. Nothing about the merge itself complains.
    """
    deltas = pd.DataFrame({
        "patient_id": ["Pelka_2021_Cell.C110", "Pelka_2021_Cell.C111"],
        "gene": TARGET, "n_normal": 262, "n_tumour": 262,
        "detect_normal": 0.032, "detect_tumour": 0.011,
    })
    deltas["short_id"] = deltas["patient_id"].astype(str).str.split(".").str[-1]
    strata = pd.DataFrame({
        "short_id": ["C110", "C111"], "atlas_methylation": "meth",
        "patient_id": ["C110", "C111"],
        "mlh1_stratum": "mlh1_methylated", "matched": True,
    })

    naive = deltas.merge(strata, on="short_id", how="left")
    assert "patient_id" not in naive.columns, (
        "if this ever passes, pandas changed its suffixing and the guard in "
        "the job is no longer describing a real failure mode"
    )

    fixed = deltas.merge(strata.drop(columns=["patient_id"]),
                         on="short_id", how="left")
    assert "patient_id" in fixed.columns
    row = arm_reading(fixed, gene=TARGET, arm=PRIMARY_STRATUM, seed=1)
    assert row["n_patients"] == 2
