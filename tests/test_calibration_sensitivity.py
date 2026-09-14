from __future__ import annotations

import itertools

import pandas as pd

from src.harness.calibration_sensitivity import (
    holdout_pairs,
    influence_summary,
    retained_pair_rows,
)


def test_holdout_pairs_enumerate_each_pair_once():
    patients = [f"P{i}" for i in range(10)]
    pairs = holdout_pairs(patients)
    assert len(pairs) == 45
    assert set(pairs) == set(itertools.combinations(sorted(patients), 2))
    assert all(sum(patient in pair for pair in pairs) == 9 for patient in patients)


def test_influence_summary_counts_and_names_changed_patients():
    influence = pd.DataFrame(
        {
            "cohort": ["smc"] * 3,
            "criterion": ["coverage_only"] * 3,
            "omitted_patient": ["P1", "P2", "P3"],
            "conclusion_changed": [True, False, True],
        }
    )
    summary = influence_summary(influence).iloc[0]
    assert summary["n_patients"] == 3
    assert summary["n_changing_conclusion"] == 2
    assert summary["changing_patients"] == "P1,P3"


def test_patient_omission_drops_only_pairs_containing_that_patient():
    sweep = pd.DataFrame(
        {
            "held_out_pair": ["P1|P2", "P1|P3", "P2|P3"],
            "seed": [11, 12, 13],
        }
    )
    retained = retained_pair_rows(sweep, "P1")
    assert retained.to_dict("records") == [{"held_out_pair": "P2|P3", "seed": 13}]
