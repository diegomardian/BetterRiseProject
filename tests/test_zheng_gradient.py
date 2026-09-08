from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix

from src.reference.jobs.zheng_gradient import (
    EXPECTED_PATIENTS,
    MIN_EPITHELIAL_PER_STAGE,
    PANEL,
    STAGES,
    STUDY_ID,
    ZhengGradientError,
    eligible_gradient_rows,
    panel_indices,
    summarise_patient,
    validate_specification,
)


def _obs(*, patients: int = EXPECTED_PATIENTS) -> pd.DataFrame:
    rows = []
    stages = {"normal": "adjacent normal", "polyp": "polyp", "carcinoma": "primary tumor"}
    for patient in range(patients):
        for stage, sample_type in stages.items():
            for _cell in range(MIN_EPITHELIAL_PER_STAGE):
                rows.append({
                    "study_id": STUDY_ID,
                    "enrichment_cell_types": "naive",
                    "sample_type": sample_type,
                    "atlas_cell_type_coarse": "Epithelial cell",
                    "patient_id": f"P{patient}",
                    "sample_id": f"P{patient}-{stage}",
                })
    return pd.DataFrame(rows)


def test_pre_read_provenance_guard_accepts_annotation_defined_population():
    assert validate_specification() == ()


def test_three_stage_cohort_is_fixed_at_exactly_three_patients():
    rows, patients = eligible_gradient_rows(_obs())
    assert len(rows) == EXPECTED_PATIENTS * len(STAGES) * MIN_EPITHELIAL_PER_STAGE
    assert len(patients) == EXPECTED_PATIENTS
    with pytest.raises(ZhengGradientError, match="expected 3"):
        eligible_gradient_rows(_obs(patients=EXPECTED_PATIENTS + 1))


def test_frozen_panel_must_be_present_exactly_once():
    with pytest.raises(ZhengGradientError, match="missing symbols"):
        panel_indices(PANEL[:-1])
    with pytest.raises(ZhengGradientError, match="duplicated symbols"):
        panel_indices(list(PANEL) + [PANEL[0]])


def test_stage_with_no_qc_cells_is_retained_as_unmeasurable_not_zero():
    symbols = list(PANEL)
    matrix = csr_matrix(np.ones((3, len(symbols)), dtype=int))
    metrics = pd.DataFrame({
        "batch": ["normal", "polyp", "carcinoma"],
        "n_counts": [6, 6, 6],
        "n_genes": [6, 6, 6],
    })
    trajectories, availability = summarise_patient(
        patient_id="P0",
        stage=np.asarray(STAGES),
        compartment=np.asarray(["epithelial"] * 3),
        matrix=matrix,
        metrics=metrics,
        keep=np.asarray([True, True, False]),
        indices=panel_indices(symbols),
    )
    table = pd.DataFrame(trajectories)
    carcinoma = table[table["stage"] == "carcinoma"]
    assert carcinoma["detection"].isna().all()
    assert not carcinoma["stage_measurable_after_qc"].any()
    availability_table = pd.DataFrame(availability).set_index("stage")
    assert availability_table.loc["carcinoma", "n_epithelial_after_qc"] == 0
