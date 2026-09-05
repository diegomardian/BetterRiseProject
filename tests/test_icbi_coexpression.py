"""The ICBI coexpression path, and the bar it must clear on Pelka.

The adaptation's whole risk is that it diverges from the committed GSE178341
result for a reason that has nothing to do with the science -- a different QC
population, a different batch key, a sorted fraction left in. So the validation
is mechanical, and these tests check that it can FAIL, not only that it passes.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.reference.icbi_slice import SliceError
from src.reference.jobs.icbi_coexpression import (
    BATCH_KEY,
    GUCA2A_DELTA_TOLERANCE,
    MIN_EPITHELIAL_PER_ARM,
    NAIVE,
    VALIDATION,
    eligible_patients,
)


def _obs(n_per_arm: int = 200, enrichment: str = NAIVE, n_patients: int = 3):
    rows = []
    for p in range(n_patients):
        for tissue, sample_type in (("normal", "adjacent normal"),
                                    ("tumour", "primary tumor")):
            for k in range(n_per_arm):
                rows.append({
                    "study_id": "S", "patient_id": f"P{p}",
                    "sample_id": f"P{p}-{tissue}",
                    "sample_type": sample_type,
                    "atlas_cell_type_coarse": "Epithelial cell" if k % 2 else "T cell",
                    "enrichment_cell_types": enrichment,
                })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Patient selection


def test_only_naive_cells_survive_the_sorted_fraction_filter():
    """A CD45-sorted fraction is immune-enriched by construction; leaving it in
    makes every compartment fraction a statement about the sort."""
    mixed = pd.concat([_obs(), _obs(enrichment="CD45+")], ignore_index=True)
    rows, patients = eligible_patients(mixed, "S")
    assert (rows["enrichment_cell_types"] == NAIVE).all()
    assert len(rows) == len(_obs())


def test_a_study_with_no_naive_cells_yields_no_patients():
    rows, patients = eligible_patients(_obs(enrichment="CD45+"), "S")
    assert patients == []


def test_patients_below_the_epithelial_floor_are_excluded():
    """Half the cells are epithelial, so n_per_arm must clear twice the floor."""
    plenty = _obs(n_per_arm=MIN_EPITHELIAL_PER_ARM * 2 + 20)
    thin = _obs(n_per_arm=MIN_EPITHELIAL_PER_ARM // 2)
    assert len(eligible_patients(plenty, "S")[1]) == 3
    assert eligible_patients(thin, "S")[1] == []


def test_a_patient_with_only_one_arm_is_excluded():
    single = _obs()
    single = single[~((single["patient_id"] == "P0") & (single["sample_type"] == "primary tumor"))]
    assert "P0" not in eligible_patients(single, "S")[1]


def test_healthy_normal_does_not_count_as_the_reference_arm():
    donors = _obs()
    donors.loc[donors["sample_type"] == "adjacent normal", "sample_type"] = "healthy normal"
    assert eligible_patients(donors, "S")[1] == []


def test_an_unknown_study_is_refused():
    with pytest.raises(SliceError, match="no cells for study"):
        eligible_patients(_obs(), "NotAStudy")


def test_the_batch_key_matches_gse178341s():
    """If this ever changes, the MAD thresholds move and the numbers with them."""
    assert BATCH_KEY == "sample_id"


# ---------------------------------------------------------------------------
# The validation bar. It must be able to fail.


def _deltas(study: str, *, guca2a: float, actb: float, n: int = 8, seed: int = 0):
    """Per-patient rows in the shape `summarise` consumes."""
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(n):
        for gene, delta in (("GUCA2A", guca2a), ("ACTB", actb),
                            ("CDX2", -0.1), ("KRT8", 0.0),
                            ("EPCAM", 0.0), ("MS4A12", -0.2)):
            rows.append({
                "study_id": study, "patient_id": f"p{i}", "gene": gene,
                "role": "target" if gene == "GUCA2A" else "control",
                "delta_detect": delta + rng.normal(0, 0.02),
                "delta_given_conditioner": delta + rng.normal(0, 0.02),
                "log2_cp10k_ratio": rng.normal(0, 0.1),
                "detect_normal": 0.5, "detect_tumour": 0.5 + delta,
            })
    return pd.DataFrame(rows)


def test_the_bar_passes_a_run_that_matches_the_committed_result(tmp_path, monkeypatch):
    import src.reference.jobs.icbi_coexpression as mod

    committed = _deltas("GSE178341", guca2a=-0.40, actb=-0.01, seed=1)
    root = tmp_path / "results" / "2026-09-04_975cf5c"
    root.mkdir(parents=True)
    committed.to_parquet(root / "coexpression_silencing.parquet")
    monkeypatch.setattr(mod, "RESULTS_DIR", tmp_path / "results")

    close = _deltas("Pelka_2021_Cell", guca2a=-0.42, actb=-0.01, seed=2)
    assert mod.check_against_committed(close)["verdict"] == "PASS"


def test_the_bar_fails_a_run_whose_target_delta_has_drifted(tmp_path, monkeypatch):
    """The check that makes the smoke test worth running at all."""
    import src.reference.jobs.icbi_coexpression as mod

    committed = _deltas("GSE178341", guca2a=-0.40, actb=-0.01, seed=1)
    root = tmp_path / "results" / "2026-09-04_975cf5c"
    root.mkdir(parents=True)
    committed.to_parquet(root / "coexpression_silencing.parquet")
    monkeypatch.setattr(mod, "RESULTS_DIR", tmp_path / "results")

    drifted = _deltas("Pelka_2021_Cell", guca2a=-0.05, actb=-0.01, seed=2)
    result = mod.check_against_committed(drifted)
    assert result["verdict"] == "FAIL"
    assert result["guca2a_drift"] > GUCA2A_DELTA_TOLERANCE


def test_the_bar_skips_rather_than_passing_when_there_is_nothing_to_compare(
    tmp_path, monkeypatch
):
    """A missing baseline must not read as a pass."""
    import src.reference.jobs.icbi_coexpression as mod

    monkeypatch.setattr(mod, "RESULTS_DIR", tmp_path / "results")
    result = mod.check_against_committed(_deltas("Pelka_2021_Cell", guca2a=-0.4, actb=0.0))
    assert result["verdict"] == "SKIPPED"


def test_the_bar_names_pelka_and_the_table_it_checks_against():
    assert VALIDATION["study_id"] == "Pelka_2021_Cell"
    assert VALIDATION["committed_study_id"] == "GSE178341"
    assert "2026-09-04_975cf5c" in VALIDATION["against"]
    assert len(VALIDATION["requirements"]) == 3
