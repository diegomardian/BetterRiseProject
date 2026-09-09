"""The WES subtype route begins with a strict, outcome-free lesion manifest."""

from __future__ import annotations

import pandas as pd
import pytest

from src.reference.jobs.wes_subtype_lineage_manifest import validate_specification
from src.reference.wes_subtype import (
    EXPECTED_GENES_PER_PATIENT,
    EXPECTED_LINEAGE_PATIENTS,
    WESSubtypeCrosswalkError,
    avenue_a_lineage_patients,
    build_lineage_manifest,
)


def _decomposition() -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    for patient in range(EXPECTED_LINEAGE_PATIENTS):
        for gene in range(EXPECTED_GENES_PER_PATIENT):
            rows.append(
                {
                    "study_id": "Chen_2021_Cell",
                    "patient_id": f"Chen_2021_Cell.HTA11_{patient}",
                    "gene": f"G{gene}",
                    "granularity_rung": "lineage",
                    "reading": "adenoma",
                }
            )
    return pd.DataFrame(rows)


def _obs() -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    for patient in range(EXPECTED_LINEAGE_PATIENTS):
        base = f"HTA11_{patient}_2000001011"
        rows.append(
            {
                "study_id": "Chen_2021_Cell",
                "patient_id": f"Chen_2021_Cell.HTA11_{patient}",
                "sample_id": f"VUMC_HTAN_validation.{base}",
                "sample_type": "polyp",
                "dataset": "VUMC_HTAN_validation",
            }
        )
    return pd.DataFrame(rows)


def test_manifest_declares_independent_genotype_claim_before_inputs_are_read():
    assert validate_specification() == ()


def test_manifest_is_one_row_per_sc_rna_polyp_biospecimen_not_participant_proxy():
    obs = _obs()
    obs.loc[len(obs)] = {
        "study_id": "Chen_2021_Cell",
        "patient_id": "Chen_2021_Cell.HTA11_0",
        "sample_id": "VUMC_HTAN_validation.HTA11_0_2000001021",
        "sample_type": "polyp",
        "dataset": "VUMC_HTAN_validation",
    }
    manifest = build_lineage_manifest(_decomposition(), obs)
    assert len(manifest) == EXPECTED_LINEAGE_PATIENTS + 1
    patient = manifest[manifest["patient_id"].eq("Chen_2021_Cell.HTA11_0")]
    assert patient["scRNA_biospecimen_id"].tolist() == [
        "HTA11_0_2000001011",
        "HTA11_0_2000001021",
    ]
    assert patient["multiple_polyp_scRNA_biospecimens"].all()
    assert not manifest["specimen_exact"].any()
    assert not manifest["suffix_match_used"].any()


def test_manifest_refuses_changed_lineage_universe_instead_of_silently_replacing_it():
    changed = _decomposition().iloc[:-1]
    with pytest.raises(WESSubtypeCrosswalkError, match="exactly 6 genes"):
        avenue_a_lineage_patients(changed)


def test_manifest_refuses_a_sample_id_without_an_explicit_htan_biospecimen_token():
    obs = _obs()
    obs.loc[0, "sample_id"] = "unparseable"
    with pytest.raises(WESSubtypeCrosswalkError, match="dataset separator"):
        build_lineage_manifest(_decomposition(), obs)
