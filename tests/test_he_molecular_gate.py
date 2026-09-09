from __future__ import annotations

import pandas as pd
import pytest

from src.reference.he_molecular_gate import (
    HEMolecularGateError,
    build_inventory,
    exit_code,
)
from src.reference.jobs.he_molecular_image_headers import candidate_images
from src.reference.jobs.he_molecular_synapse_access import candidate_entities


def _file(
    filename: str,
    biospecimen: str,
    assay: str,
    level: str,
    access: str,
    file_format: str,
) -> dict[str, str]:
    return {
        "Filename": filename,
        "Biospecimen": biospecimen,
        "Assay": assay,
        "Level": level,
        "Data Access": access,
        "File Format": file_format,
        "Synapse Id": f"syn-{filename}",
        "Data File ID": f"file-{filename}",
    }


def _files(*, include_level_1: bool = False) -> pd.DataFrame:
    rows = [
        _file("A.tif", "B-A", "H&E", "Level 2", "CRDC-GC/SB-CGC (open access)", "tif"),
        _file("B.tif", "B-B", "H&E", "Level 2", "CRDC-GC/SB-CGC (open access)", "tif"),
        _file("C.tif", "B-C", "H&E", "Level 2", "CRDC-GC/SB-CGC (open access)", "tif"),
        _file("A.vcf", "B-A", "Bulk DNA", "Level 3", "Synapse", "vcf"),
        _file("B.vcf", "B-B", "Bulk DNA", "Level 3", "Synapse", "vcf"),
        _file("loose.vcf", "B-C-other", "Bulk DNA", "Level 3", "Synapse", "vcf"),
    ]
    if include_level_1:
        rows.append(_file("thumbnail.tif", "B-A", "H&E", "Level 1", "open access", "tif"))
    return pd.DataFrame(rows)


def _biospecimen(
    biospecimen_id: str,
    participant_id: str,
    tissue_type: str,
) -> dict[str, str]:
    return {
        "HTAN Biospecimen ID": biospecimen_id,
        "Participant ID": participant_id,
        "Tumor Tissue Type": tissue_type,
    }


def _biospecimens() -> pd.DataFrame:
    return pd.DataFrame(
        [
            _biospecimen("B-A", "P-A", "Premalignant"),
            _biospecimen("B-B", "P-B", "Primary"),
            _biospecimen("B-C", "P-C", "Premalignant"),
        ]
    )


def _case(participant_id: str, diagnosis: str, site: str) -> dict[str, str]:
    return {
        "HTAN Participant ID": participant_id,
        "Primary Diagnosis": diagnosis,
        "Site of Resection or Biopsy": site,
    }


def _cases() -> pd.DataFrame:
    return pd.DataFrame(
        [
            _case("P-A", "Adenoma", "VUMC"),
            _case("P-B", "Adenocarcinoma", "VUMC"),
            _case("P-C", "Adenoma", "VUMC"),
        ]
    )


def test_exact_biospecimen_join_never_promotes_a_prefix_match():
    attrition, summary = build_inventory(_files(), _biospecimens(), _cases())
    by_id = attrition.set_index("biospecimen_id")
    assert by_id.loc["B-A", "gate_status"] == "candidate_pending_endpoint_prespecification"
    assert by_id.loc["B-B", "gate_status"] == "not_premalignant_biospecimen"
    assert by_id.loc["B-C", "gate_status"] == "no_exact_level3_vcf"
    assert summary.loc[0, "n_exact_h_and_e_level3_vcf_biospecimens"] == 2
    assert summary.loc[0, "n_premalignant_h_and_e_biospecimens"] == 2
    assert summary.loc[0, "n_premalignant_exact_matches"] == 1
    assert summary.loc[0, "n_premalignant_no_exact_vcf"] == 1
    assert summary.loc[0, "n_candidate_participants"] == 1


def test_access_and_identifiers_are_reported_without_claiming_resolution():
    attrition, summary = build_inventory(_files(), _biospecimens(), _cases())
    row = attrition.set_index("biospecimen_id").loc["B-A"]
    assert row["image_synapse_id"] == "syn-A.tif"
    assert row["image_data_file_id"] == "file-A.tif"
    assert row["image_file_format"] == "tif"
    assert not row["image_resolution_metadata_available"]
    assert row["molecular_access"] == "Synapse"
    assert row["molecular_synapse_ids"] == "syn-A.vcf"
    assert row["molecular_access_pending_determination"]
    assert summary.loc[0, "images_open_access"]
    assert summary.loc[0, "molecular_access_rollup"] == "Synapse"
    assert summary.loc[0, "n_reported_candidate_sites"] == 1


def test_molecular_access_pending_is_derived_from_candidate_access():
    files = _files()
    files.loc[files["Assay"] == "Bulk DNA", "Data Access"] = "open access"
    attrition, summary = build_inventory(files, _biospecimens(), _cases())
    row = attrition.set_index("biospecimen_id").loc["B-A"]
    assert not row["molecular_access_pending_determination"]
    assert not summary.loc[0, "molecular_access_pending_determination"]
    assert "Synapse access pending" not in summary.loc[0, "verdict"]


def test_authenticated_access_artifact_clears_only_the_synapse_pending_reason():
    access = pd.DataFrame(
        [{"biospecimen_id": "B-A", "entity_id": "syn-A.vcf", "metadata_status": "readable"}]
    )
    _, summary = build_inventory(_files(), _biospecimens(), _cases(), access)
    assert summary.loc[0, "molecular_synapse_metadata_status"] == "all_readable"
    assert not summary.loc[0, "molecular_access_pending_determination"]
    assert "Synapse access pending" not in summary.loc[0, "verdict"]
    assert "image resolution not verified" in summary.loc[0, "verdict"]


def test_access_artifact_must_cover_exactly_the_current_candidate_entities():
    access = pd.DataFrame(
        [{"biospecimen_id": "B-A", "entity_id": "syn-wrong", "metadata_status": "readable"}]
    )
    with pytest.raises(HEMolecularGateError, match="does not cover exactly"):
        build_inventory(_files(), _biospecimens(), _cases(), access)


def test_synapse_probe_uses_only_exact_premalignant_candidate_ids():
    attrition, _ = build_inventory(_files(), _biospecimens(), _cases())
    entities = candidate_entities(attrition)
    assert entities.to_dict("records") == [
        {
            "biospecimen_id": "B-A",
            "participant_id": "P-A",
            "entity_id": "syn-A.vcf",
        }
    ]


def test_header_probe_uses_only_exact_premalignant_image_ids():
    attrition, _ = build_inventory(_files(), _biospecimens(), _cases())
    images = candidate_images(attrition)
    assert images.to_dict("records") == [
        {"biospecimen_id": "B-A", "participant_id": "P-A", "entity_id": "syn-A.tif"}
    ]


def test_level_1_h_and_e_is_not_silently_counted_as_level_2():
    attrition, summary = build_inventory(_files(include_level_1=True), _biospecimens(), _cases())
    assert len(attrition) == 3
    assert summary.loc[0, "n_h_and_e_files"] == 3


def test_missing_required_column_fails_instead_of_guessing_the_join():
    files = _files().drop(columns="Biospecimen")
    with pytest.raises(HEMolecularGateError, match="Biospecimen"):
        build_inventory(files, _biospecimens(), _cases())


def test_no_level_2_h_and_e_and_unknown_biospecimen_both_refuse():
    files = _files()
    files.loc[files["Assay"] == "H&E", "Level"] = "Level 1"
    with pytest.raises(HEMolecularGateError, match="Level-2"):
        build_inventory(files, _biospecimens(), _cases())

    biospecimens = _biospecimens().query("`HTAN Biospecimen ID` != 'B-C'")
    with pytest.raises(HEMolecularGateError, match="absent from Biospecimen"):
        build_inventory(_files(), biospecimens, _cases())

    files = _files()
    files.loc[files["Assay"] == "H&E", "Biospecimen"] = ""
    with pytest.raises(HEMolecularGateError, match="missing Biospecimen"):
        build_inventory(files, _biospecimens(), _cases())


def test_duplicate_biospecimen_metadata_refuses_a_nondeterministic_join():
    biospecimens = pd.concat([_biospecimens(), _biospecimens().iloc[[0]]])
    with pytest.raises(HEMolecularGateError, match="duplicates"):
        build_inventory(_files(), biospecimens, _cases())

    cases = pd.concat([_cases(), _cases().iloc[[0]]])
    with pytest.raises(HEMolecularGateError, match="Case metadata duplicates"):
        build_inventory(_files(), _biospecimens(), cases)


def test_not_reported_case_values_do_not_become_single_site_evidence():
    cases = _cases()
    cases.loc[
        cases["HTAN Participant ID"] == "P-A", ["Primary Diagnosis", "Site of Resection or Biopsy"]
    ] = "Not Reported"
    attrition, summary = build_inventory(_files(), _biospecimens(), cases)
    row = attrition.set_index("biospecimen_id").loc["B-A"]
    assert row["gate_status"] == "candidate_pending_case_diagnosis"
    assert not summary.loc[0, "case_diagnosis_complete"]
    assert not summary.loc[0, "case_site_complete"]
    assert summary.loc[0, "n_reported_candidate_sites"] == 0
    assert exit_code(summary) == 5
