"""Metadata-only H&E-to-molecular substrate gate for HTA11.

The gate deliberately reads portal metadata only.  It does not download images,
open a VCF, pick a molecular endpoint, train an encoder, or calculate a model
metric.  A same-participant join is insufficient: an H&E image is retained
only when its HTAN biospecimen identifier exactly equals the identifier on a
Level-3 Bulk-DNA VCF.
"""

from __future__ import annotations

from collections.abc import Iterable

import pandas as pd

H_AND_E_ASSAY = "H&E"
H_AND_E_LEVEL = "Level 2"
BULK_DNA_ASSAY = "Bulk DNA"
LEVEL_3 = "Level 3"
VCF_FORMAT = "vcf"
PREMALIGNANT_TISSUE_TYPES = frozenset({"Premalignant", "Atypia - hyperplasia"})
NOT_REPORTED = "Not Reported"

FILES_REQUIRED = frozenset(
    {
        "Filename",
        "Biospecimen",
        "Assay",
        "Level",
        "Data Access",
        "File Format",
        "Synapse Id",
        "Data File ID",
    }
)
BIOSPECIMENS_REQUIRED = frozenset(
    {
        "HTAN Biospecimen ID",
        "Participant ID",
        "Tumor Tissue Type",
    }
)
CASES_REQUIRED = frozenset(
    {
        "HTAN Participant ID",
        "Primary Diagnosis",
        "Site of Resection or Biopsy",
    }
)
SYNAPSE_ACCESS_REQUIRED = frozenset({"biospecimen_id", "entity_id", "metadata_status"})


class HEMolecularGateError(ValueError):
    """The portal export is not the data contract this gate was written for."""


def _require_columns(frame: pd.DataFrame, required: Iterable[str], *, name: str) -> None:
    missing = sorted(set(required) - set(frame.columns))
    if missing:
        raise HEMolecularGateError(f"{name} metadata missing required columns {missing}")


def _text(series: pd.Series) -> pd.Series:
    """Make portal nulls explicit without treating them as biological values."""
    return series.fillna("").astype(str).str.strip()


def _unique_index(frame: pd.DataFrame, *, column: str, name: str) -> pd.DataFrame:
    duplicate = _text(frame[column]).duplicated(keep=False)
    if duplicate.any():
        values = sorted(_text(frame.loc[duplicate, column]).unique())
        raise HEMolecularGateError(f"{name} metadata duplicates {column}: {values[:5]}")
    return frame.set_index(column, drop=False)


def candidate_synapse_entities(attrition: pd.DataFrame) -> pd.DataFrame:
    """Return the exact VCF entities whose access status is relevant to the gate."""
    candidate = attrition[attrition["has_exact_level3_vcf"] & attrition["candidate_premalignant"]]
    rows: list[dict[str, str]] = []
    for row in candidate.itertuples(index=False):
        entity_ids = [
            entity_id.strip()
            for entity_id in str(row.molecular_synapse_ids).split(" | ")
            if entity_id.strip()
        ]
        if not entity_ids:
            raise HEMolecularGateError(
                f"{row.biospecimen_id} passed the VCF join without a Synapse ID"
            )
        rows.extend(
            {"biospecimen_id": str(row.biospecimen_id), "entity_id": entity_id}
            for entity_id in entity_ids
        )
    expected = pd.DataFrame(rows)
    if expected.empty:
        raise HEMolecularGateError("no exact premalignant H&E–VCF candidate entities")
    if expected.duplicated().any():
        raise HEMolecularGateError("candidate VCF entity extraction emitted a duplicate row")
    return expected.sort_values(["biospecimen_id", "entity_id"], ignore_index=True)


def _synapse_access_status(
    attrition: pd.DataFrame,
    synapse_access: pd.DataFrame | None,
) -> str:
    """Validate a prior metadata-only probe; never infer a download from it."""
    if synapse_access is None:
        return "not_checked"
    _require_columns(synapse_access, SYNAPSE_ACCESS_REQUIRED, name="Synapse access")
    observed = synapse_access.loc[:, ["biospecimen_id", "entity_id", "metadata_status"]].copy()
    if observed.duplicated(["biospecimen_id", "entity_id"]).any():
        raise HEMolecularGateError("Synapse access artifact duplicates a candidate entity")
    expected = candidate_synapse_entities(attrition)
    expected_pairs = set(map(tuple, expected[["biospecimen_id", "entity_id"]].to_numpy()))
    observed_pairs = set(
        map(tuple, observed[["biospecimen_id", "entity_id"]].astype(str).to_numpy())
    )
    if observed_pairs != expected_pairs:
        raise HEMolecularGateError(
            "Synapse access artifact does not cover exactly the candidate VCF entities"
        )
    return "all_readable" if observed["metadata_status"].eq("readable").all() else "unresolved"


def build_inventory(
    files: pd.DataFrame,
    biospecimens: pd.DataFrame,
    cases: pd.DataFrame,
    synapse_access: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return H&E-file attrition and the one-row, no-model gate summary."""
    _require_columns(files, FILES_REQUIRED, name="Files")
    _require_columns(biospecimens, BIOSPECIMENS_REQUIRED, name="Biospecimen")
    _require_columns(cases, CASES_REQUIRED, name="Case")

    h_and_e = files[
        _text(files["Assay"]).eq(H_AND_E_ASSAY) & _text(files["Level"]).eq(H_AND_E_LEVEL)
    ].copy()
    if h_and_e.empty:
        raise HEMolecularGateError("Files metadata contains no H&E Level-2 records")
    h_and_e["biospecimen_id"] = _text(h_and_e["Biospecimen"])
    if h_and_e["biospecimen_id"].eq("").any():
        raise HEMolecularGateError("H&E record missing Biospecimen")

    vcf = files[
        _text(files["Assay"]).eq(BULK_DNA_ASSAY)
        & _text(files["Level"]).eq(LEVEL_3)
        & _text(files["File Format"]).str.lower().eq(VCF_FORMAT)
    ].copy()
    vcf["biospecimen_id"] = _text(vcf["Biospecimen"])
    vcf_access = vcf.groupby("biospecimen_id", observed=True)["Data Access"].agg(
        lambda values: " | ".join(sorted(set(_text(values))))
    )
    vcf_synapse_ids = vcf.groupby("biospecimen_id", observed=True)["Synapse Id"].agg(
        lambda values: " | ".join(sorted(set(_text(values))))
    )

    bio = biospecimens.copy()
    bio["biospecimen_id"] = _text(bio["HTAN Biospecimen ID"])
    bio = _unique_index(bio, column="biospecimen_id", name="Biospecimen")
    case = cases.copy()
    case["participant_id"] = _text(case["HTAN Participant ID"])
    case = _unique_index(case, column="participant_id", name="Case")

    unknown_biospecimens = sorted(set(h_and_e["biospecimen_id"]) - set(bio.index))
    if unknown_biospecimens:
        raise HEMolecularGateError(
            f"H&E biospecimens absent from Biospecimen metadata: {unknown_biospecimens[:5]}"
        )

    h_and_e["image_access"] = _text(h_and_e["Data Access"])
    h_and_e["has_exact_level3_vcf"] = h_and_e["biospecimen_id"].isin(vcf_access.index)
    h_and_e["molecular_access"] = h_and_e["biospecimen_id"].map(vcf_access).fillna("")
    h_and_e["molecular_synapse_ids"] = h_and_e["biospecimen_id"].map(vcf_synapse_ids).fillna("")
    h_and_e["participant_id"] = h_and_e["biospecimen_id"].map(bio["Participant ID"])
    h_and_e["tumor_tissue_type"] = h_and_e["biospecimen_id"].map(bio["Tumor Tissue Type"])
    h_and_e["primary_diagnosis"] = h_and_e["participant_id"].map(case["Primary Diagnosis"])
    h_and_e["acquisition_site"] = h_and_e["participant_id"].map(case["Site of Resection or Biopsy"])

    if h_and_e[["participant_id", "primary_diagnosis", "acquisition_site"]].isna().any().any():
        missing = h_and_e.loc[
            h_and_e[["participant_id", "primary_diagnosis", "acquisition_site"]].isna().any(axis=1),
            "biospecimen_id",
        ].tolist()
        raise HEMolecularGateError(
            f"candidate metadata missing participant/case rows: {missing[:5]}"
        )

    h_and_e["candidate_premalignant"] = _text(h_and_e["tumor_tissue_type"]).isin(
        PREMALIGNANT_TISSUE_TYPES
    )
    h_and_e["diagnosis_reported"] = ~_text(h_and_e["primary_diagnosis"]).eq(NOT_REPORTED)
    h_and_e["site_reported"] = ~_text(h_and_e["acquisition_site"]).eq(NOT_REPORTED)
    h_and_e["molecular_access_pending_determination"] = _text(
        h_and_e["molecular_access"]
    ).str.contains("Synapse", case=False, regex=False)

    def reason(row: pd.Series) -> str:
        if not row["has_exact_level3_vcf"]:
            return "no_exact_level3_vcf"
        if not row["candidate_premalignant"]:
            return "not_premalignant_biospecimen"
        if not row["diagnosis_reported"]:
            return "candidate_pending_case_diagnosis"
        if not row["site_reported"]:
            return "candidate_pending_acquisition_site"
        return "candidate_pending_endpoint_prespecification"

    h_and_e["gate_status"] = h_and_e.apply(reason, axis=1)
    h_and_e["image_file"] = _text(h_and_e["Filename"])
    h_and_e["image_file_format"] = _text(h_and_e["File Format"])
    h_and_e["image_synapse_id"] = _text(h_and_e["Synapse Id"])
    h_and_e["image_data_file_id"] = _text(h_and_e["Data File ID"])
    h_and_e["image_resolution_metadata_available"] = False
    attrition = (
        h_and_e.loc[
            :,
            [
                "image_file",
                "image_file_format",
                "image_synapse_id",
                "image_data_file_id",
                "image_resolution_metadata_available",
                "biospecimen_id",
                "participant_id",
                "image_access",
                "has_exact_level3_vcf",
                "molecular_access",
                "molecular_synapse_ids",
                "molecular_access_pending_determination",
                "tumor_tissue_type",
                "primary_diagnosis",
                "acquisition_site",
                "candidate_premalignant",
                "diagnosis_reported",
                "site_reported",
                "gate_status",
            ],
        ]
        .sort_values("biospecimen_id", kind="stable")
        .reset_index(drop=True)
    )

    exact = attrition["has_exact_level3_vcf"]
    premalignant = attrition["candidate_premalignant"]
    candidate = exact & premalignant
    images_open = (
        attrition["image_access"].str.contains("open access", case=False, regex=False).all()
    )
    diagnosis_complete = attrition.loc[candidate, "diagnosis_reported"].all()
    site_complete = attrition.loc[candidate, "site_reported"].all()
    molecular_access_rollup = " | ".join(
        sorted(attrition.loc[candidate, "molecular_access"].unique())
    )
    synapse_metadata_status = _synapse_access_status(attrition, synapse_access)
    molecular_access_pending = (
        "Synapse" in molecular_access_rollup.split(" | ")
        and synapse_metadata_status != "all_readable"
    )
    reasons = [
        "image resolution not verified from metadata",
        "molecular endpoint not pre-specified",
    ]
    if molecular_access_pending:
        if synapse_metadata_status == "unresolved":
            reasons.append("molecular Synapse metadata unresolved for one or more candidates")
        else:
            reasons.append("molecular Synapse access pending determination")
    if not diagnosis_complete:
        reasons.append("case diagnosis incomplete")
    if not site_complete:
        reasons.append("case acquisition site incomplete")
    summary = pd.DataFrame(
        [
            {
                "n_h_and_e_files": int(len(attrition)),
                "n_h_and_e_biospecimens": int(attrition["biospecimen_id"].nunique()),
                "n_exact_h_and_e_level3_vcf_biospecimens": int(
                    attrition.loc[exact, "biospecimen_id"].nunique()
                ),
                "n_premalignant_h_and_e_biospecimens": int(
                    attrition.loc[premalignant, "biospecimen_id"].nunique()
                ),
                "n_premalignant_exact_matches": int(
                    attrition.loc[candidate, "biospecimen_id"].nunique()
                ),
                "n_premalignant_no_exact_vcf": int(
                    attrition.loc[premalignant & ~exact, "biospecimen_id"].nunique()
                ),
                "n_nonpremalignant_exact_matches": int(
                    attrition.loc[exact & ~premalignant, "biospecimen_id"].nunique()
                ),
                "n_candidate_participants": int(
                    attrition.loc[candidate, "participant_id"].nunique()
                ),
                "n_reported_candidate_sites": int(
                    attrition.loc[
                        candidate & attrition["site_reported"], "acquisition_site"
                    ].nunique()
                ),
                "images_open_access": bool(images_open),
                "image_resolution_metadata_available": False,
                "biospecimen_exact_crosswalk": bool(candidate.any()),
                "molecular_access_rollup": molecular_access_rollup,
                "molecular_synapse_metadata_status": synapse_metadata_status,
                "molecular_access_pending_determination": molecular_access_pending,
                "case_diagnosis_complete": bool(diagnosis_complete),
                "case_site_complete": bool(site_complete),
                "molecular_endpoint_prespecified": False,
                "verdict": "NOT LICENSED — " + "; ".join(reasons),
                "inference": "none: metadata-only source gate; no image or VCF read and no model",
            }
        ]
    )
    return attrition, summary


def exit_code(summary: pd.DataFrame) -> int:
    """Return the durable-gate exit convention: 5 means not licensed."""
    return 5 if summary.loc[0, "verdict"].startswith("NOT LICENSED") else 0
