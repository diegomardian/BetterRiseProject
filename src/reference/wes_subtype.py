"""Outcome-blind source-universe helpers for the Chen WES subtype gate.

This module builds the *transcript-side* half of a prospective WES crosswalk.
It makes no molecular assignment: a participant match, a shared numeric suffix,
or a VCF listed for a participant does not prove that a VCF came from the polyp
used in the avenue-A decomposition.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

CHEN_STUDY_ID = "Chen_2021_Cell"
LINEAGE_RUNG = "lineage"
ADENOMA_READING = "adenoma"
EXPECTED_LINEAGE_PATIENTS = 44
EXPECTED_GENES_PER_PATIENT = 6

DECOMPOSITION_REQUIRED = frozenset(
    {"study_id", "patient_id", "gene", "granularity_rung", "reading"}
)
OBS_REQUIRED = frozenset(
    {"study_id", "patient_id", "sample_id", "sample_type", "dataset"}
)
MANIFEST_REQUIRED = frozenset(
    {"patient_id", "scRNA_biospecimen_id", "specimen_exact", "crosswalk_status"}
)
PROVENANCE_REQUIRED = frozenset(
    {
        "vcf_data_file_id",
        "vcf_filename",
        "vcf_entity_id",
        "vcf_parent_data_file_id",
        "htan_participant_id",
        "vcf_assayed_biospecimen_id",
        "vcf_originating_biospecimen_id",
        "biospecimen_path",
    }
)


class WESSubtypeCrosswalkError(ValueError):
    """The fixed transcript-side universe cannot support an honest crosswalk."""


def _require_columns(frame: pd.DataFrame, required: frozenset[str], *, name: str) -> None:
    missing = sorted(required - set(frame.columns))
    if missing:
        raise WESSubtypeCrosswalkError(f"{name} is missing required columns {missing}")


def _text(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def _sc_rna_biospecimen_id(sample_id: str) -> str:
    """Extract the deposited biospecimen token; never interpret its suffix."""
    value = str(sample_id).strip()
    if "." not in value:
        raise WESSubtypeCrosswalkError(
            f"scRNA sample ID lacks its dataset separator: {value!r}"
        )
    token = value.rsplit(".", maxsplit=1)[1]
    if not token.startswith("HTA11_"):
        raise WESSubtypeCrosswalkError(
            f"scRNA sample ID does not end in an HTA11 biospecimen token: {value!r}"
        )
    return token


def avenue_a_lineage_patients(decomposition: pd.DataFrame) -> tuple[str, ...]:
    """Return exactly the committed 44-patient lineage universe, not a proxy."""
    _require_columns(decomposition, DECOMPOSITION_REQUIRED, name="decomposition")
    rows = decomposition[
        _text(decomposition["study_id"]).eq(CHEN_STUDY_ID)
        & _text(decomposition["granularity_rung"]).eq(LINEAGE_RUNG)
        & _text(decomposition["reading"]).eq(ADENOMA_READING)
    ].copy()
    if rows.empty:
        raise WESSubtypeCrosswalkError("committed decomposition has no Chen lineage adenoma rows")
    per_patient = rows.groupby("patient_id", observed=True)["gene"].nunique()
    bad = per_patient[per_patient.ne(EXPECTED_GENES_PER_PATIENT)]
    if not bad.empty:
        raise WESSubtypeCrosswalkError(
            "committed lineage universe must have exactly "
            f"{EXPECTED_GENES_PER_PATIENT} genes per patient; got {bad.to_dict()}"
        )
    patients = tuple(sorted(_text(per_patient.index.to_series())))
    if len(patients) != EXPECTED_LINEAGE_PATIENTS:
        raise WESSubtypeCrosswalkError(
            f"committed lineage universe has {len(patients)} patients, expected "
            f"{EXPECTED_LINEAGE_PATIENTS}"
        )
    prefix = f"{CHEN_STUDY_ID}."
    malformed = [patient for patient in patients if not patient.startswith(prefix)]
    if malformed:
        raise WESSubtypeCrosswalkError(
            f"Chen patient IDs lack expected prefix {prefix!r}: {malformed[:5]}"
        )
    return patients


def build_lineage_manifest(decomposition: pd.DataFrame, obs: pd.DataFrame) -> pd.DataFrame:
    """Emit one row per deposited polyp scRNA biospecimen in the 44-patient set.

    ``scRNA_biospecimen_id`` is a token to be compared for *exact equality* to
    a provenance-exported assayed/originating biospecimen.  It is not parsed or
    matched by suffix.  Multiple polyp samples remain multiple candidate
    lesions; they are not collapsed to a participant-level WES assignment.
    """
    patients = avenue_a_lineage_patients(decomposition)
    _require_columns(obs, OBS_REQUIRED, name="ICBI obs")
    selected = obs[
        _text(obs["study_id"]).eq(CHEN_STUDY_ID)
        & _text(obs["patient_id"]).isin(patients)
        & _text(obs["sample_type"]).eq("polyp")
    ].loc[:, ["patient_id", "dataset", "sample_id", "sample_type"]].copy()
    if selected.empty:
        raise WESSubtypeCrosswalkError("ICBI obs has no polyp samples for the lineage universe")
    selected = selected.drop_duplicates()
    if selected.duplicated(["patient_id", "sample_id"]).any():
        raise WESSubtypeCrosswalkError("ICBI obs duplicates a Chen patient/sample pair")
    selected["scRNA_biospecimen_id"] = selected["sample_id"].map(_sc_rna_biospecimen_id)
    if selected["scRNA_biospecimen_id"].duplicated().any():
        duplicate = selected.loc[
            selected["scRNA_biospecimen_id"].duplicated(keep=False),
            "scRNA_biospecimen_id",
        ].tolist()
        raise WESSubtypeCrosswalkError(
            f"one scRNA biospecimen token maps to multiple rows: {duplicate[:5]}"
        )
    count = selected.groupby("patient_id", observed=True)["sample_id"].transform("size")
    selected["n_polyp_scRNA_biospecimens_for_patient"] = count.astype(int)
    selected["multiple_polyp_scRNA_biospecimens"] = count.gt(1)
    selected["crosswalk_rule"] = (
        "exact equality to a VCF assayed or originating biospecimen ID only"
    )
    selected["suffix_match_used"] = False
    selected["specimen_exact"] = False
    selected["crosswalk_status"] = "awaiting_provenance_export"
    return selected.sort_values(["patient_id", "sample_id"], kind="stable").reset_index(drop=True)


def read_decomposition(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"adenoma decomposition not found: {path}")
    return pd.read_parquet(path)


def read_obs(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"ICBI obs cache not found: {path}")
    return pd.read_parquet(path, columns=sorted(OBS_REQUIRED))


def _htan_participant_id(patient_id: str) -> str:
    prefix = f"{CHEN_STUDY_ID}."
    value = str(patient_id).strip()
    if not value.startswith(prefix):
        raise WESSubtypeCrosswalkError(
            f"Chen patient ID lacks expected prefix {prefix!r}: {value!r}"
        )
    return value.removeprefix(prefix)


def _join(values: pd.Series) -> str:
    return " | ".join(sorted({str(value).strip() for value in values if str(value).strip()}))


def build_provenance_crosswalk(
    manifest: pd.DataFrame, provenance: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Match only full scRNA and WES biospecimen IDs through HTAN provenance.

    The caller must provide the normalized BigQuery export documented in
    ``docs/wes_subtype_plan.md``.  A participant hit without equality to either
    documented WES biospecimen field is retained as a failure, not promoted to
    an exact lesion match.  A pair of exact VCF records is likewise ambiguous:
    this source gate does not choose one after seeing filenames.
    """
    _require_columns(manifest, MANIFEST_REQUIRED, name="lineage manifest")
    _require_columns(provenance, PROVENANCE_REQUIRED, name="HTAN provenance export")
    if manifest.empty:
        raise WESSubtypeCrosswalkError("lineage manifest is empty")
    if manifest["scRNA_biospecimen_id"].isna().any() or _text(
        manifest["scRNA_biospecimen_id"]
    ).eq("").any():
        raise WESSubtypeCrosswalkError("lineage manifest contains a missing scRNA biospecimen ID")
    if manifest["scRNA_biospecimen_id"].duplicated().any():
        raise WESSubtypeCrosswalkError("lineage manifest duplicates an scRNA biospecimen ID")
    if provenance.empty:
        raise WESSubtypeCrosswalkError("HTAN provenance export is empty")
    prov = provenance.loc[:, sorted(PROVENANCE_REQUIRED)].copy()
    for column in prov.columns:
        prov[column] = _text(prov[column])
    if prov["vcf_data_file_id"].eq("").any() or prov["htan_participant_id"].eq("").any():
        raise WESSubtypeCrosswalkError(
            "HTAN provenance export has a VCF row missing data-file or participant ID"
        )
    if prov["vcf_data_file_id"].duplicated().any():
        raise WESSubtypeCrosswalkError("HTAN provenance export duplicates a VCF data-file ID")
    if (~prov["vcf_filename"].str.lower().str.contains(r"\.vcf(?:\.gz)?$", regex=True)).any():
        bad = prov.loc[
            ~prov["vcf_filename"].str.lower().str.contains(r"\.vcf(?:\.gz)?$", regex=True),
            "vcf_filename",
        ].tolist()
        raise WESSubtypeCrosswalkError(
            f"HTAN provenance export includes a non-VCF filename: {bad[:5]}"
        )

    rows: list[dict[str, object]] = []
    for lesion in manifest.itertuples(index=False):
        lesion_row = lesion._asdict()
        patient = _htan_participant_id(lesion_row["patient_id"])
        biospecimen = str(lesion_row["scRNA_biospecimen_id"]).strip()
        participant_vcf = prov.loc[prov["htan_participant_id"].eq(patient)].copy()
        exact = participant_vcf[
            participant_vcf["vcf_assayed_biospecimen_id"].eq(biospecimen)
            | participant_vcf["vcf_originating_biospecimen_id"].eq(biospecimen)
        ].copy()
        missing_fields = participant_vcf[
            participant_vcf["vcf_entity_id"].eq("")
            | (
                participant_vcf["vcf_assayed_biospecimen_id"].eq("")
                & participant_vcf["vcf_originating_biospecimen_id"].eq("")
            )
        ]
        if participant_vcf.empty:
            status = "no_participant_vcf"
        elif exact.empty and not missing_fields.empty:
            status = "missing_provenance_field"
        elif exact.empty:
            status = "participant_vcf_not_specimen_exact"
        elif len(exact) > 1:
            status = "ambiguous_provenance"
        elif exact["vcf_entity_id"].eq("").any():
            status = "missing_provenance_field"
        else:
            status = "specimen_exact_unique"
        rows.append(
            {
                **lesion_row,
                "htan_participant_id": patient,
                "n_participant_vcf_records": int(len(participant_vcf)),
                "n_exact_vcf_records": int(len(exact)),
                "matched_vcf_data_file_ids": _join(exact["vcf_data_file_id"]),
                "matched_vcf_entity_ids": _join(exact["vcf_entity_id"]),
                "matched_vcf_assayed_biospecimen_ids": _join(
                    exact["vcf_assayed_biospecimen_id"]
                ),
                "matched_vcf_originating_biospecimen_ids": _join(
                    exact["vcf_originating_biospecimen_id"]
                ),
                "specimen_exact": status == "specimen_exact_unique",
                "crosswalk_status": status,
            }
        )
    table = pd.DataFrame(rows).sort_values(["patient_id", "sample_id"], kind="stable")
    exact = table["specimen_exact"]
    summary = pd.DataFrame(
        [
            {
                "n_scRNA_polyp_biospecimens": int(len(table)),
                "n_lineage_patients": int(table["patient_id"].nunique()),
                "n_participant_vcf_biospecimens": int(
                    table["n_participant_vcf_records"].gt(0).sum()
                ),
                "n_specimen_exact_unique_biospecimens": int(exact.sum()),
                "n_specimen_exact_unique_patients": int(table.loc[exact, "patient_id"].nunique()),
                "n_ambiguous_provenance": int(
                    table["crosswalk_status"].eq("ambiguous_provenance").sum()
                ),
                "n_participant_only": int(
                    table["crosswalk_status"].eq("participant_vcf_not_specimen_exact").sum()
                ),
                "n_missing_provenance_field": int(
                    table["crosswalk_status"].eq("missing_provenance_field").sum()
                ),
                "n_no_participant_vcf": int(
                    table["crosswalk_status"].eq("no_participant_vcf").sum()
                ),
                "suffix_match_used": False,
                "vcf_content_read": False,
                "verdict": (
                    "SPECIMEN-EXACT WES CANDIDATES — access/header gate next"
                    if exact.any()
                    else "NO SUBSTRATE AT THIS RESOLUTION — no unique specimen-exact WES link"
                ),
            }
        ]
    )
    return table.reset_index(drop=True), summary


def read_provenance_export(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"HTAN provenance CSV not found: {path}")
    return pd.read_csv(path, dtype=str, keep_default_na=False)
