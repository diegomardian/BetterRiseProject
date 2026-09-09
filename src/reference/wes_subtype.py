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
