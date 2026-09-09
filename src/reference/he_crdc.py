"""Small, side-effect-free helpers for a bounded CRDC image-header probe."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd


class CrcdManifestError(ValueError):
    """A CRDC manifest cannot support the requested bounded probe."""


class CrcdCrosswalkError(ValueError):
    """A CRDC crosswalk is not the fixed, specimen-exact candidate set."""


_CROSSWALK_COLUMNS = (
    "biospecimen_id",
    "HTAN_Data_File_ID",
    "HTAN_Assayed_Biospecimen_ID",
    "HTAN_Originating_Biospecimen_ID",
    "drs_uri",
)
N_EXACT_HE_VCF_CANDIDATES = 18


def single_object_id(manifest_path: Path) -> str:
    """Read exactly one CRDC object identifier from a portal-generated manifest."""
    if not manifest_path.is_file():
        raise FileNotFoundError(f"CRDC manifest not found: {manifest_path}")
    try:
        payload: Any = json.loads(manifest_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise CrcdManifestError(f"CRDC manifest is not valid JSON: {manifest_path}") from exc
    if not isinstance(payload, list) or len(payload) != 1 or not isinstance(payload[0], dict):
        raise CrcdManifestError("bounded header probe requires a manifest with exactly one object")
    object_id = payload[0].get("object_id")
    if not isinstance(object_id, str) or not object_id.strip():
        raise CrcdManifestError("CRDC manifest object lacks a non-empty object_id")
    return object_id


def read_exact_he_crosswalk(path: Path) -> pd.DataFrame:
    """Read the fixed 18-row, exact H&E--VCF--CRDC BigQuery crosswalk.

    The source query must retain the two provenance columns. Their equality to
    ``biospecimen_id`` prevents an image that merely shares a participant with
    a molecular file from being inspected as an exact match.
    """
    if not path.is_file():
        raise FileNotFoundError(f"CRDC crosswalk CSV not found: {path}")
    table = pd.read_csv(path, dtype="string")
    missing = sorted(set(_CROSSWALK_COLUMNS) - set(table.columns))
    if missing:
        raise CrcdCrosswalkError(f"CRDC crosswalk lacks required columns: {missing}")
    selected = table.loc[:, _CROSSWALK_COLUMNS].copy()
    if len(selected) != N_EXACT_HE_VCF_CANDIDATES:
        raise CrcdCrosswalkError(
            "CRDC crosswalk must contain exactly "
            f"{N_EXACT_HE_VCF_CANDIDATES} rows, found {len(selected)}"
        )
    if selected.isna().any().any() or selected.eq("").any().any():
        raise CrcdCrosswalkError("CRDC crosswalk has a missing required value")
    if selected["biospecimen_id"].duplicated().any():
        raise CrcdCrosswalkError("CRDC crosswalk has duplicate biospecimen IDs")
    if selected["HTAN_Data_File_ID"].duplicated().any():
        raise CrcdCrosswalkError("CRDC crosswalk has duplicate H&E data-file IDs")
    if selected["drs_uri"].duplicated().any():
        raise CrcdCrosswalkError("CRDC crosswalk has duplicate CRDC object IDs")
    exact = selected["biospecimen_id"].eq(selected["HTAN_Assayed_Biospecimen_ID"]) & selected[
        "biospecimen_id"
    ].eq(selected["HTAN_Originating_Biospecimen_ID"])
    if not exact.all():
        raise CrcdCrosswalkError("CRDC crosswalk contains a non-exact biospecimen match")
    if not selected["drs_uri"].str.startswith("dg.").all():
        raise CrcdCrosswalkError("CRDC crosswalk has a non-DRS CRDC object ID")
    return selected.sort_values("biospecimen_id", kind="stable").reset_index(drop=True)


def signed_url_from_response(response: object) -> str:
    """Extract one HTTPS URL from Gen3's documented signed-URL response."""
    if isinstance(response, str):
        candidates = [response]
    elif isinstance(response, dict):
        candidates = [value for value in response.values() if isinstance(value, str)]
    else:
        candidates = []
    urls = [value for value in candidates if value.startswith("https://")]
    if len(urls) != 1:
        keys = sorted(response) if isinstance(response, dict) else []
        raise CrcdManifestError(
            "CRDC signed-URL response did not contain exactly one HTTPS URL "
            f"(response keys: {keys})"
        )
    return urls[0]
