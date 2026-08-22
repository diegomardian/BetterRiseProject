"""Clinical annotation from the GDC. W3.4's inputs; W3.5 will extend this.

W3.4 has to test technical factors for confounding against stage, MMR/MSI and
tumour site, so it needs those three before W3.5 curates the full survival
table. This module pulls exactly those and nothing else. **The authoritative
survival table is W3.5's job and comes from TCGA-CDR**, not from here — when it
lands, this should become its site/MSI supplier rather than a second source of
truth for anything they both carry.

THREE THINGS THE GDC WILL SILENTLY GET WRONG IF YOU LET IT
----------------------------------------------------------
1. **Cases have more than one diagnosis.** 182 of 633 COAD/READ cases do, and
   the extras include prostate, liver, breast and skin — second primaries and
   metastases. Taking ``diagnoses[0]`` assigns some patients a breast-cancer
   site. :func:`select_colorectal_diagnosis` picks by organ, explicitly.
2. **The demographic sex field is ``sex_at_birth``.** There is no
   ``demographic.gender`` in the schema; asking for it returns nothing at all
   rather than erroring, so the column silently arrives empty.
3. **MSI calls disagree with themselves.** MSI status is file-level, derived per
   aligned-reads file, and 17 cases carry both an MSI and an MSS call across
   their WGS and WXS files. Those are marked ``conflicting``, never resolved by
   a coin flip.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

GDC_API = "https://api.gdc.cancer.gov"
PROJECTS = ("TCGA-COAD", "TCGA-READ")

#: Sites that count as colorectal. Anything else on a case is a second primary
#: or a metastasis and must not supply the site or the stage.
COLORECTAL_SITES: dict[str, str] = {
    # Right / proximal — cecum through transverse colon.
    "Cecum": "right_colon",
    "Ascending colon": "right_colon",
    "Hepatic flexure of colon": "right_colon",
    "Transverse colon": "right_colon",
    # Left / distal — splenic flexure through sigmoid.
    "Splenic flexure of colon": "left_colon",
    "Descending colon": "left_colon",
    "Sigmoid colon": "left_colon",
    # Rectum.
    "Rectosigmoid junction": "rectum",
    "Rectum, NOS": "rectum",
    # Colon with no subsite recorded. Its own category — 110 cases, too many to
    # drop and impossible to assign without inventing the answer.
    "Colon, NOS": "colon_unspecified",
}

#: The split point is the splenic flexure, which matches the embryological
#: midgut/hindgut boundary, so transverse colon counts as right. Stated because
#: the other convention exists and changes the groups.
SITE_GROUPS = ("right_colon", "left_colon", "rectum", "colon_unspecified")

STAGE_ORDER = ("I", "II", "III", "IV")


class ClinicalError(RuntimeError):
    """Clinical annotation could not be assembled."""


# ---------------------------------------------------------------------------
# Fetch
# ---------------------------------------------------------------------------


def _post(endpoint: str, payload: dict[str, Any]) -> list[dict[str, Any]]:
    import requests

    response = requests.post(f"{GDC_API}/{endpoint}", json=payload, timeout=300)
    response.raise_for_status()
    data = response.json()["data"]
    hits = data["hits"]
    total = data["pagination"]["total"]
    if len(hits) < total:
        raise ClinicalError(
            f"{endpoint}: got {len(hits)} of {total} records. Raise the page size "
            f"rather than analysing a truncated cohort."
        )
    return hits


def fetch_cases(
    projects: tuple[str, ...] = PROJECTS, *, cache: str | Path | None = None
) -> list[dict[str, Any]]:
    """Case-level clinical records. Cached so the analysis is offline-repeatable."""
    if cache is not None and Path(cache).exists():
        return json.loads(Path(cache).read_text(encoding="utf-8"))
    hits = _post(
        "cases",
        {
            "filters": {
                "op": "in",
                "content": {"field": "project.project_id", "value": list(projects)},
            },
            "fields": ",".join(
                [
                    "submitter_id",
                    "project.project_id",
                    "demographic.sex_at_birth",
                    "demographic.age_at_index",
                    "diagnoses.ajcc_pathologic_stage",
                    "diagnoses.tissue_or_organ_of_origin",
                    "diagnoses.primary_diagnosis",
                ]
            ),
            "size": 5000,
            "format": "JSON",
        },
    )
    if cache is not None:
        Path(cache).parent.mkdir(parents=True, exist_ok=True)
        Path(cache).write_text(json.dumps(hits), encoding="utf-8")
    return hits


def fetch_msi(
    projects: tuple[str, ...] = PROJECTS, *, cache: str | Path | None = None
) -> list[dict[str, Any]]:
    """File-level MSI calls. One case can have several, and they can disagree."""
    if cache is not None and Path(cache).exists():
        return json.loads(Path(cache).read_text(encoding="utf-8"))
    hits = _post(
        "files",
        {
            "filters": {
                "op": "and",
                "content": [
                    {
                        "op": "in",
                        "content": {
                            "field": "cases.project.project_id",
                            "value": list(projects),
                        },
                    },
                    {"op": "not", "content": {"field": "msi_status", "value": ["MISSING"]}},
                ],
            },
            "fields": "msi_status,msi_score,cases.submitter_id,experimental_strategy",
            "size": 5000,
            "format": "JSON",
        },
    )
    if cache is not None:
        Path(cache).parent.mkdir(parents=True, exist_ok=True)
        Path(cache).write_text(json.dumps(hits), encoding="utf-8")
    return hits


# ---------------------------------------------------------------------------
# Harmonisation
# ---------------------------------------------------------------------------


def harmonise_stage(raw: str | None) -> str | None:
    """``Stage IIIB`` -> ``III``. Sub-stage letters dropped, editions merged.

    TCGA mixes AJCC editions and the sub-stage letters are not comparable across
    them, so collapsing to the four major stages is what makes stage usable as a
    covariate at this sample size. Returns None for missing or unstageable —
    never a string like "Unknown", which would sort into the middle of an
    ordered factor.
    """
    if not raw or not isinstance(raw, str):
        return None
    text = raw.strip()
    if not text.lower().startswith("stage"):
        return None
    body = text[len("stage") :].strip().upper()
    for stage in ("IV", "III", "II", "I"):
        if body.startswith(stage):
            remainder = body[len(stage) :]
            if remainder in ("", "A", "B", "C"):
                return stage
    return None


def select_colorectal_diagnosis(diagnoses: list[dict[str, Any]]) -> dict[str, Any] | None:
    """Pick the colorectal diagnosis from a case that may have several.

    Preference: a colorectal site that also carries a stage, then any colorectal
    site, then nothing. A case whose only diagnoses are non-colorectal returns
    None and is reported as unknown rather than being handed another organ's
    stage.
    """
    colorectal = [
        d for d in (diagnoses or []) if d.get("tissue_or_organ_of_origin") in COLORECTAL_SITES
    ]
    if not colorectal:
        return None
    staged = [d for d in colorectal if harmonise_stage(d.get("ajcc_pathologic_stage"))]
    return staged[0] if staged else colorectal[0]


def msi_by_patient(msi_hits: list[dict[str, Any]]) -> pd.DataFrame:
    """Collapse file-level MSI calls to one row per patient.

    Disagreement between a case's own files is recorded as ``conflicting``, not
    resolved. Picking a winner would put a fabricated label into the project's
    single pre-registered subgroup variable.
    """
    per_patient: dict[str, set[str]] = {}
    for hit in msi_hits:
        status = hit.get("msi_status")
        if not status:
            continue
        for case in hit.get("cases") or []:
            per_patient.setdefault(case["submitter_id"], set()).add(status)

    rows = [
        {
            "patient_id": patient,
            "msi_status": next(iter(statuses)) if len(statuses) == 1 else "conflicting",
            "n_distinct_calls": len(statuses),
        }
        for patient, statuses in per_patient.items()
    ]
    if not rows:
        return pd.DataFrame(columns=["patient_id", "msi_status", "n_distinct_calls"])
    return pd.DataFrame(rows).sort_values("patient_id").reset_index(drop=True)


def build_clinical_table(
    case_hits: list[dict[str, Any]], msi_hits: list[dict[str, Any]]
) -> pd.DataFrame:
    """One row per patient: project, stage, site, MSI, sex, age."""
    rows = []
    for hit in case_hits:
        diagnosis = select_colorectal_diagnosis(hit.get("diagnoses") or [])
        demographic = hit.get("demographic") or {}
        site_raw = diagnosis.get("tissue_or_organ_of_origin") if diagnosis else None
        rows.append(
            {
                "patient_id": hit["submitter_id"],
                "project": (hit.get("project") or {}).get("project_id"),
                "stage": harmonise_stage(diagnosis.get("ajcc_pathologic_stage"))
                if diagnosis
                else None,
                "site_raw": site_raw,
                "site": COLORECTAL_SITES.get(site_raw) if site_raw else None,
                "sex": demographic.get("sex_at_birth"),
                "age": demographic.get("age_at_index"),
                "n_diagnoses": len(hit.get("diagnoses") or []),
            }
        )
    clinical = pd.DataFrame(rows)
    merged = clinical.merge(msi_by_patient(msi_hits), on="patient_id", how="left")
    return merged.sort_values("patient_id").reset_index(drop=True)


def coverage_report(clinical: pd.DataFrame, patients: list[str]) -> pd.DataFrame:
    """How complete each variable is for the patients we have expression for.

    The brief singles MSI out: it is the project's one pre-registered subgroup
    variable, and if coverage is thin the variable has to change — a team
    decision taken before anyone looks at results.
    """
    sub = clinical.loc[clinical["patient_id"].isin(patients)]
    rows = []
    for column in ("stage", "site", "msi_status", "sex", "age"):
        present = int(sub[column].notna().sum())
        rows.append(
            {
                "variable": column,
                "n_patients": len(patients),
                "n_annotated": present,
                "coverage": round(present / len(patients), 4) if patients else float("nan"),
                "n_levels": int(sub[column].nunique(dropna=True)),
            }
        )
    return pd.DataFrame(rows)
