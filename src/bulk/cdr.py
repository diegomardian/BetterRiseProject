"""W3.5 — the curated clinical table, from TCGA-CDR (Liu et al. 2018).

CLAUDE.md invariant 9: **DSS and PFI are the primary endpoints. OS is secondary**
— COAD overall survival is contaminated by non-cancer death in an elderly
population.

THE CDR DISAGREES WITH INVARIANT 9, AND THAT IS WORTH KNOWING
--------------------------------------------------------------
The brief also says to honour the CDR's own recommended-use flags. Its notes
sheet says, verbatim:

    "For clinical outcome endpoints, we recommend the use of PFI for
     progression-free interval, and OS for overall survival. [...] Given the
     relatively short follow-up time, PFI is preferred over OS."

    "DSS is relatively accurate for CESC, PAAD, and UVM, and is approximated
     for other tumor types."

COAD and READ are "other tumor types". So the CDR recommends PFI and OS, and
calls DSS approximated for exactly our diseases.

These are reconcilable, and the reconciliation matters. The CDR derives DSS as
*dead **and** with tumour*, and says of it:

    "This is not a 100% accurate definition but is the best we could do with
     this dataset. Technically a patient could be with tumor but died of a car
     accident and therefore incorrectly considered as an event."

So DSS reduces the contamination invariant 9 objects to, without eliminating it
— it misclassifies in the same direction, just less often. **PFI is the endpoint
the project and the CDR agree on**, and is the one to lead with. Invariant 9
stands as written; this module implements it and records the tension rather than
quietly resolving it. Changing the invariant needs a PR and two approvals.

EXCLUSIONS ARE PER ENDPOINT, NOT PER PATIENT
--------------------------------------------
A patient with no ``tumor_status`` has no DSS but a perfectly good PFI. Dropping
them from the whole table would shrink every analysis to satisfy the strictest
one. So this emits ``usable_<endpoint>`` booleans and an exclusion reason per
endpoint, and drops nothing except redacted cases.

That is the same principle as invariant 1: "not estimable here" is a state to be
recorded, not a row to be deleted.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

#: PanCanAtlas-hosted copy of Supplemental Table S1.
CDR_URL = "https://api.gdc.cancer.gov/data/1b5f413e-a8d1-4d10-92eb-7c4ae739ed81"
CDR_FILENAME = "TCGA-CDR-SupplementalTableS1.xlsx"

#: TCGA writes missing data as bracketed sentinel strings. Left alone they
#: become categorical levels — a "[Not Available]" stage group in a Cox model.
SENTINELS = (
    "[Not Available]",
    "[Unknown]",
    "[Discrepancy]",
    "[Not Applicable]",
    "[Not Evaluated]",
    "[Completed]",
)

#: endpoint -> (role, note). Roles are the project's, per invariant 9; the notes
#: are the CDR authors' own assessment.
ENDPOINTS: dict[str, tuple[str, str]] = {
    "DSS": (
        "primary",
        "CDR calls DSS approximated for COAD/READ; derived as dead AND with "
        "tumour, so a with-tumour patient dying of another cause still counts.",
    ),
    "PFI": (
        "primary",
        "Recommended by the CDR and preferred over OS given short follow-up. "
        "The endpoint the project and the CDR agree on.",
    ),
    "OS": (
        "secondary",
        "Invariant 9: COAD OS is contaminated by non-cancer death. Reported, "
        "never led with.",
    ),
    "DFI": (
        "not_used",
        "391 of 629 COAD/READ patients have no DFI. Stage IV patients are "
        "excluded by construction. Too sparse to model.",
    ),
}


class CDRError(RuntimeError):
    """The CDR table could not be loaded or curated."""


# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------


def fetch_cdr(dest_dir: str | Path) -> Path:
    """Download the CDR workbook. Idempotent."""
    import requests

    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    target = dest_dir / CDR_FILENAME
    if target.exists() and target.stat().st_size > 0:
        return target
    response = requests.get(CDR_URL, timeout=600)
    response.raise_for_status()
    target.write_bytes(response.content)
    return target


def clean_sentinels(series: pd.Series) -> pd.Series:
    """Replace TCGA's bracketed sentinels with NaN."""
    if series.dtype == object:
        return series.replace(list(SENTINELS), np.nan)
    return series


def load_cdr(path: str | Path, types: tuple[str, ...] = ("COAD", "READ")) -> pd.DataFrame:
    """Read the CDR sheet, restricted to the given cancer types, sentinels cleaned."""
    raw = pd.read_excel(path, sheet_name="TCGA-CDR")
    if "bcr_patient_barcode" not in raw.columns:
        raise CDRError(f"{path} does not look like the TCGA-CDR sheet")
    sub = raw.loc[raw["type"].isin(types)].copy()
    if sub.empty:
        raise CDRError(f"no rows for types {types}")
    for column in sub.columns:
        sub[column] = clean_sentinels(sub[column])
    duplicated = sub["bcr_patient_barcode"].duplicated()
    if duplicated.any():
        raise CDRError(f"{int(duplicated.sum())} duplicate patient barcodes in the CDR subset")
    return sub.rename(columns={"bcr_patient_barcode": "patient_id"}).reset_index(drop=True)


# ---------------------------------------------------------------------------
# Exclusions — written out, per endpoint
# ---------------------------------------------------------------------------


def endpoint_exclusion_reason(row: pd.Series, endpoint: str) -> str | None:
    """Why this patient cannot contribute to this endpoint, or None if they can.

    The four rules, in the order they are applied:

    1. **Redacted.** The CDR's own ``Redaction`` column. Excluded from
       everything; these are cases TCGA withdrew.
    2. **No event indicator.** For DSS this means ``tumor_status`` was missing
       or discrepant, so "died of disease" is undefined — not zero.
    3. **No follow-up time.**
    4. **Non-positive follow-up time.** 22 COAD/READ patients have exactly 0
       days. A zero-length interval contributes no information to a Cox model
       and is dropped by most implementations silently; dropping it here makes
       the count visible.
    """
    if isinstance(row.get("Redaction"), str) and row["Redaction"].strip():
        return "redacted"
    event, time = row.get(endpoint), row.get(f"{endpoint}.time")
    if pd.isna(event):
        return "no_event_indicator"
    if pd.isna(time):
        return "no_followup_time"
    if float(time) <= 0:
        return "nonpositive_followup_time"
    return None


def add_usability_flags(cdr: pd.DataFrame) -> pd.DataFrame:
    """Add ``usable_<endpoint>`` and ``exclusion_<endpoint>`` for every endpoint."""
    out = cdr.copy()
    for endpoint in ENDPOINTS:
        reasons = out.apply(lambda r, e=endpoint: endpoint_exclusion_reason(r, e), axis=1)
        out[f"exclusion_{endpoint}"] = reasons
        out[f"usable_{endpoint}"] = reasons.isna()
    return out


def reconciliation(cdr: pd.DataFrame) -> pd.DataFrame:
    """How many patients were dropped from each endpoint and why.

    The brief's "done when": a reconciliation of how many patients were dropped
    and why at every step. One row per (endpoint, reason).
    """
    rows = []
    for endpoint, (role, _) in ENDPOINTS.items():
        reasons = cdr[f"exclusion_{endpoint}"]
        usable = int(reasons.isna().sum())
        rows.append(
            {
                "endpoint": endpoint,
                "role": role,
                "reason": "usable",
                "n": usable,
            }
        )
        for reason, count in reasons.dropna().value_counts().items():
            rows.append(
                {"endpoint": endpoint, "role": role, "reason": str(reason), "n": int(count)}
            )
    table = pd.DataFrame(rows)
    table["n_total"] = len(cdr)
    return table


def event_summary(cdr: pd.DataFrame) -> pd.DataFrame:
    """Events, censored and median follow-up per endpoint, on usable patients only."""
    rows = []
    for endpoint, (role, note) in ENDPOINTS.items():
        sub = cdr.loc[cdr[f"usable_{endpoint}"]]
        events = int((sub[endpoint] == 1).sum())
        rows.append(
            {
                "endpoint": endpoint,
                "role": role,
                "n_usable": int(len(sub)),
                "n_events": events,
                "n_censored": int((sub[endpoint] == 0).sum()),
                "event_rate": round(events / len(sub), 4) if len(sub) else float("nan"),
                "median_followup_days": round(float(sub[f"{endpoint}.time"].median()), 1)
                if len(sub)
                else float("nan"),
                "note": note,
            }
        )
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# The curated table
# ---------------------------------------------------------------------------


def build_curated_table(
    cdr: pd.DataFrame,
    gdc_clinical: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """CDR plus the covariates the brief lists, one row per patient.

    Stage, age, sex and the endpoints come from the CDR, which is the curated
    source. **Tumour site and MMR/MSI come from ``src/bulk/clinical.py``**, the
    GDC pull built for W3.4 — the CDR carries neither, and re-deriving them here
    would create a second source of truth for variables that already have one.

    Where both sources carry stage, the GDC value is kept as
    ``stage_gdc`` for comparison rather than merged over. Silent disagreement
    between two clinical sources is the kind of thing that surfaces in week 12.
    """
    from src.bulk.clinical import harmonise_stage

    curated = pd.DataFrame(
        {
            "patient_id": cdr["patient_id"],
            "project": "TCGA-" + cdr["type"],
            "age": cdr["age_at_initial_pathologic_diagnosis"],
            "sex": cdr["gender"].str.capitalize(),
            "stage": cdr["ajcc_pathologic_tumor_stage"].map(harmonise_stage),
            "stage_raw": cdr["ajcc_pathologic_tumor_stage"],
            "treatment_outcome_first_course": cdr["treatment_outcome_first_course"],
            "vital_status": cdr["vital_status"],
            "tumor_status": cdr["tumor_status"],
            "redacted": cdr["Redaction"].notna(),
        }
    )
    for endpoint in ENDPOINTS:
        curated[endpoint] = cdr[endpoint]
        curated[f"{endpoint}.time"] = cdr[f"{endpoint}.time"]
        curated[f"usable_{endpoint}"] = cdr[f"usable_{endpoint}"]
        curated[f"exclusion_{endpoint}"] = cdr[f"exclusion_{endpoint}"]

    if gdc_clinical is not None:
        keep = ["patient_id", "site", "msi_status", "n_distinct_calls"]
        available = [c for c in keep if c in gdc_clinical.columns]
        merged = curated.merge(
            gdc_clinical[available].rename(columns={"n_distinct_calls": "msi_n_calls"}),
            on="patient_id",
            how="left",
        )
        if "stage" in gdc_clinical.columns:
            merged = merged.merge(
                gdc_clinical[["patient_id", "stage"]].rename(columns={"stage": "stage_gdc"}),
                on="patient_id",
                how="left",
            )
        curated = merged

    return curated.sort_values("patient_id").reset_index(drop=True)


def stage_disagreement(curated: pd.DataFrame) -> pd.DataFrame:
    """Patients where the CDR and the GDC disagree about stage. Usually empty."""
    if "stage_gdc" not in curated.columns:
        return pd.DataFrame(columns=["patient_id", "stage", "stage_gdc"])
    both = curated.dropna(subset=["stage", "stage_gdc"])
    return both.loc[both["stage"] != both["stage_gdc"], ["patient_id", "stage", "stage_gdc"]]


def cohort_reconciliation(
    curated: pd.DataFrame, expression_patients: list[str]
) -> pd.DataFrame:
    """Where the CDR cohort and the expression cohort differ.

    Both counts matter. A patient with expression but no clinical record cannot
    enter a survival model; a patient with clinical data but no expression is
    irrelevant to this project but tells you the CDR is not the constraint.
    """
    cdr_ids = set(curated["patient_id"])
    expr_ids = set(expression_patients)
    return pd.DataFrame(
        [
            {"set": "CDR COAD/READ patients", "n": len(cdr_ids)},
            {"set": "patients with expression", "n": len(expr_ids)},
            {"set": "in both", "n": len(cdr_ids & expr_ids)},
            {"set": "expression but no CDR record", "n": len(expr_ids - cdr_ids)},
            {"set": "CDR but no expression", "n": len(cdr_ids - expr_ids)},
        ]
    )
