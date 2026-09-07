"""D2's locked design assembly, before the survival runner reads outcomes.

The survival specification is in ``docs/prereg_d2_survival.md``.  This module
does not choose covariates: it turns the already-locked expression-model
context into the exact D2 design, adding the pre-specified plate stratum and
continuous GUCA2A predictor.  The runner belongs in a separate change so this
module can be tested without writing or interpreting a survival result.
"""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from src.bulk.covariates import (
    PRIMARY_TUMOUR,
    build_design,
    covariate_names,
    purity_column,
    require_locked,
    total_df,
)
from src.bulk.d2_feasibility_gate import GENE, primary_tumour_values
from src.common.paths import REPO_ROOT

D2_CONTEXT = "expression_models"
D2_DESCRIPTIVE_CONTEXT = "d2_clinical_adjusted"
D2_ENDPOINTS = ("PFI", "DSS", "OS")
D2_PREREG = REPO_ROOT / "docs" / "prereg_d2_survival.md"
LEAD_ENDPOINT = "PFI"


class D2SurvivalError(RuntimeError):
    """The locked D2 design cannot be assembled honestly."""


def require_d2_locked(path: Path = D2_PREREG) -> None:
    """Refuse execution while the separate D2 contract is still a draft."""
    if not path.exists():
        raise D2SurvivalError(f"D2 pre-registration is missing: {path}")
    header = path.read_text(encoding="utf-8").splitlines()[:8]
    if not any(line.startswith("**Status:** locked") for line in header):
        raise D2SurvivalError(
            f"{path} is not locked. D2 may not read or fit survival outcomes."
        )


def d2_spec(spec: dict[str, Any]) -> dict[str, Any]:
    """Validate the inherited lock and add D2's pre-specified plate stratum."""
    require_locked(spec)
    require_d2_locked()
    endpoints = spec["endpoints"]
    required_roles = {"PFI": "primary", "DSS": "primary", "OS": "secondary"}
    for endpoint, role in required_roles.items():
        if endpoints.get(endpoint, {}).get("role") != role:
            raise D2SurvivalError(f"{endpoint} must remain {role!r} for D2")
    if not endpoints[LEAD_ENDPOINT].get("lead"):
        raise D2SurvivalError("PFI must remain D2's lead endpoint")
    if purity_column(spec) != "absolute":
        raise D2SurvivalError("D2 primary purity must be ABSOLUTE")
    if purity_column(spec, sensitivity=True) != "estimate_affy_extrapolated":
        raise D2SurvivalError("D2 sensitivity purity must be ESTIMATE")
    if "purity" not in covariate_names(spec, endpoint="PFI", context=D2_CONTEXT):
        raise D2SurvivalError("expression-model context must include purity")
    additional = spec.get("contexts", {}).get(D2_CONTEXT, {}).get("additional_required", [])
    if not any(item.get("name") == "plate" for item in additional):
        raise D2SurvivalError("expression-model context must require plate")

    out = deepcopy(spec)
    strata = list(out["model"].get("strata") or [])
    if "plate" not in strata:
        strata.append("plate")
    out["model"]["strata"] = strata
    return out


def d2_model_spec(spec: dict[str, Any], *, descriptive: bool = False) -> tuple[dict[str, Any], str]:
    """Return D2's locked primary or clinical-adjusted model context.

    The descriptive model omits purity by pre-registration §4, but it keeps
    the expression predictor and plate stratum. It is not allowed to reuse
    ``clinical_baseline``: that context is a different W3 estimand and does
    not require plate. The primary lock is validated before this derivative
    context is constructed, so the omission cannot be used to bypass D2's
    purity requirements.
    """
    locked = d2_spec(spec)
    if not descriptive:
        return locked, D2_CONTEXT

    out = deepcopy(locked)
    context = deepcopy(out["contexts"][D2_CONTEXT])
    context["exclude"] = sorted(set(context.get("exclude", [])) | {"purity"})
    overrides = deepcopy(context.get("endpoint_overrides", {}))
    for _endpoint, override in overrides.items():
        override["exclude"] = sorted(set(override.get("exclude", [])) | {"purity"})
    context["endpoint_overrides"] = overrides
    context["description"] = (
        "D2 pre-specified descriptive clinical-adjusted model: omits purity "
        "to show compositional confounding; it cannot establish the D2 claim."
    )
    out["contexts"][D2_DESCRIPTIVE_CONTEXT] = context
    return out, D2_DESCRIPTIVE_CONTEXT


def primary_tumour_annotations(manifest: pd.DataFrame) -> pd.DataFrame:
    """One pre-selected primary aliquot and its plate per participant."""
    required = {"barcode", "patient_id", "sample_type", "plate"}
    missing = sorted(required - set(manifest.columns))
    if missing:
        raise D2SurvivalError(f"sample manifest is missing {missing}")
    tumour = manifest.loc[manifest["sample_type"] == PRIMARY_TUMOUR,
                           ["patient_id", "barcode", "plate"]].copy()
    if tumour.empty:
        raise D2SurvivalError("sample manifest has no primary-tumour rows")
    if tumour["patient_id"].duplicated().any():
        raise D2SurvivalError(
            "primary-tumour manifest has duplicate participants; ingest must "
            "deduplicate aliquots before D2"
        )
    # Missing plate is a complete-case exclusion, not a malformed manifest.
    # It must reach ``build_design`` so the attrition table names and counts it.
    return tumour


def _recompute_attrition_drops(attrition: pd.DataFrame) -> pd.DataFrame:
    """Recompute drops after a caller prepends an upstream eligibility step."""
    out = attrition.copy()
    out["dropped"] = out["n"].shift(1).sub(out["n"]).fillna(0).astype(int)
    return out


def prepare_design(
    clinical: pd.DataFrame,
    purity: pd.DataFrame,
    manifest: pd.DataFrame,
    expression: pd.DataFrame,
    gene_id: str,
    spec: dict[str, Any],
    *,
    endpoint: str,
    sensitivity: bool = False,
    descriptive: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build one endpoint/purity-source D2 design and record expression loss."""
    if endpoint not in D2_ENDPOINTS:
        raise D2SurvivalError(f"D2 endpoint {endpoint!r} is not pre-specified")
    locked, context = d2_model_spec(spec, descriptive=descriptive)
    annotations = primary_tumour_annotations(manifest)
    values = primary_tumour_values(expression, manifest, gene_id).rename(GENE)
    if values.index.has_duplicates:
        raise D2SurvivalError("GUCA2A values are not one row per participant")

    if clinical["patient_id"].duplicated().any():
        raise D2SurvivalError("clinical table has duplicate participants")
    # A primary specimen is an eligibility condition separate from its plate:
    # retaining the barcode lets the attrition record distinguish no RNA sample
    # from a primary sample whose required plate annotation is missing.
    clinical_with_sample = clinical.merge(
        annotations[["patient_id", "barcode", "plate"]], on="patient_id", how="left",
        validate="one_to_one",
    )
    start = pd.DataFrame([{
        "step": "clinical table",
        "n": int(len(clinical_with_sample)),
        "n_events": int(clinical_with_sample[endpoint].sum()),
        "endpoint": endpoint,
        "context": context,
        "purity_method": purity_column(locked, sensitivity=sensitivity),
    }])
    clinical_with_plate = clinical_with_sample.loc[
        clinical_with_sample["barcode"].notna()
    ].drop(columns="barcode")
    design, attrition = build_design(
        clinical_with_plate, purity, locked, endpoint=endpoint,
        context=context, sensitivity=sensitivity,
    )
    # ``build_design``'s first row is this already-filtered frame. Rename it
    # rather than presenting it as the original clinical table, then preserve
    # the complete source table as the actual first row.
    attrition = attrition.copy()
    attrition.loc[attrition.index[0], "step"] = "has primary-tumour RNA sample"
    attrition = _recompute_attrition_drops(pd.concat([start, attrition], ignore_index=True))
    before = len(design)
    design = design.merge(values.rename(GENE), left_on="patient_id", right_index=True,
                          how="inner", validate="one_to_one")
    if design[GENE].isna().any() or not np.isfinite(design[GENE].to_numpy()).all():
        raise D2SurvivalError("final D2 design has non-finite GUCA2A")
    method = purity_column(locked, sensitivity=sensitivity)
    attrition = pd.concat([
        attrition,
        pd.DataFrame([{
            "step": "has mapped GUCA2A expression",
            "n": int(len(design)),
            "n_events": int(design[endpoint].sum()),
            "endpoint": endpoint,
            "context": context,
            "purity_method": method,
            "dropped": int(before - len(design)),
        }]),
    ], ignore_index=True)
    attrition = _recompute_attrition_drops(attrition)
    return design.reset_index(drop=True), attrition


def d2_events_per_df(
    design: pd.DataFrame, spec: dict[str, Any], *, endpoint: str, descriptive: bool = False
) -> dict[str, Any]:
    """The pre-specified 10-events-per-df lead-endpoint gate, including GUCA2A."""
    locked, context = d2_model_spec(spec, descriptive=descriptive)
    names = covariate_names(locked, endpoint=endpoint, context=context)
    df = total_df(locked, names) + 1  # continuous GUCA2A predictor
    events = int(design[endpoint].sum())
    ratio = events / df if df else float("nan")
    strata = list(locked["model"].get("strata") or [])
    if strata:
        n_strata = int(design.groupby(strata, dropna=False).ngroups)
        event_strata = design.loc[design[endpoint] == 1]
        n_event_strata = int(event_strata.groupby(strata, dropna=False).ngroups)
    else:
        n_strata = 1
        n_event_strata = int(events > 0)
    return {
        "endpoint": endpoint,
        "context": context,
        "n": int(len(design)),
        "n_events": events,
        "non_stratum_df": df,
        "events_per_df": round(ratio, 2),
        "n_nonempty_strata": n_strata,
        "n_event_contributing_strata": n_event_strata,
        "meets_lead_floor": bool(ratio >= 10) if endpoint == LEAD_ENDPOINT else None,
    }


def retained_vs_dropped(
    clinical: pd.DataFrame,
    values: pd.Series,
    design: pd.DataFrame,
    *,
    endpoint: str,
    purity_method: str,
) -> pd.DataFrame:
    """Required §6a outcome-adjacent attrition audit for GUCA2A distribution."""
    usable = clinical.loc[clinical[f"usable_{endpoint}"], ["patient_id"]]
    base = usable.merge(values.rename(GENE), left_on="patient_id", right_index=True,
                        how="inner", validate="one_to_one")
    base["retained"] = base["patient_id"].isin(set(design["patient_id"]))
    rows = []
    for retained, part in base.groupby("retained", sort=True):
        x = part[GENE].to_numpy(dtype=float)
        rows.append({
            "endpoint": endpoint,
            "purity_method": purity_method,
            "membership": "retained" if retained else "dropped",
            "n_participants": int(len(x)),
            "minimum": float(np.min(x)) if len(x) else np.nan,
            "q25": float(np.percentile(x, 25)) if len(x) else np.nan,
            "median": float(np.median(x)) if len(x) else np.nan,
            "q75": float(np.percentile(x, 75)) if len(x) else np.nan,
            "maximum": float(np.max(x)) if len(x) else np.nan,
            "iqr": float(np.percentile(x, 75) - np.percentile(x, 25)) if len(x) else np.nan,
        })
    return pd.DataFrame(rows)
