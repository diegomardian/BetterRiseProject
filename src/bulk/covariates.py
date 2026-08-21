"""W3.6 — load and enforce the pre-specified covariate set.

The point of a covariate lock is that it binds. A config nobody checks is a
document, not a lock, so this module is the enforcement:

- :func:`load_covariate_set` reads ``config/covariate_set.yaml``.
- :func:`require_locked` **refuses to let a survival model run against a config
  still marked ``proposed``.** The brief says the set is written down before any
  model is run; this is what makes that true rather than aspirational.
- :func:`build_design` assembles the modelling frame and reports what each
  filtering step cost, so the reconciliation is a by-product rather than
  something reconstructed afterwards.

W3.7 imports from here. It does not hard-code a formula.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from src.common.paths import CONFIG_DIR

COVARIATE_SET_PATH = CONFIG_DIR / "covariate_set.yaml"

#: Values in the clinical table that are labels for "we do not know", not levels.
NOT_A_LEVEL = ("conflicting",)

#: TCGA sample-type code for primary tumour. A patient's purity covariate has to
#: come from their tumour: ESTIMATE scores normal-adjacent tissue too, and gives
#: it a median "purity" of 0.69 (W3.3), which is not a number about that
#: patient's tumour at all.
PRIMARY_TUMOUR = "01"


class CovariateError(RuntimeError):
    """The covariate set is missing, malformed, or not yet locked."""


@lru_cache(maxsize=4)
def load_covariate_set(path: str | Path | None = None) -> dict[str, Any]:
    """Read the covariate config. Cached; pass a path in tests."""
    target = Path(path) if path is not None else COVARIATE_SET_PATH
    if not target.exists():
        raise CovariateError(f"{target} does not exist — W3.6 has not been done")
    spec = yaml.safe_load(target.read_text(encoding="utf-8"))
    for key in ("version", "status", "covariates", "endpoints", "model"):
        if key not in spec:
            raise CovariateError(f"{target} is missing required key {key!r}")
    if spec["status"] not in ("proposed", "locked"):
        raise CovariateError(f"status must be 'proposed' or 'locked', got {spec['status']!r}")
    return spec


def require_locked(spec: dict[str, Any]) -> None:
    """Raise unless the set is locked. Call this before fitting anything.

    A survival model fitted against a covariate set that is still under
    discussion is a model whose specification was chosen with the results
    already in view — which is exactly what pre-specification prevents.
    """
    if spec["status"] != "locked":
        raise CovariateError(
            f"covariate set {spec['version']} is status={spec['status']!r}, not 'locked'. "
            f"No survival model may run against a proposed set. Confirm it, flip "
            f"`status` to 'locked' in its own commit with a stated reason, then re-run."
        )


def covariate_names(
    spec: dict[str, Any],
    *,
    endpoint: str | None = None,
    context: str = "clinical_baseline",
) -> list[str]:
    """Covariates for an endpoint in a context, after the context's exclusions.

    Context matters because purity confounds analyses whose predictor is
    expression, and the clinical baseline has none. Asking for the covariates
    without saying which model you are building is how purity ends up costing
    65 patients in a model it does nothing for.
    """
    full = [c["name"] for c in spec["covariates"]]
    if endpoint is not None and endpoint not in spec["endpoints"]:
        raise CovariateError(
            f"unknown endpoint {endpoint!r}; known: {sorted(spec['endpoints'])}"
        )
    contexts = spec.get("contexts") or {}
    if context not in contexts:
        raise CovariateError(f"unknown context {context!r}; known: {sorted(contexts)}")

    settings = contexts[context]
    excluded = set(settings.get("exclude") or [])
    override = (settings.get("endpoint_overrides") or {}).get(endpoint or "", {})
    if override:
        excluded = set(override.get("exclude") or excluded)
    return [n for n in full if n not in excluded]


def total_df(spec: dict[str, Any], names: list[str]) -> int:
    """Degrees of freedom the named covariates cost. The number that binds."""
    by_name = {c["name"]: c for c in spec["covariates"]}
    missing = [n for n in names if n not in by_name]
    if missing:
        raise CovariateError(f"unknown covariate(s): {missing}")
    return sum(int(by_name[n]["df"]) for n in names)


def purity_column(spec: dict[str, Any], *, sensitivity: bool = False) -> str:
    """Which purity method the config names, primary or sensitivity."""
    for covariate in spec["covariates"]:
        if covariate["name"] == "purity":
            key = "sensitivity_source" if sensitivity else "source"
            if key not in covariate:
                raise CovariateError(f"purity covariate has no {key!r}")
            return str(covariate[key])
    raise CovariateError("no purity covariate in the set")


def build_design(
    clinical: pd.DataFrame,
    purity: pd.DataFrame,
    spec: dict[str, Any],
    *,
    endpoint: str,
    context: str = "clinical_baseline",
    sensitivity: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Assemble the modelling frame. Returns (design, attrition).

    ``attrition`` records n and events after each filtering step, in order. The
    brief wants a reconciliation of how many patients were dropped and why; this
    produces it as a side effect of building the frame rather than as a separate
    thing to keep in sync.
    """
    names = covariate_names(spec, endpoint=endpoint, context=context)
    method = purity_column(spec, sensitivity=sensitivity)
    strata = list(spec["model"].get("strata") or [])

    steps: list[dict[str, Any]] = []

    def record(label: str, frame: pd.DataFrame) -> None:
        steps.append(
            {
                "step": label,
                "n": int(len(frame)),
                "n_events": int((frame[endpoint] == 1).sum()) if endpoint in frame else None,
            }
        )

    frame = clinical.copy()
    record("clinical table", frame)

    frame = frame.loc[frame[f"usable_{endpoint}"]]
    record(f"usable for {endpoint}", frame)

    if "purity" in names:
        # Tumours only, BEFORE reducing to one row per patient. The purity table
        # is per (sample, method) and 50 patients have both a tumour and a
        # normal-adjacent ESTIMATE score; TCGA-AF-2689 has *only* a normal one.
        # Without this filter the patient's purity is whichever row sorts first,
        # which today is the tumour purely because "01" < "11" in the barcode.
        # That is accidental correctness, and it breaks the moment the purity
        # table is written in another order.
        if "sample_type" not in purity.columns:
            raise CovariateError(
                "the purity table has no `sample_type` column, so tumour samples "
                "cannot be told from normal-adjacent ones. ESTIMATE scores both, "
                "and a normal-adjacent score is not a covariate about the "
                "patient's tumour. Pass the table written by run_purity."
            )
        values = (
            purity.loc[
                (purity["method"] == method) & (purity["sample_type"] == PRIMARY_TUMOUR),
                ["patient_id", "purity"],
            ]
            .dropna(subset=["purity"])
            .drop_duplicates("patient_id")
        )
        if values.empty:
            raise CovariateError(
                f"no purity values for method {method!r} on primary tumours"
            )
        frame = frame.merge(values, on="patient_id", how="inner")
        record(f"has {method} purity", frame)

    # "conflicting" is a label for not-knowing, not a level to model.
    for name in names:
        if name in frame.columns and frame[name].dtype == object:
            frame[name] = frame[name].replace(list(NOT_A_LEVEL), pd.NA)

    needed = [n for n in names if n != "purity"] + strata
    for name in needed:
        if name not in frame.columns:
            raise CovariateError(f"clinical table has no column {name!r}")
        before = len(frame)
        frame = frame.loc[frame[name].notna()]
        if len(frame) != before:
            record(f"complete {name}", frame)

    keep = ["patient_id", endpoint, f"{endpoint}.time", *strata]
    keep += [n for n in names if n != "purity"]
    if "purity" in names:
        keep.append("purity")
    design = frame.loc[:, [c for c in dict.fromkeys(keep) if c in frame.columns]]
    record("final design", design)

    attrition = pd.DataFrame(steps)
    attrition["endpoint"] = endpoint
    attrition["context"] = context
    attrition["purity_method"] = method if "purity" in names else None
    attrition["dropped"] = attrition["n"].shift(1).sub(attrition["n"]).fillna(0).astype(int)
    return design.reset_index(drop=True), attrition


def events_per_df(
    design: pd.DataFrame,
    spec: dict[str, Any],
    endpoint: str,
    context: str = "clinical_baseline",
) -> dict[str, Any]:
    """Events per degree of freedom, with the conventional floor of 10 flagged."""
    names = covariate_names(spec, endpoint=endpoint, context=context)
    df = total_df(spec, names)
    events = int((design[endpoint] == 1).sum())
    ratio = events / df if df else float("nan")
    return {
        "endpoint": endpoint,
        "context": context,
        "n": int(len(design)),
        "n_events": events,
        "total_df": df,
        "events_per_df": round(ratio, 2),
        "meets_floor_of_10": bool(ratio >= 10),
    }
