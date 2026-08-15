"""Output schema and results writer. SHARED AND FROZEN.

execution_plan.md §3.3. Changing anything in this file requires a PR, two
approvals and a written reason (CLAUDE.md invariant 3).

One row per (patient, gene, granularity rung, labelling axis, study).

THE POINT OF THIS FILE
----------------------
``None`` is not ``0.0``. An unestimable intrinsic term is ``None`` with
``estimability="not_estimable"``. Reporting zero when the truthful answer is
"not estimable" is the single most likely route to a wrong conclusion in this
project, so it is enforced here by assertion, in the writer, rather than left
to review. See :func:`validate_results`.

USE
---
    from src.schema import empty_results_frame, write_results

    df = empty_results_frame()
    ...                       # fill it
    path = write_results(df, name="lee_decomposition", seed=20260815)

Never call ``df.to_parquet()`` directly for a results table. The writer is what
stamps the git sha and the seed onto every result (CLAUDE.md invariant 10).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Final

import pandas as pd

from src.common.paths import RESULTS_DIR
from src.common.provenance import provenance_record

# ---------------------------------------------------------------------------
# Controlled vocabularies
# ---------------------------------------------------------------------------

GRANULARITY_RUNGS: Final = ("epithelial", "lineage", "crypt_position", "best4")
LABELING_AXES: Final = ("stem_pole", "opposite_lineage", "chromatin", "spatial")
WEIGHTINGS: Final = ("normal", "tumour", "doubly_robust")
ESTIMABILITY: Final = ("ok", "wide_interval", "not_estimable")

# ---------------------------------------------------------------------------
# Column contract. Order is the on-disk column order.
# ---------------------------------------------------------------------------

#: column -> (pandas dtype, nullable)
COLUMNS: Final[dict[str, tuple[str, bool]]] = {
    "patient_id": ("string", False),
    "study_id": ("string", False),
    "gene": ("string", False),
    "granularity_rung": ("string", False),
    "labeling_axis": ("string", False),
    "n_cells_mature": ("Int64", False),
    "compositional": ("Float64", True),
    "intrinsic": ("Float64", True),  # None != 0.0 — this distinction is the project
    "interaction": ("Float64", True),  # reported separately, never folded in
    "weighting": ("string", False),
    "estimability": ("string", False),
    "ci_low": ("Float64", True),
    "ci_high": ("Float64", True),
}

REQUIRED_COLUMNS: Final = tuple(COLUMNS)

_ENUM_COLUMNS: Final = {
    "granularity_rung": GRANULARITY_RUNGS,
    "labeling_axis": LABELING_AXES,
    "weighting": WEIGHTINGS,
    "estimability": ESTIMABILITY,
}

#: (patient, gene, rung, axis, study) plus weighting — one row each.
KEY_COLUMNS: Final = (
    "patient_id",
    "study_id",
    "gene",
    "granularity_rung",
    "labeling_axis",
    "weighting",
)


class SchemaViolation(AssertionError):
    """Raised when a results frame breaks the frozen output contract."""


# ---------------------------------------------------------------------------
# Construction
# ---------------------------------------------------------------------------


def empty_results_frame() -> pd.DataFrame:
    """An empty frame with the frozen columns and dtypes. Start here."""
    return pd.DataFrame({col: pd.Series(dtype=dt) for col, (dt, _) in COLUMNS.items()})


def coerce_results(df: pd.DataFrame) -> pd.DataFrame:
    """Reorder to the contract and cast to the contract dtypes.

    Raises SchemaViolation on missing or unexpected columns rather than
    silently adding or dropping them.
    """
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise SchemaViolation(f"missing required columns: {missing}")
    extra = [c for c in df.columns if c not in REQUIRED_COLUMNS]
    if extra:
        raise SchemaViolation(
            f"unexpected columns: {extra}. The schema is frozen — adding a field "
            f"needs a PR with two approvals (CLAUDE.md invariant 3)."
        )
    out = df.loc[:, list(REQUIRED_COLUMNS)].copy()
    for col, (dtype, _) in COLUMNS.items():
        out[col] = out[col].astype(dtype)
    return out


# ---------------------------------------------------------------------------
# Validation — invariant 1 lives here
# ---------------------------------------------------------------------------


def validate_results(df: pd.DataFrame) -> None:
    """Assert the frozen contract. Raises SchemaViolation with a specific message.

    Checks, in order of how much damage they prevent:

    1. ``not_estimable`` rows carry ``intrinsic is None``, never ``0.0``.
       CLAUDE.md invariant 1. This is the one that matters.
    2. Enum columns contain only allowed values.
    3. Non-nullable columns are fully populated.
    4. Keys are unique.
    5. CIs are ordered and present together.
    """
    if list(df.columns) != list(REQUIRED_COLUMNS):
        raise SchemaViolation("frame is not coerced — call coerce_results() first")

    # --- 1. THE invariant -------------------------------------------------
    not_est = df["estimability"] == "not_estimable"
    offenders = df.loc[not_est & df["intrinsic"].notna()]
    if len(offenders):
        sample = offenders.loc[:, ["patient_id", "gene", "intrinsic"]].head(5)
        raise SchemaViolation(
            f"{len(offenders)} row(s) marked not_estimable carry a non-null "
            f"intrinsic value. None is not 0.0 — an unestimable intrinsic term "
            f"must be None (CLAUDE.md invariant 1).\n{sample.to_string(index=False)}"
        )
    # And the converse: an estimable row that quietly lost its estimate.
    ok_rows = df["estimability"] == "ok"
    lost = df.loc[ok_rows & df["intrinsic"].isna()]
    if len(lost):
        raise SchemaViolation(
            f"{len(lost)} row(s) marked estimability='ok' have a null intrinsic "
            f"term. Either the estimate exists or the row is not 'ok'."
        )

    # --- 2. enums ---------------------------------------------------------
    for col, allowed in _ENUM_COLUMNS.items():
        bad = set(df[col].dropna().unique()) - set(allowed)
        if bad:
            raise SchemaViolation(f"{col} contains values outside the frozen set: {sorted(bad)}")

    # --- 3. required values ----------------------------------------------
    for col, (_, nullable) in COLUMNS.items():
        if not nullable and df[col].isna().any():
            n = int(df[col].isna().sum())
            raise SchemaViolation(f"{col} is not nullable but has {n} null value(s)")

    if (df["n_cells_mature"] < 0).any():
        raise SchemaViolation("n_cells_mature has negative values")

    # --- 4. keys ----------------------------------------------------------
    dupes = df.duplicated(subset=list(KEY_COLUMNS), keep=False)
    if dupes.any():
        example = df.loc[dupes, list(KEY_COLUMNS)].head(4)
        raise SchemaViolation(
            f"{int(dupes.sum())} duplicate key rows. One row per (patient, study, "
            f"gene, rung, axis, weighting).\n{example.to_string(index=False)}"
        )

    # --- 5. intervals -----------------------------------------------------
    one_sided = df["ci_low"].isna() != df["ci_high"].isna()
    if one_sided.any():
        raise SchemaViolation(f"{int(one_sided.sum())} row(s) have only one CI bound")
    both = df["ci_low"].notna() & df["ci_high"].notna()
    inverted = both & (df["ci_low"] > df["ci_high"])
    if inverted.any():
        raise SchemaViolation(f"{int(inverted.sum())} row(s) have ci_low > ci_high")


# ---------------------------------------------------------------------------
# Writer — invariant 10 lives here
# ---------------------------------------------------------------------------


def write_results(
    df: pd.DataFrame,
    name: str,
    *,
    seed: int,
    results_dir: Path | None = None,
    notes: str | None = None,
    allow_dirty: bool = False,
) -> Path:
    """Validate, then write versioned parquet plus a provenance sidecar.

    Lands at ``results/<YYYY-MM-DD>_<sha7>/<name>.parquet`` with a matching
    ``<name>.meta.json`` recording git sha, dirty state, seed, row count and
    package versions.

    Parameters
    ----------
    seed:
        The random seed that produced this table. Required, not defaulted —
        every result carries a fixed seed (CLAUDE.md invariant 10). Pass the
        seed even for deterministic tables; ``0`` is a fine answer, silence
        is not.
    allow_dirty:
        Refuse to write from a dirty working tree unless True. A result whose
        sha does not reproduce it is worse than no result. Exploratory runs may
        pass True; anything shown to the team must not.
    """
    df = coerce_results(df)
    validate_results(df)

    prov = provenance_record(seed=seed, notes=notes)
    if prov["git_dirty"] and not allow_dirty:
        raise SchemaViolation(
            "working tree is dirty — the recorded sha would not reproduce this "
            "table. Commit first, or pass allow_dirty=True for a scratch run."
        )

    base = Path(results_dir) if results_dir is not None else RESULTS_DIR
    out_dir = base / f"{prov['date']}_{prov['git_sha_short']}"
    out_dir.mkdir(parents=True, exist_ok=True)

    parquet_path = out_dir / f"{name}.parquet"
    df.to_parquet(parquet_path, index=False)

    prov |= {
        "name": name,
        "n_rows": int(len(df)),
        "n_patients": int(df["patient_id"].nunique()),
        "n_not_estimable": int((df["estimability"] == "not_estimable").sum()),
        "schema_columns": list(REQUIRED_COLUMNS),
    }
    (out_dir / f"{name}.meta.json").write_text(json.dumps(prov, indent=2), encoding="utf-8")
    return parquet_path


def read_results(path: str | Path) -> pd.DataFrame:
    """Read a results parquet back, coerced and validated."""
    df = coerce_results(pd.read_parquet(path))
    validate_results(df)
    return df
