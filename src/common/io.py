"""Versioned table writer for tables that are NOT decomposition results.

``src/schema.py`` owns the frozen per-(patient, gene, rung, axis, study) output
contract and validates against it. Harness tables — a bake-off ranking, an
attenuation sweep, a calibration report — have entirely different shapes and do
not fit that schema. Bending them to fit would mean editing frozen code for no
benefit.

So: same ``results/{date}_{sha7}/`` convention, same provenance sidecar, no
schema validation. Callers validate their own shapes.

This duplicates about fifteen lines of directory-naming logic from
``schema.write_results``. That is deliberate — unifying them means touching
frozen code, and it is only worth a ``shared/`` PR if a third writer appears.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from src.common.paths import RESULTS_DIR
from src.common.provenance import provenance_record


class DirtyTreeError(RuntimeError):
    """Refused to write a result whose recorded sha would not reproduce it."""


def versioned_dir(seed: int, results_dir: Path | None = None) -> Path:
    """``results/{date}_{sha7}/``, created if absent."""
    prov = provenance_record(seed=seed)
    base = Path(results_dir) if results_dir is not None else RESULTS_DIR
    out = base / f"{prov['date']}_{prov['git_sha_short']}"
    out.mkdir(parents=True, exist_ok=True)
    return out


def write_versioned_table(
    df: pd.DataFrame,
    name: str,
    *,
    seed: int,
    results_dir: Path | None = None,
    notes: str | None = None,
    extra_meta: dict | None = None,
    allow_dirty: bool = False,
) -> Path:
    """Write ``df`` to versioned parquet with a provenance sidecar.

    Parameters
    ----------
    seed:
        Required, not defaulted. Every result carries a fixed seed
        (CLAUDE.md invariant 10), and a sweep without one is not a result.
    extra_meta:
        Table-specific fields folded into the sidecar — grid definition,
        replicate count, which deconvolution methods actually ran. Use it. The
        sidecar is where a reader finds out that CIBERSORTx was skipped.
    allow_dirty:
        Refuse to write from a dirty working tree unless True.
    """
    if not isinstance(df, pd.DataFrame):
        raise TypeError(f"expected a DataFrame, got {type(df).__name__}")

    prov = provenance_record(seed=seed, notes=notes)
    if prov["git_dirty"] and not allow_dirty:
        raise DirtyTreeError(
            "working tree is dirty — the recorded sha would not reproduce this "
            "table. Commit first, or pass allow_dirty=True for a scratch run."
        )

    base = Path(results_dir) if results_dir is not None else RESULTS_DIR
    out_dir = base / f"{prov['date']}_{prov['git_sha_short']}"
    out_dir.mkdir(parents=True, exist_ok=True)

    path = out_dir / f"{name}.parquet"
    df.to_parquet(path, index=False)

    prov |= {
        "name": name,
        "n_rows": int(len(df)),
        "columns": list(df.columns),
        "table_kind": "harness",
    }
    if extra_meta:
        prov |= dict(extra_meta)
    (out_dir / f"{name}.meta.json").write_text(json.dumps(prov, indent=2), encoding="utf-8")
    return path


def read_versioned_table(path: str | Path) -> tuple[pd.DataFrame, dict]:
    """Read a harness table back with its sidecar.

    Returns ``(df, meta)``. The sidecar matters as much as the table — it is
    where the seed, the sha and any skipped methods are recorded.
    """
    path = Path(path)
    df = pd.read_parquet(path)
    meta_path = path.with_suffix(".meta.json")
    meta = json.loads(meta_path.read_text(encoding="utf-8")) if meta_path.exists() else {}
    return df, meta
