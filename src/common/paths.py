"""Canonical paths. Import these instead of writing ``../../data`` anywhere.

Four people on three operating systems. Relative paths from notebooks are the
most common way a repo like this stops working on someone else's machine.
"""

from __future__ import annotations

import os
from pathlib import Path

# src/common/paths.py -> src/common -> src -> repo root
REPO_ROOT: Path = Path(__file__).resolve().parents[2]

CONFIG_DIR: Path = REPO_ROOT / "config"
RESULTS_DIR: Path = REPO_ROOT / "results"
ENV_DIR: Path = REPO_ROOT / "env"

# Data lives outside git. Point BRP_DATA_DIR at a scratch disk or a shared
# drive if you do not want 371k cells inside your repo folder.
DATA_DIR: Path = Path(os.environ.get("BRP_DATA_DIR", REPO_ROOT / "data"))
RAW_DIR: Path = DATA_DIR / "raw"
INTERIM_DIR: Path = DATA_DIR / "interim"
PROCESSED_DIR: Path = DATA_DIR / "processed"
MANIFEST_PATH: Path = REPO_ROOT / "data" / "manifest.csv"  # always in-repo


def ensure_data_dirs() -> None:
    """Create the data tree if it is missing. Safe to call repeatedly."""
    for d in (RAW_DIR, INTERIM_DIR, PROCESSED_DIR):
        d.mkdir(parents=True, exist_ok=True)


def s_matrix_path(rung: str, version: str) -> Path:
    """The W1 -> W2/W3 handoff artifact.

    ``S_matrix_{rung}_{version}.parquet`` on the fixed gene index. W3 emits bulk
    on the same index, so integration is a join and not a negotiation
    (CLAUDE.md, repo layout).
    """
    return PROCESSED_DIR / "reference" / f"S_matrix_{rung}_{version}.parquet"


def gene_index_path(version: str) -> Path:
    """The fixed gene index both arms are built on. See config/gene_index/README.md."""
    return CONFIG_DIR / "gene_index" / f"gene_index_{version}.txt"
