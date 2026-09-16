"""Locate the result tables pinned for this submission.

The manifest is deliberately version-specific. A later experiment must not
silently change a submission figure merely because its directory sorts later.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
MANIFEST = Path(__file__).with_name("results_manifest.json")


def pinned(name: str) -> Path:
    """Return the submission-pinned table named by the manifest."""
    entries = json.loads(MANIFEST.read_text(encoding="utf-8"))
    if name not in entries:
        print(
            f"{name!r} is not pinned in {MANIFEST.relative_to(REPO_ROOT)}",
            file=sys.stderr,
        )
        raise SystemExit(1)
    chosen = REPO_ROOT / entries[name]
    if not chosen.is_file():
        print(f"pinned result is missing: {chosen.relative_to(REPO_ROOT)}", file=sys.stderr)
        raise SystemExit(1)
    print(f"reading {chosen}")
    return chosen
