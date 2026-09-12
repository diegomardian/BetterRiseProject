"""The workshop-overlap and mathematical-claims matrix.

Week 1, deliverable 3. Answers reviewer objection 3 in
``docs/DECISION_2026-09-11_scope_and_pivot.md`` §3 — "this is familiar testing
advice, an internal case study, and overlapping workshop material" — by
declaring the overlap precisely rather than minimising it.

Two guards, both of which can fail:

* a row claiming ``new`` must name the evidence that makes it new. A
  contribution asserted with nothing behind it is the thing reviewers are
  looking for;
* a row that is **not** ``new`` must name where in the workshop paper it already
  appears. "Overlaps somewhere" is not a declaration.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import yaml

from src.common.paths import REPO_ROOT

OVERLAP_PATH = REPO_ROOT / "config" / "workshop_overlap.yaml"

STATUSES: frozenset[str] = frozenset(
    {"already_published", "extended", "new", "contradicted"}
)
KINDS: frozenset[str] = frozenset({"overlap", "math_claim"})

FIELDS: tuple[str, ...] = ("id", "kind", "claim", "wmhs_location", "status", "note")


class OverlapError(ValueError):
    """The overlap declaration does not declare what it must."""


def load_overlap(path: Path | None = None) -> dict:
    with open(path or OVERLAP_PATH) as handle:
        return yaml.safe_load(handle)


def build_overlap(doc: dict) -> pd.DataFrame:
    rows = doc["rows"]
    if not rows:
        raise OverlapError("the overlap matrix declares nothing")

    seen: set[str] = set()
    for row in rows:
        missing = [f for f in FIELDS if f not in row]
        if missing:
            raise OverlapError(f"row {row.get('id', '?')!r} is missing {missing}")
        if row["id"] in seen:
            raise OverlapError(f"duplicate row id {row['id']!r}")
        seen.add(row["id"])
        if row["kind"] not in KINDS:
            raise OverlapError(f"{row['id']}: kind {row['kind']!r} is not in {sorted(KINDS)}")
        if row["status"] not in STATUSES:
            raise OverlapError(
                f"{row['id']}: status {row['status']!r} is not in {sorted(STATUSES)}"
            )
        if row["status"] == "new" and not row.get("evidence"):
            raise OverlapError(
                f"{row['id']}: claims to be new and names no evidence. A "
                "contribution with nothing behind it is what objection 3 is about."
            )
        if row["status"] != "new" and not row["wmhs_location"]:
            raise OverlapError(
                f"{row['id']}: status is {row['status']!r} and wmhs_location is "
                "empty. 'Overlaps somewhere' is not a declaration."
            )

    frame = pd.DataFrame(rows)
    for column in ("evidence",):
        if column not in frame.columns:
            frame[column] = None
    frame = frame[[*FIELDS, "evidence"]].copy()
    frame["is_contribution"] = frame["status"] == "new"
    frame["needs_paper_edit"] = frame["status"] == "contradicted"
    return frame.sort_values("id").reset_index(drop=True)


def overlap_summary(frame: pd.DataFrame) -> dict:
    return {
        "n_matrix_rows": int(len(frame)),
        "n_overlap_rows": int((frame["kind"] == "overlap").sum()),
        "n_math_claims": int((frame["kind"] == "math_claim").sum()),
        "n_already_published": int((frame["status"] == "already_published").sum()),
        "n_extended": int((frame["status"] == "extended").sum()),
        "n_new": int((frame["status"] == "new").sum()),
        "n_contradicted": int((frame["status"] == "contradicted").sum()),
        "novelty_share_of_overlap_rows": round(
            float((frame.loc[frame["kind"] == "overlap", "status"] == "new").mean()), 3
        ),
    }
