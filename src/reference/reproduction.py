"""The audit's deliverables as a machine-readable reproduction index.

    from src.reference.reproduction import DELIVERABLES, build_manifest

WHY A TABLE AND NOT PROSE. The decision record's week 4 asks to *"package the
audit for another researcher"* and to *"generate every reported number from
tables"*. A reproduction section written in prose is one more document that can
drift from the code; a table that names, for each deliverable, the exact module
to run and the versioned result it produced can be checked. ``build_manifest``
resolves each table to its most recently **written** run (not its
lexicographically last name -- see ``table_resolution``), reads the sidecar, and
records whether the table is committed. A table that is absent, or a deliverable
that cannot be used for its purpose, carries a ``status`` and a reason rather
than being dropped.

THREE STATES, NOT A BOOLEAN. A challenge whose every case saw the ledger is
*committed but blocked*: the file exists and the deliverable is not done. A
boolean ``present=True`` for that row would read as success for the one thing
the submission is gated on -- the same shape as the defects this paper is
about. So ``status`` is ``present`` (committed and usable), ``blocked``
(cannot be used, whatever the file state) or ``absent`` (no committed result),
with ``table_committed`` keeping the raw file fact separate.

PATHS ARE REPO-RELATIVE. A reproduction index that records
``/Users/someone/.../BetterRiseProject-audit/results/...`` resolves for nobody
else, and week 4's "another researcher in a fresh environment" fails as
recorded. Every path here is ``results/<date>_<sha>/<table>.parquet``.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from src.reference.table_resolution import newest_by_time


@dataclass(frozen=True)
class Deliverable:
    """One committed output, its producer, and how to rebuild it."""

    #: Which week of DECISION_2026-09-11 §4 the deliverable belongs to.
    week: int
    area: str
    name: str
    table: str
    command: str
    note: str = ""
    #: True when the deliverable cannot be used for its purpose regardless of
    #: whether a table is committed -- a blocked role, not a missing file.
    blocked: bool = False


#: Every result the audit plan names, plus the two companion tables and the
#: resolver finding. Kept as data so a new deliverable is one tuple entry, not a
#: paragraph in a document.
DELIVERABLES: tuple[Deliverable, ...] = (
    Deliverable(
        1, "ledger", "classified claim-check inventory", "claim_check_inventory",
        "python -m src.reference.jobs.claim_check_inventory",
        "one worker's classification; inter-rater agreement NOT established",
    ),
    Deliverable(
        1, "challenge", "blinded guard challenge", "blinded_guard_challenge",
        "python -m src.reference.jobs.guard_challenge",
        "all cases saw_ledger=true, so is_held_out_evaluation=false and the "
        "week-one stop rule is not readable from this table. Needs the "
        "unassigned blinded challenger.",
        blocked=True,
    ),
    Deliverable(
        1, "overlap", "workshop overlap matrix", "workshop_overlap",
        "python -m src.reference.jobs.workshop_overlap",
    ),
    Deliverable(
        2, "interval", "interval stress calibration", "interval_stress_calibration",
        "python -m src.reference.jobs.interval_stress_calibration",
    ),
    Deliverable(
        2, "interval", "interval stress summary",
        "interval_stress_calibration_summary",
        "python -m src.reference.jobs.interval_stress_calibration",
    ),
    Deliverable(
        2, "cutpoint", "cutpoint crossing brackets", "cutpoint_crossing_brackets",
        "python -m src.reference.jobs.cutpoint_crossing_brackets",
    ),
    Deliverable(
        2, "cutpoint", "cutpoint crossing summary",
        "cutpoint_crossing_brackets_summary",
        "python -m src.reference.jobs.cutpoint_crossing_brackets",
    ),
    Deliverable(
        2, "challenge", "guard challenge results", "guard_challenge_results",
        "NOT RUN",
        "requires the blinded challenger; the stop rule that gates submission "
        "is not readable without it.",
        blocked=True,
    ),
    Deliverable(
        3, "adenoma", "adenoma claim sensitivity", "adenoma_claim_sensitivity",
        "python -m src.reference.jobs.adenoma_claim_sensitivity",
        "reproduction from derived inputs, not raw-data replication",
    ),
    Deliverable(
        3, "adenoma", "adenoma claim sensitivity summary",
        "adenoma_claim_sensitivity_summary",
        "python -m src.reference.jobs.adenoma_claim_sensitivity",
    ),
    Deliverable(
        3, "adenoma", "estimability attrition", "estimability_attrition",
        "python -m src.reference.jobs.adenoma_claim_sensitivity",
    ),
    Deliverable(
        3, "crowell", "Crowell summary reproduction",
        "crowell_summary_reproduction",
        "python -m src.reference.jobs.crowell_summary_reproduction",
        "reproduction from derived inputs, not raw-data replication",
    ),
    Deliverable(
        3, "disval", "DIS/VAL reproduction", "disval_reproduction",
        "python -m src.reference.jobs.disval_reproduction",
        "reproduction from derived inputs, not raw-data replication",
    ),
    Deliverable(
        4, "audit", "same-date resolver ambiguity", "table_resolution_ambiguity",
        "python -m src.reference.jobs.table_resolution_audit",
        "found by this audit; 23 ambiguous pairs, 9 with different content",
    ),
)

#: Status values, most to least usable. ``blocked`` outranks ``absent``: a
#: deliverable gated on an unfilled role is blocked whether or not a partial
#: table exists.
STATUSES: tuple[str, ...] = ("present", "blocked", "absent")


def _project_relative(path: Path, results_dir: Path) -> str:
    """``results/<date>_<sha>/<table>.parquet``, never an absolute path."""
    try:
        return f"{results_dir.name}/{path.relative_to(results_dir).as_posix()}"
    except ValueError:
        return path.name


def build_manifest(results_dir: Path) -> pd.DataFrame:
    """Resolve each deliverable to its newest run and record its provenance."""
    rows = []
    for item in DELIVERABLES:
        path = newest_by_time(results_dir, item.table)
        committed = path is not None
        if item.blocked:
            status = "blocked"
        elif committed:
            status = "present"
        else:
            status = "absent"
        row: dict[str, object] = {
            "week": item.week,
            "area": item.area,
            "deliverable": item.name,
            "table": item.table,
            "command": item.command,
            "note": item.note,
            "status": status,
            "table_committed": committed,
            "path": _project_relative(path, results_dir) if committed else None,
            "git_sha": None, "git_dirty": None, "seed": None,
            "n_rows": None, "utc_timestamp": None,
        }
        if committed:
            sidecar = path.parent / f"{item.table}.meta.json"
            if sidecar.exists():
                meta = json.loads(sidecar.read_text())
                row |= {
                    key: meta.get(key)
                    for key in ("git_sha", "git_dirty", "seed", "n_rows",
                                "utc_timestamp")
                }
        rows.append(row)
    return pd.DataFrame(rows)


def unavailable_deliverables(manifest: pd.DataFrame) -> pd.DataFrame:
    """Rows a reader must not count as done: blocked or absent."""
    if manifest.empty:
        return manifest
    return manifest[manifest["status"] != "present"]
