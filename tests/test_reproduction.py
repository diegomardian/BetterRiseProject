"""Forcing inputs for the audit reproduction index.

The index must not silently drop a deliverable: an absent table is present=False
with its note, a blocked one stays blocked, and every command is executable.
"""

from __future__ import annotations

import json

import pandas as pd

from src.reference.reproduction import (
    DELIVERABLES,
    build_manifest,
    missing_deliverables,
)


def _write_result(tmp_path, dirname, table, **meta):
    d = tmp_path / dirname
    d.mkdir()
    pd.DataFrame({"a": [1, 2]}).to_parquet(d / f"{table}.parquet")
    base = {"utc_timestamp": "2026-01-01T00:00:00+00:00", "git_sha": "abc1234",
            "git_dirty": False, "seed": 7, "n_rows": 2}
    base.update(meta)
    (d / f"{table}.meta.json").write_text(json.dumps(base))
    return d


def test_the_deliverable_table_names_are_unique():
    tables = [d.table for d in DELIVERABLES]
    assert len(tables) == len(set(tables)), "a duplicated table double-counts"


def test_every_command_is_executable_or_an_explicit_not_run():
    for item in DELIVERABLES:
        assert item.command.startswith("python -m src.") or item.command == "NOT RUN"
        assert item.note or item.command != "NOT RUN", (
            f"{item.table}: a NOT RUN deliverable must say why"
        )


def test_a_present_deliverable_carries_its_provenance(tmp_path):
    present = DELIVERABLES[0].table
    _write_result(tmp_path, "2026-01-01_aaaaaaa", present)
    manifest = build_manifest(tmp_path)
    row = manifest[manifest.table == present].iloc[0]
    assert row["present"]
    assert row["git_sha"] == "abc1234"
    assert row["seed"] == 7
    assert row["n_rows"] == 2


def test_an_absent_deliverable_is_recorded_not_dropped(tmp_path):
    manifest = build_manifest(tmp_path)
    assert len(manifest) == len(DELIVERABLES)
    assert not manifest["present"].any()
    missing = missing_deliverables(manifest)
    assert len(missing) == len(DELIVERABLES)


def test_the_blinded_challenge_stays_blocked(tmp_path):
    """The index must not let a blocked deliverable read as done."""
    manifest = build_manifest(tmp_path).set_index("table")
    assert not manifest.loc["blinded_guard_challenge", "present"]
    assert "BLOCKED" in manifest.loc["blinded_guard_challenge", "note"]
    assert not manifest.loc["guard_challenge_results", "present"]
