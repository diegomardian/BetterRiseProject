"""Forcing inputs for the shared result-table resolver.

The audit found 23 table names where name order and time order disagree within a
date. `newest` is kept because callers rely on it, but it must not be mistaken
for `newest_by_time`, and `ambiguous_same_date` must surface the disagreement
rather than resolve it silently.
"""

from __future__ import annotations

import json

import pandas as pd

from src.reference.table_resolution import (
    ambiguity_audit,
    ambiguous_same_date,
    candidate_content,
    newest,
    newest_by_time,
)


def _result(tmp_path, dirname, stamp):
    d = tmp_path / dirname
    d.mkdir()
    (d / "thing.parquet").write_bytes(b"")
    (d / "thing.meta.json").write_text(json.dumps({"utc_timestamp": stamp}))
    return d


def _pair(tmp_path):
    # e1b1cc4 sorts after da1f46b but was written earlier.
    _result(tmp_path, "2026-09-07_e1b1cc4", "2026-09-07T19:49:01+00:00")
    _result(tmp_path, "2026-09-07_da1f46b", "2026-09-07T22:17:00+00:00")


def test_name_order_and_time_order_disagree(tmp_path):
    _pair(tmp_path)
    assert newest(tmp_path, "thing").parent.name == "2026-09-07_e1b1cc4"
    assert newest_by_time(tmp_path, "thing").parent.name == "2026-09-07_da1f46b"


def test_ambiguous_same_date_reports_the_pair(tmp_path):
    _pair(tmp_path)
    rows = ambiguous_same_date(tmp_path)
    assert len(rows) == 1
    assert rows[0]["table"] == "thing"
    assert rows[0]["name_order_picks"] == "2026-09-07_e1b1cc4"
    assert rows[0]["time_order_picks"] == "2026-09-07_da1f46b"


def test_a_single_run_per_date_is_not_ambiguous(tmp_path):
    _result(tmp_path, "2026-09-06_aaaaaaa", "2026-09-06T01:00:00+00:00")
    _result(tmp_path, "2026-09-07_bbbbbbb", "2026-09-07T01:00:00+00:00")
    assert ambiguous_same_date(tmp_path) == []


def test_resolvers_return_none_when_absent(tmp_path):
    assert newest(tmp_path, "thing") is None
    assert newest_by_time(tmp_path, "thing") is None


def test_newest_by_time_falls_back_to_mtime_without_sidecars(tmp_path):
    import os
    import time

    older = tmp_path / "2026-09-07_aaaaaaa"
    newer = tmp_path / "2026-01-01_zzzzzzz"
    for d in (older, newer):
        d.mkdir()
        (d / "thing.parquet").write_bytes(b"")
    os.utime(older / "thing.parquet", (1_000_000, 1_000_000))
    now = time.time()
    os.utime(newer / "thing.parquet", (now, now))
    # No sidecars anywhere, so mtime decides and the newer file wins despite its
    # name sorting first.
    assert newest_by_time(tmp_path, "thing").parent.name == "2026-01-01_zzzzzzz"


def _result_with_table(tmp_path, dirname, stamp, values):
    d = tmp_path / dirname
    d.mkdir()
    pd.DataFrame({"a": values, "b": [1.0] * len(values)}).to_parquet(
        d / "thing.parquet"
    )
    (d / "thing.meta.json").write_text(json.dumps({"utc_timestamp": stamp}))
    return d


def test_candidate_content_separates_harmless_from_consequential(tmp_path):
    _result_with_table(tmp_path, "2026-09-07_e1b1cc4", "2026-09-07T19:49:01+00:00",
                       [1, 2, 3])
    _result_with_table(tmp_path, "2026-09-07_da1f46b", "2026-09-07T22:17:00+00:00",
                       [1, 2, 3])
    same = candidate_content(
        tmp_path, "thing", "2026-09-07_e1b1cc4", "2026-09-07_da1f46b"
    )
    assert same["frames_identical"]
    assert same["same_columns"]
    assert same["name_pick_rows"] == same["time_pick_rows"] == 3


def test_candidate_content_flags_a_superseded_run(tmp_path):
    _result_with_table(tmp_path, "2026-09-07_e1b1cc4", "2026-09-07T19:49:01+00:00",
                       [1, 2, 3])
    _result_with_table(tmp_path, "2026-09-07_da1f46b", "2026-09-07T22:17:00+00:00",
                       [9, 9, 9])
    out = candidate_content(
        tmp_path, "thing", "2026-09-07_e1b1cc4", "2026-09-07_da1f46b"
    )
    assert not out["frames_identical"]


def test_ambiguity_audit_marks_consequence(tmp_path):
    _result_with_table(tmp_path, "2026-09-07_e1b1cc4", "2026-09-07T19:49:01+00:00",
                       [1, 2, 3])
    _result_with_table(tmp_path, "2026-09-07_da1f46b", "2026-09-07T22:17:00+00:00",
                       [9, 9, 9])
    rows = ambiguity_audit(tmp_path)
    assert len(rows) == 1
    assert rows[0]["disposition"] == "differs"
