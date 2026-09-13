"""Forcing inputs for the shared result-table resolver.

The audit found 23 table names where name order and time order disagree within a
date. `newest` is kept because callers rely on it, but it must not be mistaken
for `newest_by_time`, and `ambiguous_same_date` must surface the disagreement
rather than resolve it silently.
"""

from __future__ import annotations

import json

from src.reference.table_resolution import (
    ambiguous_same_date,
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
