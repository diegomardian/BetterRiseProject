"""Resolving a result table to the run that should be read.

WHY THIS IS SHARED, AND WHY IT HAS TWO FUNCTIONS. A result directory is named
``{date}_{sha7}``. Sorting those names orders runs by **commit sha, not by
time**, and within a single day that is arbitrary: the repo committed six
Crowell multisection runs on 2026-09-07, and the lexicographically last
(``e1b1cc4``) is the second-earliest, not the canonical ``da1f46b`` written
2h28m later. Every job that does ``sorted(RESULTS_DIR.glob(...))[-1]`` carries
this silently.

Two functions, deliberately not one:

``newest``          the existing name-order resolver. Kept because many callers
                    want it and because changing it under them would be a
                    silent behavioural change, but it is only safe when a table
                    has at most one run per date.
``newest_by_time``  resolve by the sidecar's ``utc_timestamp``, falling back to
                    mtime only when no sidecar is readable. This is the one a
                    reproduction or a comparison should use.

The audit found 23 table names with ambiguous same-date runs where the two
functions disagree. This module does not pretend one rule is universally right;
it makes the choice explicit at each call site.
"""

from __future__ import annotations

import json
from pathlib import Path


def newest(results_dir: Path, name: str) -> Path | None:
    """Newest ``name`` table by directory name. Ambiguous within a date."""
    matches = sorted(results_dir.glob(f"*/{name}.parquet"))
    return matches[-1] if matches else None


def newest_by_time(results_dir: Path, name: str) -> Path | None:
    """Most recently *written* ``name`` table, by its sidecar timestamp.

    Falls back to mtime only when a sidecar is absent or unreadable, so a table
    without provenance is not silently dropped.
    """
    matches = sorted(results_dir.glob(f"*/{name}.parquet"))
    if not matches:
        return None

    def written_at(path: Path) -> tuple[str, float]:
        sidecar = path.with_suffix("").with_suffix(".meta.json")
        if not sidecar.exists():
            sidecar = path.parent / f"{name}.meta.json"
        stamp = ""
        try:
            stamp = str(json.loads(sidecar.read_text()).get("utc_timestamp") or "")
        except (OSError, ValueError):
            stamp = ""
        return stamp, path.stat().st_mtime

    return max(matches, key=written_at)


def ambiguous_same_date(results_dir: Path) -> list[dict[str, object]]:
    """Every ``(table, date)`` where name order and time order disagree.

    The audit's finding, as a function rather than a paragraph: a resolver that
    cannot say which of two runs it read is a check that cannot fail one level
    out. Returns one row per disagreement with both candidate directories, so a
    caller can report the ambiguity instead of resolving it silently.
    """
    from collections import defaultdict

    by_name: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for meta in results_dir.glob("*/*.meta.json"):
        name = meta.name[: -len(".meta.json")]
        try:
            stamp = str(json.loads(meta.read_text()).get("utc_timestamp") or "")
        except (OSError, ValueError):
            stamp = ""
        by_name[name].append((meta.parent.name, stamp))

    rows: list[dict[str, object]] = []
    for name, runs in by_name.items():
        by_date: dict[str, list[tuple[str, str]]] = defaultdict(list)
        for directory, stamp in runs:
            by_date[directory.split("_")[0]].append((directory, stamp))
        for date, group in by_date.items():
            if len(group) < 2:
                continue
            lex_last = max(directory for directory, _ in group)
            time_last = max(group, key=lambda r: (r[1], r[0]))[0]
            if lex_last != time_last:
                rows.append({
                    "table": name, "date": date, "n_same_date_runs": len(group),
                    "name_order_picks": lex_last,
                    "time_order_picks": time_last,
                })
    return sorted(rows, key=lambda r: (str(r["table"]), str(r["date"])))
