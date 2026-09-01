"""`docs/open_decisions.md` numbers have to identify a decision.

Ten numbers currently mean two different things each. That is not a cosmetic
problem: "see #13" appears in commit messages and notes across all four
workstreams, and on `main` it resolves to either *"W1 and W4 label cells
differently"* or *"Bulk GUCA2A is continuous, not bimodal"* depending on which
one the reader scrolls to first.

WHY THIS IS A RATCHET AND NOT A FIX
-----------------------------------
Renumbering rewrites cross-references in every workstream's notes and commit
messages. That is a team decision, not a tidy-up one person does inside their
own PR, so this file does not renumber anything and does not assert the
collisions are gone.

What it does assert is that **no new one appears**. `KNOWN_COLLISIONS` is a
baseline, not an endorsement: the check is a subset test, so fixing a collision
always passes and adding one always fails. The collision count went 6 → 9 → 10
in eight days, twice while a PR proposing a number was open, which is why
stopping the bleeding is worth doing before the renumbering is agreed.

The durable fix is proposed at the top of `docs/open_decisions.md`:
per-workstream number ranges, so two people picking concurrently cannot collide
without coordinating. This test is what makes that proposal enforceable once
someone adopts it.
"""

from __future__ import annotations

import re

import pytest

from src.common.paths import REPO_ROOT

DECISIONS_PATH = REPO_ROOT / "docs" / "open_decisions.md"

#: `## <number> · <title>` — the format every decision heading uses.
HEADING = re.compile(r"^## (\d+) · (.+)$")

#: Numbers that already mean two things on `main`, recorded 2026-08-30.
#: A BASELINE, NOT A TARGET. Shrinking it is always allowed; growing it is not.
#: #23 leaves this set when PR #61 renumbers its W3 entry to #26.
KNOWN_COLLISIONS = frozenset({9, 10, 11, 12, 13, 14, 15, 16, 17, 23})


def _headings() -> list[tuple[int, int, str]]:
    """(line number, decision number, title) for every decision heading."""
    out = []
    for line_no, line in enumerate(
        DECISIONS_PATH.read_text(encoding="utf-8").split("\n"), start=1
    ):
        match = HEADING.match(line)
        if match:
            out.append((line_no, int(match.group(1)), match.group(2).strip()))
    return out


def _collisions() -> dict[int, list[str]]:
    by_number: dict[int, list[str]] = {}
    for _line, number, title in _headings():
        by_number.setdefault(number, []).append(title)
    return {n: t for n, t in by_number.items() if len(t) > 1}


def test_the_decisions_file_has_decisions_in_it():
    """Guard against the parser silently matching nothing and passing."""
    assert len(_headings()) >= 20


def test_no_new_decision_number_collisions():
    """The ratchet. Fixing a collision passes; adding one fails.

    If this fails on a number you just used, someone else already has it —
    take the next free number in your workstream's range (see the top of
    docs/open_decisions.md) rather than the next number overall.
    """
    collisions = set(_collisions())
    new = collisions - KNOWN_COLLISIONS
    assert not new, (
        f"decision number(s) {sorted(new)} now mean more than one thing. "
        f"A number that does not identify a decision makes every 'see #N' in "
        f"the repo ambiguous. Pick an unused number — and if you are adding "
        f"decisions regularly, adopt the per-workstream ranges proposed at the "
        f"top of docs/open_decisions.md."
    )


def test_the_known_collisions_are_still_real():
    """Keeps the baseline honest. If a number is fixed, it should leave the
    list rather than sitting there implying a problem that no longer exists."""
    stale = KNOWN_COLLISIONS - set(_collisions())
    if stale:
        pytest.skip(
            f"#{sorted(stale)} no longer collide — remove them from "
            f"KNOWN_COLLISIONS. Skipping rather than failing so that fixing a "
            f"collision never turns main red."
        )


def test_every_decision_heading_carries_a_title():
    for line_no, number, title in _headings():
        assert title, f"decision #{number} at line {line_no} has an empty title"


def test_the_disambiguation_table_covers_every_collision():
    """While the collisions exist, a reader has to be able to resolve `#13`.

    The table at the top of the file is what makes that possible without
    renumbering, so it must not drift from the actual collisions.
    """
    text = DECISIONS_PATH.read_text(encoding="utf-8")
    _, _, after = text.partition("WHICH #N DO YOU MEAN")
    assert after, "the disambiguation table is missing from docs/open_decisions.md"
    table = after.partition("<!-- END DISAMBIGUATION -->")[0]
    for number in sorted(_collisions()):
        assert f"#{number}" in table, (
            f"#{number} collides but is not in the disambiguation table"
        )
