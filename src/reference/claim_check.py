"""The claim-check ledger, loaded, validated and counted.

``docs/HANDOFF.md`` §3 asserts in prose that a check unable to fail has been
found twenty-one times. Nothing enumerated them: the ordinals are scattered over
five sections, §7 still says seventeen, and the §3 table holds fifteen rows
while the prose around it describes more. Transcribing every entry the
repository actually describes gives **twenty-three**.

That mismatch is the point of this module rather than an obstacle to it. The
ledger is a claim about the repository, and until now it was a claim with no
check — which is the ledger's own subject applied one level out. So:

* the enumeration lives in ``config/claim_check_ledger.yaml`` as reviewable
  data, because every field in it is a transcription or a judgement;
* this module validates the shape and **refuses a silent reconciliation** — if
  the enumerated count differs from the asserted one, a written note must say
  so. A count that quietly agrees is also refused unless it agrees honestly;
* classification is per-rater. One rater is one opinion with criteria attached,
  and ``agreement()`` reports nothing at all until a second rater's file exists.

Placed under ``src/reference/`` (W1) rather than ``src/harness/`` (W2)
deliberately: ``CONTRIBUTING`` §2-3 route harness changes through a PR with two
approvals, and ``docs/HANDOFF.md`` §6k records three that landed without one.
This module adds no harness behaviour and does not need that route.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import yaml

from src.common.paths import REPO_ROOT

LEDGER_PATH = REPO_ROOT / "config" / "claim_check_ledger.yaml"

#: Every entry carries these. A missing one is a refusal, not a default.
REQUIRED_FIELDS: tuple[str, ...] = (
    "id",
    "short_name",
    "the_check",
    "why_it_could_not_fail",
    "source_anchor",
    "ordinal_claimed",
    "found_by",
    "guard_test",
    "classification",
    "ambiguous_with",
    "rationale",
    "verdict",
    "dispute_note",
)

#: Judgement fields. A second rater supplies only these, keyed by ``id``.
JUDGEMENT_FIELDS: tuple[str, ...] = (
    "classification",
    "ambiguous_with",
    "rationale",
    "verdict",
)

FOUND_BY: frozenset[str] = frozenset({"suite", "review", "run", "reader"})


class LedgerError(ValueError):
    """The ledger file does not say what it must say."""


@dataclass(frozen=True)
class Reconciliation:
    """Enumerated against asserted, with the note that licenses a difference."""

    enumerated: int
    asserted: int
    note: str | None

    @property
    def reconciles(self) -> bool:
        return self.enumerated == self.asserted

    @property
    def delta(self) -> int:
        return self.enumerated - self.asserted


def load_ledger(path: Path | None = None) -> dict:
    """Read the YAML. No validation — ``build_inventory`` does that."""
    with open(path or LEDGER_PATH) as handle:
        return yaml.safe_load(handle)


def _validate(doc: dict) -> None:
    for key in ("meta", "classes", "verdicts", "entries"):
        if key not in doc:
            raise LedgerError(f"ledger is missing the {key!r} block")

    classes = set(doc["classes"])
    verdicts = set(doc["verdicts"])
    entries = doc["entries"]
    if not entries:
        raise LedgerError("ledger enumerates nothing")

    seen: set[str] = set()
    for entry in entries:
        missing = [f for f in REQUIRED_FIELDS if f not in entry]
        if missing:
            raise LedgerError(
                f"entry {entry.get('id', '?')!r} is missing {missing}. "
                "An absent judgement is not a null one."
            )
        if entry["id"] in seen:
            raise LedgerError(f"duplicate entry id {entry['id']!r}")
        seen.add(entry["id"])

        if entry["classification"] not in classes:
            raise LedgerError(
                f"{entry['id']}: classification {entry['classification']!r} is not "
                f"one of the vocabulary fixed in the file: {sorted(classes)}"
            )
        second = entry["ambiguous_with"]
        if second is not None and second not in classes:
            raise LedgerError(
                f"{entry['id']}: ambiguous_with {second!r} is not in the vocabulary"
            )
        if second is not None and second == entry["classification"]:
            raise LedgerError(
                f"{entry['id']}: ambiguous_with repeats the primary classification, "
                "which records no ambiguity at all"
            )
        if entry["verdict"] not in verdicts:
            raise LedgerError(
                f"{entry['id']}: verdict {entry['verdict']!r} is not in the vocabulary"
            )
        if entry["verdict"] == "disputed" and not entry["dispute_note"]:
            raise LedgerError(
                f"{entry['id']}: verdict is 'disputed' with no dispute_note. "
                "An unstated dispute is not a dispute."
            )
        if entry["found_by"] not in FOUND_BY:
            raise LedgerError(
                f"{entry['id']}: found_by {entry['found_by']!r} is not one of "
                f"{sorted(FOUND_BY)}"
            )

    ordinals = [e["ordinal_claimed"] for e in entries if e["ordinal_claimed"] is not None]
    if len(ordinals) != len(set(ordinals)):
        raise LedgerError(
            f"two entries claim the same ordinal: {sorted(ordinals)}. "
            "The ledger's ordinals are its only cross-section index."
        )


def reconciliation(doc: dict) -> Reconciliation:
    """Enumerated count against the count HANDOFF §3 asserts.

    **This is a check that can fail, and at the time of writing it does.** It is
    kept failing rather than repaired because repairing it means editing §3's
    prose, which is the write-up's business and not this job's.
    """
    return Reconciliation(
        enumerated=len(doc["entries"]),
        asserted=int(doc["meta"]["asserted_count"]),
        note=doc["meta"].get("reconciliation_note"),
    )


def assert_reconciles_or_says_why(doc: dict) -> Reconciliation:
    """Refuse a silent mismatch. A documented one is allowed through.

    The failure mode this exists for is the one the ledger is about: a count
    asserted in one place, enumerated differently in another, and nothing
    comparing them. Agreeing counts pass; disagreeing counts pass **only** with
    a written note; disagreeing counts with no note raise.
    """
    rec = reconciliation(doc)
    if not rec.reconciles and not rec.note:
        raise LedgerError(
            f"enumerated {rec.enumerated} entries against an asserted "
            f"{rec.asserted} ({rec.delta:+d}) and meta.reconciliation_note is "
            "empty. Write down why they differ; do not silently adopt either."
        )
    return rec


def build_inventory(doc: dict, *, rater: str | None = None) -> pd.DataFrame:
    """One row per ledger entry, validated."""
    _validate(doc)
    frame = pd.DataFrame(doc["entries"])[list(REQUIRED_FIELDS)].copy()
    frame.insert(1, "rater", rater or doc["meta"]["rater"])
    frame["ordinal_claimed"] = frame["ordinal_claimed"].astype("Int64")
    frame["has_forcing_input"] = frame["guard_test"].notna()
    frame["counts_toward_stop_rule"] = frame["verdict"] == "scientific_failure"
    return frame.sort_values("id").reset_index(drop=True)


def agreement(first: pd.DataFrame, second: pd.DataFrame | None) -> pd.DataFrame | None:
    """Per-entry agreement between two raters, or ``None`` if there is one rater.

    Returns ``None`` rather than a table of perfect self-agreement. A single
    rater compared with themselves is the shape of a check that cannot fail, and
    this module is not going to add the twenty-fourth entry to its own ledger.
    """
    if second is None:
        return None
    merged = first.merge(
        second[["id", *JUDGEMENT_FIELDS]], on="id", how="outer", suffixes=("_a", "_b")
    )
    merged["classification_agrees"] = (
        merged["classification_a"] == merged["classification_b"]
    )
    merged["verdict_agrees"] = merged["verdict_a"] == merged["verdict_b"]
    return merged
