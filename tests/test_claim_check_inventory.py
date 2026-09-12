"""The claim-check inventory, and the inputs that force each of its guards to fail.

Written to `tests/test_checks_can_fail.py`'s rule rather than beside it: that
file is the paper's Appendix A material and is W2-adjacent, and this is week-1
work on a separate branch. Every guard in `src/reference/claim_check.py` has a
failing input here. If a guard is added there without one here, it is untested.

The pass-path tests are the ordinary half. The `_forces_` tests are the half
that can detect a validator returning clean on everything.
"""

from __future__ import annotations

import copy

import pytest

from src.reference.claim_check import (
    JUDGEMENT_FIELDS,
    REQUIRED_FIELDS,
    LedgerError,
    agreement,
    assert_reconciles_or_says_why,
    build_inventory,
    load_ledger,
    reconciliation,
)


@pytest.fixture(scope="module")
def doc() -> dict:
    return load_ledger()


@pytest.fixture
def mutable(doc) -> dict:
    return copy.deepcopy(doc)


# ---------------------------------------------------------------------------
# The pass path
# ---------------------------------------------------------------------------


def test_the_committed_ledger_loads_and_validates(doc):
    frame = build_inventory(doc)
    assert len(frame) == len(doc["entries"])
    for field in REQUIRED_FIELDS:
        assert field in frame.columns
    assert frame["id"].is_unique


def test_every_entry_names_a_source_anchor_in_the_repository(doc):
    """A classification with no traceable source is an opinion about nothing."""
    for entry in doc["entries"]:
        anchor = entry["source_anchor"]
        assert anchor and "docs/" in anchor, entry["id"]


def test_the_classification_vocabulary_is_fixed_in_the_file_not_the_code(doc):
    """The criteria must be readable beside the judgements they produced."""
    assert set(doc["classes"]) == {
        "logical_impossibility",
        "low_power",
        "implementation_error",
        "provenance_reporting_error",
    }
    for name, definition in doc["classes"].items():
        assert len(definition) > 40, name


# ---------------------------------------------------------------------------
# The reconciliation, which is the ledger's own missing check
# ---------------------------------------------------------------------------


def test_the_enumerated_count_does_not_match_the_asserted_one(doc):
    """Pinned deliberately. If someone repairs HANDOFF §3, this test tells them
    to update the note rather than letting the two drift apart again."""
    rec = reconciliation(doc)
    assert rec.enumerated == 23
    assert rec.asserted == 21
    assert not rec.reconciles
    assert rec.note, "a documented mismatch is allowed; a silent one is not"


def test_a_silent_mismatch_forces_the_reconciliation_guard_to_fail(mutable):
    mutable["meta"]["reconciliation_note"] = None
    with pytest.raises(LedgerError, match="do not silently adopt"):
        assert_reconciles_or_says_why(mutable)


def test_an_agreeing_count_passes_without_a_note(mutable):
    """The guard must not fire on the case it is not about."""
    mutable["meta"]["asserted_count"] = len(mutable["entries"])
    mutable["meta"]["reconciliation_note"] = None
    rec = assert_reconciles_or_says_why(mutable)
    assert rec.reconciles


# ---------------------------------------------------------------------------
# One forcing input per validation rule
# ---------------------------------------------------------------------------


def test_a_missing_judgement_field_forces_a_refusal(mutable):
    del mutable["entries"][0]["classification"]
    with pytest.raises(LedgerError, match="not a null one"):
        build_inventory(mutable)


def test_a_duplicate_id_forces_a_refusal(mutable):
    mutable["entries"][1]["id"] = mutable["entries"][0]["id"]
    with pytest.raises(LedgerError, match="duplicate entry id"):
        build_inventory(mutable)


def test_a_classification_outside_the_vocabulary_forces_a_refusal(mutable):
    mutable["entries"][0]["classification"] = "probably_fine"
    with pytest.raises(LedgerError, match="not .*one of the vocabulary"):
        build_inventory(mutable)


def test_an_ambiguity_that_repeats_the_primary_class_forces_a_refusal(mutable):
    entry = mutable["entries"][0]
    entry["ambiguous_with"] = entry["classification"]
    with pytest.raises(LedgerError, match="records no ambiguity at all"):
        build_inventory(mutable)


def test_a_disputed_verdict_without_a_note_forces_a_refusal(mutable):
    entry = mutable["entries"][1]
    entry["verdict"] = "disputed"
    entry["dispute_note"] = None
    with pytest.raises(LedgerError, match="unstated dispute is not a dispute"):
        build_inventory(mutable)


def test_an_unknown_discovery_route_forces_a_refusal(mutable):
    mutable["entries"][0]["found_by"] = "vibes"
    with pytest.raises(LedgerError, match="found_by"):
        build_inventory(mutable)


def test_two_entries_claiming_one_ordinal_force_a_refusal(mutable):
    """The ordinals are the ledger's only index across five document sections."""
    with_ordinals = [e for e in mutable["entries"] if e["ordinal_claimed"] is not None]
    with_ordinals[1]["ordinal_claimed"] = with_ordinals[0]["ordinal_claimed"]
    with pytest.raises(LedgerError, match="same ordinal"):
        build_inventory(mutable)


def test_an_empty_ledger_forces_a_refusal(mutable):
    mutable["entries"] = []
    with pytest.raises(LedgerError, match="enumerates nothing"):
        build_inventory(mutable)


# ---------------------------------------------------------------------------
# The agreement path, which is where this module could have added entry 24
# ---------------------------------------------------------------------------


def test_one_rater_reports_no_agreement_rather_than_perfect_agreement(doc):
    """A rater compared with themselves is a check that cannot fail."""
    assert agreement(build_inventory(doc), None) is None


def test_two_raters_disagreeing_is_visible_rather_than_averaged(mutable, doc):
    first = build_inventory(doc)
    other = mutable["entries"][0]
    other["classification"] = (
        "low_power" if other["classification"] != "low_power" else "implementation_error"
    )
    mutable["meta"]["rater"] = "rater_b"
    second = build_inventory(mutable)
    pairs = agreement(first, second)
    assert pairs is not None
    assert not bool(pairs.loc[pairs["id"] == other["id"], "classification_agrees"].iloc[0])
    assert int(pairs["classification_agrees"].sum()) == len(first) - 1


def test_a_second_rater_supplies_only_judgement_fields():
    """Transcription is not re-rated; disagreement must be about the judgement."""
    for field in JUDGEMENT_FIELDS:
        assert field in REQUIRED_FIELDS
    assert "source_anchor" not in JUDGEMENT_FIELDS
    assert "why_it_could_not_fail" not in JUDGEMENT_FIELDS
