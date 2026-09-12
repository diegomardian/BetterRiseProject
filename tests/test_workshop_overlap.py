"""The overlap matrix, and the inputs that force its guards to fail."""

from __future__ import annotations

import copy

import pytest

from src.reference.workshop_overlap import (
    OverlapError,
    build_overlap,
    load_overlap,
    overlap_summary,
)


@pytest.fixture(scope="module")
def doc() -> dict:
    return load_overlap()


@pytest.fixture
def mutable(doc) -> dict:
    return copy.deepcopy(doc)


def test_the_committed_matrix_loads(doc):
    frame = build_overlap(doc)
    assert len(frame) == len(doc["rows"])
    assert frame["id"].is_unique


def test_most_of_the_overlap_is_declared_rather_than_minimised(doc):
    """If this ever reports that nearly everything is new, something is wrong
    with the declaration, not with the paper."""
    summary = overlap_summary(build_overlap(doc))
    assert summary["n_already_published"] >= 4
    assert summary["novelty_share_of_overlap_rows"] < 0.5


def test_the_mathematical_audit_found_edits_the_paper_owes(doc):
    frame = build_overlap(doc)
    maths = frame[frame["kind"] == "math_claim"]
    assert len(maths) >= 3
    assert maths["needs_paper_edit"].sum() >= 3


def test_a_new_claim_with_no_evidence_forces_a_refusal(mutable):
    row = next(r for r in mutable["rows"] if r["status"] == "new")
    row["evidence"] = None
    with pytest.raises(OverlapError, match="names no evidence"):
        build_overlap(mutable)


def test_an_overlapping_claim_with_no_location_forces_a_refusal(mutable):
    row = next(r for r in mutable["rows"] if r["status"] == "already_published")
    row["wmhs_location"] = None
    with pytest.raises(OverlapError, match="not a declaration"):
        build_overlap(mutable)


def test_an_unknown_status_forces_a_refusal(mutable):
    mutable["rows"][0]["status"] = "sort of new"
    with pytest.raises(OverlapError, match="status"):
        build_overlap(mutable)


def test_a_duplicate_row_id_forces_a_refusal(mutable):
    mutable["rows"][1]["id"] = mutable["rows"][0]["id"]
    with pytest.raises(OverlapError, match="duplicate row id"):
        build_overlap(mutable)


def test_an_empty_matrix_forces_a_refusal(mutable):
    mutable["rows"] = []
    with pytest.raises(OverlapError, match="declares nothing"):
        build_overlap(mutable)
