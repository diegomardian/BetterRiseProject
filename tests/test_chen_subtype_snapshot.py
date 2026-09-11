"""The snapshot is a snapshot: fixed endpoints, checksummed, and label-only."""

from __future__ import annotations

import hashlib
import json

import pytest

from src.reference.jobs.chen_subtype_snapshot import (
    CLAIM_PROVENANCE,
    LABEL_PROVENANCE,
    STUDY_ID,
    SnapshotError,
    manifest_rows,
    snapshot_table,
)


def _responses() -> dict[str, bytes]:
    return {
        "study.json": json.dumps({"allSampleCount": 61}).encode(),
        "clinical_sample.json": json.dumps(
            [
                {"sampleId": "HTA11_1_2000001011", "patientId": "HTA11_1"},
                {"sampleId": "HTA11_1_2000001011", "patientId": "HTA11_1"},
                {"sampleId": "HTA11_2_2000001011", "patientId": "HTA11_2"},
            ]
        ).encode(),
        "mutations.json": json.dumps([{"sampleId": "HTA11_1_2000001011"}]).encode(),
    }


def test_the_checksum_is_of_the_bytes_that_were_written():
    responses = _responses()
    frame = snapshot_table(responses, written=True)
    row = frame.set_index("file").loc["study.json"]
    assert row["sha256"] == hashlib.sha256(responses["study.json"]).hexdigest()
    assert row["bytes"] == len(responses["study.json"])


def test_distinct_counts_not_row_counts_are_recorded():
    """Three clinical rows over two samples is two samples, not three."""
    frame = snapshot_table(_responses(), written=True).set_index("file")
    assert frame.loc["clinical_sample.json", "n_records"] == 3
    assert frame.loc["clinical_sample.json", "n_samples"] == 2
    assert frame.loc["clinical_sample.json", "n_patients"] == 2


def test_an_empty_response_is_refused_rather_than_pinned():
    responses = _responses()
    responses["study.json"] = b""
    with pytest.raises(SnapshotError):
        snapshot_table(responses, written=True)


def test_manifest_rows_carry_the_accession_and_a_source_for_every_file():
    frame = snapshot_table(_responses(), written=True)
    rows = manifest_rows(frame, downloaded_on="2026-09-10")
    assert len(rows) == len(frame)
    assert (rows["accession"] == STUDY_ID).all()
    assert rows["source_url"].notna().all()
    assert not rows["source_url"].astype(str).str.strip().eq("").any()


def test_the_snapshot_declares_morphology_labels_against_a_transcript_claim():
    """Gate 3 freezes the population, so invariant 11 binds here."""
    assert LABEL_PROVENANCE.modality == "morphology"
    assert CLAIM_PROVENANCE.modality == "transcript"
    assert LABEL_PROVENANCE.genes == ()
