"""Invariant 11's pass path: what a legitimate declaration looks like.

The inputs that must be REFUSED live in ``tests/test_checks_can_fail.py``,
beside every other guard's forcing input. This file is the other half — the
declarations the guard has to let through, including the one the proposal named:
an independently measured protein endpoint against a transcript label.
"""

from __future__ import annotations

import pytest

from src.common.label_provenance import (
    MODALITIES,
    IndependentEndpoint,
    Measurement,
    ProvenanceDeclarationError,
    check_no_circular_claim,
    check_spec_declares_provenance,
    provenance_meta,
)

# The two live analyses this invariant was written from. Both already satisfy it
# in practice -- Becker's RESULT and the Crowell multisection Amendment 1 -- and
# encoding them here is what turns "we were careful" into something that fails
# if someone stops being careful.
CROWELL_LABELS = Measurement(
    modality="morphology",
    assay="Crowell histopathology domain annotation",
    genes=(),
)
CROWELL_CLAIM = Measurement(
    modality="transcript",
    assay="10x Xenium WTx",
    genes=("GUCA2A", "MS4A12", "CDX2", "EPCAM", "KRT8"),
)
BECKER_LABELS = Measurement(
    modality="transcript",
    assay="Becker snRNA-seq mature-colonocyte marker call",
    genes=("CA1", "CA2", "AQP8", "SLC26A3", "KRT20", "CEACAM7"),
)
BECKER_CLAIM = Measurement(
    modality="transcript",
    assay="Becker snRNA-seq",
    genes=("GUCA2A", "MS4A12"),
)


def test_a_morphological_label_may_claim_any_transcript():
    assert check_no_circular_claim(labels=CROWELL_LABELS, claim=CROWELL_CLAIM) == ()


def test_a_marker_label_disjoint_from_the_targets_passes():
    """Becker's mature label and its endpoint share no gene. Invariant 2's case."""
    assert check_no_circular_claim(labels=BECKER_LABELS, claim=BECKER_CLAIM) == ()


def test_the_proposal_s_passing_fixture_an_independently_measured_endpoint():
    """CDX2 defines the label in transcript; CDX2 protein is measured by MxIF."""
    overlap = check_no_circular_claim(
        labels=Measurement(modality="transcript", assay="Chen_2021 scRNA-seq",
                           genes=("CDX2",)),
        claim=Measurement(modality="protein", assay="MxIF CDX2", genes=("CDX2",)),
        independent_endpoint=IndependentEndpoint(
            assay="MxIF CDX2",
            reason="protein is measured on the section, not read off the "
                   "transcripts that assigned the label",
        ),
    )
    assert overlap == ("CDX2",)


def test_case_does_not_decide_whether_two_genes_are_the_same_gene():
    with pytest.raises(Exception) as excinfo:
        check_no_circular_claim(
            labels=Measurement(modality="transcript", assay="A", genes=("Cdx2",)),
            claim=Measurement(modality="transcript", assay="B", genes=("CDX2",)),
        )
    assert "CDX2" in str(excinfo.value)


def test_declaring_no_genes_is_a_statement_and_is_accepted():
    spec = {
        "label_provenance": {"modality": "sample_annotation",
                             "assay": "GEO series matrix arm", "genes": []},
        "claim_provenance": {"modality": "transcript", "assay": "snRNA-seq",
                             "genes": ["GUCA2A"]},
    }
    assert check_spec_declares_provenance(spec) == ()


def test_the_sidecar_record_names_the_population_definition():
    meta = provenance_meta(CROWELL_LABELS, CROWELL_CLAIM)
    assert meta["label_provenance"]["assay"] == CROWELL_LABELS.assay
    assert meta["label_provenance"]["genes"] == ""
    assert meta["claim_provenance"]["genes"].startswith("GUCA2A;")
    assert meta["overlapping_genes"] == ""
    assert meta["independent_endpoint"] is None


def test_a_modality_outside_the_vocabulary_is_refused_rather_than_invented():
    with pytest.raises(ProvenanceDeclarationError, match="not one of"):
        Measurement(modality="vibes", assay="an instrument", genes=())
    assert "transcript" in MODALITIES and "protein" in MODALITIES
