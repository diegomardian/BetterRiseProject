"""The perturbed arm is a genotype label, and an expression label is refused."""

from __future__ import annotations

import pytest

from src.common.label_provenance import CircularClaimError, Measurement, check_no_circular_claim
from src.reference.replogle_labels import CLAIM, GUIDE_ASSIGNMENT, validate


def test_the_guide_label_passes_and_overlaps_no_gene():
    assert validate() == ()


def test_the_label_is_a_dna_assignment_not_a_transcript_reading():
    assert GUIDE_ASSIGNMENT.modality == "genotype"
    assert CLAIM.modality == "transcript"
    assert GUIDE_ASSIGNMENT.genes == ()


def test_defining_the_arm_by_the_targeted_transcript_is_refused():
    """The circularity §2 exists to forbid, and it would manufacture power."""
    by_expression = Measurement(
        modality="transcript",
        assay="Replogle 2022 Perturb-seq single-cell UMI counts",
        genes=("MLH1",),
    )
    claim_on_same = Measurement(
        modality="transcript",
        assay="Replogle 2022 Perturb-seq single-cell UMI counts",
        genes=("MLH1",),
    )
    with pytest.raises(CircularClaimError):
        check_no_circular_claim(labels=by_expression, claim=claim_on_same)
