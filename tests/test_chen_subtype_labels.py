"""Both label families pass invariant 11, and for different reasons."""

from __future__ import annotations

import pytest

from src.common.label_provenance import CircularClaimError, Measurement, check_no_circular_claim
from src.reference.chen_subtype_labels import (
    CLAIM,
    FAMILIES,
    FAMILY_G,
    FAMILY_P,
    validate_families,
)


def test_both_families_pass_the_guard_with_no_overlapping_genes():
    assert validate_families() == {"P": (), "G": ()}


def test_the_two_families_differ_in_modality_from_the_claim_and_each_other():
    """Two labels that shared a modality would share its failure modes."""
    assert FAMILY_P.modality != CLAIM.modality
    assert FAMILY_G.modality != CLAIM.modality
    assert FAMILY_P.modality != FAMILY_G.modality


def test_every_family_is_declared_against_the_same_fixed_claim():
    """A family declared against a different endpoint is not a second check."""
    assert len(FAMILIES) == 2
    assert all(m is not CLAIM for m in FAMILIES.values())


def test_the_guard_refuses_a_same_modality_label_once_genes_overlap():
    claim = Measurement(modality="transcript", assay="snRNA decomposition",
                        genes=("GUCA2A",))
    labels = Measurement(modality="transcript", assay="snRNA clustering",
                         genes=("GUCA2A",))
    with pytest.raises(CircularClaimError):
        check_no_circular_claim(labels=labels, claim=claim)


def test_the_guard_cannot_catch_a_gene_free_transcript_derived_label():
    """Why gate 2 had to read the source methods rather than trust the guard.

    Chen et al. reclassified histologically unconfirmed specimens from the
    transcriptomes. Such a label is circular, but it declares no genes, and the
    guard short-circuits on an empty overlap — so it would pass. The protection
    is that those specimens carry POLYP_TYPE 'Unknown' in the deposit and are
    excluded by rule 3, not that the guard would refuse them.
    """
    reclassified = Measurement(modality="transcript", assay=CLAIM.assay, genes=())
    assert check_no_circular_claim(labels=reclassified, claim=CLAIM) == ()


def test_family_g_names_its_own_genes_and_they_are_not_panel_genes():
    """Invariant 2: a panel gene must never appear in a label."""
    from src.common.panel import panel_genes  # noqa: PLC0415

    assert FAMILY_G.genes == ("APC", "BRAF")
    assert not set(FAMILY_G.genes) & set(panel_genes())
