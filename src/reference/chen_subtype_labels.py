"""Gate 5: the two label families of `prereg_chen_lesion_subtype`, declared once.

Invariant 11 binds per analysis, and this pre-registration runs two of them
against the same decomposition. Declaring each family beside the other is the
point: the question an approving reader has to answer is not "is this one
circular" but "are these two independent of the claim in *different* ways", and
that is not visible when the declarations live in separate job files.

Family P is a pathologist's reading of tissue morphology. Family G is a somatic
genotype from FFPE whole-exome sequencing. Neither is computed from the
transcripts that produce the decomposition, which is what invariant 11 asks.

What the invariant does **not** ask, and what §4 of the pre-registration
handles instead: whether a label is *independent of the endpoint* as opposed to
independent of the assay. Family P is a reading of tissue architecture, and the
compositional term measures tissue architecture in another modality. That is
not circular provenance and the guard here will not catch it, which is exactly
why §4 gives the two terms different standing.
"""

from __future__ import annotations

from typing import Final

from src.common.label_provenance import Measurement, check_no_circular_claim

#: The claim both families are used to qualify. Fixed once so a family cannot
#: be declared against a quietly different endpoint.
CLAIM: Final[Measurement] = Measurement(
    modality="transcript",
    assay="Chen 2021 snRNA-seq adenoma decomposition (lineage rung)",
    genes=(),
)

#: Family P — histopathology. Two pathologists categorised each polyp from
#: H&E-stained sections (Chen et al., Cell 2021). Every lesion entering an arm
#: carries a specific histologic subtype, which a transcriptomic reclassification
#: does not produce; see `docs/chen_subtype_data_gate.md`.
FAMILY_P: Final[Measurement] = Measurement(
    modality="morphology",
    assay="Chen 2021 two-pathologist H&E polyp diagnosis (cBioPortal POLYP_TYPE)",
    genes=(),
)

#: Family G — somatic genotype. Truncating APC versus BRAF V600E from the FFPE
#: whole-exome MAF. The genes named here are the *label's* genes, and they are
#: deliberately not panel targets: invariant 2 keeps the panel out of any label,
#: and APC and BRAF are not in it.
FAMILY_G: Final[Measurement] = Measurement(
    modality="genotype",
    assay="Chen 2021 FFPE whole-exome somatic calls (cBioPortal MAF, hg19)",
    genes=("APC", "BRAF"),
)

FAMILIES: Final[dict[str, Measurement]] = {"P": FAMILY_P, "G": FAMILY_G}


def validate_families() -> dict[str, tuple[str, ...]]:
    """Run the guard for both families. Called before anything is read.

    Returns the overlapping genes per family, which must be empty: the claim
    declares no genes, so any overlap would mean a label gene had been written
    into the endpoint.
    """
    return {
        name: check_no_circular_claim(labels=measurement, claim=CLAIM)
        for name, measurement in FAMILIES.items()
    }
