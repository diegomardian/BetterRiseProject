"""Gate 4: the Replogle sensitivity study's one label family, declared.

`docs/prereg_replogle_sensitivity.md` §2. The declaration is the design
decision, not a formality: a perturbed cell is defined by **which CRISPRi guide
it carries**, a DNA-level assignment, and never by whether the targeted
transcript reads low.

Defining the perturbed arm by low expression of the targeted gene and then
asking whether that gene reads low is invariant 11's circularity in its purest
form. It would also manufacture near-perfect power, which is the specific way
this study could produce a confident and meaningless curve.
"""

from __future__ import annotations

from typing import Final

from src.common.label_provenance import Measurement, check_no_circular_claim

#: The claim: can this pipeline detect that a gene is off, in cells retained by
#: construction. Measured on transcripts.
CLAIM: Final[Measurement] = Measurement(
    modality="transcript",
    assay="Replogle 2022 Perturb-seq single-cell UMI counts",
    genes=(),
)

#: The label: which guide the cell carries. Read from guide capture, not from
#: the expression matrix. Genes are empty because the guide identity is not a
#: transcript measurement of any gene -- it is a DNA barcode.
GUIDE_ASSIGNMENT: Final[Measurement] = Measurement(
    modality="genotype",
    assay="Replogle 2022 CRISPRi guide capture (per-cell guide identity)",
    genes=(),
)


def validate() -> tuple[str, ...]:
    """Run the guard. Called before any count is read."""
    return check_no_circular_claim(labels=GUIDE_ASSIGNMENT, claim=CLAIM)
