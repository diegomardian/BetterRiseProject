"""Invariant 11: no transcript-derived label may claim its own programme.

    from src.common.label_provenance import Measurement, check_no_circular_claim

    check_no_circular_claim(
        labels=Measurement(modality="morphology", assay="H&E domain call", genes=()),
        claim=Measurement(modality="transcript", assay="Xenium 5K",
                          genes=("GUCA2A", "MS4A12")),
    )

THE FAILURE THIS PREVENTS
-------------------------
Select cells by the transcripts that define maturity, then report that those
cells retain — or have lost — those same transcripts. The answer was fixed by
the selection rule before any tissue was measured, and it looks exactly like a
biological finding: a tight interval, a clean direction, a plausible story.
Invariant 2 already forbids the single worst case (targets in the labels). This
is the general form, and it binds wherever a population is defined.

WHERE IT LIVES, AND WHERE IT DOES NOT
-------------------------------------
At **analysis-specification validation** — the point where a job says what it is
about to measure — and not in the generic result writer. A writer sees a table
of numbers, by which time the population is already chosen and the columns no
longer say what defined it. A guard there could only check that two fields were
filled in.

DECLARATION IS MANDATORY AND OMISSION IS NOT A PASS
---------------------------------------------------
There is no default. An analysis whose labels use no genes must say so with an
explicit empty tuple; leaving the declaration out is refused with the same force
as declaring a circular one. This is invariant 1's rule in another place —
``None`` is not ``0.0``, and "unstated" is not "none".

THE EXCEPTION IS STRUCTURAL, NOT A NAME
---------------------------------------
An overlap is permitted only when the endpoint was **independently measured** —
CDX2 protein by MxIF against a CDX2-transcript label, say. Supplying an assay
name does not establish that. :func:`check_no_circular_claim` requires four
things that a name cannot fake: the declaration must name the assay that
actually produced the endpoint, that assay must differ from the one that
produced the labels, the modality must differ too, and the labels must not
descend from the endpoint's own assay. The stated reason is recorded beside the
result and is never, on its own, sufficient.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any, Final

#: How a measurement was made. A modality is a physical channel, not a vendor.
#: Two measurements sharing a modality share its failure modes, which is why an
#: overlap across the same modality is refused however it is justified.
MODALITIES: Final[tuple[str, ...]] = (
    "transcript",
    "protein",
    "morphology",
    "geometry",
    "spatial_position",
    "genotype",
    "sample_annotation",
)

#: Declaration keys an analysis specification must carry. Both, always.
REQUIRED_SPEC_KEYS: Final[tuple[str, ...]] = ("label_provenance", "claim_provenance")


class CircularClaimError(RuntimeError):
    """The analysis would claim the state of its own defining programme."""


class ProvenanceDeclarationError(RuntimeError):
    """The analysis did not say what defined its population, or what it claims."""


@dataclass(frozen=True)
class Measurement:
    """One measurement channel: how it was made, by what, and over which genes.

    ``genes`` is the empty tuple when the channel uses none — a morphological
    domain call, a sample-annotation arm label. It must be written out. See the
    module docstring on why omission is not a pass.

    ``derived_from`` names the assays this measurement was computed from, when
    it is downstream of another. A label derived from the endpoint's own assay
    is circular no matter what modality it is reported in.
    """

    modality: str
    assay: str
    genes: tuple[str, ...] = ()
    derived_from: tuple[str, ...] = field(default=())

    def __post_init__(self) -> None:
        if self.modality not in MODALITIES:
            raise ProvenanceDeclarationError(
                f"modality {self.modality!r} is not one of {list(MODALITIES)}. "
                f"A modality is the physical channel the measurement came "
                f"through; if none of these fits, the vocabulary is what needs "
                f"the change, in its own PR."
            )
        if not isinstance(self.assay, str) or not self.assay.strip():
            raise ProvenanceDeclarationError(
                f"a {self.modality!r} measurement was declared with no assay. "
                f"Name the instrument or annotation that produced it."
            )
        if isinstance(self.genes, str):
            raise ProvenanceDeclarationError(
                f"genes must be a tuple, got the string {self.genes!r} — which "
                f"would be read one character at a time."
            )
        for gene in self.genes:
            if not isinstance(gene, str) or not gene.strip():
                raise ProvenanceDeclarationError(f"empty gene symbol in {self.genes!r}")

    @property
    def symbols(self) -> tuple[str, ...]:
        """Gene symbols, upper-cased, so ``Cdx2`` and ``CDX2`` cannot differ."""
        return tuple(gene.strip().upper() for gene in self.genes)

    def as_row(self) -> dict[str, Any]:
        """Flat record for the result sidecar."""
        return {
            "modality": self.modality,
            "assay": self.assay,
            "genes": ";".join(self.symbols),
            "derived_from": ";".join(self.derived_from),
        }


@dataclass(frozen=True)
class IndependentEndpoint:
    """The claim that an overlapping endpoint was measured independently.

    ``assay`` must be the assay that produced the endpoint — it is checked
    against the claim's own declaration, so this cannot be a different name
    supplied to satisfy the guard. ``reason`` is recorded and never evaluated:
    it is there for the reader, and it does not carry the exception.
    """

    assay: str
    reason: str

    def __post_init__(self) -> None:
        for name in ("assay", "reason"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise ProvenanceDeclarationError(
                    f"an independent-endpoint declaration needs a {name}"
                )


def _require_measurement(value: Any, *, role: str) -> Measurement:
    if value is None:
        raise ProvenanceDeclarationError(
            f"no {role} declared. Invariant 11 requires every analysis to say "
            f"what defined its population and what it claims. If the {role} "
            f"uses no genes, declare it with `genes=()` — leaving it out is not "
            f"the same statement and is not accepted as one."
        )
    if not isinstance(value, Measurement):
        raise ProvenanceDeclarationError(
            f"{role} must be a Measurement, got {type(value).__name__}"
        )
    return value


def check_no_circular_claim(
    *,
    labels: Measurement | None,
    claim: Measurement | None,
    independent_endpoint: IndependentEndpoint | None = None,
) -> tuple[str, ...]:
    """Refuse an analysis that claims the programme its labels are made of.

    Returns the overlapping symbols, which is ``()`` in the ordinary case and
    non-empty only when a structural exception was accepted. Callers record the
    return value rather than re-deriving it.
    """
    labels = _require_measurement(labels, role="label provenance")
    claim = _require_measurement(claim, role="claim provenance")

    label_symbols = set(labels.symbols)
    overlap = tuple(symbol for symbol in claim.symbols if symbol in label_symbols)
    if not overlap:
        return ()

    named = ", ".join(sorted(set(overlap)))
    circular = (
        f"CIRCULAR CLAIM (invariant 11): {named} define the population "
        f"({labels.modality} by {labels.assay}) and are also the endpoint "
        f"({claim.modality} by {claim.assay})."
    )

    if independent_endpoint is None:
        raise CircularClaimError(
            f"{circular} Selecting on these genes fixes the answer about them "
            f"before the tissue is measured. Either change the population "
            f"definition, or declare an independently measured endpoint — an "
            f"assay name alone will not do it."
        )

    if independent_endpoint.assay != claim.assay:
        raise CircularClaimError(
            f"{circular} The independent-endpoint declaration names assay "
            f"{independent_endpoint.assay!r}, but the endpoint was declared as "
            f"measured by {claim.assay!r}. The exception must be about the "
            f"measurement that was actually made."
        )
    if claim.assay == labels.assay:
        raise CircularClaimError(
            f"{circular} Both come off {claim.assay!r}. One assay measuring a "
            f"gene twice is not an independent measurement of it."
        )
    if claim.modality == labels.modality:
        raise CircularClaimError(
            f"{circular} Both are {claim.modality!r} measurements. A second "
            f"instrument in the same channel shares the channel's failure "
            f"modes, so it cannot break the circle."
        )
    if claim.assay in labels.derived_from:
        raise CircularClaimError(
            f"{circular} The labels were derived from {claim.assay!r}, the "
            f"endpoint's own assay. A label downstream of the endpoint is "
            f"circular whichever modality it is reported in."
        )
    return overlap


def measurement_from_mapping(value: Any, *, role: str) -> Measurement:
    """Build a :class:`Measurement` from a config mapping, refusing omissions."""
    if not isinstance(value, Mapping):
        raise ProvenanceDeclarationError(
            f"{role} must be a mapping with `modality`, `assay` and `genes`, "
            f"got {type(value).__name__}"
        )
    missing = [key for key in ("modality", "assay", "genes") if key not in value]
    if missing:
        raise ProvenanceDeclarationError(f"{role} is missing {missing}")
    return Measurement(
        modality=str(value["modality"]),
        assay=str(value["assay"]),
        genes=tuple(value["genes"] or ()),
        derived_from=tuple(value.get("derived_from") or ()),
    )


def check_spec_declares_provenance(spec: Mapping[str, Any]) -> tuple[str, ...]:
    """Validate one analysis specification. This is the wiring point.

    Call it where a job validates its specification — before the population is
    built, not after the numbers exist.
    """
    if not isinstance(spec, Mapping):
        raise ProvenanceDeclarationError(
            f"an analysis specification must be a mapping, got {type(spec).__name__}"
        )
    missing = [key for key in REQUIRED_SPEC_KEYS if key not in spec]
    if missing:
        raise ProvenanceDeclarationError(
            f"analysis specification is missing {missing}. Invariant 11 has no "
            f"default: an analysis that does not say what defined its "
            f"population does not run."
        )
    labels = measurement_from_mapping(spec["label_provenance"], role="label_provenance")
    claim = measurement_from_mapping(spec["claim_provenance"], role="claim_provenance")

    declared = spec.get("independent_endpoint")
    independent = None
    if declared is not None:
        if not isinstance(declared, Mapping):
            raise ProvenanceDeclarationError(
                "independent_endpoint must be a mapping with `assay` and `reason`"
            )
        missing_keys = [key for key in ("assay", "reason") if key not in declared]
        if missing_keys:
            raise ProvenanceDeclarationError(
                f"independent_endpoint is missing {missing_keys}"
            )
        independent = IndependentEndpoint(
            assay=str(declared["assay"]), reason=str(declared["reason"])
        )
    return check_no_circular_claim(
        labels=labels, claim=claim, independent_endpoint=independent
    )


def provenance_meta(
    labels: Measurement,
    claim: Measurement,
    *,
    overlap: tuple[str, ...] = (),
    independent_endpoint: IndependentEndpoint | None = None,
) -> dict[str, Any]:
    """Sidecar record: what defined the population, and what was claimed.

    Written beside the result so a reader can check the population definition
    without reading the job, which is where every one of these arguments
    actually gets settled.
    """
    return {
        "invariant_11": "no transcript-derived label may claim its own programme",
        "label_provenance": labels.as_row(),
        "claim_provenance": claim.as_row(),
        "overlapping_genes": ";".join(overlap),
        "independent_endpoint": (
            None if independent_endpoint is None
            else {"assay": independent_endpoint.assay,
                  "reason": independent_endpoint.reason}
        ),
    }
