"""Write the outcome-blind scRNA side of the proposed Chen WES subtype gate.

This job reads only (1) patient IDs from the committed avenue-A lineage table
and (2) cached ICBI sample metadata.  It does not access HTAN, Synapse, or a
VCF; it does not select a molecular endpoint or calculate subgroup outcomes.

    python -m src.reference.jobs.wes_subtype_lineage_manifest --no-write

The output is the exact list of deposited polyp scRNA biospecimen IDs that a
future HTAN provenance export must prove related to a VCF.  Numeric suffixes
are retained as opaque identifiers and never used as a matching rule.
"""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from pathlib import Path

from src.common.io import write_versioned_table
from src.common.label_provenance import Measurement, check_no_circular_claim, provenance_meta
from src.common.paths import INTERIM_DIR, RESULTS_DIR
from src.common.provenance import DEFAULT_SEED
from src.reference.wes_subtype import (
    EXPECTED_LINEAGE_PATIENTS,
    build_lineage_manifest,
    read_decomposition,
    read_obs,
)

DEFAULT_DECOMPOSITION = RESULTS_DIR / "2026-09-06_6adcedc" / "icbi_adenoma.parquet"

LABEL_PROVENANCE = Measurement(
    modality="sample_annotation",
    assay="committed Chen avenue-A lineage patient and scRNA sample metadata",
    genes=(),
)
CLAIM_PROVENANCE = Measurement(
    modality="genotype",
    assay="future HTAN provenance-linked WES biospecimen identity",
    genes=(),
)


def validate_specification() -> tuple[str, ...]:
    """Declare a metadata-only, gene-free crosswalk before either input is read."""
    return check_no_circular_claim(labels=LABEL_PROVENANCE, claim=CLAIM_PROVENANCE)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--decomposition", type=Path, default=DEFAULT_DECOMPOSITION)
    parser.add_argument("--cache", type=Path, default=INTERIM_DIR / "icbi_obs.parquet")
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--no-write", action="store_true")
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)

    validate_specification()
    manifest = build_lineage_manifest(
        read_decomposition(args.decomposition), read_obs(args.cache)
    )
    print(manifest.to_string(index=False))
    print(
        "METADATA ONLY — no HTAN provenance export or VCF was read. "
        f"{manifest['patient_id'].nunique()}/{EXPECTED_LINEAGE_PATIENTS} lineage patients; "
        f"{len(manifest)} polyp scRNA biospecimens await exact provenance matching."
    )
    if args.no_write:
        return 0
    path = write_versioned_table(
        manifest,
        "wes_subtype_lineage_manifest",
        seed=args.seed,
        results_dir=args.results_dir,
        allow_dirty=args.allow_dirty,
        notes="outcome-blind Chen scRNA biospecimen manifest; no HTAN or VCF access",
        extra_meta={
            "plan": "docs/wes_subtype_plan.md",
            "n_lineage_patients": int(manifest["patient_id"].nunique()),
            "n_polyp_scRNA_biospecimens": int(len(manifest)),
            "n_patients_with_multiple_polyp_scRNA_biospecimens": int(
                manifest.loc[manifest["multiple_polyp_scRNA_biospecimens"], "patient_id"].nunique()
            ),
            "vcf_accessed": False,
            "molecular_endpoint_selected": False,
            "subgroup_outcome_calculated": False,
            "what_this_licenses": "only an exact HTAN provenance-export crosswalk",
            "what_this_does_not_license": (
                "participant-level linkage, VCF access, subtype assignment, or "
                "a stratified transcriptomic analysis"
            ),
            **provenance_meta(LABEL_PROVENANCE, CLAIM_PROVENANCE),
        },
    )
    print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
