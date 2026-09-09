"""Crosswalk Chen lineage polyp samples to HTAN WES provenance, without VCF reads.

Run only after exporting the normalized, VCF-only BigQuery CSV specified in
``docs/wes_subtype_plan.md``. The job reads that CSV plus the committed Step-1
manifest; it never authenticates to Synapse, downloads a VCF, inspects a VCF
header, or assigns BRAF/APC status.

    python -m src.reference.jobs.wes_subtype_provenance_crosswalk \
      --provenance-csv "$HOME/Downloads/htan_vanderbilt_vcf_provenance_r7.csv" \
      --no-write
"""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from pathlib import Path

import pandas as pd

from src.common.io import read_versioned_table, write_versioned_table
from src.common.label_provenance import Measurement, check_no_circular_claim, provenance_meta
from src.common.paths import RESULTS_DIR
from src.common.provenance import DEFAULT_SEED
from src.reference.wes_subtype import build_provenance_crosswalk, read_provenance_export

DEFAULT_MANIFEST = (
    RESULTS_DIR / "2026-09-09_7da9f79" / "wes_subtype_lineage_manifest.parquet"
)

LABEL_PROVENANCE = Measurement(
    modality="sample_annotation",
    assay="committed Chen avenue-A scRNA polyp biospecimen manifest",
    genes=(),
)
CLAIM_PROVENANCE = Measurement(
    modality="genotype",
    assay="HTAN Release-7 WES file-to-biospecimen provenance metadata",
    genes=(),
)


def validate_specification() -> tuple[str, ...]:
    return check_no_circular_claim(labels=LABEL_PROVENANCE, claim=CLAIM_PROVENANCE)


def _read_manifest(path: Path) -> pd.DataFrame:
    table, meta = read_versioned_table(path)
    if (
        meta.get("vcf_accessed") is not False
        or meta.get("molecular_endpoint_selected") is not False
    ):
        raise ValueError("lineage manifest is not the expected outcome-blind Step-1 artifact")
    return table


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--provenance-csv", required=True, type=Path)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--no-write", action="store_true")
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)

    validate_specification()
    table, summary = build_provenance_crosswalk(
        _read_manifest(args.manifest), read_provenance_export(args.provenance_csv)
    )
    print(summary.to_string(index=False))
    print("METADATA ONLY — no Synapse access, VCF header, or variant record was read.")
    if args.no_write:
        return 0
    meta = {
        "plan": "docs/wes_subtype_plan.md",
        "crosswalk_rule": "exact equality to assayed or originating biospecimen ID",
        "suffix_match_used": False,
        "vcf_content_read": False,
        "molecular_endpoint_selected": False,
        "subgroup_outcome_calculated": False,
        **provenance_meta(LABEL_PROVENANCE, CLAIM_PROVENANCE),
    }
    for result, name in (
        (table, "wes_subtype_provenance_crosswalk"),
        (summary, "wes_subtype_provenance_summary"),
    ):
        path = write_versioned_table(
            result,
            name,
            seed=args.seed,
            results_dir=args.results_dir,
            allow_dirty=args.allow_dirty,
            extra_meta=meta,
        )
        print(f"wrote {path}")
    return 0 if summary.loc[0, "n_specimen_exact_unique_biospecimens"] else 5


if __name__ == "__main__":
    sys.exit(main())
