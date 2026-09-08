"""Run the H&E molecular-prediction metadata-only substrate gate.

This reads HTA11 Files, Biospecimen, and Case exports. It does not download an
image, open a VCF, choose a molecular endpoint, or fit a model.

    python -m src.reference.jobs.he_molecular_gate \
      --files /path/to/files.tsv \
      --biospecimens /path/to/biospecimens.tsv \
      --cases /path/to/cases.tsv

When it writes a durable `NOT LICENSED` result, it exits 5 after the tables
are written. `--no-write` is only a local inspection and exits 0.
"""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from pathlib import Path

import pandas as pd

from src.common.io import write_versioned_table
from src.common.provenance import DEFAULT_SEED
from src.reference.he_molecular_gate import build_inventory, exit_code


def _read(path: Path, *, name: str) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"{name} metadata not found: {path}")
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--files", required=True, type=Path)
    parser.add_argument("--biospecimens", required=True, type=Path)
    parser.add_argument("--cases", required=True, type=Path)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--no-write", action="store_true")
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)

    attrition, summary = build_inventory(
        _read(args.files, name="Files"),
        _read(args.biospecimens, name="Biospecimen"),
        _read(args.cases, name="Case"),
    )
    print(summary.to_string(index=False))
    print("METADATA ONLY — no image or VCF read; no endpoint or model selected.")
    if args.no_write:
        return 0

    meta = {
        "prereg": "docs/he_molecular_data_gate.md",
        "assay_pair": {"image": "H&E Level 2", "molecular": "Bulk DNA Level 3 VCF"},
        "crosswalk_rule": "exact HTAN biospecimen identifier equality",
        "endpoint_selected": False,
        "image_or_vcf_read": False,
        "model_fit": False,
    }
    for table, name in (
        (attrition, "he_molecular_gate_attrition"),
        (summary, "he_molecular_gate_summary"),
    ):
        path = write_versioned_table(
            table,
            name,
            seed=args.seed,
            results_dir=args.results_dir,
            allow_dirty=args.allow_dirty,
            extra_meta=meta,
        )
        print(f"wrote {path}")
    return exit_code(summary)


if __name__ == "__main__":
    sys.exit(main())
