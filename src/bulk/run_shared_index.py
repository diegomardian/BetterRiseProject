#!/usr/bin/env python
"""Emit `config/gene_index/gene_index_1.0.0.*` — the shared index, decision #2.

    python -m src.bulk.run_shared_index --dry-run
    python -m src.bulk.run_shared_index

Reads the single-cell feature list through W1's own
``read_gse178341_index`` / ``build_gene_index`` rather than reimplementing it,
so the 43,113 here is the same 43,113 their `check_gene_index.py` reports. The
deposit is recorded in `data/manifest.csv`; features are read without loading
the matrix.

`--dry-run` prints the reconciliation and writes nothing, which is the mode to
run first — `write_gene_index` refuses to overwrite an existing version, so a
mistaken emit costs a version number.
"""

from __future__ import annotations

import argparse
import sys

from src.bulk.gene_index import GeneIndexError
from src.bulk.shared_index import (
    BULK_VERSION,
    SHARED_VERSION,
    build_shared_index,
    emit,
)
from src.common.paths import RAW_DIR

DEPOSIT = RAW_DIR / "GSE178341" / "GSE178341_crc10x_full_c295v4_submit.h5"


def reference_identifiers() -> list[str]:
    """The deposit's unversioned Ensembl ids, via W1's reader."""
    from src.reference.gene_index import build_gene_index
    from src.reference.ingest import read_gse178341_index

    if not DEPOSIT.exists():
        raise GeneIndexError(
            f"missing {DEPOSIT}.\n"
            f"It is W1's deposit and it is in data/manifest.csv with a sha256. "
            f"Fetch it there, or point BRP_DATA_DIR at a copy."
        )
    _obs, var = read_gse178341_index(DEPOSIT)
    _index, mapping, report = build_gene_index(var)
    if report.get("genomes"):
        print(f"deposit genome tag: {report['genomes']}")
    return list(mapping["ensembl_id"])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", default=SHARED_VERSION)
    parser.add_argument("--bulk-version", default=BULK_VERSION)
    parser.add_argument("--dry-run", action="store_true", help="report, write nothing")
    args = parser.parse_args()

    ids = reference_identifiers()
    print(f"reference features: {len(ids):,}")

    if args.dry_run:
        shared, report = build_shared_index(ids, bulk_version=args.bulk_version)
        print(report.to_string(index=False))
        print(f"\nwould write {len(shared):,} genes as {args.version} — nothing written")
        return 0

    index_path, map_path, report = emit(
        ids, version=args.version, bulk_version=args.bulk_version
    )
    print(report.to_string(index=False))
    print(f"\nwrote {index_path}\nwrote {map_path}")
    print(
        "\nNEXT — this index alone does not unblock the S matrices:\n"
        "  1. run_full_reference.py:208 passes gene_SYMBOLS as gene_names while "
        "this index\n     is Ensembl-keyed, so build_signature_sparse raises "
        "'no gene in the matrix\n     appears on the shared gene index'. "
        "run_pilot.py used var['ensembl_id'] and\n     was right; the "
        "full-scale job regressed.\n"
        "  2. assert_no_target_leakage compares target SYMBOLS against a gene-ID "
        "list, so\n     every invariant-2 guard is inert in Ensembl space. "
        "GUCA2A is already in all\n     four committed pilot S matrices.\n"
        "Both are W1's to fix and are raised as an issue."
    )
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except GeneIndexError as exc:
        print(f"!! {exc}", file=sys.stderr)
        sys.exit(1)
