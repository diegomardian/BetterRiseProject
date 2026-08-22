#!/usr/bin/env python
"""Emit the shared gene index. W1, week 1. **W3 is blocked on this.**

    python src/reference/jobs/emit_gene_index.py --version 1.0.0

Writes config/gene_index/gene_index_{version}.txt and .map.tsv, then reports
what W3 needs to know: the genome build, the collision count, and which frozen
panel and axis genes fail to resolve.

Read-only apart from the two files it writes, and it refuses to overwrite an
existing version — results reference the index by version, so editing one in
place makes every earlier result unreproducible without saying so.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

from src.reference.gene_index import (
    build_gene_index,
    check_panel_coverage,
    write_gene_index,
)
from src.reference.ingest import read_gse178341_index


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", default="1.0.0")
    parser.add_argument("--h5", default=None)
    parser.add_argument("--dry-run", action="store_true", help="report, write nothing")
    args = parser.parse_args()

    h5 = Path(
        args.h5
        or Path(os.environ.get("BRP_DATA_DIR", "data"))
        / "raw" / "GSE178341" / "GSE178341_crc10x_full_c295v4_submit.h5"
    )
    if not h5.exists():
        raise SystemExit(f"missing {h5}")

    _, var = read_gse178341_index(h5)
    index, mapping, report = build_gene_index(var)

    print("=" * 70)
    print("SHARED GENE INDEX")
    print("=" * 70)
    print(f"  features in deposit      {report['n_features_in']:,}")
    print(f"  genes on the index       {report['n_genes_out']:,}")
    print(f"  collapsed by dedup       {report['n_collapsed']:,} "
          f"({report['n_duplicated_ids']:,} ids had >1 feature)")
    print(f"  genome build             {', '.join(report['genomes'])}")
    print()
    print("  FOR W3: this is an hg19 liftover; TCGA STAR counts are GRCh38.")
    print("  Ensembl IDs are mostly stable across builds so the join works, but")
    print("  genes retired or reassigned between builds will not join. Strip the")
    print("  version suffix from TCGA ids before joining (ENSG00000141510.16 ->")
    print("  ENSG00000141510).")

    coverage = check_panel_coverage(mapping)
    print("\n--- frozen panel and axis coverage ---")
    summary = coverage.groupby("source", observed=True)["found"].agg(["sum", "size"])
    for source, row in summary.iterrows():
        print(f"  {source:<24} {int(row['sum'])}/{int(row['size'])} resolve")
    unmapped = coverage[~coverage["found"]]
    if len(unmapped):
        print("\n  !! UNMAPPED — a panel gene absent from the index cannot be")
        print("     tested at all, and an axis marker absent from it silently")
        print("     weakens the labels. This is a week-1 finding.")
        print(unmapped[["source", "gene_symbol"]].to_string(index=False))
    else:
        print("  all frozen panel and axis genes resolve.")

    if args.dry_run:
        print("\n--dry-run: nothing written")
        return 0

    index_path, map_path = write_gene_index(index, mapping, version=args.version)
    print(f"\nwritten:\n  {index_path}\n  {map_path}")
    print("\nCommit both, and tell W3 the version. They reindex onto it; do not")
    print("let both arms build an index (open decision #2).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
