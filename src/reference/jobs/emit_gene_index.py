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
    intersect_gene_index,
    read_shared_index,
    write_gene_index,
)
from src.reference.ingest import read_gse178341_index


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", default="1.0.0")
    # THE OPERATIVE INDEX IS THE INTERSECTION (open decision #2, corrected).
    # Without this, the emitted index is this deposit's features alone — which
    # is the right input for 0.x provisional work and the WRONG file to call
    # 1.0.0, because W3 would then hold 21,380 genes with no single-cell counts.
    parser.add_argument(
        "--intersect-with", metavar="INDEX_TXT", default=None,
        help="another arm's emitted index; 1.0.0 requires it",
    )
    parser.add_argument("--h5", default=None)
    parser.add_argument("--dry-run", action="store_true", help="report, write nothing")
    args = parser.parse_args()

    # 1.0.0 IS the intersection. Emitting the single-cell feature set under
    # that name would put a version string on every downstream result asserting
    # an agreement the file does not honour, and nothing later would catch it.
    if args.intersect_with is None and not args.version.startswith("0."):
        raise SystemExit(
            f"refusing to emit {args.version} without --intersect-with.\n\n"
            f"Open decision #2, corrected: the operative index for deconvolution "
            f"is the INTERSECTION of this deposit's features and the shared "
            f"index, because a gene present in only one matrix cannot be "
            f"deconvolved from the other. Emitting this deposit's 43,113 "
            f"features as 1.0.0 would hand W3 an index carrying 21,380 genes "
            f"with no single-cell counts.\n\n"
            f"  python {sys.argv[0]} --version {args.version} \\\n"
            f"      --intersect-with config/gene_index/gene_index_0.9.0.txt\n\n"
            f"Use a 0.x version for a provisional deposit-only index."
        )


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

    if args.intersect_with is not None:
        shared = read_shared_index(args.intersect_with)
        index, mapping, ireport = intersect_gene_index(mapping, shared)
        print()
        print("=" * 70)
        print("INTERSECTION — the operative index (open decision #2, corrected)")
        print("=" * 70)
        print(f"  this deposit             {ireport['n_ours']:,}")
        print(f"  shared index             {ireport['n_shared']:,}  "
              f"({args.intersect_with})")
        print(f"  IN BOTH                  {ireport['n_intersection']:,}")
        print(f"  deposit only             {ireport['n_ours_only']:,}  "
              f"(no bulk row)")
        print(f"  shared only              {ireport['n_shared_only']:,}  "
              f"(no single-cell counts)")
        print(f"  panel coverage           {ireport['panel_found']}/"
              f"{ireport['panel_total']}")
        print()
        print("  Sorted ascending, so this file is byte-identical whichever arm")
        print("  runs it. Verify with sha256 rather than by trusting either.")

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
