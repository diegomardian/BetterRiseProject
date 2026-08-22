"""Reconcile this deposit's features against W3's shared gene index. W1, week 2.

**Run this before anyone reindexes anything.** Open decision #2 said the thing
that must not happen is two indexes existing, and the documented fallback fired:
W1's week-1 index never reached the repo, so W3 built `gene_index_0.9.0` from the
GDC gene model. That file is correct — unversioned Ensembl key, symbols mapped,
exactly decision #3 — and this job asks the only question left, which is whether
GSE178341 can actually live on it.

The question is not rhetorical. Deconvolution needs each gene present in BOTH
matrices: a row of the S matrix for a gene the single-cell data never measured is
undefined, and a bulk gene with no single-cell counterpart contributes nothing.
So the operative index is the **intersection**, and its size is a number someone
has to look at rather than assume.

The assembly difference is the reason to check rather than trust. This deposit is
GRCh37_liftover_v28 (hg19); TCGA is GRCh38. Unversioned ENSG identifiers are
*supposed* to be stable across assemblies — that is why decision #3 chose them —
so the expected loss is small and comes from GENCODE release differences rather
than from the assembly itself. "Expected" is not "measured". This measures it,
and reports the panel genes separately because losing one of those is not a
rounding error.

    python src/reference/jobs/check_gene_index.py
    python src/reference/jobs/check_gene_index.py --version 0.9.0
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from src.common.panel import panel_genes  # noqa: E402
from src.common.paths import CONFIG_DIR  # noqa: E402
from src.reference.gene_index import GENOME_NOTE, build_gene_index  # noqa: E402
from src.reference.ingest import read_gse178341_index  # noqa: E402


def _load_shared(version: str) -> pd.DataFrame:
    directory = CONFIG_DIR / "gene_index"
    mapping = directory / f"gene_index_{version}.map.tsv"
    if not mapping.exists():
        raise SystemExit(
            f"missing {mapping}.\n"
            f"W3's index lands there when their branch merges. Until then:\n"
            f"  git show origin/w3/covariate-lock:config/gene_index/"
            f"gene_index_{version}.map.tsv > {mapping}\n"
            f"Do NOT commit it from here — it is W3's file and copying it is how "
            f"two indexes start existing (open decision #2)."
        )
    return pd.read_csv(mapping, sep="\t", dtype=str)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--version", default="0.9.0")
    args = parser.parse_args()

    shared = _load_shared(args.version)
    print(f"W3 shared index {args.version}: {len(shared):,} genes")

    data = Path(os.environ.get("BRP_DATA_DIR", "data")) / "raw" / "GSE178341"
    h5 = data / "GSE178341_crc10x_full_c295v4_submit.h5"
    if not h5.exists():
        raise SystemExit(f"missing {h5}")

    _obs, var = read_gse178341_index(h5)
    _, mine, report = build_gene_index(var)
    if report["genomes"]:
        print(f"deposit genome tag: {report['genomes']}")
    print(f"GSE178341 features:        {len(mine):,} genes")
    print(f"\n{GENOME_NOTE}\n")

    theirs = set(shared["ensembl_id"])
    ours = set(mine["ensembl_id"])
    both = theirs & ours

    print("overlap on unversioned Ensembl ID — the join that has to work:")
    print(f"  in both                  {len(both):>7,}")
    print(f"  bulk only (no sc counts) {len(theirs - ours):>7,}  "
          f"({len(theirs - ours) / len(theirs):.1%} of the shared index)")
    print(f"  sc only (no bulk row)    {len(ours - theirs):>7,}")
    print(f"\n  -> the operative index for deconvolution is the intersection: "
          f"{len(both):,} genes.")

    # Panel coverage is the part that cannot be waved through. A panel gene
    # missing from the intersection is untestable in the integrated analysis,
    # and the panel is frozen precisely so that is not negotiable after the fact.
    panel = set(panel_genes())
    symbol_col = "gene_symbol"
    shared_symbols = set(shared.loc[shared["ensembl_id"].isin(both), symbol_col])
    our_symbols = set(mine.loc[mine["ensembl_id"].isin(both), symbol_col])
    covered = panel & (shared_symbols | our_symbols)
    lost = sorted(panel - covered)

    print(f"\npanel coverage in the intersection: {len(covered)}/{len(panel)}")
    if lost:
        print(f"  !! MISSING: {lost}")
        print("     A frozen panel gene absent from the join is untestable in "
              "the integrated\n     analysis. Take this to the team before "
              "anyone reindexes — it is a\n     reason to revisit the index, not "
              "a reason to drop the gene.")
    else:
        print("  every panel gene survives the join.")

    print(
        "\nWhat to do with this number:\n"
        "  If the intersection is close to the single-cell feature count, adopt "
        "W3's\n  0.9.0 as-is and promote it to 1.0.0 — W1 does not need to emit a "
        "competing\n  index, and open decision #2 closes with W3 owning the file.\n"
        "  If it is materially smaller, 1.0.0 should BE the intersection, emitted "
        "once,\n  by whoever the team names. Either way: one index, agreed in the "
        "meeting."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
