"""W3.8 — run the independent replication in GSE39582.

    python -m src.bulk.run_replication fetch
    python -m src.bulk.run_replication run

Per invariant 4 the cohorts are never pooled. Each is estimated separately with
the same test code, and the comparison is reported side by side.
"""

from __future__ import annotations

import argparse
import sys

import pandas as pd

from src.bulk.premise import assess
from src.bulk.replication import (
    GEO_PLATFORM_URL,
    GEO_SERIES_URL,
    collapse_probes_to_symbols,
    fetch,
    fold_change_vs_normal,
    parse_platform_table,
    parse_series_matrix,
    strata,
)
from src.common.io import write_versioned_table
from src.common.paths import PROCESSED_DIR, RAW_DIR
from src.reference.ingest import append_manifest_row

SEED = 20260818

#: Set from --allow-dirty. Default False. These jobs used to pass
#: allow_dirty=True unconditionally, so the bulk arm could not write a
#: clean provenance stamp even from a spotless tree -- which is why every
#: committed bulk table records git_dirty: true.
ALLOW_DIRTY = False
GEO_DIR = RAW_DIR / "gse39582"
SERIES = GEO_DIR / "GSE39582_series_matrix.txt.gz"
PLATFORM = GEO_DIR / "GPL570_table.txt"
BULK = PROCESSED_DIR / "bulk"

GENES = ("GUCA2A", "CDX2", "MS4A12")


def step_fetch(downloaded_by: str = "jeremy749") -> None:
    fetch(GEO_SERIES_URL, SERIES)
    append_manifest_row(
        SERIES, source_url=GEO_SERIES_URL, accession="GSE39582",
        downloaded_by=downloaded_by, workstream="W3",
        notes="Marisa et al. CIT cohort, 566 tumours + 19 normal mucosa, GPL570",
    )
    fetch(GEO_PLATFORM_URL, PLATFORM)
    append_manifest_row(
        PLATFORM, source_url=GEO_PLATFORM_URL, accession="GPL570",
        downloaded_by=downloaded_by, workstream="W3",
        notes="Affymetrix HG-U133 Plus 2.0 probe -> gene symbol annotation",
    )
    print(f"wrote {SERIES}\nwrote {PLATFORM}\nmanifest rows added")


def step_run() -> None:
    if not SERIES.exists():
        raise SystemExit(f"{SERIES} missing — run `fetch` first.")

    print("parsing series matrix ...")
    expression, metadata = parse_series_matrix(SERIES)
    print(f"  {expression.shape[0]} probes x {expression.shape[1]} samples")

    genes, counts = collapse_probes_to_symbols(expression, parse_platform_table(PLATFORM))
    for key, value in counts.items():
        print(f"  {key:<28} {value:>7,}")

    masks = strata(metadata)
    print("\nstrata:")
    for name, mask in masks.items():
        print(f"  {name:<20} n={int(mask.sum()):>4}")

    rows = []
    for gene in GENES:
        if gene not in genes.index:
            print(f"WARNING: {gene} has no unique probe on this platform — skipped")
            continue
        for name, mask in masks.items():
            values = genes.loc[gene][mask]
            for a in assess(values, gene=gene, stratum=name, seed=SEED):
                rows.append({**a.to_row(), "cohort": "GSE39582", "platform": "GPL570"})

    results = pd.DataFrame(rows)
    show = results[~results["stratum"].str.endswith("|nonzero")]
    print("\n=== bimodality in GSE39582 (dip null = unimodal, so large p = one group) ===")
    with pd.option_context("display.width", 210):
        print(
            show[["gene", "stratum", "n", "median", "dip_pvalue", "bic_delta", "verdict"]]
            .to_string(index=False)
        )

    changes = pd.DataFrame([fold_change_vs_normal(genes, masks, g) for g in GENES])
    print("\n=== tumour vs normal mucosa ===")
    print(changes.to_string(index=False))

    print("\n=== side by side with TCGA (independent estimates, never pooled) ===")
    comparison = _comparison(show, changes)
    print(comparison.to_string(index=False))

    BULK.mkdir(parents=True, exist_ok=True)
    genes.loc[list(GENES)].to_csv(BULK / "gse39582_panel_expression.tsv", sep="\t")
    metadata.to_csv(BULK / "gse39582_metadata.tsv", sep="\t")
    for frame, name, note in (
        (results, "gse39582_premise_bimodality", "W3.8 replication. Same test code as W3.2."),
        (changes, "gse39582_fold_change", "W3.8 tumour vs normal mucosa"),
        (comparison, "replication_tcga_vs_gse39582", "W3.8 side-by-side. NOT pooled."),
    ):
        write_versioned_table(
            frame, name=name, seed=SEED, notes=note,
            extra_meta={
                "cohort": "GSE39582",
                "geo_platform": "GPL570",
                "n_samples": int(len(metadata)),
            },
            allow_dirty=ALLOW_DIRTY,
        )
    print("\nwrote three results tables")


def _comparison(show: pd.DataFrame, changes: pd.DataFrame) -> pd.DataFrame:
    """TCGA numbers are read from the committed W3.2 results, not retyped."""
    import glob

    tcga_files = sorted(glob.glob("results/*/tcga_premise_bimodality.parquet"))
    rows = []
    tcga = pd.read_parquet(tcga_files[-1]) if tcga_files else pd.DataFrame()
    for gene in GENES:
        tcga_row = tcga.loc[
            (tcga["gene"] == gene) & (tcga["stratum"] == "COAD+READ_tumour")
        ] if len(tcga) else pd.DataFrame()
        geo_row = show.loc[(show["gene"] == gene) & (show["stratum"] == "tumour")]
        fold = changes.loc[changes["gene"] == gene, "fold_change"]
        rows.append(
            {
                "gene": gene,
                "tcga_n": int(tcga_row["n"].iloc[0]) if len(tcga_row) else None,
                "tcga_dip_p": float(tcga_row["dip_pvalue"].iloc[0]) if len(tcga_row) else None,
                "tcga_verdict": tcga_row["verdict"].iloc[0] if len(tcga_row) else None,
                "gse39582_n": int(geo_row["n"].iloc[0]) if len(geo_row) else None,
                "gse39582_dip_p": float(geo_row["dip_pvalue"].iloc[0]) if len(geo_row) else None,
                "gse39582_verdict": geo_row["verdict"].iloc[0] if len(geo_row) else None,
                "gse39582_fold_change_vs_normal": float(fold.iloc[0]) if len(fold) else None,
            }
        )
    return pd.DataFrame(rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m src.bulk.run_replication")
    sub = parser.add_subparsers(dest="command", required=True)
    f = sub.add_parser("fetch", help="download GSE39582 and the GPL570 annotation")
    f.add_argument("--downloaded-by", default="jeremy749")
    sub.add_parser("run", help="replicate the premise check")

    parser.add_argument(
        "--allow-dirty", action="store_true",
        help="write from a dirty tree; the recorded sha will not reproduce it",
    )
    args = parser.parse_args(argv)

    global ALLOW_DIRTY

    ALLOW_DIRTY = args.allow_dirty
    if args.command == "fetch":
        step_fetch(args.downloaded_by)
    else:
        step_run()
    return 0


if __name__ == "__main__":
    sys.exit(main())
