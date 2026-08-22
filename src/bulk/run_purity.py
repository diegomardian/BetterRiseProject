"""W3.3 — run ESTIMATE and join ABSOLUTE. CLI, idempotent.

    python -m src.bulk.run_purity fetch     # ESTIMATE package + Aran table
    python -m src.bulk.run_purity run       # scores, join, agreement

Both inputs are external data: they land in ``data/raw/`` and get a
``data/manifest.csv`` row, like the GDC counts. Neither is committed — the
ESTIMATE package is GPL-2 and the Aran table is a publisher's supplementary
file.
"""

from __future__ import annotations

import argparse
import sys

import pandas as pd

from src.bulk.gdc import read_manifest
from src.bulk.gene_index import load_gene_index_map
from src.bulk.purity import (
    agreement,
    assemble_purity_table,
    estimate_scores,
    fetch_aran_table,
    fetch_estimate_package,
    load_aran_purity,
    load_common_genes,
    load_signatures,
    to_symbol_matrix,
)
from src.common.io import write_versioned_table
from src.common.paths import PROCESSED_DIR, RAW_DIR
from src.reference.ingest import append_manifest_row

SEED = 20260818
ESTIMATE_DIR = RAW_DIR / "estimate"
ARAN_DIR = RAW_DIR / "aran2015"
BULK = PROCESSED_DIR / "bulk"
INDEX_VERSION = "0.9.0"


def step_fetch(downloaded_by: str = "jeremy749") -> None:
    tarball = fetch_estimate_package(ESTIMATE_DIR)
    append_manifest_row(
        tarball,
        source_url="http://r-forge.r-project.org/src/contrib/estimate_1.0.13.tar.gz",
        accession="estimate_1.0.13",
        downloaded_by=downloaded_by,
        workstream="W3",
        notes="ESTIMATE R package (GPL-2); signatures + common genes + reference fixture",
    )
    aran = fetch_aran_table(ARAN_DIR)
    append_manifest_row(
        aran,
        source_url="https://doi.org/10.1038/ncomms9971",
        accession="aran2015_supplementary_data_1",
        downloaded_by=downloaded_by,
        workstream="W3",
        notes="Aran et al. 2015 pan-cancer purity: ESTIMATE/ABSOLUTE/LUMP/IHC/CPE",
    )
    print(f"wrote {tarball}\nwrote {aran}\nmanifest rows added")


def step_run() -> None:
    tarball = ESTIMATE_DIR / "estimate_1.0.13.tar.gz"
    aran_path = ARAN_DIR / "aran2015_ncomms9971_supplementary_data_1.xlsx"
    if not tarball.exists():
        raise SystemExit(f"{tarball} missing — run `fetch` first.")

    expression = pd.read_parquet(BULK / f"tcga_tpm_{INDEX_VERSION}.parquet")
    manifest = read_manifest(BULK / "sample_manifest.tsv").set_index("barcode")

    signatures = load_signatures(tarball)
    common = load_common_genes(tarball)
    symbol_matrix, counts = to_symbol_matrix(expression, load_gene_index_map(), common)

    print("gene mapping into ESTIMATE's space:")
    for k, v in counts.items():
        print(f"  {k:<28} {v:>7,}")
    matched = counts["unique_symbols"] / counts["common_genes_total"]
    print(f"  {'coverage of common genes':<28} {matched:>7.1%}")

    print(f"\nscoring {symbol_matrix.shape[1]} samples x {symbol_matrix.shape[0]} genes ...")
    scores = estimate_scores(symbol_matrix, signatures)

    absolute = load_aran_purity(aran_path) if aran_path.exists() else None
    if absolute is None:
        print("WARNING: Aran table absent — emitting ESTIMATE only.")
    table = assemble_purity_table(scores, absolute, manifest=manifest)

    stats = agreement(table)
    print("\nagreement:")
    for k, v in stats.items():
        print(f"  {k:<22} {v}")

    # External check: Aran published their own ESTIMATE purity for these samples.
    if absolute is not None:
        from src.bulk.purity import affymetrix_purity, sample_key

        ours = affymetrix_purity(scores["ESTIMATEScore"])
        keyed = absolute.set_index("sample_key")["aran_estimate"]
        joined = pd.DataFrame(
            {"ours": ours.to_numpy(), "aran": [
                keyed.get(sample_key(b), float("nan")) for b in scores.index
            ]},
            index=scores.index,
        ).dropna()
        if len(joined) >= 3:
            print(
                f"\nour ESTIMATE purity vs Aran's published ESTIMATE purity: "
                f"n={len(joined)}, r={joined['ours'].corr(joined['aran']):.4f}, "
                f"median |diff|={float((joined['ours'] - joined['aran']).abs().median()):.4f}"
            )

    BULK.mkdir(parents=True, exist_ok=True)
    table.to_parquet(BULK / f"tcga_purity_{INDEX_VERSION}.parquet", index=False)
    write_versioned_table(
        table,
        name="tcga_purity",
        seed=SEED,
        notes="W3.3 purity. One row per (sample, method); never coalesced.",
        extra_meta={"gene_mapping": counts, "agreement": stats, "index_version": INDEX_VERSION},
        allow_dirty=True,
    )
    print(f"\nwrote {BULK / f'tcga_purity_{INDEX_VERSION}.parquet'}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m src.bulk.run_purity")
    sub = parser.add_subparsers(dest="command", required=True)
    f = sub.add_parser("fetch", help="download ESTIMATE and the Aran table")
    f.add_argument("--downloaded-by", default="jeremy749")
    sub.add_parser("run", help="score, join and report agreement")

    args = parser.parse_args(argv)
    if args.command == "fetch":
        step_fetch(args.downloaded_by)
    else:
        step_run()
    return 0


if __name__ == "__main__":
    sys.exit(main())
