"""W3.1 — GDC ingest and normalisation. Runnable, idempotent, standalone.

    python -m src.bulk.ingest query        # what exists on the GDC
    python -m src.bulk.ingest download     # fetch STAR counts (slow, resumable)
    python -m src.bulk.ingest gene-index   # provisional shared index, decision #2 fallback
    python -m src.bulk.ingest build        # assemble matrices + sample manifest

Each step writes to ``data/`` and is safe to re-run: ``download`` skips files it
already has, ``gene-index`` refuses to overwrite a committed index, and
``build`` is a pure function of what is on disk.

WHAT LANDS WHERE
----------------
``data/processed/bulk/`` holds the matrices — too large for git, and
``data/README.md`` names ``processed/`` as the cross-workstream handoff.
``results/`` gets the small tables that are *findings*: the gene-loss report,
the panel resolution report, the portal reconciliation. Those are committed,
through ``write_versioned_table`` so they carry a sha and a seed
(CLAUDE.md invariant 10).

Nothing here uses a panel gene to define a group, so invariant 2 is not at
risk — but note that GUCA2A and CDX2 must *survive* into the matrix, because
they are the outcome variables for the week-2 premise check. The index is built
to include them; see ``src/bulk/gene_index.py``.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from src.bulk.gdc import (
    DEFAULT_SAMPLE_TYPES,
    PROJECTS,
    build_sample_manifest,
    deduplicate_aliquots,
    download_file,
    query_star_files,
    read_manifest,
    reconcile_counts,
)
from src.bulk.gene_index import (
    PAR_Y_SUFFIX,
    PROVISIONAL_VERSION,
    build_gene_index,
    gene_model_version,
    load_gene_index,
    panel_resolution_report,
    read_star_counts,
    strip_version,
    write_gene_index,
)
from src.bulk.normalise import (
    assert_counts,
    assert_log_scale,
    assert_tpm,
    counts_to_log2_cpm,
    renormalise_tpm,
)
from src.common.io import write_versioned_table
from src.common.paths import PROCESSED_DIR, RAW_DIR

#: Fixed seed. Nothing here is stochastic, but the writer requires one and
#: "0 is a fine answer, silence is not" (src/schema.py).
SEED = 20260817

#: Set from --allow-dirty. Default False. These jobs used to pass
#: allow_dirty=True unconditionally, so the bulk arm could not write a
#: clean provenance stamp even from a spotless tree -- which is why every
#: committed bulk table records git_dirty: true.
ALLOW_DIRTY = False

TCGA_RAW = RAW_DIR / "tcga"
BULK_PROCESSED = PROCESSED_DIR / "bulk"

FILES_TABLE = TCGA_RAW / "gdc_star_files.tsv"
SAMPLE_MANIFEST = BULK_PROCESSED / "sample_manifest.tsv"


# ---------------------------------------------------------------------------
# Steps
# ---------------------------------------------------------------------------


def step_query(projects: tuple[str, ...], sample_types: tuple[str, ...]) -> pd.DataFrame:
    """Ask the GDC what exists and cache the answer."""
    files = query_star_files(projects, sample_types)
    manifest = build_sample_manifest(files)
    kept, dropped = deduplicate_aliquots(manifest)

    TCGA_RAW.mkdir(parents=True, exist_ok=True)
    kept.to_csv(FILES_TABLE, sep="\t", index=False)

    print(f"GDC returned {len(manifest)} STAR-counts files for {list(projects)}")
    print(f"  kept    {len(kept)} after one-aliquot-per-(patient, sample type)")
    print(f"  dropped {len(dropped)} duplicate aliquot(s)")
    if len(dropped):
        print("  dropped barcodes:", ", ".join(dropped["barcode"].head(10)))
    print("\nReconciliation — compare these against the GDC portal facets:")
    print(reconcile_counts(kept).to_string(index=False))
    print(f"\nwrote {FILES_TABLE}")
    return kept


def step_download(limit: int | None = None) -> list[Path]:
    """Fetch the STAR-counts files listed by ``query``. Resumable."""
    if not FILES_TABLE.exists():
        raise SystemExit(f"{FILES_TABLE} missing — run `query` first.")
    files = read_manifest(FILES_TABLE)
    if limit:
        files = files.head(limit)

    paths: list[Path] = []
    for i, row in enumerate(files.itertuples(index=False), start=1):
        target = TCGA_RAW / "star" / f"{row.barcode}.rna_seq.augmented_star_gene_counts.tsv"
        paths.append(download_file(row.file_id, target))
        if i % 50 == 0 or i == len(files):
            print(f"  {i}/{len(files)}")
    print(f"{len(paths)} file(s) in {TCGA_RAW / 'star'}")
    return paths


def step_gene_index(version: str = PROVISIONAL_VERSION) -> pd.DataFrame:
    """Build the provisional shared index from one STAR file's gene model.

    One file is enough: every STAR-counts file in a GDC release carries the
    identical gene model, which is the property that makes this a shared index
    rather than a per-sample one. The build asserts it against a second file.
    """
    star_files = sorted((TCGA_RAW / "star").glob("*.tsv"))
    if not star_files:
        raise SystemExit(f"no STAR files in {TCGA_RAW / 'star'} — run `download` first.")

    first = read_star_counts(star_files[0])
    index, report = build_gene_index(first)

    if len(star_files) > 1:
        second = read_star_counts(star_files[1])
        if not first["gene_id"].equals(second["gene_id"]):
            raise SystemExit(
                f"{star_files[0].name} and {star_files[1].name} have different gene "
                f"models. The index cannot be shared across files that disagree — "
                f"check whether the download spans two GDC data releases."
            )

    model = gene_model_version(star_files[0])
    print(f"gene model: {model or 'not recorded in the file header'}")
    print("filtering:")
    print(report.summary())

    panel = panel_resolution_report(index)
    unresolved = panel.loc[panel["status"] != "resolved"]
    print(f"\npanel resolution: {len(panel) - len(unresolved)}/{len(panel)} resolved")
    if len(unresolved):
        print("  NOT CLEANLY RESOLVED — these need a decision before the matrix is used:")
        print(unresolved.to_string(index=False))

    idx_path, map_path = write_gene_index(index, version=version)
    print(f"\nwrote {idx_path}\nwrote {map_path}")

    write_versioned_table(
        report.to_frame(), name="tcga_gene_index_filtering", seed=SEED,
        notes=f"provisional index {version} from {model}", allow_dirty=ALLOW_DIRTY,
    )
    write_versioned_table(
        panel, name="tcga_panel_resolution", seed=SEED,
        notes="panel symbol -> Ensembl on the provisional index", allow_dirty=ALLOW_DIRTY,
    )
    return index


def step_build(version: str = PROVISIONAL_VERSION) -> None:
    """Assemble both matrices on the shared index and write the sample manifest."""
    if not FILES_TABLE.exists():
        raise SystemExit(f"{FILES_TABLE} missing — run `query` first.")
    manifest = read_manifest(FILES_TABLE)
    index_ids = load_gene_index(version)

    counts, tpm = _read_all(manifest)
    counts, tpm, loss = _reindex(counts, tpm, index_ids)

    log2_cpm = counts_to_log2_cpm(counts)
    tpm = renormalise_tpm(tpm)

    # Guards at the exit, not only the entry. These are the two files everyone
    # else consumes and the scale is the thing they can least afford to guess.
    assert_counts(counts, context="the emitted counts matrix")
    assert_tpm(tpm, context="the emitted TPM matrix")
    assert_log_scale(log2_cpm, context="the emitted log2-CPM matrix")

    BULK_PROCESSED.mkdir(parents=True, exist_ok=True)
    counts.to_parquet(BULK_PROCESSED / f"tcga_counts_{version}.parquet")
    tpm.to_parquet(BULK_PROCESSED / f"tcga_tpm_{version}.parquet")
    log2_cpm.to_parquet(BULK_PROCESSED / f"tcga_log2cpm_{version}.parquet")
    manifest.to_csv(SAMPLE_MANIFEST, sep="\t", index=False)

    print(f"{counts.shape[0]} samples x {counts.shape[1]} genes")
    print("reindex loss:")
    print(loss.to_string(index=False))
    print(f"\nwrote {BULK_PROCESSED}/tcga_{{counts,tpm,log2cpm}}_{version}.parquet")
    print(f"wrote {SAMPLE_MANIFEST}")

    write_versioned_table(
        reconcile_counts(manifest), name="tcga_sample_reconciliation", seed=SEED,
        notes="compare against GDC portal facet counts", allow_dirty=ALLOW_DIRTY,
    )


# ---------------------------------------------------------------------------
# Matrix assembly
# ---------------------------------------------------------------------------


def _read_all(manifest: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read every STAR file into (counts, tpm), samples x genes, keyed by barcode."""
    count_rows: dict[str, pd.Series] = {}
    tpm_rows: dict[str, pd.Series] = {}
    for row in manifest.itertuples(index=False):
        path = TCGA_RAW / "star" / f"{row.barcode}.rna_seq.augmented_star_gene_counts.tsv"
        if not path.exists():
            raise SystemExit(f"{path} missing — run `download` first.")
        star = read_star_counts(path)
        # Same key derivation as build_gene_index, or the matrix and the index
        # disagree about what a gene is called.
        star = star.assign(_key=[strip_version(g)[0] for g in star["gene_id"].astype(str)])
        star = star.loc[~star["_key"].str.endswith(PAR_Y_SUFFIX)]
        count_rows[row.barcode] = pd.Series(
            star["unstranded"].to_numpy(), index=star["_key"].to_numpy()
        )
        tpm_rows[row.barcode] = pd.Series(
            star["tpm_unstranded"].to_numpy(), index=star["_key"].to_numpy()
        )
    return pd.DataFrame(count_rows).T, pd.DataFrame(tpm_rows).T


def _reindex(
    counts: pd.DataFrame, tpm: pd.DataFrame, index_ids: list[str]
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Put both matrices on the shared index and report what that cost."""
    present = [g for g in index_ids if g in counts.columns]
    loss = pd.DataFrame(
        [
            {"step": "genes in the matrix", "n": counts.shape[1]},
            {"step": "genes on the shared index", "n": len(index_ids)},
            {"step": "on the index and in the matrix", "n": len(present)},
            {"step": "index genes absent from the matrix", "n": len(index_ids) - len(present)},
            {"step": "matrix genes not on the index", "n": counts.shape[1] - len(present)},
        ]
    )
    return counts.loc[:, present], tpm.loc[:, present], loss


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m src.bulk.ingest", description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    q = sub.add_parser("query", help="ask the GDC which STAR files exist")
    q.add_argument("--projects", nargs="+", default=list(PROJECTS))
    q.add_argument("--sample-types", nargs="+", default=list(DEFAULT_SAMPLE_TYPES))

    d = sub.add_parser("download", help="fetch the STAR files (resumable)")
    d.add_argument("--limit", type=int, default=None, help="first N files only, for a smoke test")

    g = sub.add_parser("gene-index", help="build the provisional shared index")
    g.add_argument("--version", default=PROVISIONAL_VERSION)

    b = sub.add_parser("build", help="assemble the matrices")
    b.add_argument("--version", default=PROVISIONAL_VERSION)

    parser.add_argument(
        "--allow-dirty", action="store_true",
        help="write from a dirty tree; the recorded sha will not reproduce it",
    )
    args = parser.parse_args(argv)

    global ALLOW_DIRTY

    ALLOW_DIRTY = args.allow_dirty
    if args.command == "query":
        step_query(tuple(args.projects), tuple(args.sample_types))
    elif args.command == "download":
        step_download(args.limit)
    elif args.command == "gene-index":
        step_gene_index(args.version)
    elif args.command == "build":
        step_build(args.version)
    return 0


if __name__ == "__main__":
    sys.exit(main())
