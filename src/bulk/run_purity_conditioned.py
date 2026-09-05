"""W3.2 re-run conditioned on purity (W3.3). The driver that was never written.

    python -m src.bulk.run_purity_conditioned

`purity_conditioned_check` and `purity_association` have existed in
``src/bulk/premise.py`` since W3.3 and are covered by tests. Nothing ever called
them. Two tables reached ``results/`` from an uncommitted script and have sat
there since 2026-08-18 as the only copies of their numbers, stamped `git_dirty`
against a branch that no longer exists.

That is Appendix A item 3 in the WMHS paper -- a result and the code that
produced it travelling separately -- in its last remaining instance in this
repository. This closes it.

WHY THEY COULD NOT SIMPLY BE DELETED. The obvious housekeeping move is to drop
every dirty table that has a clean twin. Sixteen of the eighteen do. These two do
not, and they are the same two with no producer, so deleting them as part of
that sweep would have destroyed the only copies of results nothing could
regenerate. The producer comes first; the delete comes after it has run.

WHAT THE ASSOCIATION TABLE SHOWS, and why it matters beyond housekeeping.
Against copy-number ABSOLUTE, purity explains r-squared 0.042 of bulk CDX2 and
0.019 of GUCA2A. Against the expression-derived ESTIMATE score, CDX2 rises to
0.104 -- two and a half times as much, on the same tumours, for the same gene.
That gap is the circularity Stage 4's instrument gate is built to exclude:
correlate an expression-derived quantity with an expression-derived purity call
and part of what comes back is the shared derivation. Same reason
``src/bulk/instrument.py`` filters on `expression_derived` rather than on a
method name.
"""

from __future__ import annotations

import argparse
import logging
import sys
from collections.abc import Sequence
from pathlib import Path

import pandas as pd

from src.bulk.premise import purity_association, purity_conditioned_check
from src.common.io import write_versioned_table
from src.common.paths import PROCESSED_DIR, RESULTS_DIR

log = logging.getLogger(__name__)

BULK = PROCESSED_DIR / "bulk"
INDEX_VERSION = "1.0.0"

#: W3.2's two genes. The premise check is about bulk GUCA2A; CDX2 travels with
#: it as the contrast, exactly as in the original table.
GENES: tuple[str, ...] = ("GUCA2A", "CDX2")

#: The seed the committed tables were written under. Kept so a re-run is
#: comparable to them rather than merely similar -- `assess` bootstraps.
SEED = 20260818


def resolve_genes(symbols: Sequence[str]) -> dict[str, str]:
    """Symbol -> unversioned Ensembl id, refusing an ambiguous mapping."""
    from src.bulk.gene_index import load_gene_index_map, resolve_symbols

    resolved, unmapped, ambiguous = resolve_symbols(
        load_gene_index_map(INDEX_VERSION), list(symbols)
    )
    if unmapped or ambiguous:
        raise SystemExit(
            f"cannot resolve {list(symbols)}: unmapped={unmapped}, "
            f"ambiguous={ambiguous}. An ambiguous symbol is a decision, not a "
            f"tie-break -- pin it here and write down why."
        )
    return resolved


def purity_series(purity: pd.DataFrame, method: str) -> pd.Series:
    """One purity call per barcode for a named method."""
    rows = purity[purity["method"] == method].dropna(subset=["purity"])
    return rows.set_index("barcode")["purity"]


def newest_purity_table() -> Path | None:
    matches = sorted(RESULTS_DIR.glob("*/tcga_purity.parquet"))
    return matches[-1] if matches else None


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expression", type=Path,
                        default=BULK / f"tcga_log2cpm_{INDEX_VERSION}.parquet")
    parser.add_argument("--manifest", type=Path, default=BULK / "sample_manifest.tsv")
    parser.add_argument("--purity", type=Path, default=None)
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    if not args.expression.exists():
        raise SystemExit(
            f"{args.expression} not found. The TCGA matrices are gitignored and "
            f"live with the data -- on the cluster, under "
            f"$BRP_DATA_DIR/processed/bulk. Build them with\n"
            f"    python -m src.bulk.ingest build"
        )
    if not args.manifest.exists():
        raise SystemExit(f"{args.manifest} not found (written by `ingest build`).")

    purity_path = args.purity or newest_purity_table()
    if purity_path is None:
        raise SystemExit(
            "no results/*/tcga_purity.parquet. Run `python -m src.bulk.run_purity "
            "run` first -- this job conditions on its output."
        )

    from src.bulk.ingest import read_manifest
    from src.bulk.normalise import assert_log_scale

    expression = pd.read_parquet(args.expression)
    assert_log_scale(expression, context=f"{args.expression.name}")
    manifest = read_manifest(args.manifest).set_index("barcode")
    purity = pd.read_parquet(purity_path)
    gene_ids = resolve_genes(GENES)

    log.info("expression %d x %d | purity %s", *expression.shape, purity_path.name)
    log.info("genes: %s", ", ".join(f"{s} ({i})" for s, i in gene_ids.items()))

    methods = sorted(purity["method"].unique())
    conditioned, association = [], []
    for method in methods:
        series = purity_series(purity, method)
        if len(series) < 30:
            log.warning("  %s: only %d calls; skipping", method, len(series))
            continue
        log.info("  %s: %d calls", method, len(series))
        conditioned.append(purity_conditioned_check(
            expression, manifest, gene_ids, series, method=method, seed=args.seed,
        ))
        association.append(purity_association(
            expression, manifest, gene_ids, series, method=method,
        ))

    if not conditioned:
        log.error("no purity method had enough calls to condition on.")
        return 2

    conditioned_table = pd.concat(conditioned, ignore_index=True)
    association_table = pd.concat(association, ignore_index=True)

    log.info("\npurity's share of each gene's variance:")
    log.info("%s", association_table.to_string(index=False))
    expression_derived = purity[purity["expression_derived"].astype(bool)]
    if not expression_derived.empty:
        log.info(
            "\nRead the r_squared column across methods, not down it. An "
            "expression-derived purity call shares a derivation with the "
            "expression it is being correlated against, so part of what it "
            "reports is that shared derivation rather than confounding."
        )

    for frame, name, note in (
        (conditioned_table, "tcga_premise_purity_conditioned",
         "W3.2 re-run conditioned on purity (W3.3). ABSOLUTE primary."),
        (association_table, "tcga_purity_expression_association",
         "How much of each gene's variance purity explains."),
    ):
        path = write_versioned_table(
            frame, name, seed=args.seed,
            results_dir=args.results_dir, allow_dirty=args.allow_dirty,
            notes=note,
            extra_meta={
                "genes": list(GENES),
                "scale": "log2(CPM+1)",
                "purity_adjusted": True,
                "purity_methods": methods,
                "purity_source": purity_path.name,
                "supersedes": (
                    "results/2026-08-18_7c49e99/, which was written by an "
                    "uncommitted script and stamped git_dirty against a branch "
                    "that no longer exists. This job is that script, committed."
                ),
                "what_this_answers": (
                    "Whether bulk GUCA2A's bimodality survives conditioning on "
                    "tumour purity, and how much of each gene's variance purity "
                    "explains at all -- the number that says whether "
                    "conditioning was ever going to matter."
                ),
            },
        )
        log.info("wrote %s", path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
