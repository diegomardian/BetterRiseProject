"""D2's outcome-blind GUCA2A variation gate.

    python -m src.bulk.d2_feasibility_gate

This job deliberately has no clinical or survival input.  It answers only the
question fixed in ``docs/d2_feasibility_gate.md``: whether the committed
primary-tumour expression matrix contains enough GUCA2A variation to justify
writing a separate marker-to-survival pre-specification.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from src.bulk.gdc import read_manifest
from src.bulk.gene_index import resolve_symbols
from src.bulk.normalise import assert_log_scale
from src.common.io import write_versioned_table
from src.common.paths import CONFIG_DIR, PROCESSED_DIR, RESULTS_DIR
from src.common.provenance import DEFAULT_SEED

log = logging.getLogger(__name__)

BULK = PROCESSED_DIR / "bulk"
INDEX_VERSION = "1.0.0"
GENE = "GUCA2A"
PRIMARY_TUMOUR = "01"
MIN_PARTICIPANTS = 100
MIN_IQR = 0.5
MIN_TAIL_FRACTION = 0.10
MIN_TAIL_N = 25


class D2FeasibilityError(RuntimeError):
    """The outcome-blind gate cannot evaluate the committed input honestly."""


def resolve_gene(gene_map: pd.DataFrame, symbol: str = GENE) -> str:
    """Resolve one symbol exactly once; ambiguity is a stop, never a choice."""
    required = {"gene_symbol", "ensembl_id"}
    missing = required - set(gene_map.columns)
    if missing:
        raise D2FeasibilityError(f"gene map is missing required columns: {sorted(missing)}")
    resolved, unmapped, ambiguous = resolve_symbols(gene_map, [symbol])
    if unmapped or ambiguous:
        raise D2FeasibilityError(
            f"{symbol} must resolve one-to-one: unmapped={unmapped}, ambiguous={ambiguous}"
        )
    return resolved[symbol]


def primary_tumour_values(
    expression: pd.DataFrame, manifest: pd.DataFrame, gene_id: str
) -> pd.Series:
    """One primary-tumour GUCA2A value per TCGA participant.

    The manifest was already deduplicated by ingest.  A duplicate here is an
    upstream contract failure, not permission to choose the more favourable
    aliquot after looking at expression.
    """
    required = {"barcode", "patient_id", "sample_type"}
    missing = required - set(manifest.columns)
    if missing:
        raise D2FeasibilityError(f"sample manifest is missing required columns: {sorted(missing)}")
    if gene_id not in expression.columns:
        raise D2FeasibilityError(f"{GENE} ({gene_id}) is absent from the expression matrix")
    if expression.index.has_duplicates:
        raise D2FeasibilityError("expression index has duplicate barcodes")
    if manifest["barcode"].duplicated().any():
        raise D2FeasibilityError("sample manifest has duplicate barcodes")

    manifest_primary = manifest[manifest["sample_type"] == PRIMARY_TUMOUR]
    if manifest_primary.empty:
        raise D2FeasibilityError("sample manifest has no primary-tumour (sample_type '01') rows")
    if manifest_primary["patient_id"].duplicated().any():
        duplicate_ids = manifest_primary.loc[
            manifest_primary["patient_id"].duplicated(), "patient_id"
        ].tolist()
        raise D2FeasibilityError(
            "more than one primary-tumour sample for a participant; ingest must "
            f"deduplicate before D2, e.g. {duplicate_ids[:3]}"
        )
    missing_primary = manifest_primary.loc[
        ~manifest_primary["barcode"].isin(expression.index), "barcode"
    ].tolist()
    if missing_primary:
        raise D2FeasibilityError(
            f"{len(missing_primary)} primary-tumour manifest barcode(s) are absent from "
            f"expression, e.g. {missing_primary[:3]}"
        )

    annotations = manifest.set_index("barcode").reindex(expression.index)
    if annotations["patient_id"].isna().any():
        missing_barcodes = expression.index[annotations["patient_id"].isna()].tolist()
        raise D2FeasibilityError(
            f"{len(missing_barcodes)} expression barcode(s) are absent from the sample manifest, "
            f"e.g. {missing_barcodes[:3]}"
        )
    primary = annotations[annotations["sample_type"] == PRIMARY_TUMOUR]
    if primary.empty:
        raise D2FeasibilityError("no primary-tumour (sample_type '01') rows overlap expression")

    values = expression.loc[primary.index, gene_id].astype(float)
    if values.isna().any() or not np.isfinite(values.to_numpy()).all():
        raise D2FeasibilityError(f"{GENE} has missing or non-finite primary-tumour values")
    values.index = primary["patient_id"].astype(str)
    values.name = GENE
    return values


def evaluate(values: pd.Series, *, gene_id: str) -> pd.DataFrame:
    """Evaluate the three pre-committed D2 conditions as one one-row table."""
    if values.index.has_duplicates:
        raise D2FeasibilityError("D2 values are not one row per participant")
    x = values.to_numpy(dtype=float)
    if not len(x):
        raise D2FeasibilityError("D2 has no primary-tumour values")
    q25, median, q75 = np.percentile(x, [25, 50, 75])
    iqr = float(q75 - q25)
    n_below = int((x < median).sum())
    n_above = int((x > median).sum())
    n = int(len(x))
    lower_required = max(MIN_TAIL_N, int(np.ceil(MIN_TAIL_FRACTION * n)))
    upper_required = lower_required
    conditions = {
        "n": n >= MIN_PARTICIPANTS,
        "iqr": iqr >= MIN_IQR,
        "below_median": n_below >= lower_required,
        "above_median": n_above >= upper_required,
    }
    verdict = "PASS" if all(conditions.values()) else "STOP"
    failed = [name for name, passed in conditions.items() if not passed]
    detail = (
        "all pre-committed conditions passed; D2 may be pre-specified separately"
        if verdict == "PASS"
        else "failed pre-committed condition(s): " + ", ".join(failed)
    )
    return pd.DataFrame([{
        "gene": GENE,
        "ensembl_id": gene_id,
        "n_participants": n,
        "minimum": float(np.min(x)),
        "q25": float(q25),
        "median": float(median),
        "q75": float(q75),
        "maximum": float(np.max(x)),
        "iqr": iqr,
        "n_exact_zero": int((x == 0).sum()),
        "n_below_median": n_below,
        "n_above_median": n_above,
        "required_tail_n": lower_required,
        "condition_n": conditions["n"],
        "condition_iqr": conditions["iqr"],
        "condition_below_median": conditions["below_median"],
        "condition_above_median": conditions["above_median"],
        "verdict": verdict,
        "detail": detail,
    }])


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--expression", type=Path,
        default=BULK / f"tcga_log2cpm_{INDEX_VERSION}.parquet",
    )
    parser.add_argument("--manifest", type=Path, default=BULK / "sample_manifest.tsv")
    parser.add_argument(
        "--gene-index-map", type=Path,
        default=CONFIG_DIR / "gene_index" / f"gene_index_{INDEX_VERSION}.map.tsv",
    )
    parser.add_argument("--results-dir", type=Path, default=RESULTS_DIR)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    for path, producer in (
        (args.expression, "python -m src.bulk.ingest build"),
        (args.manifest, "python -m src.bulk.ingest build"),
        (args.gene_index_map, "the committed 1.0.0 gene-index map"),
    ):
        if not path.exists():
            raise SystemExit(f"{path} not found; obtain it from {producer}.")

    expression = pd.read_parquet(args.expression)
    assert_log_scale(expression, context=args.expression.name)
    manifest = read_manifest(args.manifest)
    gene_map = pd.read_csv(args.gene_index_map, sep="\t", dtype={"ensembl_id": "string"})
    gene_id = resolve_gene(gene_map)
    values = primary_tumour_values(expression, manifest, gene_id)
    table = evaluate(values, gene_id=gene_id)
    row = table.iloc[0]
    log.info(
        "%s: %s (n=%d, IQR=%.4f, below/above median=%d/%d)",
        GENE, row["verdict"], row["n_participants"], row["iqr"],
        row["n_below_median"], row["n_above_median"],
    )
    log.info("%s", row["detail"])

    path = write_versioned_table(
        table,
        "d2_guca2a_feasibility",
        seed=args.seed,
        results_dir=args.results_dir,
        allow_dirty=args.allow_dirty,
        extra_meta={
            "contract": "docs/d2_feasibility_gate.md",
            "outcomes_read": False,
            "expression_scale": "log2(CPM+1)",
            "unit": "one primary-tumour RNA sample per TCGA participant",
            "conditions": {
                "min_participants": MIN_PARTICIPANTS,
                "min_iqr_log2_cpm": MIN_IQR,
                "min_tail_fraction": MIN_TAIL_FRACTION,
                "min_tail_n": MIN_TAIL_N,
            },
        },
    )
    log.info("wrote %s", path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
