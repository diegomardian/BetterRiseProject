"""W3.4 — run the batch and technical structure analysis. CLI, idempotent.

    python -m src.bulk.run_batch

Tumours only. Normal-adjacent samples are excluded from the variance analysis
because tissue type would dominate every principal component and swamp the
technical factors this is trying to measure; the clinical variables are
tumour-level anyway.

**Nothing here corrects the expression matrix** (CLAUDE.md invariant 4). The
outputs are three tables and a coverage report.
"""

from __future__ import annotations

import sys

import pandas as pd

from src.bulk.batch import confounding_table, variance_explained
from src.bulk.clinical import (
    build_clinical_table,
    coverage_report,
    fetch_cases,
    fetch_msi,
)
from src.bulk.gdc import read_manifest
from src.common.io import write_versioned_table
from src.common.paths import PROCESSED_DIR, RAW_DIR

SEED = 20260818
BULK = PROCESSED_DIR / "bulk"
CACHE = RAW_DIR / "gdc_clinical"
INDEX_VERSION = "0.9.0"

#: From the aliquot barcode. TCGA has no explicit batch column; plate is the
#: standard processing-batch proxy and TSS the site proxy.
TECHNICAL = ("tss", "plate", "centre", "analyte", "vial")
CLINICAL = ("stage", "site", "msi_status", "project")


def main(argv: list[str] | None = None) -> int:
    del argv

    manifest = read_manifest(BULK / "sample_manifest.tsv")
    tumours = manifest.loc[manifest["sample_type"] == "01"].copy()
    print(f"{len(tumours)} tumour samples, {tumours['patient_id'].nunique()} patients")

    CACHE.mkdir(parents=True, exist_ok=True)
    clinical = build_clinical_table(
        fetch_cases(cache=CACHE / "cases.json"),
        fetch_msi(cache=CACHE / "msi.json"),
    )

    patients = sorted(tumours["patient_id"].unique())
    coverage = coverage_report(clinical, patients)
    print("\n=== clinical coverage, for the patients we have expression for ===")
    print(coverage.to_string(index=False))

    msi = clinical.loc[clinical["patient_id"].isin(patients), "msi_status"]
    print("\nMSI breakdown:", dict(msi.value_counts(dropna=False)))

    annotations = tumours.merge(clinical, on="patient_id", how="left", suffixes=("", "_clin"))
    annotations = annotations.set_index("barcode")

    print("\n=== technical factor level counts ===")
    for column in TECHNICAL:
        counts = annotations[column].value_counts()
        print(
            f"  {column:<9} {counts.size:>4} levels; "
            f"largest {int(counts.iloc[0])}, singletons {(counts == 1).sum()}"
        )

    print("\n=== confounding: technical x clinical (permutation) ===")
    confounding = confounding_table(annotations, TECHNICAL, CLINICAL, seed=SEED)
    with pd.option_context("display.width", 200):
        print(confounding.to_string(index=False))

    print("\n=== variance explained (PVCA-style) ===")
    expression = pd.read_parquet(BULK / f"tcga_log2cpm_{INDEX_VERSION}.parquet")
    expression = expression.loc[expression.index.isin(annotations.index)]
    factors = annotations.loc[expression.index, list(TECHNICAL) + list(CLINICAL)]
    variance, meta = variance_explained(expression, factors, seed=SEED)
    print(f"  {meta['n_pcs_retained']} PCs covering {meta['variance_covered']:.1%} of variance")
    with pd.option_context("display.width", 200):
        print(variance.to_string(index=False))

    for frame, name, note in (
        (coverage, "tcga_clinical_coverage", "W3.4 clinical annotation coverage"),
        (confounding, "tcga_batch_confounding", "W3.4 technical x clinical, permutation tests"),
        (variance, "tcga_batch_variance_explained", "W3.4 PVCA-style. NOT a partition."),
    ):
        write_versioned_table(
            frame, name=name, seed=SEED, notes=note,
            extra_meta=meta if name.endswith("variance_explained") else None,
            allow_dirty=True,
        )
    clinical.to_csv(BULK / "clinical_annotation_w3.4.tsv", sep="\t", index=False)
    print(f"\nwrote three results tables and {BULK / 'clinical_annotation_w3.4.tsv'}")
    print("No expression matrix was written. Invariant 4: document, do not correct.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
