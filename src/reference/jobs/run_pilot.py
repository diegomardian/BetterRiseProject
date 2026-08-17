#!/usr/bin/env python
"""Run the five-patient pilot through the stages that exist. W1, week 2.

    python src/reference/jobs/run_pilot.py

The pilot is the critical path: W2's harness is blocked on it, and §8.3 moves the
week-5 gate to week 7 if it slips past week 3. It is a thin vertical slice —
every stage runs, none is tuned.

Stages implemented here:

  load      the five patients, reading only their columns off disk
  qc        per-cell metrics, per-batch thresholds with rationale, applied
  ambient   contamination from impossible genes, per sample

Stages still to be written (weeks 2-4), so the slice is not yet complete:

  malignancy   inferCNV / CopyKAT
  labels       axes 1 and 2 across four granularity rungs
  signature    _select_markers, then the pilot S matrix for W2

Artifacts land in results/{date}_{sha}/ with a provenance sidecar.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from src.common.io import write_versioned_table
from src.common.provenance import DEFAULT_SEED, set_global_seeds
from src.reference.ambient import contamination_by_sample, soup_profile_from_cells
from src.reference.ingest import (
    read_gse178341,
    read_gse178341_clusters,
    read_gse178341_metadata,
)
from src.reference.qc import (
    apply_qc,
    cell_qc_metrics,
    differential_retention,
    qc_summary,
    qc_thresholds,
)

PILOT = ["C122", "C165", "C107", "C138", "C162"]

#: Compartment label from the authors' cluster file. Safe to use for the
#: contamination mask and for the S matrix's non-epithelial columns —
#: distinguishing epithelium from immune from stroma does not depend on the
#: differentiation markers under test. NOT safe for deciding which epithelial
#: cells are mature; W1 builds those labels from the frozen axes in weeks 3-4.
#:
#: The first run used an EPCAM>=1 gate instead and returned contamination of
#: exactly 1.0 on eight samples: ambient EPCAM is everywhere, so the "epithelial"
#: mask was effectively all cells and the estimator degenerated.
EPITHELIAL_COMPARTMENT = "Epi"

#: Below this many epithelial cells a per-sample contamination estimate is noise.
MIN_CELLS_FOR_CONTAMINATION = 20


def main() -> int:
    set_global_seeds(DEFAULT_SEED)
    data = Path(os.environ.get("BRP_DATA_DIR", "data")) / "raw" / "GSE178341"
    h5 = data / "GSE178341_crc10x_full_c295v4_submit.h5"
    meta_csv = data / "GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz"
    cluster_csv = data / "GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz"
    if not h5.exists():
        raise SystemExit(f"missing {h5}")

    print("=" * 70)
    print(f"PILOT: {', '.join(PILOT)}")
    print("=" * 70)

    print("\n--- load ---")
    adata = read_gse178341(h5, patients=PILOT)
    print(f"{adata.n_obs:,} cells x {adata.n_vars:,} genes")
    metadata = read_gse178341_metadata(meta_csv)
    shared = adata.obs.index.intersection(metadata.index)
    for column in ("PROCESSING_TYPE", "MMRStatus", "MLH1Status"):
        if column in metadata.columns:
            adata.obs[column] = metadata.loc[shared, column].reindex(adata.obs.index)
    print(adata.obs.groupby(["patient_id", "tissue"], observed=True).size().to_string())

    print("\n--- qc ---")
    metrics = cell_qc_metrics(
        adata.X,
        adata.var["gene_symbol"],
        batch=adata.obs["sample_id"],
        patient_id=adata.obs["patient_id"],
        tissue=adata.obs["tissue"],
    )
    thresholds = qc_thresholds(metrics)
    passes = apply_qc(metrics, thresholds)
    summary = qc_summary(metrics, passes)
    print(f"retained {int(passes.sum()):,} of {len(passes):,} cells "
          f"({passes.mean():.1%})")
    print(summary.to_string(index=False))
    print("\nthresholds (per batch, with rationale — the week-1 deliverable):")
    print(thresholds[["batch", "metric", "lower", "upper", "n_cells", "n_failed"]]
          .to_string(index=False))

    print("\n--- differential retention: does QC cut the two arms unequally? ---")
    retention = differential_retention(metrics, passes)
    print(retention.to_string(index=False))
    if retention["flagged"].any():
        print(
            "\n!! QC removes cells at materially different rates between tumour and\n"
            "   normal in the flagged patients. The compositional term is the\n"
            "   difference between those two arms, so this biases it directly.\n"
            "   See docs/open_decisions.md #12 before trusting any composition."
        )

    print("\n--- ambient contamination (impossible genes) ---")
    symbols = list(adata.var["gene_symbol"])
    contamination = pd.DataFrame()
    if not cluster_csv.exists():
        print(f"!! missing {cluster_csv}; cannot build the epithelial mask")
    else:
        clusters = read_gse178341_clusters(cluster_csv)
        adata.obs["clTopLevel"] = clusters["clTopLevel"].reindex(adata.obs.index)
        print("compartments present:")
        print(adata.obs["clTopLevel"].value_counts().to_string())
        epithelial = (
            adata.obs["clTopLevel"].astype(str) == EPITHELIAL_COMPARTMENT
        ).to_numpy() & np.asarray(passes, dtype=bool)
        print(f"\nepithelial cells passing QC: {int(epithelial.sum()):,}")
        contamination = contamination_by_sample(
            adata.X, symbols,
            sample_id=adata.obs["sample_id"], cell_mask=epithelial,
        )
        # Sorting matters here: in a CD45-enriched sample the few epithelial
        # cells sit in immune-dominated soup, so a high estimate is expected
        # rather than a failure. And an estimate from <20 cells is noise.
        sorting = (
            adata.obs.groupby("sample_id", observed=True)["PROCESSING_TYPE"]
            .agg(lambda s: s.astype(str).mode().iat[0] if len(s) else "")
        )
        contamination["processing_type"] = (
            contamination["sample_id"].map(sorting).fillna("")
        )
        contamination["reliable"] = contamination["n_cells"] >= MIN_CELLS_FOR_CONTAMINATION
        print(contamination.to_string(index=False))

        usable = contamination[
            contamination["reliable"] & (contamination["processing_type"] == "unsorted")
        ]
        if len(usable):
            print(
                f"\nunsorted samples with >={MIN_CELLS_FOR_CONTAMINATION} epithelial "
                f"cells (n={len(usable)}): median contamination "
                f"{usable['contamination'].median():.1%}, "
                f"max {usable['contamination'].max():.1%}"
            )

        soup = soup_profile_from_cells(adata.X, symbols).sort_values(ascending=False)
        print("\ntop 10 soup genes (pooled across the pilot):")
        print(soup.head(10).to_string())

    print("\n--- pct_mito distribution, for open decision #12 ---")
    print("Pick the cap from this, not from convention. 20% is a lymphocyte number;")
    print("colonic epithelium runs higher, and much of it here is ambient.")
    if "clTopLevel" in adata.obs.columns:
        mito = metrics.copy()
        mito["compartment"] = adata.obs["clTopLevel"].astype(str).to_numpy()
        table = mito.groupby(["compartment", "tissue"], observed=True)["pct_mito"].describe(
            percentiles=[0.5, 0.75, 0.9, 0.95]
        )
        print(table[["count", "50%", "75%", "90%", "95%", "max"]].to_string())
        epi = mito[mito["compartment"] == EPITHELIAL_COMPARTMENT]
        for cap in (20, 30, 40, 50):
            kept = (epi["pct_mito"] <= cap).mean()
            by_tissue = epi.groupby("tissue", observed=True)["pct_mito"].apply(
                lambda s, c=cap: (s <= c).mean()
            )
            gap = abs(by_tissue.get("tumour", np.nan) - by_tissue.get("normal", np.nan))
            print(f"  cap {cap:>3}%: epithelium kept {kept:.1%}, tumour/normal gap {gap:.1%}")

    print("\n--- compartments available for the S matrix ---")
    print("§2.1 error 3 requires stromal, immune AND endothelial columns.")
    if "clMidwayPr" in locals().get("clusters", pd.DataFrame()).columns:
        midway = clusters["clMidwayPr"].reindex(adata.obs.index)
        print(midway.value_counts().to_string())

    print("\n--- what is in the cluster file (for the real labels) ---")
    if cluster_csv.exists():
        head = pd.read_csv(cluster_csv, nrows=3)
        print("columns:", list(head.columns))
        print(head.to_string())
    else:
        print(f"missing {cluster_csv}")

    print("\n--- writing artifacts ---")
    for frame, name in (
        (thresholds, "pilot_qc_thresholds"),
        (summary, "pilot_qc_summary"),
        (retention, "pilot_differential_retention"),
        (contamination, "pilot_contamination"),
    ):
        if len(frame):
            path = write_versioned_table(
                frame, name, seed=DEFAULT_SEED, allow_dirty=True,
                notes=f"five-patient pilot: {', '.join(PILOT)}",
            )
            print(f"  {path}")

    print("\nNEXT: malignancy calls, then labels, then the pilot S matrix for W2.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
