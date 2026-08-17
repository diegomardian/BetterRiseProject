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
    read_gse178341_metadata,
)
from src.reference.qc import apply_qc, cell_qc_metrics, qc_summary, qc_thresholds

PILOT = ["C122", "C165", "C107", "C138", "C162"]

#: Provisional epithelial call for the contamination estimate only. Real labels
#: are weeks 3-4 and must come from the frozen axes; this is a crude EPCAM+
#: gate so the ambient measurement can run on the pilot at all. EPCAM is not a
#: panel gene, so it does not violate invariant 2 — but nothing downstream may
#: use this mask.
EPITHELIAL_MARKER = "EPCAM"
EPITHELIAL_MIN_COUNTS = 1


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

    print("\n--- ambient contamination (impossible genes) ---")
    symbols = list(adata.var["gene_symbol"])
    if EPITHELIAL_MARKER not in symbols:
        print(f"!! {EPITHELIAL_MARKER} not in the matrix; skipping contamination")
        contamination = pd.DataFrame()
    else:
        column = symbols.index(EPITHELIAL_MARKER)
        epithelial = np.asarray(
            adata.X[:, column].todense()
        ).ravel() >= EPITHELIAL_MIN_COUNTS
        epithelial &= np.asarray(passes, dtype=bool)
        print(f"provisional EPCAM+ epithelial cells: {int(epithelial.sum()):,} "
              f"(crude gate, for this measurement only)")
        contamination = contamination_by_sample(
            adata.X, symbols,
            sample_id=adata.obs["sample_id"], cell_mask=epithelial,
        )
        print(contamination.to_string(index=False))
        soup = soup_profile_from_cells(adata.X, symbols).sort_values(ascending=False)
        print("\ntop 10 soup genes (pooled across the pilot):")
        print(soup.head(10).to_string())

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
