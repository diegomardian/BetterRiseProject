"""SoupX and DecontX on ONE sample, end to end. W1, week 2.

**Run this before running anything on 84 samples.** Neither correction method
has touched real data — the inferCNV stage cost four cluster runs to wrong
assumptions about an R API, and the cheapest place to find the next one is a
single sample.

Checks, in order of what they tell you:

1. Does the R side execute at all?
2. Does each method return a retention table with sane values — bounded by 1,
   not all identical, not all NaN?
3. **Do the two methods agree on which genes are soup?** That is the week-2
   deliverable, and it is the first point where the answer is scientific rather
   than mechanical.
4. Does either method's contamination estimate resemble the impossible-gene
   estimate already computed for this sample by an unrelated route?

    python src/reference/jobs/run_ambient_pilot.py --patient C122
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from src.common.provenance import DEFAULT_SEED, set_global_seeds  # noqa: E402
from src.reference.ambient import (  # noqa: E402
    compare_retention,
    contamination_by_sample,
    retention_agreement,
    run_decontx,
    run_soupx,
    soup_profile_from_cells,
)
from src.reference.ingest import (  # noqa: E402
    assign_compartments,
    read_gse178341,
    read_gse178341_clusters,
)


def main() -> int:
    set_global_seeds(DEFAULT_SEED)
    parser = argparse.ArgumentParser()
    parser.add_argument("--patient", default="C122")
    parser.add_argument("--sample", default=None,
                        help="one sample_id; default is the patient's first")
    args = parser.parse_args()

    data = Path(os.environ.get("BRP_DATA_DIR", "data")) / "raw" / "GSE178341"
    out_root = Path(os.environ.get("BRP_DATA_DIR", "data")) / "interim" / "ambient"

    adata = read_gse178341(
        data / "GSE178341_crc10x_full_c295v4_submit.h5", patients=[args.patient]
    )
    clusters = read_gse178341_clusters(
        data / "GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz"
    )
    labels = clusters["cl295v11SubFull"].reindex(adata.obs.index)
    compartment = assign_compartments(clusters).reindex(adata.obs.index)

    sample = args.sample or sorted(adata.obs["sample_id"].unique())[0]
    here = (adata.obs["sample_id"] == sample).to_numpy()
    if not here.any():
        raise SystemExit(
            f"{sample} not in {args.patient}: "
            f"{sorted(adata.obs['sample_id'].unique())}"
        )
    print(f"{sample} — {int(here.sum()):,} cells")

    matrix = adata.X[here]
    genes = adata.var["gene_symbol"]
    barcodes = adata.obs.index[here]
    cell_clusters = labels.to_numpy()[here]
    if pd.isna(cell_clusters).any():
        keep = ~pd.isna(cell_clusters)
        print(f"  dropping {int((~keep).sum()):,} cells with no cluster label")
        matrix, barcodes = matrix[keep], barcodes[keep]
        cell_clusters = cell_clusters[keep]
        here_idx = np.flatnonzero(here)[keep]
    else:
        here_idx = np.flatnonzero(here)
    print(f"  {len(set(map(str, cell_clusters)))} clusters")

    # The soup profile SoupX needs, computed from cells because no empty
    # droplets exist (decision #8).
    profile = soup_profile_from_cells(matrix, genes)

    # The independent estimate, for step 4's cross-check.
    epithelial = (compartment.to_numpy()[here_idx] == "epithelial")
    independent = float("nan")
    if epithelial.sum() >= 20:
        est = contamination_by_sample(
            matrix, genes, sample_id=[sample] * matrix.shape[0],
            cell_mask=epithelial,
        )
        independent = float(est["contamination"].iloc[0])
    print(f"  impossible-gene estimate: {independent:.1%}")

    results = {}
    for name, fn in (("soupx", run_soupx), ("decontx", run_decontx)):
        out_dir = out_root / sample / name
        print(f"\n=== {name} ===")
        kwargs = dict(barcodes=barcodes, clusters=cell_clusters, out_dir=out_dir)
        if name == "soupx":
            kwargs["soup_profile"] = profile
        try:
            fn(matrix, genes, **kwargs)
        except Exception as exc:  # noqa: BLE001 — report, do not abort the other
            print(f"!! {name} failed: {exc}")
            continue
        path = out_dir / f"{name}_retention.csv"
        if not path.exists():
            print(f"!! {name} produced no retention table at {path}")
            continue
        frame = pd.read_csv(path)
        results[name] = frame
        finite = frame["retention"].dropna()
        print(f"  retention: {len(finite):,} genes, median "
              f"{finite.median():.3f}, min {finite.min():.3f}")
        if finite.nunique() <= 1:
            print("  !! every gene has identical retention — the correction "
                  "did nothing,\n     or the matrix came back unchanged")

    if len(results) < 2:
        raise SystemExit(
            "\nonly one method produced a table — nothing to compare. Fix that "
            "before running the cohort."
        )

    comparison = compare_retention(
        results["soupx"], results["decontx"], sample_id=sample
    )
    agreement = retention_agreement(comparison)
    print("\n" + "=" * 60)
    print("DO THE TWO METHODS AGREE ON WHICH GENES ARE SOUP?")
    print("=" * 60)
    for key, value in agreement.items():
        print(f"  {key:<28} {value}")
    print(
        f"\n  independent (impossible-gene) estimate: {independent:.1%}\n"
        "  Three routes to one quantity. Agreement is reassuring; disagreement "
        "is a\n  result, not something to reconcile — and decision #16 reports "
        "rather than\n  corrects precisely because no one of them is "
        "authoritative."
    )
    print("\n  hardest-stripped genes by SoupX:")
    print(comparison.nsmallest(10, "retention_soupx")[
        ["gene", "retention_soupx", "retention_decontx"]
    ].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
