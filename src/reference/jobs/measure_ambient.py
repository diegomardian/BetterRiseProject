"""Ambient contamination across the whole cohort. W1, week 2.

**Step one of the ambient stage, and it needs no R.** The correction methods
(SoupX, DecontX) are expensive and their output has to be interpreted; this
measures the quantity they would correct, on all 62 patients, using the
impossible-gene estimator that needs neither empty droplets nor a clustering.

The number it produces decides how much the rest of the stage matters. On the
five-patient pilot the median unsorted sample sat at **1.6%** contamination. If
that holds cohort-wide, correction is removing very little — and the decision
recorded for this stage is to **measure and report rather than correct**,
because DecontX defines contamination as counts resembling other clusters and
can therefore absorb genuine low-level expression of a marker in a rare
population, which is precisely this project's signal.

If instead some samples come back at 10-20%, those samples are the finding, and
they should be excluded or flagged rather than corrected into looking fine.

Reads only the barcodes, features and per-sample columns it needs, so it is
laptop-scale despite the 9 GB deposit.

    python src/reference/jobs/measure_ambient.py
    python src/reference/jobs/measure_ambient.py --patients C122 C165
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

import pandas as pd  # noqa: E402

from src.common.io import write_versioned_table  # noqa: E402
from src.common.panel import tier_genes  # noqa: E402
from src.common.provenance import DEFAULT_SEED, set_global_seeds  # noqa: E402
from src.reference.ambient import contamination_by_sample  # noqa: E402
from src.reference.ingest import (  # noqa: E402
    assign_compartments,
    read_gse178341,
    read_gse178341_clusters,
    read_gse178341_index,
    read_gse178341_metadata,
)

#: Only unsorted samples carry an interpretable contamination estimate. A
#: CD45-enriched sample has had its epithelial fraction physically reduced, so
#: the impossible-gene ratio is measuring the enrichment, not the soup
#: (open decision #11).
UNSORTED = "unsorted"

#: Below this many epithelial cells the per-sample ratio is too noisy to read.
MIN_EPITHELIAL = 20


def main() -> int:
    set_global_seeds(DEFAULT_SEED)
    parser = argparse.ArgumentParser()
    parser.add_argument("--patients", nargs="*", default=None)
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args()

    data = Path(os.environ.get("BRP_DATA_DIR", "data")) / "raw" / "GSE178341"
    h5 = data / "GSE178341_crc10x_full_c295v4_submit.h5"
    clusters = read_gse178341_clusters(
        data / "GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz"
    )
    metadata = read_gse178341_metadata(
        data / "GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz"
    )

    obs, _var = read_gse178341_index(h5)
    patients = args.patients or sorted(obs["patient_id"].unique())
    print(f"{len(patients)} patients")

    # Per patient, not all at once: the full matrix is 9 GB and the estimator
    # only ever needs one sample's cells.
    frames = []
    for i, patient in enumerate(patients, 1):
        adata = read_gse178341(h5, patients=[patient])
        compartment = assign_compartments(clusters).reindex(adata.obs.index)
        epithelial = (compartment == "epithelial").to_numpy()
        if not epithelial.any():
            print(f"[{i}/{len(patients)}] {patient} — no epithelium, skipped")
            continue

        frame = contamination_by_sample(
            adata.X, adata.var["gene_symbol"],
            sample_id=adata.obs["sample_id"],
            cell_mask=epithelial,
        )
        frame["patient_id"] = patient
        frames.append(frame)
        usable = frame[frame["n_cells"] >= MIN_EPITHELIAL]["contamination"]
        median = float(usable.median()) if len(usable) else float("nan")
        print(f"[{i}/{len(patients)}] {patient} — {len(frame)} samples, "
              f"median {median:.1%}")

    if not frames:
        raise SystemExit("no patient produced an estimate")
    out = pd.concat(frames, ignore_index=True)

    # PROCESSING_TYPE decides which rows are interpretable at all.
    processing = (
        metadata[["PROCESSING_TYPE"]].reset_index().drop_duplicates()
        if "PROCESSING_TYPE" in metadata.columns else None
    )
    if processing is not None:
        key = processing.columns[0]
        out = out.merge(
            processing.rename(columns={key: "sample_id"}), on="sample_id", how="left"
        )

    reliable = out[
        (out["n_cells"] >= MIN_EPITHELIAL)
        & (out.get("PROCESSING_TYPE", UNSORTED) == UNSORTED)
        & out["contamination"].notna()
    ]

    print("\n" + "=" * 60)
    print("CONTAMINATION, UNSORTED SAMPLES WITH ENOUGH EPITHELIUM")
    print("=" * 60)
    if len(reliable):
        print(reliable["contamination"].describe(
            percentiles=[0.25, 0.5, 0.75, 0.9, 0.95]
        ).to_string())
        worst = reliable.nlargest(5, "contamination")[
            ["sample_id", "n_cells", "contamination"]
        ]
        print("\nworst five:")
        print(worst.to_string(index=False))
        high = reliable[reliable["contamination"] > 0.10]
        if len(high):
            print(
                f"\n!! {len(high)} sample(s) above 10% contamination. Those are "
                f"the finding —\n   exclude or flag them rather than correcting "
                f"them into looking fine. At that\n   level the soup is a "
                f"material share of every epithelial cell's counts."
            )
        else:
            print(
                "\n   No sample above 10%. If the median is also low, correction "
                "removes very\n   little, and the decision to MEASURE AND REPORT "
                "rather than correct holds:\n   DecontX can absorb genuine "
                "low-level marker expression in a rare population,\n   which is "
                "this project's signal."
            )
    else:
        print("no reliable samples — check PROCESSING_TYPE joined correctly")

    print(f"\ntargets checked against (invariant 2): {sorted(tier_genes('A'))}")
    path = write_versioned_table(
        out, "ambient_contamination", seed=DEFAULT_SEED,
        allow_dirty=args.allow_dirty,
        notes=(
            "Per-sample ambient contamination from impossible genes, all "
            "patients. Needs no empty droplets and no clustering. Only unsorted "
            "samples with >=20 epithelial cells are interpretable (#11); the "
            "rest are emitted with their counts so the exclusion is visible."
        ),
        extra_meta={
            "n_samples": int(len(out)),
            "n_reliable": int(len(reliable)),
            "median_contamination": (
                float(reliable["contamination"].median()) if len(reliable)
                else None
            ),
        },
    )
    print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
