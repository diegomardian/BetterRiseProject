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
from src.reference.ambient import (  # noqa: E402
    contamination_by_sample,
    differential_contamination,
)
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


def _report(out):
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

    # Third instance of one question: does this affect the two arms unequally?
    # Delta(mature fraction) IS the difference between them, so anything
    # asymmetric moves it directly.
    asym = differential_contamination(reliable if len(reliable) else out)
    paired = asym[asym["both_arms"]]
    print("\n" + "=" * 60)
    print("IS THE SOUP DIRTIER IN TUMOUR THAN IN MATCHED NORMAL?")
    print("=" * 60)
    if len(paired):
        print(paired.reindex(
            paired["difference"].abs().sort_values(ascending=False).index
        ).head(10).to_string(index=False))
        flagged = paired[paired["flagged"]]
        direction = float(paired["difference"].median())
        print(f"\n{len(flagged)} of {len(paired)} patients differ by more than "
              f"5 points; median tumour-minus-normal {direction:+.1%}")
        if len(flagged) > len(paired) / 4:
            print(
                "   Widespread asymmetry. Open decision #16's exclusion rule "
                "should then be\n   PER-PATIENT ON THE GAP rather than "
                "per-sample on the level: two equally\n   dirty arms are far "
                "less compromised than two arms ten points apart."
            )
    else:
        print("no patient has both arms among the interpretable samples")
    return reliable



def attach_processing_type(frame, obs, metadata):
    """Join PROCESSING_TYPE onto a per-sample frame.

    **Via obs, not directly.** ``read_gse178341_metadata`` is indexed by
    *barcode*, not by sample — a first version renamed that index to
    `sample_id`, matched nothing, and reported every sample as uninterpretable
    while the underlying numbers were fine. `obs` carries both, so joining
    through it is the only version that cannot silently miss.
    """
    if "PROCESSING_TYPE" not in metadata.columns:
        print("note: metadata has no PROCESSING_TYPE — cannot separate sorted "
              "from unsorted samples")
        return frame
    joined = obs.join(metadata[["PROCESSING_TYPE"]], how="left")
    per_sample = (
        joined.groupby("sample_id", observed=True)["PROCESSING_TYPE"]
        .agg(lambda x: x.dropna().iloc[0] if x.notna().any() else None)
        .rename("PROCESSING_TYPE")
        .reset_index()
    )
    out = frame.merge(per_sample, on="sample_id", how="left")
    matched = int(out["PROCESSING_TYPE"].notna().sum())
    print(f"PROCESSING_TYPE joined for {matched} of {len(out)} samples")
    if matched == 0:
        raise SystemExit(
            "PROCESSING_TYPE matched no sample. Every contamination estimate "
            "would be treated as uninterpretable, so stopping rather than "
            "reporting 'no reliable samples' over numbers that are fine."
        )
    return out


def main() -> int:
    set_global_seeds(DEFAULT_SEED)
    parser = argparse.ArgumentParser()
    parser.add_argument("--patients", nargs="*", default=None)
    parser.add_argument("--allow-dirty", action="store_true")
    # Re-summarise an existing table instead of re-reading the 9 GB matrix.
    # The per-sample estimates do not change; only the sorted/unsorted split
    # and the printed summary do.
    parser.add_argument("--summarise", metavar="PARQUET", default=None)
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

    if args.summarise:
        out = pd.read_parquet(args.summarise)
        out = out.drop(columns=["PROCESSING_TYPE"], errors="ignore")
        out = attach_processing_type(out, obs, metadata)
        _report(out)
        return 0

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

    out = attach_processing_type(out, obs, metadata)

    reliable = _report(out)

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
