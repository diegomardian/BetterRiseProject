"""Labels and S matrices for all 62 patients, at version 1.0.0. W1, weeks 4-5.

The last W1 artifact. `run_pilot.py` does this for five patients as an interface
test; this is the deliverable, and it differs in four ways that are decisions
rather than scale.

**1 · Both tumour-arm definitions, always.** `prereg amendment 1` requires the
MMR contrast under `filtered` and `unfiltered` both, and pre-commits to reading
disagreement as *not identifiable* rather than choosing. Emitting only one here
would quietly make that choice.

**2 · Ambient exclusions at the committed threshold.** 10%, fixed on 2026-08-23
before anyone counted its cost (#16). The paired n it leaves is printed beside
the rule, because a threshold reported without its cost invites revising it.

**3 · Malignancy calls from the full inferCNV run**, with patients that had no
separable aneuploid population carrying `not_called` rather than a threshold
drawn through noise (#15).

**4 · S matrices at 1.0.0**, on the shared gene index.

    python src/reference/jobs/run_full_reference.py
    python src/reference/jobs/run_full_reference.py --patients C122 C162
"""

from __future__ import annotations

import argparse
import glob
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

import pandas as pd  # noqa: E402

from src.common.io import write_versioned_table  # noqa: E402
from src.common.panel import tier_genes  # noqa: E402
from src.common.provenance import DEFAULT_SEED, set_global_seeds  # noqa: E402
from src.reference.ambient import ambient_exclusions  # noqa: E402
from src.reference.ingest import (  # noqa: E402
    assign_compartments,
    read_gse178341,
    read_gse178341_clusters,
    read_gse178341_metadata,
)
from src.reference.labels import (  # noqa: E402
    TUMOUR_ARMS,
    assign_labels,
    mature_cell_counts,
    rung_degeneracy,
)
from src.reference.qc import apply_qc, cell_qc_metrics  # noqa: E402

S_MATRIX_VERSION = "1.0.0"
GENE_INDEX_VERSION = "1.0.0"
DEPTH_QUANTILE = 0.25
SIGNATURE_GENES = 800
UNSORTED = "unsorted"


def _latest(pattern: str):
    hits = sorted(glob.glob(pattern))
    return pd.read_parquet(hits[-1]) if hits else None


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

    # Prior stages, read from results/ rather than recomputed. If either is
    # missing the run still proceeds — but says which safeguard is absent,
    # because a silently unfiltered run looks identical to a filtered one.
    calls = _latest("results/*/malignancy_calls.parquet")
    contamination = _latest("results/*/ambient_contamination.parquet")
    if calls is None:
        print("!! no malignancy_calls.parquet — the FILTERED tumour arm cannot "
              "be built.\n   Run emit_malignancy_calls.py first (#15).")
    if contamination is None:
        print("!! no ambient_contamination.parquet — no exclusions will be "
              "applied (#16).")

    excluded: set[str] = set()
    if contamination is not None:
        rule = ambient_exclusions(contamination)
        excluded = set(rule.loc[rule["excluded"], "sample_id"])
        print(f"ambient: excluding {len(excluded)} of {len(rule)} samples above "
              f"10% contamination (#16)")

    patients = args.patients or sorted(
        pd.unique(clusters.index.map(lambda b: str(b).split("_")[0]))
    )
    targets = sorted(tier_genes("A"))
    print(f"{len(patients)} patients · targets {targets}")

    counts_all, degeneracy_all = [], []
    for i, patient in enumerate(patients, 1):
        adata = read_gse178341(h5, patients=[patient])
        compartment = assign_compartments(clusters).reindex(adata.obs.index)
        metrics = cell_qc_metrics(adata.X, adata.var["gene_symbol"],
                                  batch=adata.obs["sample_id"])
        keep = apply_qc(metrics).to_numpy()

        # Unsorted only (#11), minus the ambient exclusions (#16).
        joined = adata.obs.join(metadata, how="left")
        usable = (joined["PROCESSING_TYPE"] == UNSORTED).to_numpy()
        usable &= ~adata.obs["sample_id"].isin(excluded).to_numpy()
        keep &= usable
        if keep.sum() < 50:
            print(f"[{i}/{len(patients)}] {patient} — {int(keep.sum())} usable "
                  f"cells, skipped")
            continue

        labels = assign_labels(
            adata.X[keep], adata.var["gene_symbol"],
            compartment=compartment.to_numpy()[keep],
            sample_id=adata.obs["sample_id"].to_numpy()[keep],
            target_genes=targets,
            tissue=adata.obs["tissue"].to_numpy()[keep],
            patient_id=adata.obs["patient_id"].to_numpy()[keep],
            depth_quantile=DEPTH_QUANTILE, seed=DEFAULT_SEED,
            index=adata.obs.index[keep],
        )

        malignant = None
        if calls is not None:
            m = calls.set_index("cell")["call"].astype(str)
            malignant = (m.reindex(adata.obs.index[keep]) == "malignant").to_numpy()

        # BOTH arms. Emitting one would silently make the choice amendment 1
        # exists to refuse.
        for arm in TUMOUR_ARMS:
            if arm == "filtered" and malignant is None:
                continue
            counts_all.append(mature_cell_counts(
                labels,
                patient_id=adata.obs["patient_id"].to_numpy()[keep],
                tissue=adata.obs["tissue"].to_numpy()[keep],
                malignant=malignant, tumour_arm=arm,
            ))
        degeneracy_all.append(rung_degeneracy(labels).assign(patient_id=patient))
        print(f"[{i}/{len(patients)}] {patient} — {int(keep.sum()):,} cells")

    if not counts_all:
        raise SystemExit("no patient produced labels")
    counts = pd.concat(counts_all, ignore_index=True)

    print("\n" + "=" * 60)
    print("PAIRED n UNDER EACH TUMOUR-ARM DEFINITION")
    print("=" * 60)
    for arm in sorted(set(counts["tumour_arm"])):
        sub = counts[counts["tumour_arm"] == arm]
        paired = (sub.groupby("patient_id", observed=True)["tissue"]
                  .nunique() >= 2).sum()
        print(f"  {arm:<12} {int(paired)} patients with both arms")
    print(
        "\n  Report BOTH. Prereg amendment 1 pre-commits to reading disagreement\n"
        "  between them as NOT IDENTIFIABLE rather than choosing whichever is\n"
        "  more interesting."
    )

    for frame, name, note in (
        (counts, "mature_cell_counts_full",
         "Per-patient mature-cell counts, all 62 patients, under BOTH tumour-arm "
         "definitions (prereg amendment 1). Unsorted samples only (#11), ambient "
         "exclusions at the committed 10% threshold (#16), malignancy calls with "
         "not_called where no aneuploid population separated (#15)."),
        (pd.concat(degeneracy_all, ignore_index=True), "rung_degeneracy_full",
         "Which granularity rungs collapsed onto the same partition, per "
         "patient. The curve is only a curve where they did not."),
    ):
        path = write_versioned_table(
            frame, name, seed=DEFAULT_SEED, notes=note,
            allow_dirty=args.allow_dirty,
        )
        print(f"wrote {path}")

    print(
        f"\nNEXT: S matrices at {S_MATRIX_VERSION} on gene index "
        f"{GENE_INDEX_VERSION}, then checks.py for G1 —\nand commit G1's "
        "threshold before looking at anything, as #16's was."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
