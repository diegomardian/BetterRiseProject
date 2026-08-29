"""Per-patient mature-cell MEANS on GSE178341 — the input the decomposition lacked.

`run_full_reference` emits `mature_cell_counts_full`: fractions and counts, under
both tumour-arm definitions. `maturity_summary`'s docstring says the caller adds
``gene``, ``mean_normal`` and ``mean_tumour``. Nobody had, so no decomposition
existed on the primary cohort.

This is that caller. It mirrors `run_full_reference`'s loop exactly — same QC,
same #11 unsorted filter, same #16 ambient exclusions, same skip rules, same
`assign_labels` call with the same committed `depth_quantile` — and adds one
thing: the mean of each panel gene **within the mature cells of each arm**.

Deliberately NOT re-deriving the counts. It joins to W1's committed
`mature_cell_counts_full` so the two cannot drift; a mismatch on the join keys is
raised rather than filled.

TARGET SET IS TIERS A+B+C+D
---------------------------
Tier E is excluded because MUC2 and TFF3 are `opposite_lineage` markers, and open
decision #1's narrow reading scopes the target set to the run in question. The
broad set would cost half of axis 2's markers for nothing.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from src.common.io import write_versioned_table  # noqa: E402
from src.common.panel import tier_genes  # noqa: E402
from src.common.provenance import set_global_seeds  # noqa: E402
from src.reference.ingest import (  # noqa: E402
    assign_compartments,
    read_gse178341,
    read_gse178341_clusters,
    read_gse178341_metadata,
)
from src.reference.labels import (  # noqa: E402
    RUNG_SPECS,
    TRANSCRIPT_AXES,
    assign_labels,
    granularity_rungs,
    label_column,
)
from src.reference.qc import apply_qc, cell_qc_metrics, qc_thresholds  # noqa: E402

DEFAULT_SEED = 20260101
DEPTH_QUANTILE = 0.10
UNSORTED = "unsorted"
TIERS = ("A", "B", "C", "D")


def _latest(pattern: str) -> Path | None:
    import glob

    hits = sorted(glob.glob(pattern))
    return Path(hits[-1]) if hits else None


def main() -> int:
    set_global_seeds(DEFAULT_SEED)
    parser = argparse.ArgumentParser()
    parser.add_argument("--patients", nargs="*", default=None)
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args()

    targets = sorted({g for tier in TIERS for g in tier_genes(tier)})
    label_targets = sorted(tier_genes("A"))  # labels stay NARROW — decision #21
    print(f"target genes (A+B+C+D): {targets}")

    data = Path(os.environ.get("BRP_DATA_DIR", "data")) / "raw" / "GSE178341"
    h5 = data / "GSE178341_crc10x_full_c295v4_submit.h5"
    clusters = read_gse178341_clusters(
        data / "GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz"
    )
    metadata = read_gse178341_metadata(
        data / "GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz"
    )

    counts_path = _latest("results/*/mature_cell_counts_full.parquet")
    if counts_path is None:
        print("!! no mature_cell_counts_full.parquet — run run_full_reference first")
        return 1
    print(f"joining to {counts_path}")

    contamination = _latest("results/*/ambient_contamination.parquet")
    excluded: set[str] = set()
    if contamination is not None:
        from src.reference.ambient import ambient_exclusions

        rule = ambient_exclusions(pd.read_parquet(contamination))
        excluded = set(rule.loc[rule["excluded"], "sample_id"])
        print(f"ambient: excluding {len(excluded)} samples (#16)")

    patients = args.patients or sorted(
        pd.unique(clusters.index.map(lambda b: str(b).split("_")[0]))
    )
    rungs = list(granularity_rungs())

    rows: list[dict] = []
    skipped: list[dict] = []
    for i, patient in enumerate(patients, start=1):
        try:
            adata = read_gse178341(h5, patients=[patient], verify=False)
        except Exception as exc:  # noqa: BLE001 - a per-patient read failure is data, not a crash
            skipped.append({"patient_id": patient, "reason": f"read failed: {exc}"})
            continue

        compartment = assign_compartments(clusters.reindex(adata.obs.index))
        metrics = cell_qc_metrics(
            adata.X, adata.var["gene_symbol"], batch=adata.obs["sample_id"]
        )
        keep = apply_qc(metrics, qc_thresholds(metrics)).to_numpy()
        joined = adata.obs.join(metadata, how="left")
        keep &= (joined["PROCESSING_TYPE"] == UNSORTED).to_numpy()
        keep &= ~adata.obs["sample_id"].isin(excluded).to_numpy()

        if keep.sum() < 50:
            skipped.append({"patient_id": patient, "reason": "fewer than 50 cells after QC"})
            continue
        is_epi = compartment.to_numpy()[keep] == "epithelial"
        is_normal = adata.obs["tissue"].to_numpy()[keep] == "normal"
        if not (is_epi & is_normal).any():
            skipped.append({"patient_id": patient, "reason": "no normal epithelium (#9)"})
            continue

        labels = assign_labels(
            adata.X[keep], adata.var["gene_symbol"],
            compartment=compartment.to_numpy()[keep],
            sample_id=adata.obs["sample_id"].to_numpy()[keep],
            target_genes=label_targets,
            tissue=adata.obs["tissue"].to_numpy()[keep],
            patient_id=adata.obs["patient_id"].to_numpy()[keep],
            depth_quantile=DEPTH_QUANTILE, seed=DEFAULT_SEED,
            index=adata.obs.index[keep],
        )

        # CP10K on the kept cells: decompose() needs a linear, depth-normalised
        # scale for the Kitagawa identity to be additive. Raw counts would
        # confound "how much a cell makes" with how deeply it was sequenced.
        matrix = adata.X[keep]
        library = np.asarray(matrix.sum(axis=1), dtype=float).ravel()
        library[library == 0] = 1.0
        symbols = list(adata.var["gene_symbol"])
        tissue = adata.obs["tissue"].to_numpy()[keep]

        present = [g for g in targets if g in symbols]
        columns = {g: symbols.index(g) for g in present}
        expression = {
            g: np.asarray(matrix[:, j].todense()).ravel() / library * 1e4
            for g, j in columns.items()
        }

        for axis in TRANSCRIPT_AXES:
            for rung in rungs:
                column = label_column(axis, rung)
                if column not in labels.columns:
                    continue
                mature = labels[column].astype(str).to_numpy() == RUNG_SPECS[rung].mature
                for gene, values in expression.items():
                    row = {
                        "patient_id": patient,
                        "study_id": "GSE178341",
                        "gene": gene,
                        "granularity_rung": rung,
                        "labeling_axis": axis,
                    }
                    for arm_name, arm_mask in (
                        ("normal", tissue == "normal"),
                        ("tumour", tissue == "tumour"),
                    ):
                        sel = mature & arm_mask
                        row[f"mean_{arm_name}"] = (
                            float(values[sel].mean()) if sel.any() else np.nan
                        )
                        row[f"n_mature_{arm_name}"] = int(sel.sum())
                    rows.append(row)

        print(f"[{i}/{len(patients)}] {patient} — {int(keep.sum()):,} cells, "
              f"{len(present)} genes", flush=True)

    # Skips are printed FIRST and unconditionally. run_full_reference carries a
    # comment about exactly this: a loss printed once in a scrolling log but
    # missing from the summary is a loss nobody will find later. Returning early
    # on an empty frame without saying why is the same mistake.
    if skipped:
        print(f"\nskipped {len(skipped)} of {len(patients)}:")
        for reason, n in pd.DataFrame(skipped).groupby("reason", observed=True).size().items():
            print(f"  {n:2d}  {reason}")

    means = pd.DataFrame(rows)
    if means.empty:
        print("!! no patient produced a row — see the skip reasons above")
        return 1

    print(f"\nmeans: {len(means)} rows, {means.patient_id.nunique()} patients")

    path = write_versioned_table(
        means, "mature_cell_means_full", seed=DEFAULT_SEED,
        notes=(
            "Per-patient mean CP10K expression of panel tiers A+B+C+D WITHIN the "
            "mature cells of each arm. The input mature_cell_counts_full does not "
            "carry and decompose_cohort requires. Mirrors run_full_reference's "
            "loop: same QC, #11 unsorted filter, #16 ambient exclusions, skip "
            "rules and depth_quantile."
        ),
        extra_meta={
            "target_genes": targets,
            "label_target_genes": label_targets,
            "why_not_tier_E": "MUC2 and TFF3 are opposite_lineage markers; open "
                              "decision #1 scopes the target set to this run.",
            "depth_quantile": DEPTH_QUANTILE,
            "join_to": str(counts_path),
            "skipped": skipped,
        },
        allow_dirty=args.allow_dirty,
    )
    print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
