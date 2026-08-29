#!/usr/bin/env python
"""The per-gene summary the estimator decomposes. W1, §6.1.

`src.estimator.kitagawa.decompose_cohort` needs ten columns per
(patient, study, gene, rung, axis). `maturity_summary()` supplies seven; the
three it cannot are **gene**, **mean_normal** and **mean_tumour**, because they
need the target gene's expression and `maturity_summary` never sees a matrix.
This job supplies them.

WHAT THE MEANS ARE
------------------
`decompose()` splits ``total = compositional + intrinsic + interaction`` as::

    compositional = (frac_tumour - frac_normal) * mean_normal
    intrinsic     = frac_normal * (mean_tumour - mean_normal)

so ``mean_*`` is the gene's mean expression **within mature cells** of that arm.
The compositional term is "the mature cells left"; the intrinsic term is "the
mature cells that stayed turned it down".

CP10K, NOT RAW COUNTS
---------------------
`src/harness/pseudobulk.py` takes raw-count means, which is right for synthetic
cells sequenced to a common depth. **Real arms are not.** W1 measured tumour and
normal 1.64x apart in median depth, 20 of 32 patients above 1.5x
(`results/2026-08-27_e5ebdc3/`). A raw-count `mean_tumour - mean_normal` is then
partly a library-size difference, and it lands entirely in the **intrinsic**
term — the one this project exists to measure.

Counts are scaled to 10,000 per cell and **not** logged. Logging would break the
decomposition: ``total = frac * mean`` is an identity on a linear scale and not
on a log one.

WHICH GENES
-----------
Panel tiers A, B, C and D. Tier E is excluded because MUC2 and TFF3 are
`opposite_lineage` markers, and open decision #1's narrow reading scopes the
target set to the run in question — the same call #49 made. Verified at call
time: no target gene appears in either labelling axis.

    qsub src/reference/jobs/build_decomposition_summary.sh
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import sparse  # noqa: E402

from src.common.io import write_versioned_table  # noqa: E402
from src.common.panel import axis_genes, granularity_rungs, tier_genes  # noqa: E402
from src.common.provenance import DEFAULT_SEED, set_global_seeds  # noqa: E402
from src.harness.depth_confound import match_arm_depth  # noqa: E402
from src.reference.ingest import (  # noqa: E402
    assign_compartments,
    read_gse178341,
    read_gse178341_clusters,
    read_gse178341_index,
    read_gse178341_metadata,
)
from src.reference.labels import (  # noqa: E402
    NON_EPITHELIAL,
    TRANSCRIPT_AXES,
    UNRESOLVED,
    _positions,
    assign_labels,
    cell_type_vector,
    label_column,
)
from src.reference.qc import apply_qc, cell_qc_metrics, qc_thresholds  # noqa: E402

STUDY_ID = "GSE178341"
UNSORTED = "unsorted"
DEPTH_QUANTILE = 0.25
MATURE = "mature_colonocyte"
CP10K = 1e4

#: Tiers A-D. Not E — see the module docstring.
TARGET_TIERS = ("A", "B", "C", "D")


def target_genes() -> list[str]:
    genes = sorted({g for tier in TARGET_TIERS for g in tier_genes(tier)})
    # Invariant 2, checked here rather than assumed: a gene that is both a
    # target and a label marker makes a silenced cell unreadable from an absent
    # one at the rung it is the control for.
    axis_markers = {g for a in TRANSCRIPT_AXES for g in axis_genes(a)}
    clash = sorted(set(genes) & axis_markers)
    if clash:
        raise SystemExit(
            f"target gene(s) {clash} are also labelling markers. Invariant 2: a "
            f"target may not define the labels it is measured against."
        )
    return genes


def main() -> int:
    set_global_seeds(DEFAULT_SEED)
    parser = argparse.ArgumentParser()
    parser.add_argument("--patients", nargs="*", default=None)
    parser.add_argument("--allow-dirty", action="store_true")
    # DECISION #24.1: match_arm_depth is the PRIMARY read for lineage and
    # crypt_position. A depth floor removes the mechanism by which depth reaches
    # the maturity call; it does not make the arms comparable, and W1 measured
    # them 1.64x apart in median depth after the floor. This equalises the
    # distributions by construction — no threshold to choose.
    #
    # A separate flag rather than the default, because the unmatched read is
    # what results/2026-08-28_8965a6f was built from and the two must stay
    # comparable. The name says which one a table came from.
    parser.add_argument("--match-depth", action="store_true")
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
    genes = target_genes()
    print(f"{len(patients)} patients · {len(genes)} target genes: {genes}")

    rows, skipped = [], []
    for i, patient in enumerate(patients, 1):
        adata = read_gse178341(h5, patients=[patient])
        compartment = assign_compartments(clusters).reindex(adata.obs.index)
        metrics = cell_qc_metrics(adata.X, adata.var["gene_symbol"],
                                  batch=adata.obs["sample_id"])
        keep = apply_qc(metrics, qc_thresholds(metrics)).to_numpy()
        joined = adata.obs.join(metadata, how="left")
        keep &= (joined["PROCESSING_TYPE"] == UNSORTED).to_numpy()
        if keep.sum() < 50:
            skipped.append((patient, "fewer than 50 cells survived QC"))
            print(f"[{i}/{len(patients)}] {patient} — too few cells")
            continue

        tissue = adata.obs["tissue"].to_numpy()[keep]
        comp = compartment.to_numpy()[keep]
        if not ((comp == "epithelial") & (tissue == "normal")).any():
            skipped.append((patient, "no normal epithelium (#9/#11/#16)"))
            print(f"[{i}/{len(patients)}] {patient} — no reference arm")
            continue

        names = adata.var["gene_symbol"]
        positions, found = _positions(names, genes)
        if len(found) < len(genes):
            skipped.append((patient, f"missing target genes {sorted(set(genes)-set(found))}"))
            print(f"[{i}/{len(patients)}] {patient} — missing targets")
            continue

        block = adata.X[keep]
        labels = assign_labels(
            block, names, compartment=comp,
            sample_id=adata.obs["sample_id"].to_numpy()[keep],
            target_genes=sorted(tier_genes("A")), tissue=tissue,
            patient_id=adata.obs["patient_id"].to_numpy()[keep],
            depth_quantile=DEPTH_QUANTILE, seed=DEFAULT_SEED,
            index=adata.obs.index[keep],
        )

        # CP10K on the target columns only. Scaling is per cell, so slicing
        # first is exact — the denominator is the cell's FULL library size.
        totals = np.asarray(block.sum(axis=1), dtype=float).ravel()
        totals[totals == 0] = np.nan
        sub = block[:, positions]
        sub = sub.toarray() if sparse.issparse(sub) else np.asarray(sub)
        expr = sub / totals[:, None] * CP10K

        is_normal, is_tumour = tissue == "normal", tissue == "tumour"

        # Matched WITHIN patient, and after labelling rather than before: the
        # labels are a property of the cell, the matching is a property of the
        # comparison. Subsampling before labelling would move the per-sample cut
        # points too, which is a different intervention.
        if args.match_depth:
            depth_all = np.asarray(block.sum(axis=1), dtype=float).ravel()
            matched = match_arm_depth(depth_all, tissue, seed=DEFAULT_SEED)
        else:
            matched = np.ones(len(tissue), dtype=bool)

        for axis in TRANSCRIPT_AXES:
            for rung in granularity_rungs():
                call = cell_type_vector(labels, axis, rung)
                raw = labels[label_column(axis, rung)].astype(str).to_numpy()
                scored = (raw != UNRESOLVED) & (raw != NON_EPITHELIAL) & matched
                mature = (call == MATURE) & scored

                # Fractions over the RESOLVED epithelium, matching
                # mature_cell_counts: a cell that could not be measured is not a
                # cell measured to be immature (#14).
                den_n, den_t = (scored & is_normal).sum(), (scored & is_tumour).sum()
                if den_n == 0 or den_t == 0:
                    continue
                mat_n, mat_t = mature & is_normal, mature & is_tumour
                frac_n, frac_t = mat_n.sum() / den_n, mat_t.sum() / den_t

                for j, gene in enumerate(found):
                    # None, not 0.0, when an arm has no mature cell to average
                    # over. An absent mean is not a mean of zero, and writing
                    # zero here would put a fabricated number straight into the
                    # intrinsic term (invariant 1).
                    mn = float(np.nanmean(expr[mat_n, j])) if mat_n.any() else None
                    mt = float(np.nanmean(expr[mat_t, j])) if mat_t.any() else None
                    rows.append({
                        "patient_id": patient, "study_id": STUDY_ID, "gene": gene,
                        "granularity_rung": rung, "labeling_axis": axis,
                        "frac_mature_normal": float(frac_n),
                        "frac_mature_tumour": float(frac_t),
                        "mean_normal": mn, "mean_tumour": mt,
                        # Estimability reads the TUMOUR arm: it is the depleted
                        # one, so it is the binding constraint on whether an
                        # intrinsic term can be asked for at all.
                        "n_cells_mature": int(mat_t.sum()),
                        "n_mature_normal": int(mat_n.sum()),
                        "n_scored_normal": int(den_n), "n_scored_tumour": int(den_t),
                        # #24.1: n after matching must travel with any matched
                        # number, because the cost is the point.
                        "depth_matched": bool(args.match_depth),
                        "n_cells_before_matching": int(len(tissue)),
                        "n_cells_after_matching": int(matched.sum()),
                    })
        print(f"[{i}/{len(patients)}] {patient} — "
              f"{len([r for r in rows if r['patient_id']==patient])} rows")

    if not rows:
        raise SystemExit("no patient produced a summary")
    out = pd.DataFrame(rows)

    print("\n" + "=" * 64)
    print("SUMMARY FOR decompose_cohort")
    print("=" * 64)
    print(f"  rows                {len(out):,}")
    print(f"  patients            {out.patient_id.nunique()}")
    print(f"  genes x rungs x axes {out.gene.nunique()} x "
          f"{out.granularity_rung.nunique()} x {out.labeling_axis.nunique()}")
    print(f"  mean_normal absent  {int(out.mean_normal.isna().sum()):,}")
    print(f"  mean_tumour absent  {int(out.mean_tumour.isna().sum()):,}")
    print("\n  n_cells_mature (tumour arm), by rung:")
    print(out.groupby("granularity_rung")["n_cells_mature"]
          .describe()[["50%", "max"]].round(1).to_string())
    if skipped:
        print(f"\n  {len(skipped)} patient(s) skipped:")
        for p, why in skipped[:8]:
            print(f"    {p}  {why}")

    print("\n  NOT a decomposition. This is the input decompose_cohort consumes;")
    print("  it applies no cutpoint and nulls no intrinsic term.")

    path = write_versioned_table(
        out,
        "decomposition_summary_matched" if args.match_depth
        else "decomposition_summary",
        seed=DEFAULT_SEED,
        allow_dirty=args.allow_dirty,
        notes=(
            "Per-gene summary for src.estimator.kitagawa.decompose_cohort, "
            "GSE178341. mean_normal/mean_tumour are CP10K expression means "
            "WITHIN mature cells of each arm — not raw counts, because W1 "
            "measured the arms 1.64x apart in median depth and a raw-count "
            "difference lands in the intrinsic term. Not logged: total = frac * "
            "mean is an identity on a linear scale only. Fractions are over "
            "resolved epithelium (#14). n_cells_mature is the tumour arm's "
            "count, the binding constraint on estimability. Tiers A-D; E "
            "excluded because MUC2/TFF3 are opposite_lineage markers (#1)."
            + (
                " DEPTH-MATCHED (#24.1): cells subsampled per patient so both "
                "arms share a depth distribution, applied AFTER labelling since "
                "the labels are a property of the cell and the matching is a "
                "property of the comparison. n before and after matching travel "
                "with every row."
                if args.match_depth else
                " UNMATCHED — the depth floor only. #24.1 makes the matched read "
                "primary for lineage and crypt_position; run with --match-depth."
            )
        ),
        extra_meta={
            "n_patients": int(out.patient_id.nunique()),
            "n_rows": int(len(out)),
            "genes": genes,
            "depth_quantile": DEPTH_QUANTILE,
            "normalisation": "CP10K, not logged",
            "depth_matched": bool(args.match_depth),
            "cells_retained_by_matching": (
                float(out.n_cells_after_matching.sum() / out.n_cells_before_matching.sum())
                if args.match_depth else 1.0
            ),
        },
    )
    print(f"\nwrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
