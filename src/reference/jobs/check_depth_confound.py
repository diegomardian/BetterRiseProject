#!/usr/bin/env python
"""Is W1's compositional term separable from sequencing depth? W1.

**Why.** W2's first real-data run ([PR #45](https://github.com/diegomardian/BetterRiseProject/pull/45))
found that on Lee/SMC the maturity call is not separable from depth: 32% of
epithelial cells sampled zero stem markers, scored at the top of an inverted
axis, and were 4.7x shallower. The arms were not depth-matched either — normal
4,519 UMI against tumour 19,244 — which converts that sensitivity into a
**46-point apparent compositional loss in the hypothesised direction, out of
dropout alone**.

W2 notes W1 already handles half of this: `assign_labels` thins marker counts to
a common depth and marks the shallow ones `unresolved_depth`, *"not scored, not
counted as immature"*. The population W4 calls most mature is the one W1 refuses
to score.

**But their diagnostic has two conditions and thinning only answers one.**

1. ``maturity_tracks_depth`` — within an arm, do deeper cells get a different
   call? Thinning is aimed at this.
2. ``arms_are_depth_matched`` — were the two arms sequenced comparably? **Never
   checked on GSE178341**, and decision #14 already flagged at pilot scale that
   "the depth floor cuts the two arms unequally, and it can flip the sign of the
   compositional term."

Either alone is a caveat. Both together and the compositional term is not
separable from the depth imbalance, and nothing built on it should be quoted
without saying so.

This runs W2's diagnostic, unmodified, over W1's own labels — same QC, same
exclusions, same depth thinning as `run_full_reference.py`, so the answer is
about the labels actually shipped rather than a reconstruction of them.

    qsub src/reference/jobs/check_depth_confound.sh
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from src.common.io import write_versioned_table  # noqa: E402
from src.common.panel import granularity_rungs, tier_genes  # noqa: E402
from src.common.provenance import DEFAULT_SEED, set_global_seeds  # noqa: E402
from src.reference.ingest import (  # noqa: E402
    assign_compartments,
    read_gse178341,
    read_gse178341_clusters,
    read_gse178341_index,
    read_gse178341_metadata,
)
from src.reference.labels import (  # noqa: E402
    TRANSCRIPT_AXES,
    assign_labels,
    cell_type_vector,
)
from src.reference.qc import apply_qc, cell_qc_metrics, qc_thresholds  # noqa: E402

try:
    from src.harness.depth_confound import depth_confound_report
except ImportError as exc:  # pragma: no cover
    raise SystemExit(
        "src.harness.depth_confound is not available. It lands with W2's PR #45; "
        "until that merges, run this from a worktree of that branch:\n\n"
        "  git worktree add ../w2check origin/w2/real-data-depth-confound\n"
        "  cd ../w2check && python src/reference/jobs/check_depth_confound.py\n\n"
        f"({exc})"
    ) from exc

UNSORTED = "unsorted"
DEPTH_QUANTILE = 0.25
MATURE = "mature_colonocyte"


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

    rows = []
    for i, patient in enumerate(patients, 1):
        adata = read_gse178341(h5, patients=[patient])
        compartment = assign_compartments(clusters).reindex(adata.obs.index)
        metrics = cell_qc_metrics(adata.X, adata.var["gene_symbol"],
                                  batch=adata.obs["sample_id"])
        keep = apply_qc(metrics, qc_thresholds(metrics)).to_numpy()
        joined = adata.obs.join(metadata, how="left")
        keep &= (joined["PROCESSING_TYPE"] == UNSORTED).to_numpy()
        if keep.sum() < 50:
            print(f"[{i}/{len(patients)}] {patient} — too few usable cells")
            continue

        tissue = adata.obs["tissue"].to_numpy()[keep]
        if len(set(tissue.tolist())) < 2:
            print(f"[{i}/{len(patients)}] {patient} — one arm only, skipped")
            continue

        labels = assign_labels(
            adata.X[keep], adata.var["gene_symbol"],
            compartment=compartment.to_numpy()[keep],
            sample_id=adata.obs["sample_id"].to_numpy()[keep],
            target_genes=sorted(tier_genes("A")),
            tissue=tissue,
            patient_id=adata.obs["patient_id"].to_numpy()[keep],
            depth_quantile=DEPTH_QUANTILE, seed=DEFAULT_SEED,
            index=adata.obs.index[keep],
        )

        # Depth is total UMIs per cell, on the same cells the labels describe.
        depth = np.asarray(adata.X[keep].sum(axis=1)).ravel()

        for axis in TRANSCRIPT_AXES:
            for rung in granularity_rungs():
                # cell_type_vector, not a hand-mapping: which bin counts as
                # mature is recorded in RUNG_SPECS and duplicating that mapping
                # is how the two drift apart (labels.py says so).
                is_mature = cell_type_vector(labels, axis, rung) == MATURE
                report = depth_confound_report(depth, is_mature, tissue)
                rows.append({
                    "patient_id": patient, "labeling_axis": axis,
                    "granularity_rung": rung,
                    "confounded": bool(report["confounded"]),
                    "maturity_tracks_depth": bool(report["maturity_tracks_depth"]),
                    "arms_are_depth_matched": bool(report["arms_are_depth_matched"]),
                    # Key names are depth_ratio_between_arms and
                    # worst_within_arm_rho — read them off the module, do not
                    # guess. A first version guessed and died on patient one.
                    "depth_ratio": float(report["depth_ratio_between_arms"]),
                    "worst_rho": float(report["worst_within_arm_rho"]),
                    "n_cells": int(keep.sum()),
                })
        mine = [r for r in rows if r["patient_id"] == patient]
        finite = [abs(r["worst_rho"]) for r in mine if np.isfinite(r["worst_rho"])]
        n_nan = len(mine) - len(finite)
        worst = f"{max(finite):.2f}" if finite else "n/a"
        print(f"[{i}/{len(patients)}] {patient} — depth ratio "
              f"{mine[-1]['depth_ratio']:.2f}, worst |rho| {worst}"
              + (f", {n_nan}/{len(mine)} rungs UNCOMPUTABLE" if n_nan else ""))

    if not rows:
        raise SystemExit("no patient produced a paired diagnostic")
    out = pd.DataFrame(rows)

    print("\n" + "=" * 64)
    print("IS THE COMPOSITIONAL TERM SEPARABLE FROM DEPTH?")
    print("=" * 64)
    print(out.groupby(["labeling_axis", "granularity_rung"])[
        ["confounded", "maturity_tracks_depth", "arms_are_depth_matched"]
    ].mean().round(3).to_string())
    print(f"\ndepth ratio (tumour:normal): median "
          f"{out['depth_ratio'].median():.2f}, max {out['depth_ratio'].max():.2f}")
    n_conf = int(out["confounded"].sum())
    if n_conf:
        print(f"\n!! {n_conf} of {len(out)} patient-axis-rung combinations are "
              f"CONFOUNDED.\n   Both conditions hold: the call tracks depth AND "
              f"the arms are not\n   depth-matched. A decomposition from those "
              f"must not be quoted without\n   saying so — W2's PR #45, and "
              f"decision #14 predicted it at pilot scale.")
    else:
        print("\n   No combination is confounded on both conditions. W1's depth "
              "thinning\n   is doing what decision #14 hoped, and the arms are "
              "comparably sequenced.")

    path = write_versioned_table(
        out, "depth_confound_reference", seed=DEFAULT_SEED,
        allow_dirty=args.allow_dirty,
        notes=(
            "W2's depth_confound_report (PR #45) run unmodified over W1's own "
            "labels on GSE178341, same QC and depth thinning as "
            "run_full_reference. Answers the half of W2's diagnostic W1 had not "
            "checked: whether the two arms are depth-matched. Decision #14 "
            "flagged at pilot scale that the depth floor cuts the arms unequally."
        ),
        extra_meta={
            "n_patients": int(out["patient_id"].nunique()),
            "n_confounded": n_conf,
            "median_depth_ratio": float(out["depth_ratio"].median()),
        },
    )
    print(f"\nwrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
