#!/usr/bin/env python
"""Why does the maturity score tie? W1, answering the open half of issue #42.

**The finding this investigates.** `crypt_position` collapses to a two-bin split
on ~90% of patients: the 33rd and 67th percentiles of the maturity score land on
the same value, so `labels.py` falls back to `(crypt_bottom, crypt_top)` and the
rung stops being independent of `lineage`. Measured at full scale in decision
#14 — median Jaccard 1.000, identical in 37 of 60 patient-axis pairs.

**What is NOT known is why**, and that decides which fix is even possible:

- If the score ties because it is built from **too few markers**, a richer axis
  separates tertiles and the fix is an axis change — frozen, so PR plus two
  approvals, but a concrete proposal with a number behind it.
- If it ties because a large block of cells expresses **none of the markers**,
  no marker count of the same kind helps. Those cells are identical by
  construction, and the granularity curve honestly has three points.

The second is the live hypothesis, and it has a mechanism. `stem_pole` is
[LGR5, ASCL2, MKI67, OLFM4, SMOC2] — stem markers. A *mature* cell expresses
none of them, so every mature cell scores identically at the top of an inverted
axis. That is exactly what W2 found on Lee (issue #44): 32.0% of epithelial
cells with no stem marker, all scoring maximally mature. W1's depth thinning
removes the cells that are shallow, but a cell with adequate depth that
genuinely expresses no stem marker still ties with every other such cell.

**This measures which it is, and does not propose a fix.** Proposing one would
mean choosing an axis after seeing which choice separates, and the axes are
frozen precisely so that cannot happen quietly (CLAUDE.md invariant 3).

    qsub src/reference/jobs/measure_score_ties.sh
"""

from __future__ import annotations

import argparse
import itertools
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from src.common.io import write_versioned_table  # noqa: E402
from src.common.panel import axis_genes, tier_genes  # noqa: E402
from src.common.provenance import DEFAULT_SEED, set_global_seeds  # noqa: E402
from src.reference.ingest import (  # noqa: E402
    assign_compartments,
    read_gse178341,
    read_gse178341_clusters,
    read_gse178341_index,
    read_gse178341_metadata,
)
from src.reference.labels import TRANSCRIPT_AXES, score_markers  # noqa: E402
from src.reference.qc import apply_qc, cell_qc_metrics, qc_thresholds  # noqa: E402

UNSORTED = "unsorted"
DEPTH_QUANTILE = 0.25

#: The two cuts `crypt_position` needs. If they land on the same value the rung
#: degenerates to two bins — labels.py:378-383.
TERTILES = (1 / 3, 2 / 3)


def _tie_structure(score: np.ndarray) -> dict:
    """How much of the score is a single repeated value, and does it eat a cut."""
    values, counts = np.unique(score, return_counts=True)
    biggest = int(counts.max())
    mode_value = float(values[int(counts.argmax())])
    q33, q67 = (float(np.quantile(score, q)) for q in TERTILES)
    return {
        "n_cells": int(score.size),
        "n_distinct": int(values.size),
        "largest_tie": biggest,
        "largest_tie_frac": biggest / score.size,
        "tie_is_the_minimum": bool(mode_value == float(values.min())),
        "q33": q33,
        "q67": q67,
        # THE failure: both tertile cuts on one value, so the middle bin is empty.
        "tertiles_collide": bool(q33 == q67),
    }


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
        epithelial = keep & (compartment.to_numpy() == "epithelial")
        if epithelial.sum() < 50:
            print(f"[{i}/{len(patients)}] {patient} — too few epithelial cells")
            continue

        matrix = adata.X[epithelial]
        names = adata.var["gene_symbol"]
        totals = np.asarray(matrix.sum(axis=1)).ravel().astype(float)
        target = float(np.quantile(totals, DEPTH_QUANTILE))
        # Same population assign_labels scores: those clearing the depth target.
        scored = totals >= target
        if scored.sum() < 50:
            print(f"[{i}/{len(patients)}] {patient} — too few above depth target")
            continue
        sub = matrix[scored]

        for axis in TRANSCRIPT_AXES:
            markers = list(axis_genes(axis))
            present = [m for m in markers if m in set(names)]
            if len(present) < 2:
                continue

            full = score_markers(sub, names, present, context=f"{axis}/{patient}",
                                 target_genes=sorted(tier_genes("A")),
                                 depth_target=target, seed=DEFAULT_SEED)
            row = {"patient_id": patient, "labeling_axis": axis,
                   "n_markers": len(present), "subset_size": len(present)}
            row.update(_tie_structure(full))

            # HOW MANY CELLS EXPRESS NO MARKER AT ALL? If the tie block is these
            # cells, marker COUNT is not the lever — they are identical by
            # construction and always will be.
            positions = [list(names).index(m) for m in present]
            detected = np.asarray((sub[:, positions] > 0).sum(axis=1)).ravel()
            row["frac_zero_markers"] = float((detected == 0).mean())
            rows.append(row)

            # THE SWEEP. Score on every subset of size k, k = 1..K-1, and ask
            # whether the tertiles still collide. If collision falls steadily
            # with k, more markers is the lever. If it is flat, it is not.
            for k in range(1, len(present)):
                collide, ties = [], []
                for combo in itertools.combinations(present, k):
                    s = score_markers(
                        sub, names, list(combo), context=f"{axis}/{patient}/k{k}",
                        target_genes=sorted(tier_genes("A")),
                        depth_target=target, seed=DEFAULT_SEED,
                    )
                    st = _tie_structure(s)
                    collide.append(st["tertiles_collide"])
                    ties.append(st["largest_tie_frac"])
                rows.append({
                    "patient_id": patient, "labeling_axis": axis,
                    "n_markers": len(present), "subset_size": k,
                    "tertiles_collide": float(np.mean(collide)),
                    "largest_tie_frac": float(np.mean(ties)),
                    "frac_zero_markers": None, "n_cells": int(sub.shape[0]),
                })

        done = [r for r in rows if r["patient_id"] == patient
                and r["subset_size"] == r["n_markers"]]
        summary = ", ".join(
            f"{r['labeling_axis'][:4]} tie {r['largest_tie_frac']:.0%}"
            f"{' COLLIDE' if r['tertiles_collide'] else ''}" for r in done
        )
        print(f"[{i}/{len(patients)}] {patient} — {summary}")

    if not rows:
        raise SystemExit("no patient produced a score")
    out = pd.DataFrame(rows)
    # Mixed dtypes per column: the full-marker rows carry a bool (one
    # measurement) while the sweep rows carry a float (a rate over subsets), and
    # the sweep leaves the full-set-only fields None. pyarrow cannot type an
    # object column and the run died at the write with everything computed —
    # after printing the summary, which is the only reason the first run was not
    # a total loss. Cast explicitly rather than letting inference decide.
    for column in ("tertiles_collide", "tie_is_the_minimum", "frac_zero_markers",
                   "largest_tie_frac", "q33", "q67"):
        if column in out.columns:
            out[column] = pd.to_numeric(out[column], errors="coerce").astype(float)
    for column in ("n_distinct", "largest_tie"):
        if column in out.columns:
            out[column] = pd.to_numeric(out[column], errors="coerce").astype("Int64")

    full = out[out["subset_size"] == out["n_markers"]]

    print("\n" + "=" * 66)
    print("WHY DOES THE MATURITY SCORE TIE?")
    print("=" * 66)
    print(full.groupby("labeling_axis")[
        ["largest_tie_frac", "frac_zero_markers", "tertiles_collide",
         "tie_is_the_minimum"]
    ].mean().round(3).to_string())

    print("\n--- does adding markers help? (collision rate by subset size) ---")
    print(out.pivot_table(index="subset_size", columns="labeling_axis",
                          values="tertiles_collide", aggfunc="mean")
          .round(3).to_string())

    corr = full[["largest_tie_frac", "frac_zero_markers"]].corr().iloc[0, 1]
    print(f"\ncorrelation(largest tie block, cells with no marker detected) = {corr:.3f}")

    # TWO SEPARATE QUESTIONS, and an earlier version of this block conflated
    # them and printed the wrong conclusion.
    #
    #   WHAT is the tie block?      -> the correlation answers this
    #   Can it be made SMALLER?     -> only the sweep answers this
    #
    # A high correlation says the tied cells are the ones expressing no marker.
    # It says nothing about whether more markers would leave fewer such cells —
    # and the sweep showed collision falling monotonically with k, so it does.
    sweep = out.pivot_table(index="subset_size", columns="labeling_axis",
                            values="tertiles_collide", aggfunc="mean")
    if corr > 0.8:
        print(
            "\n  MECHANISM: the tie block IS the cells expressing no marker.\n"
            "  They score identically because they have nothing to score on."
        )
    else:
        print(
            "\n  MECHANISM: the tie block is NOT simply the undetected cells.\n"
            "  Something else is compressing the score; look at n_distinct."
        )
    for axis in sweep.columns:
        series = sweep[axis].dropna()
        if len(series) < 2:
            continue
        first, last = float(series.iloc[0]), float(series.iloc[-1])
        drop = first - last
        verdict = (
            "marker count IS a lever" if drop > 0.1 else
            "marker count is NOT a lever — collision barely moves"
        )
        plateau = (
            " and the curve has NOT plateaued, so more would help further"
            if len(series) > 2 and (float(series.iloc[-2]) - last) > 0.02 else
            " and the curve has flattened, so more of the same buys little"
        )
        print(f"  {axis}: {first:.2f} -> {last:.2f} across "
              f"{int(series.index[0])}..{int(series.index[-1])} markers — "
              f"{verdict}{plateau}.")
    print(
        "\n  NOT ANSWERED HERE: whether the remaining zero-marker cells are\n"
        "  biologically negative or merely detection-limited. A mature cell\n"
        "  expresses no stem marker however many you look for, so there is a\n"
        "  floor this measurement cannot locate. Proposing an axis change on an\n"
        "  extrapolation past that floor would be choosing markers to pass."
    )

    path = write_versioned_table(
        out, "score_tie_structure", seed=DEFAULT_SEED,
        allow_dirty=args.allow_dirty,
        notes=(
            "Why the maturity score ties, answering the open half of issue #42. "
            "Per patient per axis: tie structure of the full marker set, the "
            "fraction of cells expressing no marker at all, and a sweep over "
            "every marker subset of size k to test whether marker COUNT is the "
            "lever. Measures which fix is possible; deliberately proposes none, "
            "since choosing an axis after seeing which one separates is what the "
            "freeze exists to prevent."
        ),
        extra_meta={"n_patients": int(out["patient_id"].nunique())},
    )
    print(f"\nwrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
