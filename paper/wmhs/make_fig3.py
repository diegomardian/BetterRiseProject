#!/usr/bin/env python
"""Appendix Figure 2 — the recovery curve, and why it carries no estimator.

Every number on the axes comes from
the newest ``results/<date>_<sha>/calibration_gap_recovery.parquet`` (the 50-
replicate run). Nothing is transcribed. Run from the repo root:

    python paper/wmhs/make_fig3.py

The caption's central claim — ``max|i_hat - i_realised| = 0`` — is read out of
the table's own ``max_abs_residual_vs_realised`` column and printed, so the
figure and the sentence beside it cannot drift apart.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
from _tables import newest

matplotlib.use("Agg")
# Drawn at the width it is printed at (\linewidth = 5.5in), so point
# sizes below are the sizes that reach the page. Axis labels take theirs
# from rcParams, so pin them to match the tick labels.
matplotlib.rcParams.update({"axes.labelsize": 7.5})
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

TABLE_NAME = "calibration_gap_recovery"
OUT = Path("paper/wmhs/figures/fig3_generator_statistic.pdf")

DETECTABLE_SHIFT = 0.5
COMMITTED_OK, COMMITTED_WIDE = 50, 20
SHADE = "#ececec"

POOL_STYLE = {
    "pooled": ("#1f77b4", "draw pool as implemented"),
    "reference": ("#d62728", "a second reading of the same choice"),
}


def main() -> int:
    rec = pd.read_parquet(newest(TABLE_NAME))

    residual = float(rec["max_abs_residual_vs_realised"].max())
    n_rows = int(rec["n_replicates"].sum())
    n_undefined = int(rec["n_ratio_undefined"].sum())
    print(f"oracle sweep rows           : {n_rows:,}")
    print(f"max |i_hat - i_realised|    : {residual:.3g}")
    print(f"rows with an undefined ratio: {n_undefined:,} (shift = 1.0, the null arm)")

    at_effect = rec[(rec["shift"] == DETECTABLE_SHIFT) & (rec["grid"] == "extended")]

    fig, ax = plt.subplots(figsize=(5.5, 2.0))
    ax.axvspan(1, COMMITTED_WIDE, color=SHADE, zorder=0)
    ax.axvspan(COMMITTED_WIDE, COMMITTED_OK, color=SHADE, alpha=0.55, zorder=0)
    for edge in (COMMITTED_WIDE, COMMITTED_OK):
        ax.axvline(edge, color="#999999", lw=0.8, zorder=1)
    ax.axhline(1.0, ls="--", lw=0.9, color="#444444", zorder=2)

    for pool, (colour, label) in POOL_STYLE.items():
        rows = (
            at_effect[at_effect["pool"] == pool]
            .groupby("median_n_cells_mature")
            .agg(
                median=("ratio_median", "median"),
                lo=("ratio_q25", "median"),
                hi=("ratio_q75", "median"),
            )
            .sort_index()
        )
        rows = rows[rows.index > 0]
        ax.fill_between(rows.index, rows.lo, rows.hi, color=colour, alpha=0.15, lw=0)
        ax.plot(rows.index, rows["median"], "-o", ms=3.5, lw=1.4,
                color=colour, label=label)

    ax.set_xscale("log")
    ax.set_xlim(1, 1000)
    ax.set_ylim(0.4, 2.1)
    ax.set_xticks([1, 10, 20, 50, 90, 200, 500])
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_xlabel("cells available to the intrinsic comparison, $n$")
    ax.set_ylabel("recovered / true")
    ax.tick_params(labelsize=7.5)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.legend(loc="upper right", frameon=False, fontsize=8)
    fig.tight_layout()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
