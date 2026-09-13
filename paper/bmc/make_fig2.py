#!/usr/bin/env python
"""BMC Figure 2 --- cutpoint crossings as brackets, SMC and KUL3 separate.

Every number comes from the newest
``results/<date>_<sha>/cutpoint_crossing_brackets_summary.parquet`` and its
per-seed companion ``cutpoint_crossing_brackets.parquet``. Run from the repo
root:

    python paper/bmc/make_fig2.py

Each row is one (draw pool, criterion). The thick bar is the bracket across
seeds; the thin lines are the individual seeds, so an unstable crossing is
visible as a spread rather than reported as one number. A row with no bar is a
crossing the point rule cannot express at all. The marker is the point
threshold the bracket replaces.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams.update({"axes.labelsize": 7.5})

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _tables import table  # noqa: E402

SUMMARY = "cutpoint_crossing_brackets_summary"
PER_SEED = "cutpoint_crossing_brackets"
OUT = ROOT / "paper" / "bmc" / "figures" / "fig2_brackets.pdf"

COHORTS = (("smc", "SMC (GSE132465)"), ("kul3", "KUL3 (GSE144735)"))
ROWS = (
    ("pooled", "wide"),
    ("pooled", "ok"),
    ("reference", "wide"),
    ("reference", "ok"),
)
COLOUR = {"wide": "#1f77b4", "ok": "#d62728"}
LABEL = {"wide": "wide (coverage)", "ok": "ok (coverage + discrimination)"}


def _row(summary, per_seed, cohort, pool, criterion):
    cell = summary[(summary.cohort == cohort) & (summary.pool == pool)
                   & (summary.criterion == criterion)].iloc[0]
    seeds = per_seed[(per_seed.cohort == cohort) & (per_seed.pool == pool)
                     & (per_seed.criterion == criterion)]
    return cell, seeds


def main() -> int:
    summary = pd.read_parquet(table(SUMMARY))
    per_seed = pd.read_parquet(table(PER_SEED))

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.6), sharex=True)
    for ax, (cohort, title) in zip(axes, COHORTS, strict=True):
        for index, (pool, criterion) in enumerate(ROWS):
            y = len(ROWS) - 1 - index
            cell, seeds = _row(summary, per_seed, cohort, pool, criterion)
            colour = COLOUR[criterion]

            if cell["n_with_upper"] == 0:
                ax.text(22, y, "not estimable", fontsize=6.5, color="0.45",
                        va="center")
                continue

            # Individual seeds as thin lines, jittered so identical brackets
            # are visible as repeated rather than as one.
            for offset, (_, seed) in enumerate(seeds.iterrows()):
                low, high = seed["lower_n_cells_mature"], seed["upper_n_cells_mature"]
                jitter = 0.10 * (offset - (len(seeds) - 1) / 2)
                if np.isfinite(low):
                    ax.plot([low, high], [y + jitter, y + jitter], color=colour,
                            alpha=0.28, linewidth=0.6, zorder=2)
                else:
                    ax.plot([high], [y + jitter], marker="|", color=colour,
                            alpha=0.35, markersize=4, zorder=2)

            stable = bool(cell["bracket_stable"])
            lo, hi = cell["lower_min"], cell["upper_max"]
            if not np.isfinite(lo):
                lo = cell["upper_min"]
            ax.plot([lo, hi], [y, y], color=colour, linewidth=3.2,
                    solid_capstyle="butt", alpha=1.0 if stable else 0.45,
                    zorder=3)
            if not stable:
                ax.text(hi * 1.05, y, "unstable", fontsize=6, color=colour,
                        va="center")

            point = cell["point_min"]
            if np.isfinite(point):
                ax.plot([point], [y], marker="D", color="0.15", markersize=3.2,
                        zorder=4)

        ax.set_yticks(range(len(ROWS)))
        ax.set_yticklabels(
            [f"{pool} · {criterion}" for pool, c in ROWS], fontsize=7,
        )
        ax.set_title(title, fontsize=8)
        ax.set_xscale("log")
        ax.set_xlim(18, 1100)
        ax.set_xticks([20, 50, 100, 200, 400, 800])
        ax.set_xticklabels(["20", "50", "100", "200", "400", "800"])
        ax.tick_params(labelsize=7)
        ax.grid(axis="x", color="0.9", linewidth=0.6)
        ax.set_axisbelow(True)

    axes[0].set_ylabel("draw pool · criterion", fontsize=7.5)
    for ax in axes:
        ax.set_xlabel("mature cells")

    handles = [
        plt.Line2D([], [], color=COLOUR["wide"], linewidth=3, label=LABEL["wide"]),
        plt.Line2D([], [], color=COLOUR["ok"], linewidth=3, label=LABEL["ok"]),
        plt.Line2D([], [], color="0.15", marker="D", linestyle="none",
                   markersize=3.2, label="point threshold"),
        plt.Line2D([], [], color="0.4", linewidth=3, alpha=0.45,
                   label="unstable across seeds"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=7,
               frameon=False, bbox_to_anchor=(0.5, -0.05))
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight")
    print(f"wrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
