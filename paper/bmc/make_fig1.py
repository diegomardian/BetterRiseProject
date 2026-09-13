#!/usr/bin/env python
"""BMC Figure 1 --- false-positive rate by interval method and generator.

Every number on the axes comes from the newest
``results/<date>_<sha>/interval_stress_calibration_summary.parquet``. Nothing is
transcribed. Run from the repo root:

    python paper/bmc/make_fig1.py

Six panels, one per generator; a line per interval method; the shaded band is
the min--max across the seven cohorts at that sample size. The dashed rule is
the nominal 5% and the dotted rule is nominal plus the two-point tolerance, so
"which methods stay inside the bar" is readable off the figure directly.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams.update({"axes.labelsize": 7.5})

import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _tables import table  # noqa: E402

TABLE_NAME = "interval_stress_calibration_summary"
OUT = ROOT / "paper" / "bmc" / "figures" / "fig1_interval_stress.pdf"

NOMINAL = 5.0
TOLERANCE = 2.0  # MISCALIBRATION_TOLERANCE, in percentage points

REGIMES = (
    ("gaussian", "Gaussian"),
    ("skewed", "Skewed"),
    ("heavy_tailed", "Heavy-tailed"),
    ("spike_slab", "Spike-slab"),
    ("boundary_rare", "Boundary (near zero)"),
    ("boundary_saturated", "Boundary (saturated)"),
)
METHODS = {
    "percentile": ("#d62728", "percentile bootstrap"),
    "bca": ("#ff7f0e", "BCa"),
    "student_t": ("#1f77b4", "Student-$t$"),
}


def main() -> int:
    frame = pd.read_parquet(table(TABLE_NAME))
    print(f"cohorts: {frame['cohort'].nunique()}, cells: {len(frame)}")
    worst_student = frame[frame["method"] == "student_t"]["fpr_median"].max()
    print(f"worst Student-t cell: {100 * worst_student:.1f}%")

    fig, axes = plt.subplots(2, 3, figsize=(7.0, 3.6), sharex=True, sharey=True)
    for ax, (key, label) in zip(axes.ravel(), REGIMES, strict=True):
        block = frame[frame["regime"] == key]
        for method, (colour, _) in METHODS.items():
            rows = block[block["method"] == method]
            agg = (
                rows.groupby("n_patients")
                .agg(median=("fpr_median", "median"),
                     low=("fpr_min", "min"),
                     high=("fpr_max", "max"))
                .reset_index()
                .sort_values("n_patients")
            )
            ax.plot(agg["n_patients"], 100 * agg["median"], "-o", color=colour,
                    markersize=2.5, linewidth=1.0)
            ax.fill_between(agg["n_patients"], 100 * agg["low"],
                            100 * agg["high"], color=colour, alpha=0.12,
                            linewidth=0)
        ax.axhline(NOMINAL, color="0.3", linestyle="--", linewidth=0.7)
        ax.axhline(NOMINAL + TOLERANCE, color="0.3", linestyle=":",
                   linewidth=0.7)
        ax.set_title(label, fontsize=8)
        ax.set_xscale("log")
        ticks = sorted(frame["n_patients"].unique())
        ax.set_xticks(ticks)
        ax.set_xticklabels([str(t) for t in ticks])
        ax.set_ylim(0, 25)
        ax.tick_params(labelsize=7)

    for ax in axes[:, 0]:
        ax.set_ylabel("false-positive rate (%)")
    for ax in axes[1, :]:
        ax.set_xlabel("patients per cohort")

    handles = [
        plt.Line2D([], [], color=colour, marker="o", markersize=3, label=label)
        for colour, label in METHODS.values()
    ]
    handles.append(plt.Line2D([], [], color="0.3", linestyle="--",
                              label=f"nominal {NOMINAL:.0f}%"))
    handles.append(plt.Line2D([], [], color="0.3", linestyle=":",
                              label=f"nominal + {TOLERANCE:.0f}pp"))
    fig.legend(handles=handles, loc="lower center", ncol=5, fontsize=7,
               frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight")
    print(f"wrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
