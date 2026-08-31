#!/usr/bin/env python
"""Figure 1, read from the committed calibration table.

Every number on the axes comes from
the newest ``results/<date>_<sha>/calibration_gap_bins_r500.parquet``, and the
path it used is printed. Nothing is transcribed. Run from the repo root:

    python paper/wmhs/make_fig1.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
from _tables import newest

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

TABLE_NAME = "calibration_gap_bins_r500"
OUT = Path("paper/wmhs/figures/fig1_calibration.pdf")

COVERAGE_TARGET = 0.90
DISCRIMINATION_TARGET = 0.80
COMMITTED_OK, COMMITTED_WIDE = 50, 20

COV = "#1f77b4"
DISC = "#d62728"
SHADE = "#ececec"


def panel(ax, rows: pd.DataFrame, *, title: str, subtitle: str) -> None:
    grouped = rows.groupby("n_cells_mature").agg(
        coverage=("coverage", "median"),
        discrimination=("discrimination", "median"),
        cov_lo=("coverage", "min"),
        cov_hi=("coverage", "max"),
        disc_lo=("discrimination", "min"),
        disc_hi=("discrimination", "max"),
    )
    x = grouped.index.to_numpy()

    ax.axvspan(1, COMMITTED_WIDE, color=SHADE, zorder=0)
    ax.axvspan(COMMITTED_WIDE, COMMITTED_OK, color=SHADE, alpha=0.55, zorder=0)
    for edge in (COMMITTED_WIDE, COMMITTED_OK):
        ax.axvline(edge, color="#999999", lw=0.8, zorder=1)

    ax.axhline(COVERAGE_TARGET, ls="--", lw=0.9, color=COV, zorder=2)
    ax.axhline(DISCRIMINATION_TARGET, ls="--", lw=0.9, color=DISC, zorder=2)

    # Seed-to-seed spread, so the reader sees the crossing's precision directly.
    ax.fill_between(x, grouped.cov_lo, grouped.cov_hi, color=COV, alpha=0.15, lw=0)
    ax.fill_between(x, grouped.disc_lo, grouped.disc_hi, color=DISC, alpha=0.15, lw=0)
    ax.plot(x, grouped.coverage, "-o", ms=3.5, lw=1.4, color=COV, label="coverage")
    ax.plot(x, grouped.discrimination, "-s", ms=3.5, lw=1.4, color=DISC,
            label="discrimination")

    ax.set_xscale("log")
    ax.set_xlim(1, 1000)
    ax.set_ylim(0, 1.05)
    ax.set_xticks([1, 10, 20, 50, 90, 200, 500])
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_yticks([0, 0.5, DISCRIMINATION_TARGET, COVERAGE_TARGET, 1])
    ax.set_yticklabels(["0", ".5", ".8", ".9", "1"])
    ax.set_xlabel("cells available to the intrinsic comparison, $n$")
    ax.set_title(f"{title}\n{subtitle}", fontsize=8.5, linespacing=1.5)
    ax.tick_params(labelsize=7.5)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def main() -> int:
    bins = pd.read_parquet(newest(TABLE_NAME))
    ext = bins[bins["grid"] == "extended"]

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.3), sharey=True)

    panel(
        axes[0],
        ext[ext["pool"] == "pooled"],
        title="draw pool as implemented (both tissues pooled)",
        subtitle="no $n$ meets both targets: the routine returns no cutpoint",
    )
    axes[0].set_ylabel("probability")

    ref = ext[ext["pool"] == "reference"]
    panel(
        axes[1],
        ref,
        title="a second reading of the same choice (reference tissue only)",
        subtitle="calibrated $ok = 90$ (7 of 8 seeds); committed value is 50",
    )
    axes[1].axvline(90, color="#111111", ls=":", lw=1.2, zorder=3)
    axes[1].annotate(
        "calibrated\n$ok = 90$", xy=(90, 0.12), xytext=(150, 0.06),
        fontsize=7.5, color="#111111",
        arrowprops=dict(arrowstyle="-", lw=0.7, color="#111111"),
    )

    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=2, frameon=False,
               fontsize=8, bbox_to_anchor=(0.5, -0.06))
    axes[0].text(1.15, COVERAGE_TARGET + 0.015, "targets fixed before the sweep",
                 fontsize=7, color="#555555", va="bottom")

    fig.tight_layout(rect=(0, 0.02, 1, 1))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
