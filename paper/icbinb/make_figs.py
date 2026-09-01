#!/usr/bin/env python
"""Figures 1 and 2 of the ICBINB paper, read from committed tables.

    python paper/icbinb/make_figs.py

WHAT THIS DELIBERATELY DOES NOT PLOT. An earlier caption described the ceiling
figure as every patient-by-resolution point at "its own prevalence and
correlation". No committed table carries a per-arm correlation --- only
``worst_rho``, the maximum over the two arms --- so drawing that scatter would
mean pairing a max-over-arms statistic with one arm's prevalence, which is
precisely the mispairing the paper reports making three times. The figure shows
what is backed: the bound, the tolerance, the crossing, and where this cohort's
prevalences actually fall against it.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
from _tables import newest

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import yaml  # noqa: E402

OUT = Path("paper/icbinb/figures")
TOLERANCE = 0.20


def bound(p):
    return np.sqrt(3 * p * (1 - p))


def crossing(tol=TOLERANCE):
    """Smallest prevalence at which the tolerance is attainable."""
    from scipy.optimize import brentq

    return brentq(lambda p: bound(p) - tol, 1e-9, 0.5)


def fig_ceiling() -> None:
    counts = pd.read_parquet(newest("mature_cell_counts_full"))
    p_cross = crossing()

    fig, ax = plt.subplots(figsize=(9.0, 2.7))
    grid = np.geomspace(1e-4, 0.5, 400)
    ax.plot(grid, bound(grid), lw=1.6, color="#1f77b4",
            label=r"attainable bound $\sqrt{3p(1-p)}$")
    ax.axhline(TOLERANCE, ls="--", lw=1.1, color="#d62728",
               label=f"tolerance $|\\rho| \\geq {TOLERANCE:.2f}$")
    ax.axvline(p_cross, ls=":", lw=1.1, color="#444444")
    ax.axvspan(1e-4, p_cross, color="#ececec", zorder=0)
    ax.text(1.6e-4, 0.62, "the check cannot fire\nat any prevalence in here",
            fontsize=7.5, color="#555555", va="top")
    ax.annotate(f"$p = {p_cross*100:.4f}\\%$", xy=(p_cross, 0.02),
                xytext=(p_cross * 1.5, 0.10), fontsize=7.5, color="#111111",
                arrowprops=dict(arrowstyle="-", lw=0.7, color="#111111"))

    # Where this cohort's labels actually sit. Prevalence only: no committed
    # table carries a per-arm correlation to pair with it.
    tumour = counts[counts["tissue"] != "normal"]
    for rung, marker in zip(sorted(tumour["granularity_rung"].unique()),
                            "os^vD", strict=False):
        rows = tumour[tumour["granularity_rung"] == rung]["mature_fraction"]
        rows = rows[(rows > 0) & (rows < 1)]
        if rows.empty:
            continue
        med = float(rows.median())
        ax.plot([med], [bound(med)], marker=marker, ms=6, color="#2ca02c",
                mec="white", mew=0.8, ls="none", label=f"{rung} (median $p$)")

    ax.set_xscale("log")
    ax.set_xlim(1e-4, 0.5)
    ax.set_ylim(0, 0.95)
    ax.set_xlabel("prevalence $p$ of the labelled state (diseased arm)")
    ax.set_ylabel(r"largest attainable $|\rho|$")
    ax.tick_params(labelsize=7.5)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.legend(frameon=False, fontsize=7.5, loc="upper left")
    fig.tight_layout()
    fig.savefig(OUT / "fig_ceiling.pdf", bbox_inches="tight")
    print(f"  crossing at p = {p_cross*100:.4f}%")
    print(f"  wrote {OUT/'fig_ceiling.pdf'}")


def fig_tiers() -> None:
    df = pd.read_parquet(newest("decomposition_summary"))
    panel = yaml.safe_load(open("config/panel.yaml"))
    tier = {g: t for t, v in panel["tiers"].items() for g in v["genes"]}
    df = df.assign(tier=df["gene"].map(tier))
    df = df[df["tier"].isin(["A", "B", "D"])].copy()
    df["rel"] = df["mean_tumour"] / df["mean_normal"] - 1
    df = df[np.isfinite(df["rel"])]

    per = (df.groupby(["tier", "gene", "granularity_rung", "labeling_axis"])
             ["rel"].median().reset_index())

    label = {"A": "A · compositional control", "B": "B · intrinsic control",
             "D": "D · retained control"}
    colour = {"A": "#1f77b4", "B": "#2ca02c", "D": "#d62728"}

    fig, ax = plt.subplots(figsize=(9.0, 2.2))
    for i, t in enumerate(["A", "B", "D"]):
        rows = per[per["tier"] == t]["rel"].to_numpy()
        rows = np.clip(rows, -1.05, 0.5)
        ax.scatter(rows, np.full_like(rows, i) + np.random.default_rng(0)
                   .uniform(-0.13, 0.13, len(rows)),
                   s=22, color=colour[t], alpha=0.65, edgecolor="white", lw=0.5)
        ax.plot([np.median(rows)], [i], marker="|", ms=22, mew=2.2,
                color=colour[t])

    ax.axvline(-1.0, ls=":", lw=1.0, color="#444444")
    ax.axvline(0.0, ls="-", lw=0.8, color="#999999")
    ax.text(-0.995, 2.42, "total loss ($m_T = 0$)", fontsize=7, color="#555555")
    ax.text(0.01, 2.42, "no change", fontsize=7, color="#555555")
    ax.set_yticks(range(3))
    ax.set_yticklabels([label[t] for t in ["A", "B", "D"]], fontsize=8)
    ax.set_ylim(-0.5, 2.6)
    ax.set_xlim(-1.1, 0.35)
    ax.set_xlabel("relative per-cell change $m_T/m_N - 1$, one point per gene "
                  "$\\times$ stratum")
    ax.tick_params(labelsize=7.5)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT / "fig_tiers.pdf", bbox_inches="tight")

    for t in ["A", "B", "D"]:
        r = per[per["tier"] == t]["rel"]
        print(f"  tier {t}: {r.min():.4f} .. {r.max():.4f}")
    print(f"  wrote {OUT/'fig_tiers.pdf'}")


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    fig_ceiling()
    fig_tiers()
    return 0


if __name__ == "__main__":
    sys.exit(main())
