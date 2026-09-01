#!/usr/bin/env python
"""Figures 1 and 2 of the ICBINB paper, read from committed tables.

    python paper/icbinb/make_figs.py

The ceiling figure plots each arm's correlation against THAT ARM's prevalence,
from ``depth_confound_per_arm``. It could not before: the shipped table keeps
only ``worst_rho``, the maximum over arms, and pairing that with one arm's
prevalence is the mispairing the paper reports making three times.
``src/reference/jobs/depth_confound_per_arm.py`` persists the per-arm rows,
which is what made an honest scatter possible.
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
COLOURS = {"best4": "#d62728", "lineage": "#1f77b4",
           "crypt_position": "#2ca02c", "epithelial": "#9467bd"}


def bound(p):
    return np.sqrt(3 * p * (1 - p))


def crossing(tol: float = TOLERANCE) -> float:
    """Smallest prevalence at which the tolerance is attainable."""
    from scipy.optimize import brentq

    return brentq(lambda p: bound(p) - tol, 1e-9, 0.5)


def fig_ceiling() -> None:
    per_arm = pd.read_parquet(newest("depth_confound_per_arm"))
    p_cross = crossing()

    fig, (left, right) = plt.subplots(1, 2, figsize=(9.2, 2.9))

    grid = np.geomspace(1e-4, 1.0, 500)
    left.plot(grid, bound(grid), lw=1.6, color="#111111", zorder=3,
              label=r"bound $\sqrt{3p(1-p)}$")
    left.axhline(TOLERANCE, ls="--", lw=1.1, color="#d62728", zorder=3,
                 label=f"tolerance {TOLERANCE:.2f}")
    left.axvspan(1e-4, p_cross, color="#ececec", zorder=0)
    left.axvspan(1 - p_cross, 1.0, color="#ececec", zorder=0)

    for rung, g in per_arm.groupby("granularity_rung"):
        ok = g[np.isfinite(g["prevalence"]) & np.isfinite(g["abs_rho"])]
        if ok.empty:
            continue
        left.scatter(ok["prevalence"].clip(1.2e-4, 1.0), ok["abs_rho"], s=13,
                     alpha=0.55, color=COLOURS.get(rung, "#777777"),
                     edgecolor="white", lw=0.3, zorder=4, label=rung)

    left.set_xscale("log")
    left.set_xlim(1e-4, 1.0)
    left.set_ylim(0, 0.95)
    left.set_xlabel("prevalence $p$ of the mature call, within arm")
    left.set_ylabel(r"$|\rho|$ against depth")
    left.set_title("each arm at its own prevalence,\nagainst what that "
                   "prevalence permits", fontsize=8.5, linespacing=1.4)
    left.legend(frameon=False, fontsize=6.5, loc="upper left", ncol=2)

    order = ["epithelial", "crypt_position", "lineage", "best4"]
    present = [r for r in order if r in set(per_arm["granularity_rung"])]
    data = [per_arm.loc[per_arm["granularity_rung"] == r, "rho_vs_ceiling"]
            .replace([np.inf, -np.inf], np.nan).dropna().to_numpy()
            for r in present]
    parts = right.boxplot(data, vert=False, widths=0.6, patch_artist=True,
                          medianprops=dict(color="#111111", lw=1.4),
                          flierprops=dict(marker=".", ms=3, alpha=0.5))
    for patch, r in zip(parts["boxes"], present, strict=False):
        patch.set_facecolor(COLOURS.get(r, "#777777"))
        patch.set_alpha(0.35)
    right.set_yticklabels(present, fontsize=7.5)
    right.set_xlabel(r"$|\rho|$ as a share of its own bound")
    right.set_title("normalised by what each could reach",
                    fontsize=8.5, linespacing=1.4)

    for ax in (left, right):
        ax.tick_params(labelsize=7.5)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)

    fig.tight_layout()
    fig.savefig(OUT / "fig_ceiling.pdf", bbox_inches="tight")

    unreachable = int((~per_arm["tolerance_is_reachable"]).sum())
    print(f"  crossing at p = {p_cross*100:.4f}%")
    print(f"  {unreachable} of {len(per_arm)} arm-rows cannot reach the tolerance")
    for rung, g in per_arm.groupby("granularity_rung"):
        bad = int((~g["tolerance_is_reachable"]).sum())
        print(f"    {rung:16s} median p {g['prevalence'].median():.4f}  "
              f"unreachable {bad}/{len(g)}  "
              f"median rho/ceiling {g['rho_vs_ceiling'].median():.3f}")
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
        rows = np.clip(per[per["tier"] == t]["rel"].to_numpy(), -1.05, 0.5)
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
