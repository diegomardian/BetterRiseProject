#!/usr/bin/env python
"""Figure 3 — the curve cannot tell these apart; one column can.

Reads the newest ``results/<date>_<sha>/trial_recovery.parquet``. Run from the
repo root:

    python paper/wmhs/make_fig4.py

Left panel: the recovery curve for four estimators of the same treatment effect.
Three of the four sit near 1 and tighten with cohort size, and two of those three
are exactly blind — their curves are the same line to machine precision.

Right panel: the residual against the *realised* effect, log scale. It is
identically zero for the two degenerate estimators and non-zero for the other
two, at every cohort size. That is the whole diagnostic.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
from _tables import newest

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

TABLE_NAME = "trial_recovery"
OUT = Path("paper/wmhs/figures/fig4_trial_recovery.pdf")

STYLE = {
    "gcomp-from-generator": ("#1f77b4", "-", "o", "G-computation, generator's strata"),
    "ipw-saturated":        ("#7fb3d5", "--", "s", "IPW, saturated propensity"),
    "ipw-cross-fitted":     ("#2ca02c", "-", "^", "IPW, cross-fitted"),
    "ols-stratum-dummies":  ("#9467bd", "-", "D", "OLS, stratum dummies"),
    "unadjusted":           ("#d62728", "-", "v", "unadjusted difference"),
}
FLOOR = 1e-16  # so exact zeros are drawable on a log axis


def main() -> int:
    df = pd.read_parquet(newest(TABLE_NAME))

    fig, (left, right) = plt.subplots(1, 2, figsize=(9.4, 3.4))

    for name, (colour, ls, marker, label) in STYLE.items():
        rows = df[df["estimator"] == name].sort_values("n_patients")
        if rows.empty:
            continue
        x = rows["n_patients"].to_numpy()
        left.fill_between(x, rows.ratio_q25, rows.ratio_q75, color=colour,
                          alpha=0.12, lw=0)
        left.plot(x, rows.ratio_median, ls=ls, marker=marker, ms=3.5, lw=1.4,
                  color=colour, label=label)
        right.plot(x, rows.max_residual_vs_realised.clip(lower=FLOOR), ls=ls,
                   marker=marker, ms=3.5, lw=1.4, color=colour)

    left.axhline(1.0, ls=":", lw=0.9, color="#444444")
    left.set_xscale("log")
    left.set_ylim(0.3, 2.9)
    left.set_ylabel("recovered / requested")
    left.set_title("the recovery curve\nthree of the four look fine",
                   fontsize=8.5, linespacing=1.5)

    right.set_xscale("log")
    right.set_yscale("log")
    right.set_ylim(FLOOR / 5, 1e2)
    right.axhspan(FLOOR / 5, 1e-12, color="#ececec", zorder=0)
    right.set_ylabel(r"max $|\hat\theta - \theta_{\mathrm{realised}}|$")
    right.set_title("the one-line check\nzero means the curve saw nothing",
                    fontsize=8.5, linespacing=1.5)
    right.text(120, 2e-15, "identically zero", fontsize=7, color="#555555")

    for ax in (left, right):
        ax.set_xlabel("patients per cohort")
        ax.tick_params(labelsize=7.5)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)

    handles, labels = left.get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, frameon=False,
               fontsize=7.5, bbox_to_anchor=(0.5, -0.07))

    fig.tight_layout(rect=(0, 0.03, 1, 1))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight")
    print(f"wrote {OUT}")

    worst = df.groupby("estimator")["max_residual_vs_realised"].max()
    for name, value in worst.items():
        print(f"  {name:22s} max residual {value:.3g}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
