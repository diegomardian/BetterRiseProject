"""W3.7 — run the baseline survival models. Clinical covariates only.

    python -m src.bulk.run_survival

The deliverable is a verdict, not a hazard ratio: does stage come out
prognostic? If it does not, something is wrong upstream and the brief says to
investigate rather than report.
"""

from __future__ import annotations

import argparse
import sys

import pandas as pd

from src.bulk.covariates import build_design, events_per_df, load_covariate_set
from src.bulk.gdc import read_manifest
from src.bulk.survival import (
    fit_cox,
    kaplan_meier_by_stage,
    power_note,
    proportional_hazards_check,
    stage_sanity_check,
)
from src.common.io import versioned_dir, write_versioned_table
from src.common.paths import PROCESSED_DIR

SEED = 20260818

#: Set from --allow-dirty. Default False. These jobs used to pass
#: allow_dirty=True unconditionally, so the bulk arm could not write a
#: clean provenance stamp even from a spotless tree -- which is why every
#: committed bulk table records git_dirty: true.
ALLOW_DIRTY = False
BULK = PROCESSED_DIR / "bulk"
CONTEXT = "clinical_baseline"
ENDPOINTS = ("PFI", "DSS", "OS")


def _plot(design_by_endpoint: dict[str, pd.DataFrame], out_path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from lifelines import KaplanMeierFitter

    endpoints = list(design_by_endpoint)
    fig, axes = plt.subplots(1, len(endpoints), figsize=(5 * len(endpoints), 4.2), squeeze=False)
    for col, endpoint in enumerate(endpoints):
        ax = axes[0][col]
        design = design_by_endpoint[endpoint]
        for stage in ("I", "II", "III", "IV"):
            group = design.loc[design["stage"] == stage]
            if group.empty:
                continue
            fitter = KaplanMeierFitter()
            fitter.fit(
                group[f"{endpoint}.time"],
                group[endpoint],
                label=f"{stage} (n={len(group)})",
            )
            fitter.plot_survival_function(ax=ax, ci_show=False)
        ax.set_title(f"{endpoint} by stage")
        ax.set_xlabel("days")
        ax.set_ylabel("survival probability")
        ax.set_ylim(0, 1.02)
        ax.legend(fontsize=8, frameon=False)
    fig.suptitle(
        "W3.7 baseline — clinical covariates only, stratified by project", fontsize=11
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


def main(argv: list[str] | None = None) -> int:
    global ALLOW_DIRTY
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--allow-dirty", action="store_true",
        help="write from a dirty tree; the recorded sha will not reproduce it",
    )
    ALLOW_DIRTY = parser.parse_args(argv).allow_dirty
    spec = load_covariate_set()
    print(
        f"covariate set {spec['version']} — status={spec['status']}, "
        f"locked by {spec['locked_by']}"
    )

    clinical = pd.read_csv(BULK / "clinical_curated.tsv", sep="\t")
    purity = pd.read_parquet(BULK / "tcga_purity_0.9.0.parquet")
    manifest = read_manifest(BULK / "sample_manifest.tsv")
    patients = set(manifest.loc[manifest["sample_type"] == "01", "patient_id"])
    clinical = clinical.loc[clinical["patient_id"].isin(patients)]

    coefficients, ph_tables, km_tables, verdicts, attritions = [], [], [], [], []
    designs: dict[str, pd.DataFrame] = {}

    for endpoint in ENDPOINTS:
        design, attrition = build_design(
            clinical, purity, spec, endpoint=endpoint, context=CONTEXT
        )
        designs[endpoint] = design
        attritions.append(attrition)

        print(f"\n{'=' * 66}\n{endpoint} — {power_note(design, spec, endpoint, CONTEXT)}")

        fitter, tidy = fit_cox(design, spec, endpoint=endpoint, context=CONTEXT)
        coefficients.append(tidy)
        with pd.option_context("display.width", 160):
            print(tidy[["term", "hazard_ratio", "ci_low", "ci_high", "p"]].to_string(index=False))

        verdict = stage_sanity_check(tidy, endpoint=endpoint)
        verdicts.append(verdict)
        mark = "PASS" if verdict["passed"] else "FAIL"
        print(
            f"  [{mark}] stage IV HR = {verdict['stage_IV_hazard_ratio']} "
            f"({verdict['stage_IV_ci_low']:.2f}–{verdict['stage_IV_ci_high']:.2f}), "
            f"p = {verdict['stage_IV_p']:.3g} · {verdict['verdict']}"
        )

        km, logrank = kaplan_meier_by_stage(design, endpoint=endpoint)
        km["logrank_p"] = logrank
        km_tables.append(km)
        print(f"  log-rank across stages: p = {logrank:.3g}")
        print(
            km[["stage", "n", "n_events", "survival_at_1y", "survival_at_3y"]].to_string(
                index=False
            )
        )

        ph = proportional_hazards_check(
            fitter, design, spec, endpoint=endpoint, context=CONTEXT
        )
        ph_tables.append(ph)
        violations = sorted(set(ph.loc[ph["violates_ph"], "term"]))
        detail = (
            "violated by " + ", ".join(violations)
            if violations
            else "no term violates at p<0.05"
        )
        print(f"  proportional hazards: {detail}")

    stats = pd.DataFrame(
        [events_per_df(designs[e], spec, e, CONTEXT) for e in ENDPOINTS]
    )
    out_dir = versioned_dir(SEED)
    figure = _plot(designs, out_dir / "w3.7_km_by_stage.png")

    for frame, name, note in (
        (pd.concat(coefficients), "tcga_baseline_cox", "W3.7 Cox coefficients, clinical only"),
        (pd.concat(ph_tables), "tcga_baseline_ph_tests", "W3.7 Schoenfeld residual tests"),
        (pd.concat(km_tables), "tcga_baseline_km", "W3.7 Kaplan-Meier by stage"),
        (pd.DataFrame(verdicts), "tcga_baseline_sanity", "W3.7 stage sanity check verdicts"),
        (pd.concat(attritions), "tcga_baseline_attrition", "W3.7 cohort attrition per endpoint"),
        (stats, "tcga_baseline_power", "W3.7 events per degree of freedom"),
    ):
        write_versioned_table(frame, name=name, seed=SEED, notes=note, allow_dirty=ALLOW_DIRTY)

    print(f"\nfigure: {figure}")
    print("wrote six results tables")
    passed = all(v["passed"] for v in verdicts)
    headline = (
        "PASS — stage is prognostic on every endpoint"
        if passed
        else "FAIL — investigate upstream"
    )
    print(f"\nSANITY CHECK: {headline}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
