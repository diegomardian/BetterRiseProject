"""Run the locked D2 survival specification; never choose a model at runtime.

    python -m src.bulk.run_d2_survival

The outcome-blind feasibility gate and the locked D2 pre-registration precede
this job. This runner consumes the committed TCGA inputs, reports every
complete-case exclusion, and fits only the fixed continuous-GUCA2A Cox models.
It deliberately has no cut point, Kaplan--Meier plot, interaction, or subgroup
branch.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

from src.bulk.covariates import load_covariate_set
from src.bulk.d2_feasibility_gate import (
    GENE,
    INDEX_VERSION,
    primary_tumour_values,
    resolve_gene,
)
from src.bulk.d2_survival import (
    D2_ENDPOINTS,
    LEAD_ENDPOINT,
    d2_events_per_df,
    d2_model_spec,
    prepare_design,
    retained_vs_dropped,
)
from src.bulk.gdc import read_manifest
from src.bulk.normalise import assert_log_scale
from src.bulk.survival import fit_cox, proportional_hazards_check
from src.common.io import write_versioned_table
from src.common.paths import CONFIG_DIR, PROCESSED_DIR, RESULTS_DIR
from src.common.provenance import DEFAULT_SEED

log = logging.getLogger(__name__)

BULK = PROCESSED_DIR / "bulk"
PURITY_METHODS = (("absolute", False), ("estimate_affy_extrapolated", True))
MODEL_CLINICAL = "clinical_adjusted"
MODEL_PRIMARY = "purity_adjusted"


class D2RunError(RuntimeError):
    """The locked D2 outcome run cannot be assembled or read."""


def _add_model_columns(
    frame: pd.DataFrame,
    *,
    model: str,
    purity_method: str | None,
) -> pd.DataFrame:
    out = frame.copy()
    out["model"] = model
    out["purity_method"] = purity_method
    return out


def _term(coefficients: pd.DataFrame, endpoint: str, purity_method: str) -> pd.Series:
    hit = coefficients.loc[
        (coefficients["endpoint"] == endpoint)
        & (coefficients["model"] == MODEL_PRIMARY)
        & (coefficients["purity_method"] == purity_method)
        & (coefficients["term"] == GENE)
    ]
    if len(hit) != 1:
        raise D2RunError(
            f"expected one {GENE} coefficient for {endpoint}/{purity_method}, got {len(hit)}"
        )
    return hit.iloc[0]


def d2_verdict(
    statistics: pd.DataFrame,
    coefficients: pd.DataFrame,
    ph: pd.DataFrame,
) -> pd.DataFrame:
    """Apply the locked §5--6 decision tree; no p-value search is possible."""
    lead = statistics.loc[
        (statistics["endpoint"] == LEAD_ENDPOINT)
        & (statistics["model"] == MODEL_PRIMARY)
        & (statistics["purity_method"] == "absolute")
    ]
    if len(lead) != 1:
        raise D2RunError("D2 has no unique primary-ABSOLUTE PFI estimability row")
    if not bool(lead.iloc[0]["meets_lead_floor"]):
        return pd.DataFrame([{
            "verdict": "NOT ESTIMABLE",
            "detail": (
                "PFI has fewer than 10 events per non-stratum df; no Cox-based "
                "directional conclusion is read."
            ),
            "lead_estimable": False,
            "lead_interval_excludes_one": None,
            "ph_violated": None,
            "endpoint_discordance": None,
            "dss_interval_includes_one": None,
            "purity_source_sensitive": None,
        }])

    pfi = _term(coefficients, "PFI", "absolute")
    dss = _term(coefficients, "DSS", "absolute")
    sensitivity = _term(coefficients, "PFI", "estimate_affy_extrapolated")
    # ENDPOINT-SPECIFIC, and it was not. §5 labels "the corresponding endpoint"
    # PH VIOLATED; without this filter a violation in the OS model -- which §6
    # says is descriptive and cannot reverse the primary conclusion -- silences
    # the lead endpoint, under a detail line naming PFI.
    ph_terms = ph.loc[
        (ph["endpoint"] == LEAD_ENDPOINT)
        & (ph["model"] == MODEL_PRIMARY)
        & (ph["purity_method"] == "absolute")
        & ph["term"].isin([GENE, "purity"])
        & ph["violates_ph"]
    ]
    ph_violated = not ph_terms.empty
    # §6's rules are ORDERED: "if the PFI interval includes 1, no supported
    # association" comes BEFORE "if the hazard ratio is at least 1, direction
    # contrary". Testing `hazard_ratio >= 1` first reported a null interval with
    # a point estimate a hair above 1 as a contrary FINDING, which is the
    # direction of overstatement this document exists to prevent.
    pfi_excludes_one = bool(pfi["ci_high"] < 1 or pfi["ci_low"] > 1)
    pfi_supported = bool(pfi["hazard_ratio"] < 1 and pfi["ci_high"] < 1)
    # §6 defines discordance and then QUALIFIES it in the next sentence: "An
    # imprecise DSS interval alone is reported as such, not treated as a
    # contradiction." The sign test alone ignored that qualifier, so two null
    # endpoints whose point estimates happen to straddle 1 -- PFI 0.970 and DSS
    # 1.001, both intervals containing 1 -- were reported as a contradiction
    # between endpoints. Discordance now requires the DSS interval to exclude 1,
    # which is what makes an opposite direction a claim rather than noise.
    dss_interval_includes_one = bool(dss["ci_low"] <= 1 <= dss["ci_high"])
    dss_opposite = bool(
        (pfi["hazard_ratio"] - 1) * (dss["hazard_ratio"] - 1) < 0
        and not dss_interval_includes_one
    )
    ci_decision_changed = bool((pfi["ci_high"] < 1) != (sensitivity["ci_high"] < 1))
    direction_reversed = bool(
        (pfi["hazard_ratio"] - 1) * (sensitivity["hazard_ratio"] - 1) < 0
    )
    sensitive = ci_decision_changed or direction_reversed

    if ph_violated:
        verdict = "PH VIOLATED"
        detail = (
            "PFI primary model has a Schoenfeld violation for GUCA2A or purity; "
            "no Cox-based directional conclusion is read."
        )
    elif sensitive:
        verdict = "PURITY-SOURCE SENSITIVE"
        # The label is pre-specified and stands. The numbers go in the detail so
        # it cannot be read as a robust association whose robustness is in
        # question -- when both intervals contain 1 there is no association for
        # the purity source to be sensitive about, and a reader must be able to
        # see that from the verdict row alone.
        detail = (
            f"PFI direction or confidence-interval decision changes between "
            f"ABSOLUTE and ESTIMATE purity. ABSOLUTE HR "
            f"{float(pfi['hazard_ratio']):.4f} "
            f"[{float(pfi['ci_low']):.4f}, {float(pfi['ci_high']):.4f}]; "
            f"ESTIMATE HR {float(sensitivity['hazard_ratio']):.4f} "
            f"[{float(sensitivity['ci_low']):.4f}, "
            f"{float(sensitivity['ci_high']):.4f}]. "
            + ("BOTH INTERVALS INCLUDE 1: the lead result is null under either "
               "purity source, and the label is a caveat on a null."
               if not pfi_excludes_one
               and not bool(sensitivity["ci_high"] < 1 or sensitivity["ci_low"] > 1)
               else "")
        )
    elif dss_opposite:
        verdict = "ENDPOINT DISCORDANCE"
        detail = (
            "DSS and PFI GUCA2A point estimates have opposite directions; no "
            "general survival claim is read."
        )
    elif pfi_supported:
        verdict = "SUPPORTED ASSOCIATION"
        detail = "PFI ABSOLUTE-purity-adjusted GUCA2A HR is below 1 with a 95% CI excluding 1."
    elif not pfi_excludes_one:
        verdict = "NO SUPPORTED ASSOCIATION"
        detail = (
            f"PFI ABSOLUTE-purity-adjusted GUCA2A interval includes 1 "
            f"(HR {float(pfi['hazard_ratio']):.3f})."
        )
    else:
        verdict = "DIRECTION CONTRARY"
        detail = "PFI ABSOLUTE-purity-adjusted GUCA2A interval excludes 1 and lies above it."
    return pd.DataFrame([{
        "verdict": verdict,
        "detail": detail,
        "lead_estimable": True,
        "lead_interval_excludes_one": pfi_excludes_one,
        "ph_violated": ph_violated,
        "endpoint_discordance": dss_opposite,
        # §6: "An imprecise DSS interval alone is reported as such, not treated
        # as a contradiction." Recorded so the reader is not left to infer which
        # of the two a discordance row means.
        "dss_interval_includes_one": dss_interval_includes_one,
        "purity_source_sensitive": sensitive,
    }])


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clinical", type=Path, default=BULK / "clinical_curated.tsv")
    parser.add_argument("--purity", type=Path, default=BULK / "tcga_purity_0.9.0.parquet")
    parser.add_argument("--manifest", type=Path, default=BULK / "sample_manifest.tsv")
    parser.add_argument(
        "--expression", type=Path, default=BULK / f"tcga_log2cpm_{INDEX_VERSION}.parquet"
    )
    parser.add_argument(
        "--gene-index-map", type=Path,
        default=CONFIG_DIR / "gene_index" / f"gene_index_{INDEX_VERSION}.map.tsv",
    )
    parser.add_argument("--results-dir", type=Path, default=RESULTS_DIR)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    for path in (args.clinical, args.purity, args.manifest, args.expression, args.gene_index_map):
        if not path.exists():
            raise SystemExit(f"{path} not found")
    clinical = pd.read_csv(args.clinical, sep="\t")
    purity = pd.read_parquet(args.purity)
    manifest = read_manifest(args.manifest)
    expression = pd.read_parquet(args.expression)
    assert_log_scale(expression, context=args.expression.name)
    gene_id = resolve_gene(pd.read_csv(args.gene_index_map, sep="\t"))
    values = primary_tumour_values(expression, manifest, gene_id)
    spec = load_covariate_set()

    coefficients, ph_tables, attritions, statistics, distributions = [], [], [], [], []
    primary_designs: dict[tuple[str, str], pd.DataFrame] = {}
    for method, sensitivity in PURITY_METHODS:
        for endpoint in D2_ENDPOINTS:
            design, attrition = prepare_design(
                clinical, purity, manifest, expression, gene_id, spec,
                endpoint=endpoint, sensitivity=sensitivity,
            )
            primary_designs[(endpoint, method)] = design
            attritions.append(
                _add_model_columns(attrition, model=MODEL_PRIMARY, purity_method=method)
            )
            stats = d2_events_per_df(design, spec, endpoint=endpoint)
            statistics.append(pd.DataFrame([{
                **stats, "model": MODEL_PRIMARY, "purity_method": method,
            }]))
            distributions.append(retained_vs_dropped(
                clinical, values, design, endpoint=endpoint, purity_method=method
            ))

    # The no-purity model is descriptive and deliberately has its own complete-case population.
    for endpoint in D2_ENDPOINTS:
        design, attrition = prepare_design(
            clinical, purity, manifest, expression, gene_id, spec,
            endpoint=endpoint, descriptive=True,
        )
        attritions.append(_add_model_columns(attrition, model=MODEL_CLINICAL, purity_method=None))
        stats = d2_events_per_df(design, spec, endpoint=endpoint, descriptive=True)
        statistics.append(pd.DataFrame([{**stats, "model": MODEL_CLINICAL, "purity_method": None}]))

    stats_frame = pd.concat(statistics, ignore_index=True)
    lead_ok = bool(stats_frame.loc[
        (stats_frame["endpoint"] == LEAD_ENDPOINT)
        & (stats_frame["model"] == MODEL_PRIMARY)
        & (stats_frame["purity_method"] == "absolute"), "meets_lead_floor"
    ].item())
    if lead_ok:
        for (endpoint, method), design in primary_designs.items():
            locked, context = d2_model_spec(spec)
            fitter, tidy = fit_cox(
                design, locked, endpoint=endpoint, context=context, extra_continuous=(GENE,)
            )
            coefficients.append(_add_model_columns(tidy, model=MODEL_PRIMARY, purity_method=method))
            ph = proportional_hazards_check(
                fitter, design, locked, endpoint=endpoint, context=context,
                extra_continuous=(GENE,),
            )
            ph_tables.append(_add_model_columns(ph, model=MODEL_PRIMARY, purity_method=method))

        for endpoint in D2_ENDPOINTS:
            design, _ = prepare_design(
                clinical, purity, manifest, expression, gene_id, spec,
                endpoint=endpoint, descriptive=True,
            )
            locked, context = d2_model_spec(spec, descriptive=True)
            fitter, tidy = fit_cox(
                design, locked, endpoint=endpoint, context=context, extra_continuous=(GENE,)
            )
            coefficients.append(_add_model_columns(tidy, model=MODEL_CLINICAL, purity_method=None))
            ph = proportional_hazards_check(
                fitter, design, locked, endpoint=endpoint, context=context,
                extra_continuous=(GENE,),
            )
            ph_tables.append(_add_model_columns(ph, model=MODEL_CLINICAL, purity_method=None))

    coefficient_frame = (
        pd.concat(coefficients, ignore_index=True) if coefficients else pd.DataFrame()
    )
    ph_frame = pd.concat(ph_tables, ignore_index=True) if ph_tables else pd.DataFrame()
    verdict = (
        d2_verdict(stats_frame, coefficient_frame, ph_frame)
        if lead_ok else d2_verdict(stats_frame, pd.DataFrame(), pd.DataFrame())
    )
    log.info("D2 verdict: %s", verdict.iloc[0]["verdict"])

    common = {
        "contract": "docs/prereg_d2_survival.md",
        "predictor": "GUCA2A per +1 log2-CPM, continuous",
        "primary_purity": "absolute",
        "sensitivity_purity": "estimate_affy_extrapolated",
        "strata": ["project", "plate"],
        "lead_endpoint": LEAD_ENDPOINT,
    }
    outputs = [
        (stats_frame, "d2_survival_estimability"),
        (pd.concat(attritions, ignore_index=True), "d2_survival_attrition"),
        (pd.concat(distributions, ignore_index=True), "d2_guca2a_retained_vs_dropped"),
        (verdict, "d2_survival_verdict"),
    ]
    if not coefficient_frame.empty:
        outputs.extend([
            (coefficient_frame, "d2_survival_coefficients"),
            (ph_frame, "d2_survival_ph_tests"),
        ])
    for frame, name in outputs:
        path = write_versioned_table(
            frame, name, seed=args.seed, results_dir=args.results_dir,
            allow_dirty=args.allow_dirty, extra_meta=common,
        )
        log.info("wrote %s", path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
