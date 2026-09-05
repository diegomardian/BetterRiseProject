"""Stage 4, step 2: the instrument gate, then the pre-registered variance arm.

    python -m src.bulk.run_stage4_variance --fractions results/<date>_<sha>/stage4_fractions.parquet

Split from ``run_deconvolution`` because deconvolution is a per-sample fit over
hundreds of samples and this is cheap arithmetic on its output. Re-running the
verdict should not mean re-running the instrument.

ORDER MATTERS AND IS NOT NEGOTIABLE. The locked pre-specification says the
positive control GATES the analysis: "if the instrument fails, no R-squared is
reported", because a low R-squared means either "fraction does not explain this
gene" or "deconvolution does not work here" and those have opposite
consequences. So the gate runs first and a failure returns before any regression
is fitted -- not computed-then-suppressed, which would leave the number sitting
in a variable for someone to report later.

TWO GATES, NOT ONE. The pre-committed instrument check reads the non-epithelial
aggregate. Stage 4's predictor is the epithelial-internal split, and the
committed log-scale reference drives that to exactly 0.0 while the aggregate
still passes at r = 0.881 (see src/bulk/deconvolution.py). So the predictor
check from step 1 is re-read here and also blocks. That is not an amendment to
the locked spec: its own estimability clause already requires an unestimable
fraction to be None rather than 0.0, and a column of zeros for every patient is
that coercion at cohort scale.
"""

from __future__ import annotations

import argparse
import logging
import sys
from collections.abc import Sequence
from pathlib import Path

import pandas as pd

from src.bulk.instrument import InstrumentError, gate_verdict, run_instrument_check
from src.bulk.prespec import load_prespec, matched_null_genes, outcome_genes, require_locked_prespec
from src.bulk.variance_arm import (
    PLATE_AS,
    Attrition,
    VarianceArmError,
    benjamini_hochberg,
    build_design,
    compare_to_null,
    gene_r_squared,
    negative_control_verdict,
    primary_verdict,
    resolve_r_squared_kinds,
    secondary_verdict,
)
from src.common.io import write_versioned_table
from src.common.paths import PROCESSED_DIR, RESULTS_DIR
from src.common.provenance import DEFAULT_SEED

log = logging.getLogger(__name__)

BULK = PROCESSED_DIR / "bulk"
HOUSEKEEPING = ("ACTB", "GAPDH")


def _newest(name: str) -> Path | None:
    matches = sorted(RESULTS_DIR.glob(f"*/{name}.parquet"))
    return matches[-1] if matches else None


def resolve_outcome_ids(
    symbols: list[str], index_map: pd.DataFrame, available: pd.Index
) -> dict[str, str]:
    """Symbol -> the id it is stored under in the expression matrix.

    The pre-specification names GUCA2A and CDX2; the matrices are keyed by
    unversioned Ensembl id. A symbol that resolves to two ids is NOT tie-broken
    here -- `resolve_symbols` refuses to, deliberately, and a panel gene mapping
    to two ids is a decision someone has to write down. Symbols already present
    as columns are passed through, so a symbol-keyed matrix (tests, GSE39582)
    still works.
    """
    from src.bulk.gene_index import resolve_symbols

    out: dict[str, str] = {}
    resolved, _unmapped, ambiguous = resolve_symbols(index_map, symbols)
    for symbol in symbols:
        if symbol in available:                      # already symbol-keyed
            out[symbol] = symbol
        elif symbol in ambiguous:
            log.error("  %s maps to %s. Pin one in the caller and say why; this "
                      "is not a tie-break.", symbol, ambiguous[symbol])
        elif symbol in resolved and resolved[symbol] in available:
            out[symbol] = resolved[symbol]
    return out


def run_gate(fractions: pd.DataFrame, purity: pd.DataFrame, rung: str) -> tuple[bool, list, str]:
    """The locked positive control, per method. Returns (proceed, results, message)."""
    results = []
    for method in sorted(fractions["method"].unique()):
        try:
            results.append(run_instrument_check(fractions, purity, method=method, rung=rung))
        except InstrumentError as exc:
            log.warning("  %s: instrument check unevaluable -- %s", method, exc)
    proceed, message = gate_verdict(results)
    return proceed, results, message


def run_variance(
    expression: pd.DataFrame,
    fractions: pd.DataFrame,
    covariates: pd.DataFrame | None,
    index_map: pd.DataFrame,
    spec: dict,
    *,
    rung: str,
    method: str,
    covariate_names: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """The pre-registered arm for one (rung, method). Returns (verdicts, nulls)."""
    attrition = Attrition()
    design = build_design(
        fractions[fractions["method"] == method], covariates,
        covariate_names=covariate_names, attrition=attrition,
    )
    aligned = expression.loc[design["sample_id"]]
    design = design.reset_index(drop=True)
    aligned = aligned.reset_index(drop=True)

    genes = outcome_genes(spec)
    kinds = resolve_r_squared_kinds(spec)

    # The spec names symbols; the expression matrix is on the 1.0.0 Ensembl
    # index. Resolving here rather than assuming: a symbol looked up directly
    # would simply not be found, and "skipping GUCA2A" is a message nobody would
    # read as "the analysis did not happen".
    ids = resolve_outcome_ids([*genes, *HOUSEKEEPING], index_map, aligned.columns)

    fits, nulls, null_rows = {}, {}, []
    for gene in [*genes, *HOUSEKEEPING]:
        gene_id = ids.get(gene)
        if gene_id is None:
            log.warning("  %s is not in the expression matrix; skipping", gene)
            continue
        series = aligned[gene_id].rename(gene)
        fits[gene] = gene_r_squared(series, design, covariate_names)

    for gene in genes:
        if gene not in fits:
            continue
        null_genes = matched_null_genes(aligned, index_map, ids[gene], spec)
        nulls[gene] = [
            gene_r_squared(aligned[g], design, covariate_names)
            for g in null_genes if g in aligned.columns
        ]
        log.info("  %s (%s): %d matched null genes", gene, ids[gene], len(nulls[gene]))
        null_rows.extend(
            {"granularity_rung": rung, "method": method, "target_gene": gene,
             "null_gene": f.gene, "marginal_r2": f.marginal_r2, "partial_r2": f.partial_r2}
            for f in nulls[gene]
        )

    rows = []
    for kind in kinds:
        comparisons = {
            g: compare_to_null(fits[g], nulls[g], kind) for g in genes if g in nulls
        }
        if len(comparisons) < 2:
            log.warning("  only %d outcome gene(s) comparable; no verdict", len(comparisons))
            continue
        guca2a, cdx2 = comparisons["GUCA2A"], comparisons["CDX2"]
        primary, primary_detail = primary_verdict(guca2a, cdx2)
        secondary, secondary_detail = secondary_verdict(guca2a, cdx2)
        control, control_detail = negative_control_verdict(
            [fits[h] for h in HOUSEKEEPING if h in fits], guca2a.null_median, kind,
        )
        adjusted = benjamini_hochberg(
            {g: 1.0 - c.percentile for g, c in comparisons.items()}
        )
        for gene, comparison in comparisons.items():
            rows.append({
                "cohort": "TCGA-COAD/READ", "granularity_rung": rung, "method": method,
                "r_squared_kind": kind, "gene": gene, "n_samples": fits[gene].n,
                "r_squared": comparison.r2, "covariate_r2": fits[gene].covariate_r2,
                "n_null_genes": comparison.n_null,
                "null_median": comparison.null_median,
                "null_p05": comparison.null_p05, "null_p95": comparison.null_p95,
                "percentile": comparison.percentile, "excess": comparison.excess,
                "exceeds_null": comparison.exceeds_null,
                "p_value": 1.0 - comparison.percentile,
                "p_value_bh": adjusted[gene],
                "primary_verdict": primary, "primary_detail": primary_detail,
                "secondary_verdict": secondary, "secondary_detail": secondary_detail,
                "negative_controls": control, "negative_controls_detail": control_detail,
                "plate_entered_as": PLATE_AS,
                "covariates": ", ".join(covariate_names) or "none",
            })
        log.info("  [%s] primary=%s secondary=%s controls=%s",
                 kind, primary, secondary, control)
        log.info("      %s", primary_detail)

    return pd.DataFrame(rows), pd.DataFrame(null_rows)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fractions", type=Path, default=None,
                        help="stage4_fractions.parquet; newest committed if omitted")
    parser.add_argument("--predictor-checks", type=Path, default=None)
    parser.add_argument("--expression", type=Path,
                        default=BULK / "tcga_log2cpm_1.0.0.parquet")
    parser.add_argument("--purity", type=Path, default=None)
    parser.add_argument("--covariates", type=Path, default=None)
    parser.add_argument("--rung", default="lineage")
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    spec = load_prespec()
    require_locked_prespec(spec)      # refuses a proposed spec. It is locked.
    log.info("pre-specification %s, locked %s by %s",
             spec["version"], spec["locked_on"], spec["locked_by"])

    fractions_path = args.fractions or _newest("stage4_fractions")
    if fractions_path is None:
        raise SystemExit(
            "no stage4_fractions.parquet. Run `python -m src.bulk.run_deconvolution "
            "--rung lineage --linearise-reference` first."
        )
    fractions = pd.read_parquet(fractions_path)
    fractions = fractions[fractions["granularity_rung"] == args.rung]
    if fractions.empty:
        raise SystemExit(f"{fractions_path} carries no {args.rung} rows")
    log.info("fractions: %s (%d rows, %s rung)", fractions_path.name, len(fractions), args.rung)

    # --- Gate 2 first, because it is the cheap one and it blocks. -----------
    checks_path = args.predictor_checks or _newest("stage4_predictor_checks")
    if checks_path is not None:
        checks = pd.read_parquet(checks_path)
        checks = checks[checks["granularity_rung"] == args.rung]
        unusable = checks[checks["verdict"] != "usable"]
        if not checks.empty and len(unusable) == len(checks):
            log.error("STOP. No method produced a usable predictor at the %s rung:", args.rung)
            for _, row in unusable.iterrows():
                log.error("  %s: %s -- %s", row["method"], row["verdict"], row["detail"])
            return 3
        for _, row in unusable.iterrows():
            log.warning("  excluding %s: %s", row["method"], row["verdict"])
        fractions = fractions[~fractions["method"].isin(unusable["method"])]

    # --- Gate 1: the locked positive control. ------------------------------
    purity_path = args.purity or _newest("tcga_purity")
    if purity_path is None:
        raise SystemExit("no tcga_purity.parquet; the instrument gate cannot run.")
    purity = pd.read_parquet(purity_path)
    log.info("instrument gate against %s", purity_path.name)
    proceed, gate_results, message = run_gate(fractions, purity, args.rung)
    log.info("%s", message)

    gate_frame = pd.DataFrame([r.as_row() for r in gate_results])
    if not gate_frame.empty:
        gate_frame.insert(0, "cohort", "TCGA-COAD/READ")
        path = write_versioned_table(
            gate_frame, "stage4_instrument_gate", seed=args.seed,
            results_dir=args.results_dir, allow_dirty=args.allow_dirty,
            extra_meta={
                "prespec_version": spec["version"], "prespec_locked_on": spec["locked_on"],
                "purity_source": purity_path.name,
                "what_this_answers": (
                    "The locked prespec's positive control: does the deconvolved "
                    "non-epithelial fraction track (1 - ABSOLUTE purity) at "
                    "r >= 0.5? It does NOT certify the epithelial-internal split "
                    "Stage 4 regresses on -- see stage4_predictor_checks."
                ),
            },
        )
        log.info("wrote %s", path)

    if not proceed:
        log.error("\nSTOP per the locked pre-specification. No R-squared is reported.")
        return 4

    # --- The arm. ----------------------------------------------------------
    if not args.expression.exists():
        raise SystemExit(
            f"{args.expression} not found. The gate passed, so the run is "
            f"authorised, but the log2-CPM matrix lives with the data on the "
            f"cluster. Point --expression at tcga_log2cpm_1.0.0.parquet."
        )
    from src.bulk.gene_index import load_gene_index_map
    from src.bulk.normalise import assert_log_scale

    expression = pd.read_parquet(args.expression)
    assert_log_scale(expression, context=f"{args.expression.name} (the outcome)")
    index_map = load_gene_index_map("1.0.0")

    covariates = pd.read_parquet(args.covariates) if args.covariates else None
    covariate_names: list[str] = []
    if covariates is not None:
        from src.bulk.covariates import covariate_names as spec_covariates
        from src.bulk.covariates import load_covariate_set

        wanted = spec_covariates(load_covariate_set(), context="expression_models")
        covariate_names = [c for c in [*wanted, "plate"] if c in covariates.columns]
        log.info("covariates: %s (plate as %s)", ", ".join(covariate_names), PLATE_AS)
    else:
        log.warning("no --covariates: reporting the MARGINAL R-squared only. The "
                    "locked spec asks for the adjusted one too.")

    verdict_frames, null_frames = [], []
    for method in sorted(fractions["method"].unique()):
        log.info("%s / %s rung", method, args.rung)
        try:
            verdicts, nulls = run_variance(
                expression, fractions, covariates, index_map, spec,
                rung=args.rung, method=method, covariate_names=covariate_names,
            )
        except VarianceArmError as exc:
            log.error("  REFUSED: %s", exc)
            continue
        verdict_frames.append(verdicts)
        null_frames.append(nulls)

    if not verdict_frames or all(f.empty for f in verdict_frames):
        log.error("no verdict could be formed on any method.")
        return 5

    for frame, name in ((pd.concat(verdict_frames, ignore_index=True), "stage4_variance_verdicts"),
                        (pd.concat(null_frames, ignore_index=True), "stage4_matched_nulls")):
        path = write_versioned_table(
            frame, name, seed=args.seed,
            results_dir=args.results_dir, allow_dirty=args.allow_dirty,
            extra_meta={
                "prespec_version": spec["version"],
                "prespec_locked_on": spec["locked_on"],
                "prespec_locked_by": spec["locked_by"],
                "outcome_scale": spec["model"]["outcome_scale"],
                "outcome_genes": outcome_genes(spec),
                "matched_null_seed": spec["matched_null"]["seed"],
                "plate_entered_as": PLATE_AS,
                "r_squared_kinds": list(resolve_r_squared_kinds(spec)),
                "pooling": "never; TCGA only in this run (invariant 4)",
                "interpretation_note": spec["prediction"]["interpretation_note"],
                "what_this_answers": (
                    "The locked Stage 4 prediction. A percentile within an "
                    "abundance-matched null, never a raw R-squared and never two "
                    "compared across genes (issue #54). Both the partial and the "
                    "marginal are carried because the lock does not say which "
                    "the arms are stated on."
                ),
            },
        )
        log.info("wrote %s", path)

    verdicts = pd.concat(verdict_frames, ignore_index=True)
    disagreement = verdicts.groupby(["method", "r_squared_kind"])["primary_verdict"].first()
    if disagreement.groupby(level=0).nunique().gt(1).any():
        log.warning("\nThe partial and marginal R-squared give DIFFERENT primary "
                    "verdicts. The locked spec does not say which its arms are "
                    "stated on, so both are reported and neither is chosen.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
