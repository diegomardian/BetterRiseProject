"""A1: does the reference recipe cost the instrument its margin?

    python -m src.harness.jobs.run_pseudobulk_recovery

Builds two references from the SAME Lee/SMC cells and the SAME genes --
``mean(CP10K)`` and ``expm1(mean(log1p(CP10K)))`` -- then deconvolves pseudobulk
whose compartment fractions we set, and reports the gap between them.

    cross   SMC reference -> KUL3 pseudobulk. The leg that matters: cross-cohort,
            which is TCGA's situation.
    within  SMC reference -> held-out SMC pseudobulk. The confound check. If
            this recovers well and `cross` collapses, the loss is BATCH rather
            than recipe. Batch hits both recipes roughly equally so the gap
            should survive it, but reading a gap without knowing which regime
            you are in is how a batch effect gets reported as a recipe effect.

A DECISION EXPERIMENT, NOT A RESULT. It says nothing about colorectal cancer.
Its output decides whether the W1 linear rebuild is worth doing, and nothing
else. The absolute correlations here will beat TCGA's, because pseudobulk has no
ambient contamination and no library preparation between the cells and the
mixture. **The gap transfers; the absolutes do not.**

Why pseudobulk truth rather than the Stage 4 gate's own: the gate correlates the
deconvolved non-epithelial fraction against (1 - ABSOLUTE purity), and purity is
the malignant-cell share rather than the epithelial share, so the gate carries a
definitional ceiling nobody has measured. Fractions we set remove that confound
entirely -- any shortfall is the instrument.
"""

from __future__ import annotations

import argparse
import logging
import sys
from collections.abc import Sequence
from pathlib import Path

import pandas as pd

from src.common.io import write_versioned_table
from src.common.paths import RESULTS_DIR
from src.common.provenance import DEFAULT_SEED
from src.harness.pseudobulk import patient_holdout
from src.harness.pseudobulk_recovery import (
    RECIPES,
    RecoveryError,
    build_reference,
    make_pseudobulk,
    recipe_gap,
    recover,
)

log = logging.getLogger(__name__)

DEFAULT_S_MATRIX = RESULTS_DIR / "2026-08-26_63ead2e" / "S_matrix_lineage_1.0.0.parquet"

#: The Stage 4 gate missed by this much: 0.479 against a 0.5 threshold. A recipe
#: gap at or above it makes the linear rebuild the explanation.
GATE_SHORTFALL = 0.021


def signature_symbols(s_matrix: Path) -> list[str]:
    """The committed signature's genes, as SYMBOLS.

    The S matrix is keyed by Ensembl id and Lee's GEO matrix by symbol, so the
    two intersect in ZERO genes until this mapping happens. That is a silently
    empty deconvolution rather than an error, and it is the same identifier-space
    class the WMHS appendix records twice. An ambiguous symbol is dropped rather
    than tie-broken, matching `resolve_symbols`' own refusal to guess.
    """
    from src.bulk.gene_index import load_gene_index_map

    genes = pd.read_parquet(s_matrix).set_index("gene").index
    index_map = load_gene_index_map("1.0.0").set_index("ensembl_id")
    info = index_map.reindex(genes)
    # Read as a string, matching `matched_null_genes`' handling of `on_panel`.
    # `.fillna(False)` on an all-NaN object column downcasts, which pandas
    # warns about and which turns a missing gene into a False rather than a
    # missing one.
    ambiguous = info["symbol_ambiguous"].astype(str).str.lower().eq("true")
    symbols = info.loc[~ambiguous, "gene_symbol"].dropna()
    if ambiguous.any():
        log.info("dropped %d ambiguous symbol(s)", int(ambiguous.sum()))
    if symbols.empty:
        raise RecoveryError(
            f"{s_matrix.name}'s genes resolved to no symbols at all. The matrix "
            f"is Ensembl-keyed and Lee is symbol-keyed; without this mapping the "
            f"deconvolution runs on an empty gene set and reports nothing."
        )
    return sorted(set(symbols))


def load(cohort: str, genes: list[str]):
    from src.estimator.lee_io import load_lee_cohort

    log.info("loading Lee/%s ...", cohort.upper())
    return load_lee_cohort(cohort, target_genes=["GUCA2A"], extra_genes=genes)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--s-matrix", type=Path, default=DEFAULT_S_MATRIX)
    parser.add_argument("--n-pseudobulk", type=int, default=200)
    parser.add_argument("--n-cells", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    try:
        genes = signature_symbols(args.s_matrix)
    except RecoveryError as exc:
        log.error("%s", exc)
        return 2
    log.info("signature: %d genes, resolved to symbols", len(genes))

    smc = load("smc", genes)
    kul3 = load("kul3", genes)
    usable = [g for g in genes
              if g in smc.expression.columns and g in kul3.expression.columns]
    log.info("genes present in both cohorts: %d of %d (%.1f%%)",
             len(usable), len(genes), 100 * len(usable) / len(genes))
    if len(usable) < 200:
        log.error("only %d shared genes. nu-SVR's robustness is a "
                  "high-dimensionality property; this is too few.", len(usable))
        return 2

    smc_type = smc.cells["author_cell_type"]
    kul3_type = kul3.cells["author_cell_type"]
    log.info("SMC compartments:  %s", smc_type.value_counts().to_dict())
    log.info("KUL3 compartments: %s", kul3_type.value_counts().to_dict())

    # Split by PATIENT, using the harness's own function rather than a second
    # implementation of it. Cells within a patient are not independent draws,
    # and a reference built from the cells its own pseudobulk is drawn from
    # measures memorisation.
    patients = sorted(smc.cells["patient_id"].unique())
    train, held = patient_holdout(
        patients, n_held_out=max(1, len(patients) // 2), seed=args.seed
    )
    log.info("SMC patients: %d reference %s | %d held out %s",
             len(train), list(train), len(held), list(held))
    in_reference = smc.cells["patient_id"].isin(train).to_numpy()

    ref_expression = smc.expression.loc[in_reference, usable]
    ref_type = smc_type[in_reference]

    results = []
    legs = (
        ("cross", kul3.expression[usable], kul3_type),
        ("within", smc.expression.loc[~in_reference, usable], smc_type[~in_reference]),
    )
    for leg, expression, cell_type in legs:
        try:
            bulk, truth = make_pseudobulk(
                expression, cell_type, n_samples=args.n_pseudobulk,
                n_cells=args.n_cells, seed=args.seed,
            )
        except RecoveryError as exc:
            log.error("  %s leg: %s", leg, exc)
            continue
        log.info("\n%s leg: %d pseudobulk samples of %d cells | compartments %s",
                 leg, len(bulk), args.n_cells, list(truth.columns))
        for recipe in RECIPES:
            reference = build_reference(ref_expression, ref_type, recipe=recipe)
            got = recover(bulk, truth, reference, leg=leg, recipe=recipe)
            for r in got:
                log.info("  %-6s %-6s %-6s  r(non-epi)=%6.3f  r(epi)=%6.3f  "
                         "epi-zeros=%d/%d", leg, recipe, r.method,
                         r.r_non_epithelial, r.r_epithelial,
                         r.epithelial_exact_zero, r.n_samples)
            results.extend(got)

    if not results:
        log.error("no leg produced a result")
        return 3

    table = pd.DataFrame([r.as_row() for r in results])
    gaps = recipe_gap(results)
    log.info("\n%s", "=" * 72)
    log.info("THE NUMBER THIS EXPERIMENT IS FOR: linear minus log1p")
    log.info("%s", "=" * 72)
    log.info("%s", gaps.to_string(index=False))

    cross = gaps[gaps["leg"] == "cross"]
    if not cross.empty and "r_non_epithelial_gap" in cross.columns:
        best = float(cross["r_non_epithelial_gap"].max())
        verdict = (
            "AT OR ABOVE the shortfall: the recipe is a live explanation for the "
            "gate failure, and the W1 linear rebuild is justified."
            if best >= GATE_SHORTFALL else
            "BELOW the shortfall: the recipe is not what is costing the margin. "
            "The gate failure is a bulk-deconvolution limit and the rebuild will "
            "not reach it."
        )
        log.info(
            "\nLargest cross-cohort gap on the gate's own quantity: %+.4f\n"
            "The Stage 4 gate missed by %.3f (0.479 against 0.5).\n  -> %s\n\n"
            "This is a statement about the INSTRUMENT. It is not a result about "
            "the biology,\nand the absolute correlations here do not transfer to "
            "real bulk.",
            best, GATE_SHORTFALL, verdict,
        )

    for frame, name in ((table, "pseudobulk_recovery"),
                        (gaps, "pseudobulk_recovery_gap")):
        path = write_versioned_table(
            frame, name, seed=args.seed,
            results_dir=args.results_dir, allow_dirty=args.allow_dirty,
            extra_meta={
                "legs": {
                    "cross": "SMC reference -> KUL3 pseudobulk (cross-cohort, like TCGA)",
                    "within": "SMC reference -> held-out SMC pseudobulk (batch confound check)",
                },
                "recipes": {
                    "linear": "mean(CP10K), what a linear mixture actually sums",
                    "log1p": "expm1(mean(log1p(CP10K))), build_signature's recipe",
                },
                "ground_truth": (
                    "realised compartment fractions of the sampled cells, not the "
                    "requested Dirichlet draw"
                ),
                "n_genes": len(usable),
                "n_pseudobulk": args.n_pseudobulk,
                "n_cells_per_sample": args.n_cells,
                "smc_reference_patients": list(train),
                "smc_holdout_patients": list(held),
                "gate_shortfall": GATE_SHORTFALL,
                "what_this_answers": (
                    "Whether the reference RECIPE costs the instrument the 0.021 "
                    "the Stage 4 gate missed by. A decision about the "
                    "deconvolution leg's future."
                ),
                "what_this_does_not_answer": (
                    "Anything about colorectal cancer. And the absolute "
                    "correlations do not transfer to real bulk -- pseudobulk has "
                    "no ambient contamination and no library preparation. The GAP "
                    "transfers; the absolutes do not."
                ),
            },
        )
        log.info("wrote %s", path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
