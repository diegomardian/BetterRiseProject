"""Stage 4, step 1: bulk -> mature-colonocyte fraction, per rung, per method.

    python -m src.bulk.run_deconvolution --rung lineage
    python -m src.bulk.run_deconvolution --rung all --linearise-reference

Reads TCGA TPM (linear, samples x genes, gene index 1.0.0) and the committed
``S_matrix_{rung}_1.0.0.parquet``, and writes one row per (sample, method, rung)
with the two aggregates Stage 4 needs: the predictor
(``mature_colonocyte_fraction``) and the instrument gate's quantity
(``non_epithelial_fraction``).

WHY ``--linearise-reference`` IS NOT THE DEFAULT. The committed matrices are on
the mean-of-log1p scale and the bulk is linear, so the run is refused outright
unless you say which repair you are applying. See
``src.bulk.deconvolution.linearise``: it is an approximation, biased low by
Jensen, and the right fix is a linearly-built reference from W1. Making it the
default would hide a known misspecification behind a flag nobody reads.

The refusal is not the instrument gate doing its job. The gate cannot see this
failure at all -- it reads the non-epithelial aggregate, which survives the
mismatch that zeroes the predictor. That asymmetry is the reason this driver
checks the predictor itself before writing anything.
"""

from __future__ import annotations

import argparse
import logging
import sys
from collections.abc import Sequence
from pathlib import Path

import pandas as pd

from src.bulk.deconvolution import (
    LINEAR_CP10K,
    LOG1P_CP10K,
    DeconvolutionError,
    check_predictor,
    deconvolve_cohort,
    default_methods,
    linearise,
    load_reference,
    mature_column,
    require_usable_predictor,
)
from src.common.io import write_versioned_table
from src.common.paths import PROCESSED_DIR, RESULTS_DIR
from src.common.provenance import DEFAULT_SEED

log = logging.getLogger(__name__)

BULK = PROCESSED_DIR / "bulk"
DEFAULT_S_MATRIX_DIR = RESULTS_DIR / "2026-08-26_63ead2e"
RUNGS: tuple[str, ...] = ("epithelial", "lineage", "crypt_position", "best4")

#: Panel targets, which must never appear in the reference (invariant 2).
def _panel_targets() -> list[str]:
    from src.common.panel import panel_genes

    return list(panel_genes())


def load_bulk(path: Path) -> pd.DataFrame:
    """TCGA TPM, samples x genes, on the shared index.

    Asserted linear on the way in rather than assumed. The whole point of this
    module is that a scale nobody checked is how the predictor becomes a
    constant, so reading a matrix without asking what scale it is on would be
    the same mistake one level up.
    """
    from src.bulk.normalise import assert_linear_scale

    if not path.exists():
        raise SystemExit(
            f"{path} not found.\n"
            f"The TCGA matrices are gitignored and live where the data is -- on "
            f"the cluster, under $BRP_DATA_DIR/processed/bulk. Build them with\n"
            f"    python -m src.bulk.ingest build\n"
            f"or point --bulk at an existing tcga_tpm_1.0.0.parquet."
        )
    bulk = pd.read_parquet(path)
    assert_linear_scale(bulk, context=f"{path.name} (the deconvolution input)")
    return bulk


def run_one_rung(
    bulk: pd.DataFrame,
    rung: str,
    *,
    s_matrix_dir: Path,
    linearise_reference: bool,
) -> tuple[pd.DataFrame, list, dict[str, str], str]:
    """Fractions and predictor checks for one rung."""
    from src.bulk.deconvolution import summarise_fractions

    path = s_matrix_dir / f"S_matrix_{rung}_1.0.0.parquet"
    reference = load_reference(path, rung=rung, scale=LOG1P_CP10K, targets=_panel_targets())
    note = "committed reference used as-is"
    if linearise_reference:
        reference = linearise(reference)
        note = (
            "reference linearised with expm1 from the committed mean-of-log1p "
            "matrix. APPROXIMATE: a geometric rather than arithmetic mean of "
            "CP10K+1, biased low by Jensen, worst for the most dispersed genes. "
            "The correct fix is a linearly-built reference from W1."
        )
    log.info("  %s", reference.describe())

    long, skipped = deconvolve_cohort(
        bulk, reference, bulk_scale=LINEAR_CP10K, methods=default_methods()
    )
    summary = summarise_fractions(long, rung)
    checks = [
        check_predictor(summary, rung=rung, method=method)
        for method in sorted(summary["method"].unique())
    ]
    for check in checks:
        level = log.info if check.usable else log.warning
        level("  %-6s %-15s %s", check.method, check.verdict, check.detail)
    return summary, checks, skipped, note


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rung", default="lineage", choices=(*RUNGS, "all"))
    parser.add_argument("--bulk", type=Path, default=BULK / "tcga_tpm_1.0.0.parquet")
    parser.add_argument("--s-matrix-dir", type=Path, default=DEFAULT_S_MATRIX_DIR)
    parser.add_argument(
        "--linearise-reference", action="store_true",
        help="expm1 the committed log reference so it can meet linear bulk. "
             "Required, and deliberately not the default -- read "
             "src.bulk.deconvolution.linearise before using it.",
    )
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    bulk = load_bulk(args.bulk)
    log.info("bulk: %d samples x %d genes from %s", *bulk.shape, args.bulk.name)

    rungs = RUNGS if args.rung == "all" else (args.rung,)
    frames, all_checks, notes = [], [], {}
    skipped_any: dict[str, str] = {}
    for rung in rungs:
        log.info("%s rung (mature bin: %s)", rung, mature_column(rung))
        try:
            summary, checks, skipped, note = run_one_rung(
                bulk, rung,
                s_matrix_dir=args.s_matrix_dir,
                linearise_reference=args.linearise_reference,
            )
        except DeconvolutionError as exc:
            log.error("  REFUSED: %s", exc)
            return 2
        frames.append(summary)
        all_checks.extend(checks)
        skipped_any.update(skipped)
        notes[rung] = note

    table = pd.concat(frames, ignore_index=True)
    table.insert(0, "cohort", "TCGA-COAD/READ")

    checks_frame = pd.DataFrame([
        {
            "granularity_rung": c.rung, "method": c.method, "n_samples": c.n_samples,
            "n_exact_zero": c.n_exact_zero, "fraction_sd": c.sd,
            "is_constant": c.is_constant, "verdict": c.verdict, "detail": c.detail,
        }
        for c in all_checks
    ])

    # Write the tables BEFORE the refusal, so a failed run leaves the evidence
    # of its failure behind rather than only a traceback. A refused predictor is
    # a Stage 4 result and it needs a stamped table like any other.
    for frame, name in ((table, "stage4_fractions"), (checks_frame, "stage4_predictor_checks")):
        path = write_versioned_table(
            frame, name, seed=args.seed,
            results_dir=args.results_dir, allow_dirty=args.allow_dirty,
            extra_meta={
                "cohort": "TCGA-COAD/READ",
                "rungs": list(rungs),
                "reference_scale": LINEAR_CP10K if args.linearise_reference else LOG1P_CP10K,
                "reference_note": notes,
                "bulk_source": str(args.bulk.name),
                "bulk_scale": LINEAR_CP10K,
                "methods_skipped": skipped_any,
                "what_this_answers": (
                    "Per (sample, method, rung): the mature-colonocyte fraction "
                    "Stage 4 regresses on, and the non-epithelial fraction its "
                    "instrument gate checks against ABSOLUTE purity. Emitted "
                    "side by side because the second can be right while the "
                    "first is a constant -- see src/bulk/deconvolution.py."
                ),
                "invariant_1": (
                    "A rung with no maturity call carries mature_colonocyte_"
                    "fraction = None with estimability='not_estimable', never 0.0."
                ),
            },
        )
        log.info("wrote %s", path)

    try:
        require_usable_predictor(all_checks)
    except DeconvolutionError as exc:
        log.error("\n%s", exc)
        log.error(
            "\nSTOP. No R-squared may be computed from this run. The tables above "
            "record the refusal and are the Stage 4 result for it."
        )
        return 3

    usable = [c for c in all_checks if c.usable]
    log.info("\n%d of %d (rung, method) combinations produced a usable predictor",
             len(usable), len(all_checks))
    return 0


if __name__ == "__main__":
    sys.exit(main())
