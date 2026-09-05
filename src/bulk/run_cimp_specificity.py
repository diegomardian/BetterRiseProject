"""Is the CIMP association with GUCA2A locus-specific, or axis-level?

Registered in ``docs/prereg_cimp_specificity.md`` before this was written. Read
that first: it states the decision rule, the expected outcome, and what no result
here can establish. This module implements it and nothing else.

THE SHAPE OF THE TEST. GUCA2A is lower in the methylator phenotype. So are CDX2
and MS4A12. The question is whether GUCA2A falls *beyond* markers of the same
population, which is what locus-specific promoter silencing would predict, or
*with* them, which is what upstream suppression or a shift in differentiated
content would predict. Only the first would justify the 450k methylation leg.

WHY WITHIN-SAMPLE DIFFERENCES. The estimand pairs two genes inside one array
before comparing strata. A per-sample difference cancels loading and scanner
effects that hit both genes together, and it is lower variance than subtracting
two independently estimated stratum medians. The stratum contrast is then a
difference of those paired differences.

TWO REFERENCES, EXPECTED TO DISAGREE. CDX2 drives GUCA2A, so it is the
mechanistically relevant reference. MS4A12 is colonocyte-restricted with no
regulatory relationship, so it is the population reference. The committed medians
already put GUCA2A between them -- it falls less than CDX2 and more than MS4A12 --
and the pre-registration commits to reading that as *not specific*, because
support requires both references to agree. Picking the agreeable one afterwards
is the degree of freedom the two-reference design exists to remove.

    python -m src.bulk.run_cimp_specificity
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from src.common.io import write_versioned_table
from src.common.paths import PROCESSED_DIR
from src.common.provenance import DEFAULT_SEED

log = logging.getLogger(__name__)

#: Set from --allow-dirty. Default False.
ALLOW_DIRTY = False

BULK = PROCESSED_DIR / "bulk"
EXPRESSION = BULK / "gse39582_panel_expression.tsv"
METADATA = BULK / "gse39582_metadata.tsv"

TARGET = "GUCA2A"

#: The two references and what each one tests. Reported separately, always.
REFERENCES: dict[str, str] = {
    "CDX2": "upstream transcription factor; drives the target",
    "MS4A12": "colonocyte-restricted, no regulatory relationship to the target",
}

N_BOOTSTRAP = 10_000

#: Adjusters available without new compute. Purity is deliberately absent: the
#: pre-registration admits it only if both references favour the target, which
#: is the surprise case.
ADJUSTERS = ("mmr.status", "tumor.location")


def paired_differences(expression: pd.DataFrame, reference: str) -> pd.Series:
    """Per sample, ``log2 target − log2 reference``.

    Both genes come off the same array in the same hybridisation, so anything
    that scales a sample as a whole divides out here rather than being adjusted
    for later.
    """
    missing = [g for g in (TARGET, reference) if g not in expression.index]
    if missing:
        raise KeyError(f"{missing} absent from the panel expression table")
    return (expression.loc[TARGET] - expression.loc[reference]).astype(float)


def stratum_contrast(
    differences: pd.Series, positive: np.ndarray, negative: np.ndarray,
    *, seed: int = DEFAULT_SEED, n_boot: int = N_BOOTSTRAP,
) -> dict:
    """``mean[D | CIMP+] − mean[D | CIMP−]``, stratified bootstrap over samples.

    The two strata are resampled separately, so an interval reflects the
    91-and-405 design rather than a pooled draw that could return an implausible
    split. One sample per patient in this cohort, so resampling samples is
    resampling patients and invariant 5 is satisfied rather than sidestepped.

    Sign convention: NEGATIVE means the target fell further in CIMP+ than the
    reference did, which is the direction locus-specific silencing predicts.
    """
    a = differences.to_numpy(dtype=float)[positive]
    b = differences.to_numpy(dtype=float)[negative]
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if len(a) < 3 or len(b) < 3:
        return {"n_positive": len(a), "n_negative": len(b), "contrast": np.nan,
                "ci_low": np.nan, "ci_high": np.nan, "excludes_zero": False}

    rng = np.random.default_rng(seed)
    draws = (rng.choice(a, size=(n_boot, len(a)), replace=True).mean(axis=1)
             - rng.choice(b, size=(n_boot, len(b)), replace=True).mean(axis=1))
    lo, hi = (float(x) for x in np.percentile(draws, [2.5, 97.5]))
    return {"n_positive": int(len(a)), "n_negative": int(len(b)),
            "contrast": float(a.mean() - b.mean()),
            "ci_low": lo, "ci_high": hi, "excludes_zero": bool(lo * hi > 0)}


def verdict(rows: pd.DataFrame) -> tuple[bool, str]:
    """The pre-registered decision rule, applied without discretion.

    Support requires BOTH references to favour the target with intervals
    excluding zero. One reference agreeing is the configuration the committed
    medians already show, and calling that support would be choosing the
    reference after seeing the answer.
    """
    unadjusted = rows.loc[rows["adjustment"] == "none"]
    if len(unadjusted) < len(REFERENCES):
        return False, "UNDEFINED: not every reference was estimable"
    favourable = unadjusted.loc[
        (unadjusted["contrast"] < 0) & unadjusted["excludes_zero"]
    ]
    if len(favourable) == len(REFERENCES):
        return True, (
            f"SPECIFIC: {TARGET} falls further than both references with "
            f"intervals excluding zero. This is the pre-registered surprise and "
            f"it triggers purity adjustment and a full re-read before any claim."
        )
    named = ", ".join(
        f"{r.reference} {r.contrast:+.3f} [{r.ci_low:+.3f}, {r.ci_high:+.3f}]"
        for r in unadjusted.itertuples()
    )
    return False, (
        f"NOT SPECIFIC: support needs both references and {len(favourable)} of "
        f"{len(REFERENCES)} favour {TARGET} -- {named}. The methylator "
        f"association reads as axis-level or compositional, not locus-specific."
    )


def _adjusted_masks(metadata: pd.DataFrame, positive: np.ndarray,
                    negative: np.ndarray, column: str):
    """Within-level CIMP masks, so a level with only one arm drops out.

    Categorical adjustment by stratification rather than regression: the levels
    are few and the estimand stays a difference of means, which keeps the
    bootstrap the same object it is elsewhere in this file.
    """
    if column not in metadata:
        return
    values = metadata[column].to_numpy()
    for level in sorted({v for v in values if isinstance(v, str) and v}):
        at = values == level
        if (positive & at).sum() >= 3 and (negative & at).sum() >= 3:
            yield level, positive & at, negative & at


def run(expression: pd.DataFrame, metadata: pd.DataFrame,
        *, seed: int = DEFAULT_SEED) -> pd.DataFrame:
    """One row per (reference, adjustment level). Unadjusted rows carry 'none'."""
    from src.bulk.replication import strata

    masks = strata(metadata)
    for needed in ("tumour|CIMP+", "tumour|CIMP-"):
        if needed not in masks:
            raise KeyError(f"{needed} absent — this cohort's CIMP annotation is why it is used")
    positive, negative = masks["tumour|CIMP+"], masks["tumour|CIMP-"]

    rows = []
    for reference, role in REFERENCES.items():
        differences = paired_differences(expression, reference)
        base = {"target": TARGET, "reference": reference, "reference_role": role}
        rows.append(base | {"adjustment": "none", "level": "all"}
                    | stratum_contrast(differences, positive, negative, seed=seed))
        for column in ADJUSTERS:
            for level, pos, neg in _adjusted_masks(metadata, positive, negative, column):
                rows.append(base | {"adjustment": column, "level": level}
                            | stratum_contrast(differences, pos, neg, seed=seed))
    return pd.DataFrame(rows)


def main(argv=None) -> int:
    global ALLOW_DIRTY
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument(
        "--allow-dirty", action="store_true",
        help="write from a dirty tree; the recorded sha will not reproduce it",
    )
    args = parser.parse_args(argv)
    ALLOW_DIRTY = args.allow_dirty
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    for path in (EXPRESSION, METADATA):
        if not path.exists():
            log.error("%s missing — run `python -m src.bulk.run_replication run` first", path)
            return 1
    expression = pd.read_csv(EXPRESSION, sep="\t", index_col=0)
    metadata = pd.read_csv(METADATA, sep="\t", index_col=0)
    log.info("%d genes x %d samples", *expression.shape)

    table = run(expression, metadata, seed=args.seed)
    supported, reading = verdict(table)
    for r in table.loc[table["adjustment"] == "none"].itertuples():
        log.info("  %-8s vs %-8s  %+0.3f  [%+0.3f, %+0.3f]  n=%d/%d",
                 TARGET, r.reference, r.contrast, r.ci_low, r.ci_high,
                 r.n_positive, r.n_negative)
    log.info("%s", reading)

    path = write_versioned_table(
        table, "gse39582_cimp_specificity", seed=args.seed,
        results_dir=args.results_dir, allow_dirty=ALLOW_DIRTY,
        extra_meta={
            "preregistered_in": "docs/prereg_cimp_specificity.md",
            "question": "is the CIMP association with GUCA2A locus-specific or axis-level",
            "estimand": "mean[log2 target - log2 reference | CIMP+] - same | CIMP-",
            "decision_rule": "support requires BOTH references, intervals excluding zero",
            "expected_before_running": "not specific; references disagree in direction",
            "verdict": reading,
            "supported": supported,
            "cannot_establish": (
                "Bulk is fraction x per-cell mean, so no outcome here separates "
                "silencing from colonocyte-specific compositional loss. A null "
                "closes a route; a positive only promotes the question to a leg "
                "with the same structure."
            ),
            "n_bootstrap": N_BOOTSTRAP,
        },
    )
    log.info("wrote %s (%d rows)", path, len(table))
    return 0


if __name__ == "__main__":
    sys.exit(main())
