"""Gate 6: who the labels miss, and whether the missed patients differ.

    python -m src.reference.jobs.chen_subtype_attrition --no-write

Gate condition 6 of `docs/prereg_chen_lesion_subtype.md` §9. Family P covers 23
of the 44 avenue-A patients. The 21 it misses are not missing at random in any
way anyone has checked, and D2 §6a is the standing reason to look: there, the
dropped participants' GUCA2A ran *lower* than the retained set's, in the
direction of the outcome, and nobody would have known without comparing.

**This job computes no arm contrast.** Coverage is orthogonal to `AD` versus
`SER`: the question here is whether the 23 patients who carry any label differ
from the 21 who carry none. Running the contrast before this comparison has
been inspected is what the condition's ordering exists to prevent.

**`None` is not `0.0`.** A `not_estimable` intrinsic term is excluded from every
mean and counted beside it (invariant 1). A group whose values are missing is
not a group whose values are zero.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from src.common.io import write_versioned_table
from src.common.label_provenance import check_no_circular_claim, provenance_meta
from src.common.provenance import DEFAULT_SEED
from src.reference.chen_subtype_labels import CLAIM, FAMILY_P
from src.reference.jobs.chen_subtype_id_provenance import (
    ARMS,
    DEFAULT_CROSSWALK,
    SNAPSHOT,
    exact_join,
    read_labels,
)

DECOMPOSITION = Path("results/2026-09-06_5791c01/adenoma_decomposition.parquet")
RUNG = "lineage"
TERMS = ("compositional", "intrinsic", "interaction")

LABEL_PROVENANCE = FAMILY_P
CLAIM_PROVENANCE = CLAIM


class AttritionError(RuntimeError):
    """The attrition comparison cannot be established as the pre-registration fixes it."""


def coverage_by_patient(joined: pd.DataFrame) -> pd.DataFrame:
    """One row per avenue-A patient: does family P reach them at all?

    Coverage is "carries any `AD` or `SER` lesion", which includes the
    conflicted patient excluded from the arms. That patient is covered by the
    labels and excluded by a rule; folding them into the uncovered group would
    describe the label's reach wrongly.
    """
    labelled = joined[joined["label_matched"] & joined["POLYP_TYPE"].isin(ARMS)]
    covered = set(labelled["patient_id"])
    rows = []
    for patient, block in joined.groupby("patient_id", sort=True):
        types = sorted({str(v) for v in block["POLYP_TYPE"].dropna()})
        rows.append(
            {
                "patient_id": patient,
                "covered": patient in covered,
                "n_specimens": len(block),
                "n_specimens_with_a_cbioportal_record": int(block["label_matched"].sum()),
                "polyp_type_values": ";".join(types) if types else "",
                "uncovered_reason": (
                    ""
                    if patient in covered
                    else ("label is Unknown" if types else "no cBioPortal record")
                ),
            }
        )
    return pd.DataFrame(rows)


def compare(decomposition: pd.DataFrame, coverage: pd.DataFrame) -> pd.DataFrame:
    """Covered versus uncovered, per gene, term and weighting. Descriptive."""
    values = decomposition[decomposition["granularity_rung"] == RUNG].merge(
        coverage[["patient_id", "covered"]], on="patient_id", how="inner"
    )
    if values.empty:
        raise AttritionError(f"no {RUNG} rows joined to a coverage flag")

    rows = []
    long = values.melt(
        id_vars=["patient_id", "gene", "weighting", "covered", "estimability"],
        value_vars=list(TERMS), var_name="term", value_name="value",
    )
    for (gene, term, weighting), block in long.groupby(["gene", "term", "weighting"], sort=True):
        row: dict[str, object] = {"gene": gene, "term": term, "weighting": weighting}
        stats = {}
        for label, group in (("covered", block[block["covered"]]),
                             ("uncovered", block[~block["covered"]])):
            # Invariant 1: a missing value is dropped and counted, never zeroed.
            usable = group["value"].dropna()
            row[f"n_{label}"] = int(len(group))
            row[f"n_{label}_not_estimable"] = int(group["value"].isna().sum())
            row[f"mean_{label}"] = float(usable.mean()) if len(usable) else None
            row[f"median_{label}"] = float(usable.median()) if len(usable) else None
            stats[label] = usable
        a, b = stats["covered"], stats["uncovered"]
        row["mean_difference"] = (
            float(a.mean() - b.mean()) if len(a) and len(b) else None
        )
        # A descriptive magnitude, on a pooled SD. Not a test, and no threshold
        # is attached to it: none was pre-committed, and inventing one after
        # seeing the table is the choice this project refuses elsewhere.
        if len(a) > 1 and len(b) > 1:
            pooled = np.sqrt(
                ((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1))
                / (len(a) + len(b) - 2)
            )
            row["standardised_difference"] = (
                float((a.mean() - b.mean()) / pooled) if pooled > 0 else None
            )
        else:
            row["standardised_difference"] = None
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["term", "gene", "weighting"]).reset_index(drop=True)


def summarise(coverage: pd.DataFrame, comparison: pd.DataFrame) -> dict[str, object]:
    # to_numeric first: an all-None column is object dtype, where abs() raises.
    magnitude = pd.to_numeric(
        comparison["standardised_difference"], errors="coerce"
    ).abs()
    top = (
        comparison.assign(_magnitude=magnitude)
        .dropna(subset=["_magnitude"])
        .sort_values("_magnitude", ascending=False)
        .head(1)
    )
    return {
        "n_patients": int(len(coverage)),
        "n_covered": int(coverage["covered"].sum()),
        "n_uncovered": int((~coverage["covered"]).sum()),
        "n_uncovered_label_unknown": int(
            (coverage["uncovered_reason"] == "label is Unknown").sum()
        ),
        "n_uncovered_no_record": int(
            (coverage["uncovered_reason"] == "no cBioPortal record").sum()
        ),
        # Invariant 1, reported before the arms are ever split: if every
        # unestimable intrinsic term sits on one side of the coverage line,
        # the two groups are not the same size in the term that matters.
        "n_intrinsic_not_estimable_covered": int(
            comparison.loc[comparison["term"] == "intrinsic", "n_covered_not_estimable"].sum()
        ),
        "n_intrinsic_not_estimable_uncovered": int(
            comparison.loc[comparison["term"] == "intrinsic", "n_uncovered_not_estimable"].sum()
        ),
        "largest_abs_standardised_difference": (
            float(top["_magnitude"].iloc[0]) if len(top) else None
        ),
        "largest_difference_at": (
            f"{top['gene'].iloc[0]}/{top['term'].iloc[0]}/{top['weighting'].iloc[0]}"
            if len(top) else ""
        ),
        "systematic_difference_rule": "none pre-committed; inspected by hand",
        "arm_contrast_computed": False,
        "verdict": "ATTRITION TABLE PRODUCED — inspect before any arm contrast",
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--decomposition", type=Path, default=DECOMPOSITION)
    parser.add_argument("--crosswalk", type=Path, default=DEFAULT_CROSSWALK)
    parser.add_argument("--snapshot", type=Path, default=SNAPSHOT)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--no-write", action="store_true")
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)

    check_no_circular_claim(labels=LABEL_PROVENANCE, claim=CLAIM_PROVENANCE)

    joined = exact_join(read_labels(args.snapshot), pd.read_parquet(args.crosswalk))
    coverage = coverage_by_patient(joined)
    decomposition = pd.read_parquet(args.decomposition)
    comparison = compare(decomposition, coverage)
    outcome = summarise(coverage, comparison)

    print(coverage["uncovered_reason"].value_counts(dropna=False).to_string())
    print()
    intrinsic = comparison[
        (comparison["term"] == "intrinsic") & (comparison["weighting"] == "doubly_robust")
    ]
    print("intrinsic, doubly_robust — the term §5 makes primary:")
    print(intrinsic[
        ["gene", "n_covered", "n_uncovered", "mean_covered", "mean_uncovered",
         "mean_difference", "standardised_difference"]
    ].to_string(index=False))
    print()
    for key, value in outcome.items():
        print(f"  {key:42} {value}")
    print("\nNO ARM CONTRAST WAS COMPUTED — coverage only.")

    if args.no_write:
        return 0
    meta = {
        **provenance_meta(LABEL_PROVENANCE, CLAIM_PROVENANCE),
        "prereg": "docs/prereg_chen_lesion_subtype.md",
        "gate_condition": 6,
        "decomposition": str(args.decomposition),
        "granularity_rung": RUNG,
        **outcome,
    }
    for table, name in (
        (coverage, "chen_subtype_coverage"),
        (comparison, "chen_subtype_attrition_comparison"),
    ):
        path = write_versioned_table(
            table, name, seed=args.seed, results_dir=args.results_dir,
            allow_dirty=args.allow_dirty, extra_meta=meta,
        )
        print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
