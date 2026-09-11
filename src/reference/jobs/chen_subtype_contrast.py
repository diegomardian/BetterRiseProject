"""§5: the intrinsic term between family P arms, AD minus SER.

    python -m src.reference.jobs.chen_subtype_contrast --no-write

The locked design is `docs/prereg_chen_lesion_subtype.md`. This is the first
step in it that reads an outcome, and it is the only irreversible one, so every
choice it makes was fixed before it ran:

* the arms and their sizes (§3, verified by gate 4);
* the floor of 8 estimable patients per arm (§7), which decides whether
  anything is reported at all;
* Welch's t for the between-arm difference, and all three weightings with none
  primary (§5's two under-determinations, fixed at lock);
* the compositional term is **descriptive only** (§4) and carries no test.

`None` is not `0.0`. A `not_estimable` patient is dropped from the arm and
counted beside every estimate (invariant 1).
"""

from __future__ import annotations

import argparse
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from src.common.io import write_versioned_table
from src.common.label_provenance import check_no_circular_claim, provenance_meta
from src.common.provenance import DEFAULT_SEED
from src.reference.chen_subtype_labels import CLAIM, FAMILY_P
from src.reference.interval_calibration import NOMINAL_ALPHA
from src.reference.jobs.chen_subtype_attrition import DECOMPOSITION, RUNG
from src.reference.jobs.chen_subtype_id_provenance import (
    ARMS,
    DEFAULT_CROSSWALK,
    SNAPSHOT,
    classify_patients,
    exact_join,
    read_labels,
)

#: §7, last row. Fixed before the estimability counts were seen.
MIN_ESTIMABLE_PER_ARM = 8
TERMS = ("intrinsic", "compositional", "interaction")
#: §4. Only the intrinsic term is a test; the others are reported, not read.
TESTED_TERM = "intrinsic"

LABEL_PROVENANCE = FAMILY_P
CLAIM_PROVENANCE = CLAIM


class ContrastError(RuntimeError):
    """The contrast cannot be computed as the locked design fixes it."""


def welch_difference(a: np.ndarray, b: np.ndarray, *, alpha: float = NOMINAL_ALPHA) -> dict:
    """AD minus SER, with a Welch t interval. Unequal variances by construction."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.size < 2 or b.size < 2:
        raise ContrastError("a Welch interval needs at least 2 observations per arm")
    va, vb = a.var(ddof=1) / a.size, b.var(ddof=1) / b.size
    se = float(np.sqrt(va + vb))
    if se == 0.0:
        raise ContrastError("zero standard error; the arms carry no variation")
    df = (va + vb) ** 2 / (va**2 / (a.size - 1) + vb**2 / (b.size - 1))
    crit = float(stats.t.ppf(1 - alpha / 2, df))
    diff = float(a.mean() - b.mean())
    return {
        "mean_difference": diff,
        "ci_low": diff - crit * se,
        "ci_high": diff + crit * se,
        "welch_df": float(df),
        "interval_method": "welch_t",
    }


def arm_values(values: pd.DataFrame, arm: str, gene: str, term: str, weighting: str) -> pd.Series:
    block = values[
        (values["arm"] == arm) & (values["gene"] == gene) & (values["weighting"] == weighting)
    ]
    return block[term]


def per_gene(values: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for gene in sorted(values["gene"].unique()):
        for term in TERMS:
            for weighting in sorted(values["weighting"].unique()):
                row: dict[str, object] = {
                    "gene": gene, "term": term, "weighting": weighting,
                    "tested": term == TESTED_TERM,
                    "standing": "test" if term == TESTED_TERM else "descriptive_only",
                }
                arms = {}
                for arm in ARMS:
                    series = arm_values(values, arm, gene, term, weighting)
                    usable = series.dropna()
                    arms[arm] = np.asarray(usable, dtype=float)
                    row[f"n_{arm}"] = int(len(series))
                    row[f"n_{arm}_estimable"] = int(len(usable))
                    row[f"n_{arm}_not_estimable"] = int(series.isna().sum())
                    row[f"mean_{arm}"] = float(usable.mean()) if len(usable) else None
                below = [a for a in ARMS if len(arms[a]) < MIN_ESTIMABLE_PER_ARM]
                if below:
                    row.update({
                        "mean_difference": None, "ci_low": None, "ci_high": None,
                        "welch_df": None, "interval_method": "none",
                        "excludes_zero": None,
                        "estimability": "not_estimable",
                        "reason": (
                            f"{', '.join(below)} below the §7 floor of "
                            f"{MIN_ESTIMABLE_PER_ARM} estimable patients"
                        ),
                    })
                else:
                    out = welch_difference(arms["AD"], arms["SER"])
                    row.update(out)
                    row["excludes_zero"] = bool(
                        out["ci_low"] > 0.0 or out["ci_high"] < 0.0
                    )
                    row["estimability"] = "estimated"
                    row["reason"] = ""
                rows.append(row)
    return pd.DataFrame(rows)


def per_pair(values: pd.DataFrame) -> pd.DataFrame:
    """§5: every pair, not the targets against one housekeeping comparator.

    The within-patient gene difference first, then the arm difference of it —
    so the pair contrast is paired where the data are paired.
    """
    rows = []
    genes = sorted(values["gene"].unique())
    for weighting in sorted(values["weighting"].unique()):
        wide = values[values["weighting"] == weighting].pivot_table(
            index=["patient_id", "arm"], columns="gene", values=TESTED_TERM
        )
        for left, right in combinations(genes, 2):
            if left not in wide or right not in wide:
                continue
            delta = (wide[left] - wide[right]).dropna()
            arms = {
                arm: np.asarray(
                    delta[delta.index.get_level_values("arm") == arm], dtype=float
                )
                for arm in ARMS
            }
            row: dict[str, object] = {
                "gene": left, "other": right, "term": TESTED_TERM,
                "weighting": weighting, "contrast": f"{left} - {right}",
                "n_AD": int(len(arms["AD"])), "n_SER": int(len(arms["SER"])),
            }
            if any(len(arms[a]) < MIN_ESTIMABLE_PER_ARM for a in ARMS):
                row.update({"mean_difference": None, "ci_low": None, "ci_high": None,
                            "excludes_zero": None, "estimability": "not_estimable"})
            else:
                out = welch_difference(arms["AD"], arms["SER"])
                row.update({k: out[k] for k in ("mean_difference", "ci_low", "ci_high")})
                row["excludes_zero"] = bool(out["ci_low"] > 0.0 or out["ci_high"] < 0.0)
                row["estimability"] = "estimated"
            rows.append(row)
    return pd.DataFrame(rows)


def falsifier_branch(gene_table: pd.DataFrame, weighting: str) -> dict[str, object]:
    """§7, evaluated on one weighting. The branches are pre-committed."""
    block = gene_table[
        (gene_table["term"] == TESTED_TERM) & (gene_table["weighting"] == weighting)
    ].set_index("gene")
    if block["estimability"].eq("not_estimable").any():
        return {"weighting": weighting, "branch": "NOT ESTIMABLE",
                "consequence": "below the §7 floor; no interval reported on any statistic"}
    cdx2 = bool(block.loc["CDX2", "excludes_zero"])
    guca2a = bool(block.loc["GUCA2A", "excludes_zero"])
    if cdx2 and not guca2a:
        branch, consequence = (
            "PREDICTION REPLICATES",
            "tier-level within-study heterogeneity in one cohort; not a mechanism claim",
        )
    elif not cdx2 and not guca2a:
        branch, consequence = (
            "NOT SEPARABLE AT THIS RESOLUTION",
            "item 3 closes; 13 versus 9 patients is the pre-stated reason, "
            "and this is not a biological null",
        )
    elif guca2a and not cdx2:
        branch, consequence = (
            "CONTRADICTS THE PUBLISHED PREDICTION",
            "reported without a post-hoc mechanism; the SER arm's hyperplastic "
            "majority is inspected for the obvious artefact",
        )
    else:
        branch, consequence = (
            "BOTH EXCLUDE ZERO",
            "not a branch §7 anticipated; reported as such and read as neither "
            "confirmation nor contradiction",
        )
    return {"weighting": weighting, "branch": branch, "consequence": consequence,
            "cdx2_excludes_zero": cdx2, "guca2a_excludes_zero": guca2a}


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
    patients = classify_patients(joined)
    assigned = patients[patients["status"] == "assigned"][["patient_id", "arm"]]
    decomposition = pd.read_parquet(args.decomposition)
    values = decomposition[decomposition["granularity_rung"] == RUNG].merge(
        assigned, on="patient_id", how="inner"
    )
    if values.empty:
        raise ContrastError("no arm-assigned patient joined to a decomposition value")

    gene_table = per_gene(values)
    pair_table = per_pair(values)
    branches = pd.DataFrame(
        [falsifier_branch(gene_table, w) for w in sorted(values["weighting"].unique())]
    )

    tested = gene_table[gene_table["term"] == TESTED_TERM]
    print("intrinsic — AD minus SER, Welch t (§5); the only tested term (§4):")
    print(tested[["gene", "weighting", "n_AD_estimable", "n_SER_estimable",
                  "mean_difference", "ci_low", "ci_high", "excludes_zero"]].to_string(index=False))
    print("\n§7 branches, one per weighting, none primary:")
    print(branches.to_string(index=False))
    agree = branches["branch"].nunique() == 1
    print(f"\n  branches agree across weightings: {agree}")
    print("  compositional and interaction terms are reported, not read (§4).")

    if args.no_write:
        return 0
    meta = {
        **provenance_meta(LABEL_PROVENANCE, CLAIM_PROVENANCE),
        "prereg": "docs/prereg_chen_lesion_subtype.md",
        "prereg_locked_at": "3a6d446",
        "granularity_rung": RUNG,
        "min_estimable_per_arm": MIN_ESTIMABLE_PER_ARM,
        "interval_method": "welch_t",
        "tested_term": TESTED_TERM,
        "compositional_is_descriptive_only": True,
        "branches_agree_across_weightings": bool(agree),
        "branch": branches["branch"].iloc[0] if agree else "DISAGREES ACROSS WEIGHTINGS",
    }
    for table, name in (
        (gene_table, "chen_subtype_arm_contrast"),
        (pair_table, "chen_subtype_arm_contrast_pairs"),
        (branches, "chen_subtype_falsifier_branch"),
    ):
        path = write_versioned_table(
            table, name, seed=args.seed, results_dir=args.results_dir,
            allow_dirty=args.allow_dirty, extra_meta=meta,
        )
        print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
