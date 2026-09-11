"""POST-HOC: do the arms differ in the technical covariates, not the biology?

    python -m src.reference.jobs.chen_subtype_control_diagnostic --no-write

**This is post-hoc and it is not in the locked design.** It was written after
`results/2026-09-10_dddd34f/` showed `ACTB` and `KRT8` separating the arms more
strongly than any target, and it claims nothing about the target programme —
§8 bars the interesting readings and this does not go near them.

What it asks is narrower: are the two arms different in **sequencing depth,
cell count, or mature-cell resolution**? Those are properties of the
measurement, not of the tissue programme, and they are already committed in
`icbi_adenoma.parquet`. A yes makes the weighting-unstable result easier to
read. A no leaves the control movement unexplained, which is also worth
recording rather than leaving as an untested hunch.

No interval is reported and no threshold is applied. This is a description.
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
from src.reference.jobs.chen_subtype_attrition import RUNG
from src.reference.jobs.chen_subtype_id_provenance import (
    ARMS,
    DEFAULT_CROSSWALK,
    SNAPSHOT,
    classify_patients,
    exact_join,
    read_labels,
)

INPUTS = Path("results/2026-09-06_765eb29/icbi_adenoma.parquet")
#: Measurement properties only. Nothing here is a transcript programme.
COVARIATES = (
    "depth_normal", "depth_tumour", "depth_ratio",
    "n_normal", "n_tumour",
    "n_cells_resolved_normal", "n_cells_resolved_tumour",
    "n_cells_epithelial_normal", "n_cells_epithelial_tumour",
    "unresolved_share_normal", "unresolved_share_tumour", "unresolved_arm_gap",
    "frac_mature_normal", "frac_mature_tumour",
)

LABEL_PROVENANCE = FAMILY_P
CLAIM_PROVENANCE = CLAIM


def describe(inputs: pd.DataFrame, assigned: pd.DataFrame) -> pd.DataFrame:
    """Per covariate: each arm's median, and a standardised difference."""
    block = inputs[inputs["granularity_rung"] == RUNG].merge(
        assigned, on="patient_id", how="inner"
    )
    # One row per patient: these are patient-level properties repeated per gene.
    patients = block.groupby(["patient_id", "arm"], as_index=False).first()
    rows = []
    for covariate in COVARIATES:
        if covariate not in patients:
            continue
        groups = {
            arm: patients.loc[patients["arm"] == arm, covariate].dropna().to_numpy(dtype=float)
            for arm in ARMS
        }
        a, b = groups["AD"], groups["SER"]
        row: dict[str, object] = {
            "covariate": covariate,
            "n_AD": int(a.size), "n_SER": int(b.size),
            "median_AD": float(np.median(a)) if a.size else None,
            "median_SER": float(np.median(b)) if b.size else None,
            "mean_AD": float(a.mean()) if a.size else None,
            "mean_SER": float(b.mean()) if b.size else None,
        }
        if a.size > 1 and b.size > 1:
            pooled = np.sqrt(
                ((a.size - 1) * a.var(ddof=1) + (b.size - 1) * b.var(ddof=1))
                / (a.size + b.size - 2)
            )
            row["standardised_difference"] = (
                float((a.mean() - b.mean()) / pooled) if pooled > 0 else None
            )
        else:
            row["standardised_difference"] = None
        rows.append(row)
    frame = pd.DataFrame(rows)
    frame["post_hoc"] = True
    frame["interval_reported"] = False
    return frame


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, default=INPUTS)
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
    table = describe(pd.read_parquet(args.inputs), assigned)

    ordered = table.reindex(
        pd.to_numeric(table["standardised_difference"], errors="coerce")
        .abs().sort_values(ascending=False).index
    )
    print("POST-HOC — technical covariates only, no interval, no threshold:")
    print(ordered[["covariate", "n_AD", "n_SER", "median_AD", "median_SER",
                   "standardised_difference"]].to_string(index=False))

    if args.no_write:
        return 0
    path = write_versioned_table(
        table, "chen_subtype_control_diagnostic", seed=args.seed,
        results_dir=args.results_dir, allow_dirty=args.allow_dirty,
        extra_meta={
            **provenance_meta(LABEL_PROVENANCE, CLAIM_PROVENANCE),
            "prereg": "docs/prereg_chen_lesion_subtype.md",
            "pre_registered": False,
            "post_hoc": True,
            "claims_nothing_about_the_target_programme": True,
            "interval_reported": False,
        },
    )
    print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
