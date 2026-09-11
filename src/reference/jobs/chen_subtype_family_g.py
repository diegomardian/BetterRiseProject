"""Family G closed on the floor, and the Synapse chase closed with it.

    python -m src.reference.jobs.chen_subtype_family_g --no-write

§3's family G is truncating `APC` against `BRAF` V600E. §7's floor is 8
estimable patients per arm, fixed before any count was seen. The arms are 5 and
4, so the contrast is **NOT ESTIMABLE** and no interval is reported on any
statistic.

**The mutation calls are already open.** They are in the pinned cBioPortal MAF,
not behind Synapse. So the certification, the 403 on `syn23520239`, and the
download ACL on the Level 3 VCFs were never the binding constraint on this
analysis — the arm sizes were, and they would not have changed. This job exists
to record that as a measurement rather than leave it as an inference.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

from src.common.io import write_versioned_table
from src.common.label_provenance import check_no_circular_claim, provenance_meta
from src.common.paths import RAW_DIR
from src.common.provenance import DEFAULT_SEED
from src.reference.chen_subtype_labels import CLAIM, FAMILY_G
from src.reference.jobs.chen_subtype_attrition import DECOMPOSITION, RUNG
from src.reference.jobs.chen_subtype_contrast import MIN_ESTIMABLE_PER_ARM
from src.reference.jobs.chen_subtype_id_provenance import DEFAULT_CROSSWALK

MUTATIONS = RAW_DIR / "cbioportal_hta11" / "mutations.json"
#: §3: "truncating APC (nonsense, frameshift, splice, start/stop-loss)".
TRUNCATING = frozenset({
    "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
    "Splice_Site", "Nonstop_Mutation", "Translation_Start_Site",
})

LABEL_PROVENANCE = FAMILY_G
CLAIM_PROVENANCE = CLAIM


def read_mutations(path: Path = MUTATIONS) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"pinned MAF not found: {path}; run gate 3 first")
    return pd.DataFrame(
        [
            {
                "sampleId": row["sampleId"],
                "gene": (row.get("gene") or {}).get("hugoGeneSymbol"),
                "mutation_type": row.get("mutationType"),
                "protein_change": row.get("proteinChange"),
            }
            for row in json.loads(path.read_text())
        ]
    )


def arms(mutations: pd.DataFrame, crosswalk: pd.DataFrame) -> dict[str, set[str]]:
    """Positive-only arms on the specimen-exact candidates, as §3 fixes them."""
    exact = crosswalk[crosswalk["specimen_exact"]].drop_duplicates("scRNA_biospecimen_id")
    scoped = mutations[mutations["sampleId"].isin(set(exact["scRNA_biospecimen_id"]))]
    apc = set(
        scoped.loc[
            (scoped["gene"] == "APC") & (scoped["mutation_type"].isin(TRUNCATING)), "sampleId"
        ]
    )
    braf = set(
        scoped.loc[
            (scoped["gene"] == "BRAF") & (scoped["protein_change"] == "V600E"), "sampleId"
        ]
    )
    patient = exact.set_index("scRNA_biospecimen_id")["patient_id"]
    to_patients = lambda ids: {patient[i] for i in ids if i in patient.index}  # noqa: E731
    return {"APC_truncating": to_patients(apc - braf), "BRAF_V600E": to_patients(braf - apc)}


def floor_table(arm_patients: dict[str, set[str]], decomposition: pd.DataFrame) -> pd.DataFrame:
    lineage = decomposition[decomposition["granularity_rung"] == RUNG]
    estimable = set(lineage.loc[lineage["intrinsic"].notna(), "patient_id"])
    rows = []
    for arm, patients in sorted(arm_patients.items()):
        n_estimable = len(patients & estimable)
        rows.append({
            "arm": arm,
            "n_patients": len(patients),
            "n_estimable": n_estimable,
            "floor": MIN_ESTIMABLE_PER_ARM,
            "clears_floor": n_estimable >= MIN_ESTIMABLE_PER_ARM,
            "estimability": (
                "estimated" if n_estimable >= MIN_ESTIMABLE_PER_ARM else "not_estimable"
            ),
            "mean_difference": None,
            "ci_low": None,
            "ci_high": None,
            "label_source": "pinned cBioPortal MAF — open access, not Synapse",
        })
    return pd.DataFrame(rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mutations", type=Path, default=MUTATIONS)
    parser.add_argument("--crosswalk", type=Path, default=DEFAULT_CROSSWALK)
    parser.add_argument("--decomposition", type=Path, default=DECOMPOSITION)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--no-write", action="store_true")
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)

    check_no_circular_claim(labels=LABEL_PROVENANCE, claim=CLAIM_PROVENANCE)

    table = floor_table(
        arms(read_mutations(args.mutations), pd.read_parquet(args.crosswalk)),
        pd.read_parquet(args.decomposition),
    )
    print(table.drop(columns=["label_source"]).to_string(index=False))
    closed = not table["clears_floor"].any()
    verdict = (
        "NOT ESTIMABLE — both arms below the §7 floor" if closed
        else "AN ARM CLEARS THE FLOOR"
    )
    print(f"\n  verdict: {verdict}")
    print("  labels came from the open MAF; Synapse access was never the binding constraint.")

    if args.no_write:
        return 0
    path = write_versioned_table(
        table, "chen_subtype_family_g_floor", seed=args.seed,
        results_dir=args.results_dir, allow_dirty=args.allow_dirty,
        extra_meta={
            **provenance_meta(LABEL_PROVENANCE, CLAIM_PROVENANCE),
            "prereg": "docs/prereg_chen_lesion_subtype.md",
            "prereg_locked_at": "3a6d446",
            "min_estimable_per_arm": MIN_ESTIMABLE_PER_ARM,
            "interval_reported": False,
            "synapse_access_required": False,
            "verdict": "NOT ESTIMABLE" if closed else "AN ARM CLEARS THE FLOOR",
        },
    )
    print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
