"""Gate 4: the labelled lesions are the same physical specimens as the scRNA ones.

    python -m src.reference.jobs.chen_subtype_id_provenance --no-write

Gate condition 4 of `docs/prereg_chen_lesion_subtype.md` §9, on the standard
`wes_subtype_plan.md` §1 rule: **exact equality, suffix families not trusted**.

Two biospecimens from one participant differ only in a suffix — `HTA11_866`
carries both `..._2000001011` and `..._3004761011`, and they are different
lesions with different diagnoses. A prefix or fuzzy match would silently merge
them and inflate an arm. This job matches on the full identifier and refuses
anything else.

It reads the pinned snapshot, never the live API, and it reads no decomposition
value.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

from src.common.io import write_versioned_table
from src.common.label_provenance import (
    Measurement,
    check_no_circular_claim,
    provenance_meta,
)
from src.common.paths import RAW_DIR
from src.common.provenance import DEFAULT_SEED

SNAPSHOT = RAW_DIR / "cbioportal_hta11" / "clinical_sample.json"
DEFAULT_CROSSWALK = Path("results/2026-09-09_0bf9734/wes_subtype_provenance_crosswalk.parquet")
ARMS = ("AD", "SER")

#: Frozen in prereg §3 before any outcome was read. This job's purpose is to
#: fail loudly if the snapshot no longer reproduces them.
#:
#: The prereg's table mixes two denominators in adjacent columns: its lesion
#: counts are *labelled* lesions, which include the conflicted patient's two,
#: while its patient counts are *analysable* patients, which exclude that
#: patient. Both are right and they do not describe the same set, so both are
#: checked here under names that say which is which.
EXPECTED = {
    ("AD", "labelled_lesions"): 16, ("AD", "analysable_lesions"): 15, ("AD", "patients"): 13,
    ("SER", "labelled_lesions"): 11, ("SER", "analysable_lesions"): 10, ("SER", "patients"): 9,
}
EXPECTED_CONFLICTED = ("HTA11_6801",)
EXPECTED_MULTI_CONCORDANT = ("HTA11_6818", "HTA11_8622", "HTA11_866")

LABEL_PROVENANCE = Measurement(
    modality="morphology",
    assay="Chen 2021 pathologist polyp diagnosis via pinned cBioPortal HTA11 snapshot",
    genes=(),
)
CLAIM_PROVENANCE = Measurement(
    modality="transcript",
    assay="Chen 2021 snRNA-seq adenoma decomposition (lineage rung)",
    genes=(),
)


class IdProvenanceError(RuntimeError):
    """The identifier join cannot be established as the pre-registration fixes it."""


def read_labels(path: Path = SNAPSHOT) -> pd.DataFrame:
    """Sample-level labels from the pinned snapshot, with case normalised."""
    if not path.is_file():
        raise IdProvenanceError(f"pinned snapshot not found: {path}; run gate 3 first")
    records = json.loads(path.read_text())
    frame = pd.DataFrame(
        [
            {"sampleId": r["sampleId"], "patientId": r["patientId"],
             "attr": r["clinicalAttributeId"], "value": r["value"]}
            for r in records
        ]
    )
    wide = frame.pivot_table(
        index=["sampleId", "patientId"], columns="attr", values="value", aggfunc="first"
    ).reset_index()
    # Rule 4: case is normalised before any field is read.
    for column in ("ADVANCED", "ATYPIA"):
        if column in wide:
            wide[column] = wide[column].astype(str).str.strip().str.casefold()
    return wide


def exact_join(labels: pd.DataFrame, crosswalk: pd.DataFrame) -> pd.DataFrame:
    """Join on the full identifier only; a suffix family is not a match."""
    specimens = crosswalk.drop_duplicates("scRNA_biospecimen_id")
    joined = specimens.merge(
        labels, left_on="scRNA_biospecimen_id", right_on="sampleId", how="left",
        validate="one_to_one",
    )
    joined["label_matched"] = joined["sampleId"].notna()

    # A match that is not full-string equality is a defect, not a near miss.
    matched = joined[joined["label_matched"]]
    if not (matched["scRNA_biospecimen_id"] == matched["sampleId"]).all():
        raise IdProvenanceError("a non-identical identifier was accepted as a match")

    # Prefix families must stay distinct: same participant, different specimen.
    families = matched.groupby(
        matched["scRNA_biospecimen_id"].str.rsplit("_", n=1).str[0]
    )["scRNA_biospecimen_id"].nunique()
    joined.attrs["n_prefix_families_with_multiple_specimens"] = int((families > 1).sum())
    return joined


def classify_patients(joined: pd.DataFrame) -> pd.DataFrame:
    """Per patient: the arm, or the reason there is not one."""
    labelled = joined[joined["label_matched"] & joined["POLYP_TYPE"].isin(ARMS)]
    rows = []
    for patient, block in labelled.groupby("patient_id", sort=True):
        arms = sorted(set(block["POLYP_TYPE"]))
        rows.append(
            {
                "patient_id": patient,
                "cbio_patient_id": block["patientId"].iloc[0],
                "n_labelled_lesions": len(block),
                "labels": ";".join(sorted(block["POLYP_TYPE"])),
                # Rule 1: conflicting patients are excluded, never assigned.
                "arm": arms[0] if len(arms) == 1 else None,
                "status": "assigned" if len(arms) == 1 else "excluded_conflicting",
            }
        )
    return pd.DataFrame(rows).sort_values("patient_id").reset_index(drop=True)


def arm_summary(patients: pd.DataFrame, joined: pd.DataFrame) -> pd.DataFrame:
    assigned = patients[patients["status"] == "assigned"]
    labelled = joined[joined["label_matched"] & joined["POLYP_TYPE"].isin(ARMS)]
    rows = []
    for arm in ARMS:
        block = assigned[assigned["arm"] == arm]
        rows.append(
            {
                "arm": arm,
                # Every lesion carrying the label, conflicted patients included.
                "labelled_lesions": int((labelled["POLYP_TYPE"] == arm).sum()),
                # The lesions that actually enter the arm.
                "analysable_lesions": int(block["n_labelled_lesions"].sum()),
                "patients": int(len(block)),
                "expected_labelled_lesions": EXPECTED[(arm, "labelled_lesions")],
                "expected_analysable_lesions": EXPECTED[(arm, "analysable_lesions")],
                "expected_patients": EXPECTED[(arm, "patients")],
            }
        )
    frame = pd.DataFrame(rows)
    frame["matches_prereg"] = (
        (frame["labelled_lesions"] == frame["expected_labelled_lesions"])
        & (frame["analysable_lesions"] == frame["expected_analysable_lesions"])
        & (frame["patients"] == frame["expected_patients"])
    )
    return frame


def verdict(summary: pd.DataFrame, patients: pd.DataFrame, joined: pd.DataFrame) -> dict:
    conflicted = tuple(
        sorted(patients.loc[patients["status"] == "excluded_conflicting", "cbio_patient_id"])
    )
    multi = tuple(
        sorted(
            patients.loc[
                (patients["status"] == "assigned") & (patients["n_labelled_lesions"] > 1),
                "cbio_patient_id",
            ]
        )
    )
    reproduced = bool(summary["matches_prereg"].all())
    named = conflicted == EXPECTED_CONFLICTED and multi == EXPECTED_MULTI_CONCORDANT
    return {
        "arms_reproduce_prereg": reproduced,
        "named_patients_reproduce": named,
        "conflicted_patients": ";".join(conflicted),
        "multi_lesion_concordant_patients": ";".join(multi),
        "n_specimens_matched": int(joined["label_matched"].sum()),
        "n_specimens_unmatched": int((~joined["label_matched"]).sum()),
        "n_prefix_families_with_multiple_specimens": int(
            joined.attrs.get("n_prefix_families_with_multiple_specimens", 0)
        ),
        "suffix_family_match_used": False,
        "decomposition_value_read": False,
        "verdict": (
            "ID PROVENANCE HOLDS — arms and named patients reproduce"
            if reproduced and named
            else "ID PROVENANCE FAILS — the snapshot does not reproduce the frozen arms"
        ),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--crosswalk", type=Path, default=DEFAULT_CROSSWALK)
    parser.add_argument("--snapshot", type=Path, default=SNAPSHOT)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--no-write", action="store_true")
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)

    # Invariant 11: refuse before the first read, not at the writer.
    check_no_circular_claim(labels=LABEL_PROVENANCE, claim=CLAIM_PROVENANCE)

    labels = read_labels(args.snapshot)
    crosswalk = pd.read_parquet(args.crosswalk)
    joined = exact_join(labels, crosswalk)
    patients = classify_patients(joined)
    summary = arm_summary(patients, joined)
    outcome = verdict(summary, patients, joined)

    print(summary.to_string(index=False))
    print()
    for key, value in outcome.items():
        print(f"  {key:44} {value}")
    print("\nEXACT EQUALITY ONLY — no suffix family was matched, no decomposition value read.")

    if args.no_write:
        return 0 if outcome["arms_reproduce_prereg"] else 5
    meta = {
        **provenance_meta(LABEL_PROVENANCE, CLAIM_PROVENANCE),
        "prereg": "docs/prereg_chen_lesion_subtype.md",
        "gate_condition": 4,
        "crosswalk": str(args.crosswalk),
        "snapshot": str(args.snapshot),
        **outcome,
    }
    for table, name in (
        (summary, "chen_subtype_arm_sizes"),
        (patients, "chen_subtype_patient_assignment"),
    ):
        path = write_versioned_table(
            table, name, seed=args.seed, results_dir=args.results_dir,
            allow_dirty=args.allow_dirty, extra_meta=meta,
        )
        print(f"wrote {path}")
    return 0 if outcome["arms_reproduce_prereg"] else 5


if __name__ == "__main__":
    sys.exit(main())
