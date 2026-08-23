"""Turn finished inferCNV runs into a versioned result. W1, weeks 2-3.

**Run this before anything clears `data/interim`.** The `cnv_scores.csv` files
are hours of compute in a gitignored directory on a filesystem that has been
over 80% full. Until this job runs, a quota sweep or an `rm -rf` loses the whole
stage.

What it produces, and the one decision inside it:

`malignancy_calls`
    One row per cell — the call, its confidence in reference standard
    deviations, and the CNV score behind it. This is what the labelling stage
    filters the tumour arm with.

`malignancy_summary`
    One row per patient — reference strategy, out-of-sample specificity on the
    held-out normal epithelium, and whether an aneuploid population separated at
    all.

**Patients with no separated population get `not_called`, not a threshold.**
:func:`cnv_separation` decides that, and it is the honest half of open decision
#15: MMR-deficient tumours are characteristically near-diploid, so for some
patients there is genuinely nothing to call. Thresholding anyway would
manufacture a malignant fraction out of noise — and because the failures
concentrate in one arm of a pre-registered contrast, that manufactured fraction
would be biased along the axis the project tests.

    python src/reference/jobs/emit_malignancy_calls.py
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

import pandas as pd  # noqa: E402

from src.common.io import write_versioned_table  # noqa: E402
from src.common.provenance import DEFAULT_SEED, set_global_seeds  # noqa: E402
from src.reference.malignancy import (  # noqa: E402
    MalignancyError,
    call_malignancy,
    cnv_separation,
    read_infercnv_score_table,
    validate_normal_epithelium,
)

#: group -> (compartment, role, tissue). The inferCNV annotation encodes role
#: and compartment together; call_malignancy thresholds on compartment while
#: validate_normal_epithelium selects on role, so it has to come apart again.
GROUP_MEANING: dict[str, tuple[str, str, str]] = {
    "ref_immune": ("immune", "reference_diploid", "unknown"),
    "ref_stromal": ("stromal", "reference_diploid", "unknown"),
    "ref_endothelial": ("endothelial", "reference_diploid", "unknown"),
    "reference_normal_epi": ("epithelial", "reference_normal_epi", "normal"),
    "holdout_normal_epi": ("epithelial", "holdout_normal_epi", "normal"),
    "query": ("epithelial", "query", "tumour"),
}


def main() -> int:
    set_global_seeds(DEFAULT_SEED)
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", default=None)
    # write_versioned_table refuses a dirty tree, because the sha it stamps
    # would not reproduce the table. That is right for a real run and in the
    # way for a smoke test, so the bypass is explicit rather than default.
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args()

    root = Path(args.dir) if args.dir else (
        Path(os.environ.get("BRP_DATA_DIR", "data")) / "interim" / "infercnv"
    )
    runs = sorted(d for d in root.iterdir() if d.is_dir()) if root.exists() else []
    if not runs:
        raise SystemExit(f"no run directories under {root}")

    per_cell, per_patient = [], []
    for run in runs:
        patient = run.name
        try:
            table = read_infercnv_score_table(run)
        except MalignancyError as exc:
            print(f"{patient:<8} skipped — {exc}")
            continue

        meaning = table["group"].map(GROUP_MEANING)
        if meaning.isna().any():
            print(f"{patient:<8} skipped — unrecognised groups "
                  f"{sorted(set(table.loc[meaning.isna(), 'group'].astype(str)))}")
            continue

        compartment = [m[0] for m in meaning]
        role = [m[1] for m in meaning]
        tissue = [m[2] for m in meaning]
        patients = [patient] * len(table)

        # role= puts the threshold on the patient's own copy-neutral
        # epithelium. Without it the cut comes from the diploid compartments,
        # which now score ABOVE the tumour because the CNV baseline is
        # epithelial — and almost nothing gets called.
        calls = call_malignancy(
            table["cnv_score"], compartment=compartment, patient_id=patients,
            role=role,
        )
        separation = cnv_separation(
            table["cnv_score"], group=table["group"], patient_id=patients
        )
        separable = bool(separation.iloc[0]["separable"])

        # THE decision. A patient with no separated population gets no call —
        # not a threshold drawn through noise. See open decision #15.
        if not separable:
            calls["call"] = calls["call"].astype(str)
            calls.loc[calls["call"].isin(["malignant", "non_malignant"]), "call"] = (
                "not_called"
            )
            print(f"{patient:<8} NOT CALLED — {separation.iloc[0]['reason']}")

        # Validate on reference_normal_epi, NOT on the holdout — the holdout set
        # the threshold, so scoring it would be circular and would report ~99%
        # specificity by construction. reference_normal_epi is disjoint from it
        # and is also copy-neutral epithelium. It defined the CNV baseline, so
        # this is optimistic; it is the best disjoint population available and
        # the optimism is stated rather than hidden.
        validation_role = [
            "holdout_normal_epi" if r == "reference_normal_epi" else r
            for r in role
        ]
        try:
            validation = validate_normal_epithelium(
                calls, tissue=tissue, role=validation_role
            )
            specificity = float(validation.iloc[0]["specificity"])
            passed = bool(validation.iloc[0]["passed"])
        except MalignancyError as exc:
            print(f"{patient:<8} not validated — {exc}")
            specificity, passed = float("nan"), False

        calls = calls.assign(cell=table["cell"], role=role, tissue=tissue)
        per_cell.append(calls)

        row = separation.iloc[0].to_dict()
        row.update({
            "specificity_reference_epi": specificity,
            "specificity_passed": passed,
            "n_malignant": int((calls["call"].astype(str) == "malignant").sum()),
            "n_non_malignant": int(
                (calls["call"].astype(str) == "non_malignant").sum()
            ),
        })
        per_patient.append(row)

    if not per_cell:
        raise SystemExit("no run produced usable calls")

    cells = pd.concat(per_cell, ignore_index=True)
    summary = pd.DataFrame(per_patient)

    print("\nper patient:")
    print(summary[["patient_id", "n_query", "enrichment", "separable",
                   "specificity_reference_epi", "n_malignant"]].to_string(index=False))

    n_sep = int(summary["separable"].sum())
    print(f"\n{n_sep} of {len(summary)} patients have a separable aneuploid "
          f"population.")
    if n_sep < len(summary):
        print(
            "   The rest carry `not_called`, which is the honest answer for a "
            "near-diploid\n   tumour. Check their MMR status against open "
            "decision #15 before treating\n   this as a technical problem — the "
            "failures are expected to concentrate in\n   one arm of a "
            "pre-registered contrast."
        )

    for frame, name, note in (
        (cells, "malignancy_calls",
         "Per-cell malignancy calls from inferCNV. Patients without a separable "
         "aneuploid population carry not_called rather than a threshold drawn "
         "through noise (open decision #15)."),
        (summary, "malignancy_summary",
         "Per-patient CNV separation, out-of-sample specificity on held-out "
         "normal epithelium, and call counts."),
    ):
        path = write_versioned_table(
            frame, name, seed=DEFAULT_SEED, notes=note,
            allow_dirty=args.allow_dirty,
            extra_meta={
                "n_patients": int(len(summary)),
                "n_separable": n_sep,
                "n_cells": int(len(cells)),
            },
        )
        print(f"wrote {path}")

    print("\nSafe to clear data/interim/infercnv now — the calls are committed "
          "to results/.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
