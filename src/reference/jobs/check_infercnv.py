"""Is the inferCNV run right, not merely finished? W1, weeks 2-3.

Four questions, in increasing order of how much they tell you:

1. **Did it produce a score at all?** ``cnv_scores.csv`` exists, one row per cell.
2. **Did the gene join work?** The gene-order file is keyed on symbol, so if
   GENCODE v28 and the deposit disagree on naming, inferCNV does not error — it
   infers copy number from whatever survived. A large drop here means the
   answer came from a fraction of the genome.
3. **Do the scores order the way biology says they must?** Diploid reference
   cells lowest, tumour epithelium highest. If they do not, something upstream
   is wrong and no threshold will rescue it.
4. **Is normal epithelium misread as tumour?** The held-out 30% was never in the
   baseline, so scoring it is an out-of-sample test — the one execution_plan.md
   asks for by name ("normal epithelium not misread as tumour").

**Questions 3 and 4 are both needed, and 4 alone is not enough.** That is not
obvious, and it is the reason this script runs both.

execution_plan.md §4 names question 4 — "normal epithelium not misread as
tumour" — as the "done when" for this stage. But the threshold in
`call_malignancy` is a quantile of the *reference* cells' own scores, so if the
reference groups arrive mislabelled the threshold inverts along with them. On a
deliberately inverted fixture, tumour epithelium scored below the diploid
reference and question 4 still reported **98% specificity and passed**, while
every call was backwards.

Question 3 catches that, because it asks about the ordering rather than about a
count either side of a line that moved.

    python src/reference/jobs/check_infercnv.py
    python src/reference/jobs/check_infercnv.py --dir /path/to/infercnv
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

import pandas as pd  # noqa: E402

from src.reference.malignancy import (  # noqa: E402
    MalignancyError,
    call_malignancy,
    read_infercnv_score_table,
    validate_normal_epithelium,
)

#: The group each cell was written with, taken back apart. The name encodes role
#: and compartment together — convenient for inferCNV, and it has to be undone
#: here because call_malignancy thresholds on compartment while
#: validate_normal_epithelium selects on role.
GROUP_MEANING: dict[str, tuple[str, str, str]] = {
    # group -> (compartment, role, tissue)
    "ref_immune": ("immune", "reference_diploid", "unknown"),
    "ref_stromal": ("stromal", "reference_diploid", "unknown"),
    "ref_endothelial": ("endothelial", "reference_diploid", "unknown"),
    "reference_normal_epi": ("epithelial", "reference_normal_epi", "normal"),
    "holdout_normal_epi": ("epithelial", "holdout_normal_epi", "normal"),
    "query": ("epithelial", "query", "tumour"),
}

#: The comparator, and it is NOT the diploid compartments.
#:
#: Once the baseline is matched normal epithelium alone, immune/stromal/
#: endothelial cells are no longer references — they are other cell types being
#: scored against an epithelial baseline, so they deviate for cell-type reasons
#: that have nothing to do with copy number. Comparing tumour epithelium to them
#: measures the epithelial-vs-immune difference, not aneuploidy.
#:
#: `holdout_normal_epi` is the right comparator: the same cell type as the
#: query, from the same patient, never in the baseline. The difference between
#: it and the query is copy number and nothing else.
COMPARATOR = "holdout_normal_epi"

#: Reported for context, never used as the comparator. Expected to be HIGH
#: against an epithelial baseline; that is cell-type difference, not a problem.
NON_COMPARATOR_GROUPS = ("ref_immune", "ref_stromal", "ref_endothelial")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", default=None)
    args = parser.parse_args()

    root = Path(args.dir) if args.dir else (
        Path(os.environ.get("BRP_DATA_DIR", "data")) / "interim" / "infercnv"
    )
    if not root.exists():
        raise SystemExit(f"no {root} — has anything run?")

    runs = sorted(d for d in root.iterdir() if d.is_dir())
    if not runs:
        raise SystemExit(f"no run directories under {root}")

    rows = []
    for run in runs:
        patient = run.name
        try:
            table = read_infercnv_score_table(run)
        except Exception as exc:  # noqa: BLE001 - report, do not abort the sweep
            print(f"{patient:<8} NO SCORES — {exc}")
            continue

        if "group" not in table.columns or table["group"].isna().all():
            print(f"{patient:<8} scores present but no annotations.tsv to join")
            continue

        medians = table.groupby("group", observed=True)["cnv_score"].median()
        query = table[table["group"] == "query"]["cnv_score"]
        holdout = table[table["group"] == COMPARATOR]["cnv_score"]
        if not len(query) or not len(holdout):
            print(f"{patient:<8} missing query or {COMPARATOR} — cannot compare")
            continue

        # The MEDIAN query cell is not expected to be malignant. The tumour
        # sample contains substantial non-malignant epithelium — that is the
        # whole reason malignancy calling exists — so the median is diluted by
        # it. What has to separate is the upper tail: a subpopulation clearly
        # above the copy-neutral distribution.
        tail = float(query.quantile(0.90))
        holdout_tail = float(holdout.quantile(0.90))
        rows.append({
            "patient_id": patient,
            "n_cells": len(table),
            "holdout_median": float(holdout.median()),
            "query_median": float(query.median()),
            "holdout_p90": holdout_tail,
            "query_p90": tail,
            "tail_ratio": tail / holdout_tail if holdout_tail else float("nan"),
            "median_above": bool(query.median() > holdout.median()),
            "frac_query_above_holdout_p90": float((query > holdout_tail).mean()),
        })

        print(f"\n--- {patient} — {len(table):,} cells ---")
        print(medians.rename("median_cnv_score").to_string())

    if not rows:
        raise SystemExit("\nno run produced a readable score table")

    summary = pd.DataFrame(rows)
    print("\n" + "=" * 66)
    print("DO THE SCORES ORDER THE WAY BIOLOGY REQUIRES?")
    print("=" * 66)
    print(summary.to_string(index=False))

    # Scale first: an inverted ordering on a flattened matrix is a symptom, and
    # chasing the reference groups when the real problem is the score's scale
    # wastes a day. Scale is judged on the QUERY's tail, not on copy-neutral
    # cells — normal epithelium is supposed to sit near 1.
    summary["implied_per_gene_dev"] = summary["query_p90"] ** 0.5
    flat = summary[summary["implied_per_gene_dev"] < 0.05]
    if len(flat):
        print(
            "\n!! implied per-gene deviation is under 0.05 for "
            + ", ".join(flat["patient_id"])
            + ".\n   A real copy-number change moves the smoothed residual by "
              "0.1-0.3, so this\n   matrix has been flattened and the score is "
              "measuring residual noise.\n"
              "   Two known causes, in the order they were ruled out here:\n"
              "     - MULTIPLE reference groups. STEP 08 runs use_bounds=TRUE, "
              "which zeroes\n       observation deviation inside the range of "
              "the reference-group means.\n       Groups spanning different cell "
              "types make that range wide. Check\n       ref_group_names in "
              "run_infercnv.R — one matched-normal group is right.\n"
              "     - denoise = TRUE, which sets values within 1.5 SD of the "
              "reference mean\n       TO the mean. Already off by default; "
              "turning it off alone moved the\n       pilot from 0.017 to 0.023 "
              "and fixed nothing.\n"
              "   `grep 'fraction exactly 1' logs/` tells them apart: bounding "
              "leaves a\n   large fraction at exactly 1 even with denoise off."
        )

    broken = summary[~summary["median_above"]]
    if len(broken):
        print(
            "\n!! " + ", ".join(broken["patient_id"])
            + " score tumour epithelium at or BELOW their own held-out normal\n"
              "   epithelium — same cell type, same patient, so the only "
              "difference left is\n   copy number. Either those tumours carry "
              "little detectable CNV, or the\n   reference groups reached "
              "inferCNV mislabelled."
        )

    # The number that actually decides whether malignancy calling is possible.
    # A tumour sample is a MIXTURE, so the median says little; what matters is
    # whether a subpopulation separates from the copy-neutral distribution.
    weak = summary[summary["frac_query_above_holdout_p90"] < 0.20]
    if len(weak):
        print(
            "\n!! " + ", ".join(weak["patient_id"])
            + " have under 20% of tumour epithelium above the held-out\n"
              "   normal's 90th percentile. There is no separated population to "
              "call malignant.\n   That may be real — a near-diploid tumour, or "
              "a sample that is mostly\n   non-malignant epithelium — but those "
              "patients cannot carry a malignancy\n   call, and saying so is "
              "better than thresholding noise."
        )
    else:
        print(
            "\n   Every patient has a tumour subpopulation separated from its "
            "own copy-neutral\n   epithelium. That is what malignancy calling "
            "needs, and it is measured\n   against the cell-type-matched "
            "holdout rather than against other cell types."
        )

    # THE out-of-sample check. execution_plan.md §4 lists it as the "done when"
    # for this stage: normal epithelium must not be misread as tumour. The
    # held-out 30% never entered the baseline, so this is genuinely out of
    # sample — validating on baseline cells would be circular.
    print("\n" + "=" * 66)
    print("IS NORMAL EPITHELIUM BEING MISREAD AS TUMOUR?")
    print("(held-out cells only — they never entered the CNV baseline)")
    print("=" * 66)
    validations = []
    for run in runs:
        try:
            table = read_infercnv_score_table(run)
        except Exception:  # noqa: BLE001 - already reported above
            continue
        if "group" not in table.columns or table["group"].isna().all():
            continue

        meaning = table["group"].map(GROUP_MEANING)
        if meaning.isna().any():
            unknown = sorted(set(table.loc[meaning.isna(), "group"].astype(str)))
            print(f"{run.name:<8} unrecognised group(s) {unknown} — skipping")
            continue

        try:
            calls = call_malignancy(
                table["cnv_score"],
                compartment=[m[0] for m in meaning],
                patient_id=[run.name] * len(table),
            )
            validations.append(
                validate_normal_epithelium(
                    calls,
                    tissue=[m[2] for m in meaning],
                    role=[m[1] for m in meaning],
                )
            )
        except MalignancyError as exc:
            print(f"{run.name:<8} cannot validate — {exc}")

    if validations:
        report = pd.concat(validations, ignore_index=True)
        print(report.to_string(index=False))
        failed = report[~report["passed"]]
        if len(failed):
            print(
                "\n!! " + ", ".join(failed["patient_id"])
                + " call held-out NORMAL epithelium malignant too often.\n"
                  "   execution_plan.md §4 says stop here. Every downstream "
                  "compositional and\n   intrinsic number would be computed over "
                  "a tumour arm contaminated with\n   normal cells, or a normal "
                  "arm stripped of them."
            )
        else:
            print(
                "\n   Every patient clears the specificity floor, on cells the "
                "baseline never saw.\n   This is the check §4 names as the "
                "'done when' for malignancy calling."
            )

    print(
        "\nNEXT: write these calls to a versioned parquet under results/ before\n"
        "anything clears data/interim — cnv_scores.csv is gitignored and the run\n"
        "is hours of compute."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
