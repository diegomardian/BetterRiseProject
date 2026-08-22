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

Question 3 is the one worth waiting for. A run can finish, produce a
well-formed CSV, and have the tumour cells scoring *below* the reference —
which would mean the reference groups were mislabelled, and every downstream
number would be confidently backwards.

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

from src.reference.malignancy import read_infercnv_score_table  # noqa: E402

#: Groups that should carry the LOWEST scores — they defined the baseline.
REFERENCE_GROUPS = ("ref_immune", "ref_stromal", "ref_endothelial",
                    "reference_normal_epi")


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
        reference = table[table["group"].isin(REFERENCE_GROUPS)]["cnv_score"]
        query = table[table["group"] == "query"]["cnv_score"]
        holdout = table[table["group"] == "holdout_normal_epi"]["cnv_score"]

        ordered = (
            bool(query.median() > reference.median())
            if len(query) and len(reference) else None
        )
        rows.append({
            "patient_id": patient,
            "n_cells": len(table),
            "median_reference": float(reference.median()) if len(reference) else float("nan"),
            "median_holdout": float(holdout.median()) if len(holdout) else float("nan"),
            "median_query": float(query.median()) if len(query) else float("nan"),
            "query_above_reference": ordered,
        })

        print(f"\\n--- {patient} — {len(table):,} cells ---")
        print(medians.rename("median_cnv_score").to_string())

    if not rows:
        raise SystemExit("\\nno run produced a readable score table")

    summary = pd.DataFrame(rows)
    print("\\n" + "=" * 66)
    print("DO THE SCORES ORDER THE WAY BIOLOGY REQUIRES?")
    print("=" * 66)
    print(summary.to_string(index=False))

    broken = summary[summary["query_above_reference"] == False]  # noqa: E712
    if len(broken):
        print(
            "\\n!! " + ", ".join(broken["patient_id"])
            + " score tumour epithelium BELOW the diploid reference.\\n"
              "   That is not a threshold problem and no cutoff will fix it — it "
              "means the\\n   reference groups reached inferCNV mislabelled, or "
              "the gene order file is\\n   wrong for this build. Every downstream "
              "call from these patients would be\\n   confidently backwards. Stop "
              "and find the cause."
        )
    else:
        print(
            "\\n   Tumour epithelium scores above the diploid reference in every "
            "patient.\\n   That is the ordering aneuploidy predicts, and it is "
            "the cheapest evidence\\n   that the reference groups arrived intact."
        )

    print(
        "\\nNEXT: call_malignancy() thresholds these per patient against the\\n"
        "reference cells' own scores, then validate_normal_epithelium() scores\\n"
        "the held-out normal epithelium — the out-of-sample check that\\n"
        "execution_plan.md asks for by name."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
