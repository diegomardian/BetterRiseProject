"""Can the adenoma question be asked on the ICBI atlas? Path C, option 1.

    python -m src.reference.jobs.icbi_adenoma_feasibility

WHY THIS IS THE FIRST THING TO ASK. Every route to the mechanism question on
carcinoma data has now terminated (docs/HANDOFF.md §2). The remaining family is
different SUBSTRATE rather than more of the same, and its cheapest member is
the adenoma arm of an atlas already on disk: `sample_type` carries a `polyp`
category the carcinoma pipeline excludes, and 106,958 cells sit in it.

An adenoma reading would be worth having for two reasons the carcinoma one
cannot supply. `best4` -- the resolution where the question is actually posed --
is unestimable in carcinoma at a median of 3 mature cells, and adenomas retain
far more differentiated epithelium. And a matched normal->polyp->carcinoma
series constrains WHEN the loss happens, which no two-arm comparison can.

THE ANSWER IS YES, AND THE FIRST RUN GOT IT WRONG. Run 2026-09-05:

    Chen_2021_Cell    93,913 polyp cells, 106 patients, 100% raw counts,
                      44 with a matched normal at >= 100 epithelial cells/arm
    Zheng_2022        13,045 polyp cells, 3 patients with the full gradient

`Chen_2021_Cell` is the VUMC/HTAN polyp atlas -- `dataset` reads
VUMC_HTAN_discovery / _validation / _cohort3 / _CRC and the sample ids are
HTA11_*. It is the cohort a data hunt would have gone looking for, and it was
already on disk.

An earlier version of this job reported ZERO usable pairs for it, because its
reference samples are labelled `healthy normal` and that label was read as
"a different donor". Fifty-one of its patients carry BOTH a polyp and a
`healthy normal` under the same patient_id: it is the patient's own unaffected
mucosa. The reading only ever compares two arms of one patient and every caller
groups by patient_id first, so the per-patient grouping -- not the label -- is
what rules out cross-donor pairing. Excluding the label cost the entire cohort.

That was the fifth vocabulary error in this repository and the first one made
inside a feasibility verdict, where it would have redirected weeks of work
toward fetching data already present.

WHAT THIS STILL HANDS THE NEXT COHORT. Check `patient_id` overlap between arms,
not the arm LABELS. A summary saying "polyp and normal samples present" does not
distinguish 44 usable pairs from 0, and neither does a cell count.

This reads only the cached obs. No expression, no cluster.
"""

from __future__ import annotations

import argparse
import logging
import sys
from collections.abc import Sequence
from pathlib import Path

import pandas as pd

from src.common.io import write_versioned_table
from src.common.paths import INTERIM_DIR
from src.common.provenance import DEFAULT_SEED
from src.reference.icbi_slice import COMPARTMENT_MAP

log = logging.getLogger(__name__)

#: The arms an adenoma reading uses. The reference is EITHER normal label: what
#: makes it the patient's own is the patient_id grouping, not the string.
POLYP, TUMOUR = "polyp", "primary tumor"
REFERENCE = ("adjacent normal", "healthy normal")

#: Same floors as the carcinoma feasibility, so the two are comparable.
CELL_THRESHOLDS: tuple[int, ...] = (50, 100, 200, 500)

EPITHELIAL = tuple(k for k, v in COMPARTMENT_MAP.items() if v == "epithelial")


def polyp_studies(obs: pd.DataFrame) -> pd.DataFrame:
    """One row per study carrying polyp cells, with what it can support."""
    rows = []
    for study, block in obs[obs["sample_type"] == POLYP].groupby("study_id", observed=True):
        whole = obs[obs["study_id"] == study]
        naive = whole[whole["enrichment_cell_types"].astype(str) == "naive"]
        epithelial = naive[
            naive["atlas_cell_type_coarse"].astype(str).isin(EPITHELIAL)
        ]
        counts = (
            epithelial.groupby(["patient_id", "sample_type"], observed=True)
            .size().unstack(fill_value=0)
        )
        for arm in (POLYP, TUMOUR, *REFERENCE):
            if arm not in counts.columns:
                counts[arm] = 0
        # Either normal label serves, because both are this patient's own.
        counts["reference"] = counts[list(REFERENCE)].max(axis=1)

        row = {
            "study_id": str(study),
            "polyp_cells": int(len(block)),
            "polyp_patients": int(block["patient_id"].nunique()),
            "raw_counts_share": float(
                (block["matrix_type"].astype(str) == "raw counts").mean()
            ),
            "has_adjacent_normal": bool((whole["sample_type"] == REFERENCE[0]).any()),
            "has_healthy_normal": bool((whole["sample_type"] == REFERENCE[1]).any()),
            "patients_with_both_arms": int(
                ((counts[POLYP] > 0) & (counts["reference"] > 0)).sum()
            ),
            "median_epithelial_genes": float(epithelial["n_genes"].median())
            if len(epithelial) else float("nan"),
        }
        for threshold in CELL_THRESHOLDS:
            paired = counts[
                (counts[POLYP] >= threshold) & (counts["reference"] >= threshold)
            ]
            gradient = paired[paired[TUMOUR] >= threshold]
            row[f"paired_at_{threshold}"] = int(len(paired))
            row[f"gradient_at_{threshold}"] = int(len(gradient))
        rows.append(row)
    return pd.DataFrame(rows).sort_values("polyp_cells", ascending=False, ignore_index=True)


def verdict(table: pd.DataFrame, threshold: int = 100) -> tuple[str, str]:
    """Is the adenoma question askable here? Needs MIN_PREMISE_PATIENTS pairs."""
    from src.reference.jobs.coexpression_silencing import MIN_PREMISE_PATIENTS

    if table.empty:
        return "NOT VIABLE", "the atlas carries no polyp cells at all"
    best = int(table[f"paired_at_{threshold}"].max())
    total = int(table[f"paired_at_{threshold}"].sum())
    if best >= MIN_PREMISE_PATIENTS * 3:
        return "VIABLE", (
            f"one study carries {best} patients with matched polyp and adjacent "
            f"normal at >= {threshold} epithelial cells per arm."
        )
    if total >= MIN_PREMISE_PATIENTS:
        return "MARGINAL", (
            f"{total} matched patient(s) across {int((table[f'paired_at_{threshold}'] > 0).sum())} "
            f"study/studies, best single study {best}. At or near the "
            f"{MIN_PREMISE_PATIENTS}-patient floor, so a premise could be "
            f"computed but not resolved. Worth running only as one study in a "
            f"meta-analysis."
        )
    return "NOT VIABLE", (
        f"{total} matched patient(s) at >= {threshold} epithelial cells per arm, "
        f"below the {MIN_PREMISE_PATIENTS} floor. Check patient_id overlap "
        f"between arms before concluding this -- reading the arm LABELS instead "
        f"is what once put this cohort's usable pairs at zero."
    )


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=INTERIM_DIR / "icbi_obs.parquet")
    parser.add_argument("--threshold", type=int, default=100)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    if not args.cache.exists():
        raise SystemExit(
            f"{args.cache} not found. Run\n"
            f"    python -m src.reference.jobs.pull_icbi_metadata\n"
            f"first (it accepts a local atlas path as well as the URL)."
        )
    obs = pd.read_parquet(args.cache)
    log.info("%s cells in the cached obs", f"{len(obs):,}")

    table = polyp_studies(obs)
    if table.empty:
        log.error("no polyp cells in this atlas build")
        return 2

    log.info("\n%s", table.to_string(index=False))
    state, detail = verdict(table, args.threshold)
    log.info("\n%s\nPATH C, OPTION 1 (ICBI polyp arm): %s\n%s\n%s",
             "=" * 72, state, detail, "=" * 72)
    log.info(
        "\nCheck patient_id OVERLAP between arms, not the arm labels. A summary\n"
        "saying 'polyp and normal samples present' does not distinguish 44\n"
        "usable pairs from 0, and neither does a cell count."
    )

    path = write_versioned_table(
        table, "icbi_adenoma_feasibility", seed=DEFAULT_SEED,
        results_dir=args.results_dir, allow_dirty=args.allow_dirty,
        extra_meta={
            "verdict": state, "verdict_detail": detail,
            "threshold_epithelial_per_arm": args.threshold,
            "arms": {"polyp": POLYP, "reference": list(REFERENCE), "carcinoma": TUMOUR},
            "healthy_normal_included": (
                "`healthy normal` counts as the reference arm. What makes it "
                "the patient's own is the patient_id grouping, not the label -- "
                "51 of Chen_2021_Cell's patients carry both a polyp and a "
                "`healthy normal` under one patient_id. Reading the label "
                "instead put this cohort's usable pairs at 0 against 44."
            ),
            "what_this_answers": (
                "Whether the adenoma question can be asked on data already on "
                "disk -- the cheapest member of the different-substrate family "
                "that remains after every carcinoma route terminated."
            ),
            "what_this_does_not_answer": (
                "Whether an adenoma reading would resolve the premise on a "
                "cohort that DOES have matched normals. This is availability, "
                "not power."
            ),
        },
    )
    log.info("wrote %s", path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
