#!/usr/bin/env python
"""Pick the five pilot patients, and check the remaining batch variables. W1, week 1-2.

    python src/reference/jobs/select_pilot.py

Two jobs, both wanted before the pilot runs:

1. **Batch confounding.** Chemistry was already cleared (between-patient only).
   SOURCE_HOSPITAL, TISSUE_PROCESSING_TEAM and PROCESSING_TYPE are the other
   batch variables in the metadata and have not been checked. A variable that
   tracks tissue within patients would contaminate the paired comparison the
   same way chemistry would have.

2. **Pilot selection, written down before any results exist.** The plan is
   explicit that the five are chosen and justified up front — otherwise the
   pilot becomes a place to look for patients that behave well.

Criteria, applied in order and all recorded in the output:

  - matched tumour AND normal (an unmatched patient cannot exercise the
    decomposition at all)
  - **unsorted cells in both arms.** PROCESSING_TYPE is per-sample, and
    CD45pMACS / LiveMACS / mixUnsortCD45MACS samples have had their cell-type
    composition deliberately altered. The compositional term is
    Delta(mature epithelial fraction), so comparing a sorted sample against an
    unsorted one measures the sort, not the tumour. Counts below are therefore
    computed on unsorted cells only.
  - at least MIN_TUMOUR / MIN_NORMAL unsorted cells in each arm, so positivity
    is not the thing being tested at this stage
  - spans the tier-B strata: at least one mlh1_methylated, one
    mlh1_intact_mmrd, one mmr_proficient, so G2's controls are exercised early
  - spans cohort depth: smallest, median and largest eligible tumours, so the
    pipeline meets its range rather than its easy cases

Deterministic. Writes results/pilot_selection.csv plus the rationale.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import pandas as pd

from src.common.provenance import DEFAULT_SEED, provenance_record
from src.reference.ingest import (
    assign_mlh1_strata,
    patient_cohort_table,
    read_gse178341_index,
    read_gse178341_metadata,
)

MIN_TUMOUR = 500
MIN_NORMAL = 300
N_PILOT = 5

BATCH_VARIABLES = ("SOURCE_HOSPITAL", "TISSUE_PROCESSING_TEAM", "PROCESSING_TYPE")

#: The only PROCESSING_TYPE whose cell-type composition is unmanipulated.
#: CD45pMACS is immune enrichment, LiveMACS is viability selection, and
#: mixUnsortCD45MACS is a deliberate mixture — none can carry a compositional
#: estimate. See docs/open_decisions.md #11.
UNSORTED = "unsorted"


def batch_confounding(obs: pd.DataFrame, metadata: pd.DataFrame) -> None:
    """Report each batch variable against tissue, and within-patient consistency."""
    joined = obs.join(metadata, how="left")
    for variable in BATCH_VARIABLES:
        if variable not in joined.columns:
            print(f"\n--- {variable}: not in metadata ---")
            continue
        print(f"\n--- {variable} vs tissue ---")
        print(pd.crosstab(joined[variable], joined["tissue"]))
        per_patient = joined.groupby("patient_id", observed=True)[variable].nunique()
        mixed = per_patient[per_patient > 1]
        print(
            f"patients on a single {variable}: {int((per_patient == 1).sum())} "
            f"of {len(per_patient)}"
        )
        if len(mixed):
            print(f"  MIXED WITHIN PATIENT -> possible confounder: {list(mixed.index)}")
        else:
            print("  constant within every patient -> the paired design absorbs it")


def select(table: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    """Apply the criteria. Returns (selection, rationale lines)."""
    notes: list[str] = []

    eligible = table[
        table["matched"]
        & (table["n_tumour_unsorted"] >= MIN_TUMOUR)
        & (table["n_normal_unsorted"] >= MIN_NORMAL)
    ].copy()
    notes.append(
        f"{int(table['matched'].sum())} of {len(table)} patients are matched; "
        f"{int(table['compositionally_usable'].sum())} of those have UNSORTED cells "
        f"in both arms; {len(eligible)} also clear >={MIN_TUMOUR} tumour and "
        f">={MIN_NORMAL} normal unsorted cells."
    )
    notes.append(
        "Sorted samples (CD45pMACS, LiveMACS, mixUnsortCD45MACS) are excluded from "
        "the counts: their cell-type composition is manipulated, so they cannot "
        "carry a compositional estimate (open decision #11)."
    )
    if len(eligible) < N_PILOT:
        raise SystemExit(f"only {len(eligible)} eligible patients; loosen the thresholds")

    eligible = eligible.sort_values("n_tumour_unsorted")
    chosen: list[str] = []

    # One from each stratum first, so tier B's controls are present from day one.
    for stratum in ("mlh1_methylated", "mlh1_intact_mmrd", "mmr_proficient"):
        pool = eligible[(eligible["mlh1_stratum"] == stratum) & (~eligible.index.isin(chosen))]
        if pool.empty:
            notes.append(f"WARNING: no eligible patient in stratum {stratum}")
            continue
        # Median-depth representative of the stratum.
        pick = pool.index[len(pool) // 2]
        chosen.append(pick)
        notes.append(f"{pick}: median-depth representative of {stratum}")

    # Then the extremes of depth, so the pipeline meets its range.
    for label, pick in (
        ("smallest eligible tumour", eligible.index[0]),
        ("largest eligible tumour", eligible.index[-1]),
    ):
        if pick not in chosen and len(chosen) < N_PILOT:
            chosen.append(pick)
            notes.append(f"{pick}: {label}, to exercise the depth range")

    # Fill any shortfall with the next most balanced patients.
    while len(chosen) < N_PILOT:
        pool = eligible[~eligible.index.isin(chosen)]
        if pool.empty:
            break
        pick = (pool["n_normal"] / pool["n_tumour"]).sub(1).abs().idxmin()
        chosen.append(pick)
        notes.append(f"{pick}: most balanced tumour/normal ratio among the remainder")

    return table.loc[chosen], notes


def main() -> int:
    data = Path(os.environ.get("BRP_DATA_DIR", "data")) / "raw" / "GSE178341"
    h5 = data / "GSE178341_crc10x_full_c295v4_submit.h5"
    meta_csv = data / "GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz"
    if not h5.exists():
        raise SystemExit(f"missing {h5} — download it first")

    obs, _ = read_gse178341_index(h5)
    metadata = read_gse178341_metadata(meta_csv)

    print("=" * 70)
    print("BATCH CONFOUNDING — the variables chemistry did not cover")
    print("=" * 70)
    batch_confounding(obs, metadata)

    table = assign_mlh1_strata(patient_cohort_table(obs, metadata), metadata)

    # Compositional counts must come from unsorted samples only. PROCESSING_TYPE
    # is per-sample, so this filters cells rather than dropping patients — a
    # patient with both a CD45-sorted and an unsorted tumour keeps the latter.
    joined = obs.join(metadata, how="left")
    unsorted = patient_cohort_table(joined[joined["PROCESSING_TYPE"] == UNSORTED])
    for arm in ("n_tumour", "n_normal"):
        table[f"{arm}_unsorted"] = (
            unsorted[arm].reindex(table.index).fillna(0).astype(int)
        )
    table["compositionally_usable"] = (table["n_tumour_unsorted"] > 0) & (
        table["n_normal_unsorted"] > 0
    )

    print("\n" + "=" * 70)
    print("COMPOSITIONAL COHORT — unsorted samples only")
    print("=" * 70)
    print(f"matched:                       {int(table['matched'].sum())} of {len(table)}")
    print(f"unsorted in both arms:         {int(table['compositionally_usable'].sum())}")
    deep = (table["n_tumour_unsorted"] >= MIN_TUMOUR) & (
        table["n_normal_unsorted"] >= MIN_NORMAL
    )
    print(f"  of those, adequate depth:    {int(deep.sum())}")
    print("\nby tier-B stratum (unsorted in both arms):")
    print(table[table["compositionally_usable"]]["mlh1_stratum"].value_counts().to_string())

    print("\n" + "=" * 70)
    print("PILOT SELECTION")
    print("=" * 70)
    selection, notes = select(table)
    for note in notes:
        print(f"  - {note}")
    print()
    columns = [
        c
        for c in ("n_tumour", "n_normal", "n_tumour_unsorted", "n_normal_unsorted",
                  "n_samples", "chemistry", "MMRStatus", "MLH1Status", "mlh1_stratum")
        if c in selection.columns
    ]
    print(selection[columns].to_string())

    out = Path("results")
    out.mkdir(exist_ok=True)
    selection.to_csv(out / "pilot_selection.csv")
    provenance = provenance_record(seed=DEFAULT_SEED)
    (out / "pilot_selection_rationale.txt").write_text(
        "Five-patient pilot selection — chosen before any expression was examined.\n\n"
        + "\n".join(f"- {n}" for n in notes)
        + f"\n\ncriteria: matched, UNSORTED n_tumour>={MIN_TUMOUR} and "
        f"n_normal>={MIN_NORMAL}, strata spanned, depth range spanned. "
        f"Sorted samples excluded from all counts — open decision #11.\n"
        f"git_sha: {provenance['git_sha']}\nseed: {DEFAULT_SEED}\n",
        encoding="utf-8",
    )
    print("\nwritten: results/pilot_selection.csv and _rationale.txt")
    print(f"\npilot patients: {list(selection.index)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
