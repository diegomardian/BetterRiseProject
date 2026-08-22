"""Emit the per-patient cohort table. W1, week 2.

Does two jobs that happen to need the same table, and both are about making a
number visible **before** it becomes convenient to have it be different.

**Open decision #9 — the 26 patients with no matched normal.** The compositional
term is Delta(mature fraction) against the patient's *own* normal, so an
unmatched patient contributes to neither arm. Emitting them as `not_estimable`
rows in the results frame was the obvious move and is the wrong one: the schema
has no field for *why*, so "no normal arm" would be indistinguishable from "too
few mature cells", and `gate_g4_verdict` counts patients below the mature-cell
threshold — 26 guaranteed-zero rows would push a pre-committed gate criterion
toward failure for a reason that has nothing to do with positivity. So they are
visible here instead, in a table of their own, where nothing downstream counts
them by accident.

**Open decision #10 — the refined tier-B test.** The `mlh1_stratum` column is the
cohort composition the G2 prediction in `docs/prereg_g2_mlh1.md` is committed
against. Writing it as a versioned artifact, with a git sha and before any
expression has been examined, is what makes that pre-registration checkable
rather than a claim about what someone remembers intending.

    python src/reference/jobs/emit_cohort_table.py
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from src.common.io import write_versioned_table  # noqa: E402
from src.common.provenance import DEFAULT_SEED, set_global_seeds  # noqa: E402
from src.reference.ingest import (  # noqa: E402
    assign_mlh1_strata,
    patient_cohort_table,
    read_gse178341_index,
    read_gse178341_metadata,
)

#: Processing type that leaves composition untouched (open decision #11).
UNSORTED = "unsorted"


def main() -> int:
    set_global_seeds(DEFAULT_SEED)
    data = Path(os.environ.get("BRP_DATA_DIR", "data")) / "raw" / "GSE178341"
    h5 = data / "GSE178341_crc10x_full_c295v4_submit.h5"
    meta_csv = data / "GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz"
    for path in (h5, meta_csv):
        if not path.exists():
            raise SystemExit(f"missing {path}")

    obs, _var = read_gse178341_index(h5)
    metadata = read_gse178341_metadata(meta_csv)
    table = assign_mlh1_strata(patient_cohort_table(obs, metadata), metadata)

    # The compositional n is not the matched n: sorted samples cannot carry a
    # composition (#11), so a patient matched only through a CD45-enriched arm
    # is matched for the intrinsic term and not for the compositional one.
    # PROCESSING_TYPE is per-SAMPLE, so this filters cells rather than dropping
    # patients: someone with both a CD45-sorted and an unsorted tumour keeps the
    # unsorted one. Same index join select_pilot.py uses — deliberately copied
    # rather than reinvented, because the two tables have to agree.
    joined = obs.join(metadata, how="left")
    unsorted = patient_cohort_table(joined[joined["PROCESSING_TYPE"] == UNSORTED])
    for arm in ("n_tumour", "n_normal"):
        table[f"{arm}_unsorted"] = (
            unsorted[arm].reindex(table.index).fillna(0).astype(int)
        )
    table["compositionally_usable"] = (table["n_tumour_unsorted"] > 0) & (
        table["n_normal_unsorted"] > 0
    )

    table = table.reset_index()
    n = len(table)
    n_matched = int(table["matched"].sum())
    n_compositional = int(table["compositionally_usable"].sum())

    print(f"patients                       {n:>4}")
    print(f"matched tumour + normal        {n_matched:>4}  "
          f"({n - n_matched} contribute to neither term — decision #9)")
    print(f"matched AND unsorted in both   {n_compositional:>4}  "
          f"<- the real compositional n (decision #11)")
    print("\nmlh1_stratum x matched — the composition decision #10 commits against:")
    print(
        table.groupby(["mlh1_stratum", "matched"], observed=True)
        .size().rename("n_patients").reset_index().to_string(index=False)
    )

    path = write_versioned_table(
        table, "cohort_table", seed=DEFAULT_SEED,
        notes=(
            "Per-patient cohort composition for GSE178341. Backs open decision "
            "#9 (unmatched patients are visible here rather than as "
            "not_estimable rows, which would contaminate gate G4) and open "
            "decision #10 (the mlh1_stratum composition the G2 prediction in "
            "docs/prereg_g2_mlh1.md is committed against, written before any "
            "expression was examined)."
        ),
        extra_meta={
            "n_patients": n,
            "n_matched": n_matched,
            "n_matched_unsorted": n_compositional,
            "strata": table["mlh1_stratum"].value_counts().to_dict(),
        },
    )
    print(f"\nwrote {path}")
    print("COMMIT THIS. Its value is the timestamp — a pre-registration written "
          "after\nthe fact is not a pre-registration (CONTRIBUTING §4: results go "
          "in git).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
