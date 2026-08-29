#!/usr/bin/env python
"""The decomposition, and whether any rung may be quoted. W1 + W4, §6.1.

Runs `decompose_cohort` and `bootstrap_over_patients` over W1's per-gene summary
and writes the frozen-schema result. **The estimator is W4's and is called, not
reimplemented** — every cutpoint and every `intrinsic = None` decision happens
inside `decompose_cohort` (invariant 1, and `src/estimator/README.md` says not to
duplicate `classify_estimability`).

WHAT THIS IS FOR, BEYOND THE NUMBERS
------------------------------------
§6.1 asks for the per-patient split at every rung and axis. Producing it is the
easy half. The hard half is that **each rung is disqualified by a different
check, and all four checks are correct**:

  epithelial       every scored cell is mature, so Δfraction ≈ 0 by
                   construction — the compositional term cannot move
  lineage          depth-confounded on 14–17% of patients (#45 diagnostic,
                   `results/2026-08-27_e5ebdc3/`)
  crypt_position   depth-confounded, and a two-bin split on ~90% of patients
                   (#42), so not independent of `lineage`
  best4            clean, and never reaches the estimability cutpoint — max 45
                   mature cells against a floor of 50, which is G4's finding
                   (#48) on W1's own cohort

So this emits a `quotable` column per rung with the reason attached, and the
summary prints the matrix. **If nothing is quotable, that is the result** — a
decomposition nobody may cite is a finding about the cohort, not a failure to
compute one.

    qsub src/reference/jobs/run_decomposition.sh
"""

from __future__ import annotations

import argparse
import glob
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

import pandas as pd  # noqa: E402

from src.common.provenance import DEFAULT_SEED, set_global_seeds  # noqa: E402
from src.estimator.kitagawa import (  # noqa: E402
    attach_intrinsic_ci,
    bootstrap_over_patients,
    decompose_cohort,
)
from src.schema import coerce_results, write_results  # noqa: E402

#: Why each rung may or may not be quoted. Written here, before the numbers are
#: read, so the disqualifications are not chosen after seeing which rung is
#: interesting. Every entry cites the measurement that established it.
QUOTABILITY: dict[str, tuple[bool, str]] = {
    "epithelial": (
        False,
        "degenerate by construction: every scored cell is mature, so the "
        "compositional term cannot move (222/232 rows at mature_fraction 1.000)",
    ),
    "lineage": (
        False,
        "depth-confounded on 17.2% of patients (W2's diagnostic, PR #45); the "
        "arms are 1.64x apart in median depth on 20 of 32 patients",
    ),
    "crypt_position": (
        False,
        "depth-confounded on 14.1% of patients, AND a two-bin split on ~90% of "
        "them (#42), so not an independent point on the granularity curve",
    ),
    "best4": (
        False,
        "clean on depth (median |rho| 0.069) but never estimable: max 45 mature "
        "cells against a cutpoint of 50 — G4's finding (#48) on this cohort",
    ),
}


def main() -> int:
    set_global_seeds(DEFAULT_SEED)
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", default=None)
    parser.add_argument("--n-boot", type=int, default=2000)
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args()

    path = args.summary or sorted(
        glob.glob("results/*/decomposition_summary.parquet"), key=os.path.getmtime
    )[-1]
    summary = pd.read_parquet(path)
    print(f"summary: {path}")
    print(f"  {len(summary):,} rows · {summary.patient_id.nunique()} patients · "
          f"{summary.gene.nunique()} genes")

    # A row with no mature cell in an arm has no mean to difference. Dropped
    # rather than filled: an absent mean is not a mean of zero (invariant 1),
    # and decompose() would otherwise return a term computed from a fabrication.
    usable = summary.dropna(subset=["mean_normal", "mean_tumour"])
    dropped = len(summary) - len(usable)
    if dropped:
        print(f"  {dropped} row(s) dropped: an arm had no mature cell to average")

    print("\nrunning decompose_cohort …")
    result = decompose_cohort(usable)
    print(f"  {len(result):,} rows (one per weighting)")

    print(f"running bootstrap_over_patients, n_boot={args.n_boot} …")
    print("  patients, not cells — invariant 5")
    boot = bootstrap_over_patients(usable, n_boot=args.n_boot, seed=DEFAULT_SEED)

    # attach_intrinsic_ci, not a hand-rolled merge. The bootstrap frame is LONG
    # (one row per term) and decompose_cohort's is WIDE (three term columns), so
    # merging on the shared keys fans out 3x — a first version produced 2,916
    # rows from 972 and every third one was wrong.
    #
    # More importantly the schema has ONE ci_low/ci_high slot for three point
    # estimates, and which term fills it is decision #10, closed 2026-08-22:
    # the cohort-level INTRINSIC band, because estimability is defined for
    # intrinsic and not for the other two. This helper is where that decision
    # lives; re-deriving it here would be the fifth instance of duplicating a
    # mapping until the two drift.
    merged = attach_intrinsic_ci(result, boot)
    print(f"  CIs attached: {int(merged.ci_low.notna().sum()):,} of {len(merged):,} rows")
    print("  NOTE: a cohort-level band broadcast onto each patient row, not a")
    print("        per-patient interval (#10). Say so when presenting it.")

    print("\n" + "=" * 68)
    print("ESTIMABILITY, BY RUNG")
    print("=" * 68)
    print(pd.crosstab(merged.granularity_rung, merged.estimability).to_string())

    print("\n" + "=" * 68)
    print("MAY ANY RUNG BE QUOTED?")
    print("=" * 68)
    for rung in sorted(merged.granularity_rung.unique()):
        ok, why = QUOTABILITY.get(rung, (True, ""))
        mark = "QUOTABLE" if ok else "no"
        print(f"  {rung:<16} {mark:<9} {why}")
    if not any(ok for ok, _ in QUOTABILITY.values()):
        print(
            "\n  NO RUNG CLEARS EVERY CHECK. That is the result, not a failure to\n"
            "  produce one: each rung is disqualified by a different measurement\n"
            "  and all four measurements are correct. §5's consequence for G4 —\n"
            "  'non-identifiability with diagnostics becomes the headline result,\n"
            "  not a caveat' — is the honest reading of this table."
        )

    out = coerce_results(merged)
    written = write_results(
        out, "decomposition_gse178341", seed=DEFAULT_SEED,
        allow_dirty=args.allow_dirty,
        notes=(
            "Kitagawa decomposition, GSE178341, all tiers A-D, four rungs, both "
            "live axes, with patient-level bootstrap CIs (invariant 5). Estimator "
            "is W4's decompose_cohort, called not reimplemented. Every rung "
            "carries `quotable` and the measurement that disqualified it: "
            "epithelial degenerate by construction, lineage and crypt_position "
            "depth-confounded (PR #45), best4 never estimable (G4, #48). Rows "
            "where an arm had no mature cell are dropped, not filled — an absent "
            "mean is not a mean of zero.\n\nQUOTABILITY, per rung — the schema is "
            "frozen so this travels in the sidecar rather than as a column:\n"
            + "\n".join(
                f"  {rung}: {'quotable' if ok else 'NOT quotable'} — {why}"
                for rung, (ok, why) in sorted(QUOTABILITY.items())
            )
        ),
    )
    print(f"\nwrote {written}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
