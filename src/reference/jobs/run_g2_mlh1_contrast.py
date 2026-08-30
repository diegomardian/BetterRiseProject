#!/usr/bin/env python
"""The refined tier-B test for G2, pre-registered 2026-08-22 and never run.

`docs/prereg_g2_mlh1.md`, written **before any target-gene expression was
examined**. This executes it and changes nothing about it.

WHY IT IS WORTH RUNNING AFTER G2 HAS ALREADY FAILED
---------------------------------------------------
G2's tier-B arm failed because MLH1 pooled across all patients shows nothing.
**MLH1 silencing is a property of one stratum, not of MMR deficiency**, so
pooling dilutes a stratum-specific effect toward zero — which is what was
observed. The pre-registration predicts exactly that and says where to look:

    mlh1_methylated       promoter hypermethylation silences transcription -> HIGH
    mlh1_intact_mmrd      MMRd via MSH2/MSH6/PMS2, MLH1 untouched          -> NEAR ZERO

And it is a **stronger validation than tier separation**, in the prereg's words:

    Tier separation can be satisfied by an estimator that merely tracks
    expression level. This cannot.

The negative control is mechanistic. `mlh1_intact_mmrd` patients reach the same
MSI-H phenotype through a different gene — same disease, same selective pressure,
MLH1 specifically spared. **If the intrinsic term fires there anyway, the term is
measuring something else.**

This does not rescue G2. G2 failed as pre-registered and stays failed. §6 of the
prereg: *supporting evidence for G2, not its primary basis.*

WHAT THIS REFUSES TO DO
-----------------------
- **No p-value.** §4: "there is no reading of this cohort on which a p-value from
  it means anything." 12 against 7 is the ceiling.
- **No pooling of weightings.** All three reported separately.
- **No silent dropping.** `not_estimable` patients are excluded and *counted* —
  dropping them quietly biases the arm that runs out of mature cells first, which
  is the arm the prediction is about.
- **No point estimate below 3 per arm.** §5: the test then returns no evidence,
  neither for nor against, and says so.

    python src/reference/jobs/run_g2_mlh1_contrast.py
"""

from __future__ import annotations

import argparse
import glob
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from src.common.io import write_versioned_table  # noqa: E402
from src.common.provenance import DEFAULT_SEED, set_global_seeds  # noqa: E402

GENE = "MLH1"
TREATMENT = "mlh1_methylated"
CONTROL = "mlh1_intact_mmrd"

#: §5: below this per arm the test returns no evidence, neither way.
MIN_PER_ARM = 3

#: Patient-level bootstrap for the interval (invariant 5). An interval, not a
#: test — §4 is explicit that no p-value from this cohort means anything.
N_BOOT = 2000


def _difference_ci(treat: np.ndarray, control: np.ndarray, *, seed: int) -> tuple:
    """Median difference, treatment minus control, with a patient bootstrap."""
    rng = np.random.default_rng(seed)
    point = float(np.median(treat) - np.median(control))
    draws = np.empty(N_BOOT)
    for i in range(N_BOOT):
        t = rng.choice(treat, size=treat.size, replace=True)
        c = rng.choice(control, size=control.size, replace=True)
        draws[i] = np.median(t) - np.median(c)
    return point, float(np.quantile(draws, 0.025)), float(np.quantile(draws, 0.975))


def main() -> int:
    set_global_seeds(DEFAULT_SEED)
    parser = argparse.ArgumentParser()
    parser.add_argument("--decomposition", default=None)
    parser.add_argument("--cohort", default=None)
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args()

    # PREFER THE MATCHED READ, and say which one was used rather than letting a
    # fallback decide silently. The first run of this job fell through to the
    # unmatched decomposition because the matched one lived on another branch,
    # and reported a result computed on the depth-confounded read without
    # anything in the output saying so.
    matched = sorted(glob.glob("results/*/decomposition_gse178341_matched.parquet"),
                     key=os.path.getmtime)
    unmatched = sorted(glob.glob("results/*/decomposition_gse178341.parquet"),
                       key=os.path.getmtime)
    if args.decomposition:
        dpath = args.decomposition
    elif matched:
        dpath = matched[-1]
    elif unmatched:
        dpath = unmatched[-1]
        print("!! no matched decomposition found — falling back to the UNMATCHED "
              "read.\n   lineage and crypt_position carry a depth confound "
              "there (#24.1, #45), so\n   any contrast from this run inherits "
              "it. Merge w1/matched-decomposition.")
    else:
        raise SystemExit("no decomposition found under results/")
    cpath = args.cohort or sorted(
        glob.glob("results/*/cohort_table.parquet"), key=os.path.getmtime
    )[-1]
    dec = pd.read_parquet(dpath)
    cohort = pd.read_parquet(cpath)[["patient_id", "mlh1_stratum"]]
    is_matched = dpath.endswith("_matched.parquet")
    print(f"decomposition: {dpath}")
    print(f"  read: {'DEPTH-MATCHED (#24.1)' if is_matched else 'UNMATCHED'}")
    print(f"cohort table : {cpath}")

    mlh1 = dec[dec["gene"] == GENE].merge(cohort, on="patient_id", how="left")
    if mlh1.empty:
        raise SystemExit(f"no {GENE} rows in {dpath}")

    # Excluded AND counted. §4: dropping them silently biases the arm that runs
    # out of mature cells first, which is the arm the prediction is about.
    excluded = mlh1[mlh1["estimability"] == "not_estimable"]
    usable = mlh1[mlh1["estimability"] != "not_estimable"].dropna(subset=["intrinsic"])

    print(f"\n{GENE} rows: {len(mlh1):,} · usable {len(usable):,} · "
          f"not_estimable {len(excluded):,}")
    print("\npatients per stratum, among usable rows:")
    print(usable.groupby("mlh1_stratum")["patient_id"].nunique().to_string())
    print("\nnot_estimable patients per stratum — the exclusion is not neutral:")
    print(excluded.groupby("mlh1_stratum")["patient_id"].nunique().to_string()
          or "  none")

    rows = []
    for (rung, axis, weighting), block in usable.groupby(
        ["granularity_rung", "labeling_axis", "weighting"], observed=True
    ):
        t = block[block["mlh1_stratum"] == TREATMENT].groupby("patient_id")["intrinsic"].median()
        c = block[block["mlh1_stratum"] == CONTROL].groupby("patient_id")["intrinsic"].median()
        row = {
            "granularity_rung": rung, "labeling_axis": axis, "weighting": weighting,
            "n_methylated": int(t.size), "n_intact_mmrd": int(c.size),
            "median_methylated": float(t.median()) if t.size else None,
            "median_intact_mmrd": float(c.median()) if c.size else None,
        }
        # §5: fewer than 3 per arm and this returns NO EVIDENCE. The point
        # estimate is withheld rather than reported beside a wide interval,
        # because a number on the page is read whatever the caveat says.
        if t.size < MIN_PER_ARM or c.size < MIN_PER_ARM:
            row |= {"difference": None, "ci_low": None, "ci_high": None,
                    "verdict": "no_evidence",
                    "reason": f"{t.size} vs {c.size} patients, below the "
                              f"pre-registered floor of {MIN_PER_ARM} per arm"}
        else:
            d, lo, hi = _difference_ci(t.to_numpy(), c.to_numpy(), seed=DEFAULT_SEED)
            # Direction is the whole test (§1). More negative = more intrinsic
            # loss, so the prediction is difference < 0.
            row |= {"difference": d, "ci_low": lo, "ci_high": hi,
                    "verdict": "as_predicted" if d < 0 else "falsified",
                    "reason": ("methylated shows the larger intrinsic loss"
                               if d < 0 else
                               "intact_mmrd shows loss comparable to or greater "
                               "than methylated — §5's first falsifier")}
        rows.append(row)

    out = pd.DataFrame(rows)

    print("\n" + "=" * 70)
    print("REFINED TIER-B TEST — MLH1, methylated vs intact-MMRd")
    print("=" * 70)
    show = ["granularity_rung", "labeling_axis", "weighting", "n_methylated",
            "n_intact_mmrd", "difference", "ci_low", "ci_high", "verdict"]
    print(out[show].to_string(index=False))

    evaluable = out[out["verdict"] != "no_evidence"]
    print("\n" + "=" * 70)
    print("AGAINST THE PRE-REGISTERED FALSIFIERS")
    print("=" * 70)
    if evaluable.empty:
        print(
            "  NO EVIDENCE, either way. Every combination fell below the "
            f"{MIN_PER_ARM}-per-arm floor\n  §5 fixed in advance. The cohort "
            "cannot run this test; that is a statement\n  about the cohort, not "
            "about MLH1, and the point estimates are withheld\n  rather than "
            "printed beside an interval nobody could use."
        )
    else:
        n_pred = int((evaluable["verdict"] == "as_predicted").sum())
        print(f"  evaluable combinations: {len(evaluable)}")
        print(f"  as predicted (methylated loses more): {n_pred}")
        print(f"  falsified: {len(evaluable) - n_pred}")
        # §5's second falsifier: a reversal ACROSS rungs means labelling
        # artifact rather than biology.
        by_rung = evaluable.groupby("granularity_rung")["difference"].median()
        signs = set(np.sign(by_rung.dropna()))
        print(f"\n  direction by rung: "
              f"{', '.join(f'{r} {v:+.3f}' for r, v in by_rung.items())}")
        if len(signs) > 1:
            print(
                "\n  !! THE DIRECTION REVERSES ACROSS RUNGS. §5's second "
                "falsifier: that\n     means a labelling artifact rather than "
                "biology, whatever the\n     individual intervals say."
            )
        else:
            print("\n  Direction is consistent across rungs — §5's second "
                  "falsifier does not fire.")

    print("\n  NOT a hypothesis test and no p-value is computed (§4). "
          "Supporting\n  evidence for G2, not its primary basis (§6) — G2 failed "
          "as\n  pre-registered and this does not change that.")

    path = write_versioned_table(
        out, "g2_mlh1_stratified_contrast", seed=DEFAULT_SEED,
        allow_dirty=args.allow_dirty,
        notes=(
            "The refined tier-B test pre-registered in docs/prereg_g2_mlh1.md on "
            "2026-08-22, before any target-gene expression was examined. "
            "mlh1_methylated against mlh1_intact_mmrd, matched patients, all "
            "three weightings reported separately, not_estimable excluded AND "
            "counted. Reported as a difference with a patient-level bootstrap "
            "interval, never as a hypothesis test — §4 states no p-value from "
            "this cohort would mean anything. Below 3 patients per arm the "
            "point estimate is WITHHELD, per §5. Supporting evidence for G2, "
            "not its primary basis; G2 failed as pre-registered and this does "
            "not change that."
            + (
                " Computed on the DEPTH-MATCHED decomposition (#24.1)."
                if is_matched else
                " WARNING: computed on the UNMATCHED decomposition, where "
                "lineage and crypt_position carry a depth confound."
            )
        ),
        extra_meta={
            "decomposition": dpath, "cohort_table": cpath,
            "depth_matched_read": bool(is_matched),
            "n_not_estimable_rows": int(len(excluded)),
            "min_per_arm": MIN_PER_ARM, "n_boot": N_BOOT,
            "prereg": "docs/prereg_g2_mlh1.md, 2026-08-22",
        },
    )
    print(f"\nwrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
