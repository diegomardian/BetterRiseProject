"""Can the ICBI atlas resolve the coexpression premise? Sizing, before compute.

    python -m src.reference.jobs.icbi_premise_feasibility

The premise behind the silencing reading came back UNRESOLVED on all three
cohorts: each one's control interval straddles the 0.5 log2 tolerance, so
whether the arms are comparable is undecided at that many patients. That is a
POWER problem, not a structural one -- unlike the decomposition's algebraic
collapse and the survivorship confound, which no number of extra patients
touches. So it is the one blocker more data can move, and this job measures
whether the ICBI atlas actually carries the data.

It answers the question in the shape invariant 4 requires. The useful quantity
is NOT the pooled patient count -- pooling is forbidden, and per-study intervals
do not tighten just because other studies exist. It is **how many independent
study-level estimates there are to meta-analyse**. We have three. If the atlas
carries fourteen, that is the difference between three undecided readings and a
meta-analysis.

Reads the cached obs from `pull_icbi_metadata` and touches no expression data.
"""

from __future__ import annotations

import argparse
import logging
import sys
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd

from src.common.io import write_versioned_table
from src.common.paths import INTERIM_DIR
from src.common.provenance import DEFAULT_SEED
from src.reference.icbi import (
    NORMAL_SAMPLE_TYPES,
    TUMOUR_SAMPLE_TYPES,
    ICBIError,
    assert_vocabulary_matches,
)

log = logging.getLogger(__name__)

#: Coarse labels that are epithelium. `Cancer cell` is included deliberately:
#: in the diseased arm the mature colonocytes we are asking about are labelled
#: Cancer cell, and excluding it drops most tumour-arm epithelium. Counting only
#: `Epithelial cell` gives 7 usable patients against 136 -- a nineteen-fold
#: difference produced entirely by a label choice, so it is named here rather
#: than buried in a filter.
EPITHELIAL_LABELS: tuple[str, ...] = ("Epithelial cell", "Cancer cell")

#: Cells per arm below which the premise check has nothing to measure. Not a
#: positivity cutpoint -- the calibrated one governs the decomposition, and this
#: is a coarser question -- but the same idea, and reported at several values so
#: the number is visible rather than assumed.
CELL_THRESHOLDS: tuple[int, ...] = (50, 100, 200, 500, 1000)

#: The studies already read. Their contribution is subtracted, because "the
#: atlas has 136 paired patients" is not the same claim as "the atlas adds 85".
ALREADY_READ: tuple[str, ...] = ("Lee_2020_Nat_Genet", "Pelka_2021_Cell")


def paired_epithelial_counts(obs: pd.DataFrame) -> pd.DataFrame:
    """Epithelial cells per (study, patient, arm), for patients holding both."""
    assert_vocabulary_matches(
        obs["sample_type"], TUMOUR_SAMPLE_TYPES + NORMAL_SAMPLE_TYPES,
        field="sample_type",
    )
    assert_vocabulary_matches(
        obs["atlas_cell_type_coarse"], EPITHELIAL_LABELS, field="atlas_cell_type_coarse",
    )
    arms = obs[obs["sample_type"].isin(TUMOUR_SAMPLE_TYPES + NORMAL_SAMPLE_TYPES)].copy()
    arms["arm"] = np.where(arms["sample_type"].isin(TUMOUR_SAMPLE_TYPES), "tumour", "normal")

    epithelial = arms[arms["atlas_cell_type_coarse"].astype(str).isin(EPITHELIAL_LABELS)]
    counts = (
        epithelial.groupby(["study_id", "patient_id", "arm"], observed=True)
        .size().unstack(fill_value=0)
    )
    for arm in ("tumour", "normal"):
        if arm not in counts.columns:
            counts[arm] = 0
    both = counts[(counts["tumour"] > 0) & (counts["normal"] > 0)].copy()
    both["min_arm"] = both[["tumour", "normal"]].min(axis=1)
    return both.reset_index()


def feasibility(counts: pd.DataFrame, depth: pd.Series) -> pd.DataFrame:
    """One row per threshold: patients, studies, and how many are NEW."""
    rows = []
    for threshold in CELL_THRESHOLDS:
        usable = counts[counts["min_arm"] >= threshold]
        new = usable[~usable["study_id"].isin(ALREADY_READ)]
        rows.append({
            "min_cells_per_arm": threshold,
            "n_patients": len(usable),
            "n_studies": usable["study_id"].nunique(),
            "n_patients_new": len(new),
            "n_studies_new": new["study_id"].nunique(),
            "n_patients_already_read": len(usable) - len(new),
            "median_epithelial_genes": float(
                depth.reindex(usable["study_id"].unique()).median()
            ) if len(usable) else float("nan"),
        })
    return pd.DataFrame(rows)


def per_study(counts: pd.DataFrame, depth: pd.Series, threshold: int = 100) -> pd.DataFrame:
    """The table invariant 4 actually needs: one row per candidate study."""
    usable = counts[counts["min_arm"] >= threshold]
    out = (
        usable.groupby("study_id", observed=True)
        .agg(n_paired_patients=("patient_id", "nunique"),
             median_min_arm=("min_arm", "median"),
             max_min_arm=("min_arm", "max"))
        .reset_index()
    )
    out["median_epithelial_genes"] = out["study_id"].map(depth).astype(float)
    out["already_read"] = out["study_id"].isin(ALREADY_READ)
    out["min_cells_per_arm"] = threshold
    return out.sort_values("n_paired_patients", ascending=False, ignore_index=True)


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
            f"first -- it reads /obs over range requests and caches it."
        )
    obs = pd.read_parquet(args.cache)
    log.info("%s cells in the cached obs", f"{len(obs):,}")

    try:
        counts = paired_epithelial_counts(obs)
    except ICBIError as exc:
        log.error("%s", exc)
        return 2

    epithelial = obs[obs["atlas_cell_type_coarse"].astype(str).isin(EPITHELIAL_LABELS)]
    depth = epithelial.groupby("study_id", observed=True)["n_genes"].median()

    table = feasibility(counts, depth)
    studies = per_study(counts, depth, args.threshold)

    log.info("\npaired patients with epithelium in BOTH arms: %d", len(counts))
    log.info("%s", table.to_string(index=False))
    log.info("\ncandidate studies at >= %d cells per arm:", args.threshold)
    log.info("%s", studies.to_string(index=False))

    at_threshold = table[table["min_cells_per_arm"] == args.threshold].iloc[0]
    log.info(
        "\nWHAT THIS MEANS. We currently have THREE study-level estimates, all "
        "UNRESOLVED.\nThe atlas carries %d studies at this threshold, of which "
        "%d are new (%d new patients).\nThe gain is the number of independent "
        "estimates to meta-analyse, not the patient\ncount -- per-study "
        "intervals do not tighten because other studies exist (invariant 4).",
        int(at_threshold["n_studies"]), int(at_threshold["n_studies_new"]),
        int(at_threshold["n_patients_new"]),
    )

    # Two facts that decide whether any of this is usable, checked rather than
    # assumed. Integrated expression would void invariant 4 outright.
    raw_share = float((obs["matrix_type"].astype(str) == "raw counts").mean())
    genomes = obs["reference_genome"].astype(str).nunique()
    log.info("\nraw counts: %.1f%% of cells (integration would void invariant 4)",
             100 * raw_share)
    log.info("reference genomes across studies: %d -- the shared gene index has "
             "to absorb that", genomes)

    for frame, name in ((table, "icbi_premise_feasibility"),
                        (studies, "icbi_premise_candidate_studies")):
        path = write_versioned_table(
            frame, name, seed=DEFAULT_SEED,
            results_dir=args.results_dir, allow_dirty=args.allow_dirty,
            extra_meta={
                "source": "ICBI CRC atlas /obs, read over HTTP range requests",
                "epithelial_labels": list(EPITHELIAL_LABELS),
                "epithelial_label_note": (
                    "Cancer cell is counted as epithelium. Excluding it gives 7 "
                    "usable patients against 136 -- a nineteen-fold difference "
                    "from a label choice, so it is a stated decision."
                ),
                "already_read": list(ALREADY_READ),
                "raw_counts_share": raw_share,
                "n_reference_genomes": genomes,
                "what_this_answers": (
                    "Whether the ICBI atlas can supply enough INDEPENDENT "
                    "study-level estimates to meta-analyse the coexpression "
                    "premise, which is the one blocker that is a power problem "
                    "rather than a structural one. Sizing only -- no expression "
                    "data is read."
                ),
                "what_this_does_not_answer": (
                    "Whether the premise then resolves. More studies narrow a "
                    "meta-analytic interval; they do not touch the "
                    "decomposition's algebraic collapse or the survivorship "
                    "confound, and neither is a power problem."
                ),
            },
        )
        log.info("wrote %s", path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
