"""Is the adenoma result driven by one specimen collection? A split, not a replication.

    python -m src.reference.jobs.disval_stability

Pre-registered in ``docs/prereg_disval_stability.md``, committed before this ran.
**Read §1 before reading any number here**, because the thing this is most
likely to be misquoted as is the thing it explicitly is not.

WHAT IT IS NOT. `Chen_2021_Cell`'s pooled 44-patient answer is already known and
committed. **A split of analysed data cannot confirm that data.** Same lab, same
platform, same COLON MAP population, same cells. The "single-cohort" qualifier
on the adenoma result stays, and nothing produced here licenses removing it.

WHAT IT IS. A test of exactly one failure mode: that the finding is an artefact
of a single collection batch. `Chen_2021_Cell` is four specimen collections, and
discovery and validation are independent collections about a year apart whose
patient sets are disjoint but for one. Batch-drivenness is a real and common way
for a single-cohort result to be wrong, it is the only such way this data can
address, and no other available analysis addresses it.

WHY IT EXISTS NOW. `docs/prereg_becker_replication.md` was the designated
closing path for avenue A's post-hoc-statistic problem, and its Amendment 2
records that Becker's paired cohort is **four donors**. At n=4 this project's
own interval fires 18.8% under a true null. So the statistic gap needed a
closing path that does not depend on Becker; this is the strongest one the
committed data supports, and it needs no cluster and no download.

THE ASYMMETRY THAT GOVERNS HOW THIS IS READ, from §5 and §6. At n=15 the
Student-t half-width is **1.80x** what it is at n=43, so a contrast within a
factor of ~1.8 of zero can miss for width alone. **A half-miss is therefore weak
evidence and a sign reversal is strong** — width does not flip signs. That is
fixed in the pre-registration so a half-failure cannot later be read as a
refutation, nor a half-success as independent confirmation.
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
from src.common.paths import INTERIM_DIR, RESULTS_DIR
from src.common.provenance import DEFAULT_SEED
from src.reference.interval_calibration import student_t_interval
from src.reference.jobs.adenoma_decomposition_scales import (
    LOAD_BEARING,
    MEDIAN_STATISTICS,
    STATISTICS,
    TARGET_BLOCK,
    _median_interval,
    scale_free,
)

log = logging.getLogger(__name__)

STUDY = "Chen_2021_Cell"

#: The specimen collections, from the atlas's `dataset` column. `_CRC` scores no
#: adenoma patient and is listed so its absence is a recorded fact.
COLLECTIONS: tuple[str, ...] = (
    "VUMC_HTAN_discovery", "VUMC_HTAN_validation", "VUMC_HTAN_cohort3",
)
CARCINOMA_COLLECTION = "VUMC_HTAN_CRC"

#: The two collections §1 calls independent — disjoint patient sets, about a
#: year apart. `cohort3` is carried as a third set, never merged into either.
PRIMARY_HALVES: tuple[str, str] = ("VUMC_HTAN_discovery", "VUMC_HTAN_validation")

#: The claim under test: GUCA2A against every member of the other block.
#: Pre-registered §4 as the four contrasts that survived the agreement rule.
PRIMARY_GENE = "GUCA2A"

RUNG = "lineage"
WEIGHTING = "doubly_robust"


class StabilityError(ValueError):
    """The split cannot be formed as pre-registered."""


def newest(name: str) -> Path | None:
    matches = sorted(RESULTS_DIR.glob(f"*/{name}.parquet"))
    return matches[-1] if matches else None


def collection_of(obs_cache: Path) -> pd.DataFrame:
    """Patient -> specimen collection, from the cached obs.

    A patient appearing in two collections is returned with both, joined by
    ``|``. That is not tidied here: §3 excludes such a patient from both halves
    and the caller has to see it to do that.

    The carcinoma collection is dropped BEFORE grouping. Three patients span
    `_CRC` and `cohort3`, and folding that in would make them look like
    cross-collection patients when they are not.
    """
    obs = pd.read_parquet(obs_cache, columns=["study_id", "patient_id", "dataset"])
    rows = obs[obs["study_id"] == STUDY].copy()
    if rows.empty:
        raise StabilityError(f"no {STUDY} cells in {obs_cache}")
    rows["dataset"] = rows["dataset"].astype(str)
    rows = rows[rows["dataset"] != CARCINOMA_COLLECTION]

    grouped = (rows.groupby("patient_id")["dataset"]
               .agg(lambda v: "|".join(sorted(set(v)))))
    unknown = sorted({d for value in grouped for d in value.split("|")}
                     - set(COLLECTIONS))
    if unknown:
        raise StabilityError(
            f"unknown collection(s) {unknown}. The split is defined on "
            f"{list(COLLECTIONS)}; a value nobody has seen means the atlas is "
            f"not the one this was verified against."
        )
    return grouped.rename("collection").reset_index()


def assign_halves(split: pd.DataFrame, scored: Sequence[str]) -> pd.DataFrame:
    """One row per scored patient with its collection, and the exclusions named.

    **A patient in more than one collection is excluded from every half**, per
    §3: assigning it to either would put the same patient on both sides of a
    comparison the whole design requires to be disjoint.
    """
    scored = pd.Index([str(p) for p in scored])
    lookup = dict(zip(split["patient_id"].astype(str), split["collection"],
                      strict=True))
    rows = []
    for patient in scored:
        collection = lookup.get(patient) or lookup.get(f"{STUDY}.{patient}")
        rows.append({
            "patient_id": patient,
            "collection": collection,
            "shared": bool(collection and "|" in collection),
            "found": collection is not None,
        })
    frame = pd.DataFrame(rows)
    missing = frame.loc[~frame["found"], "patient_id"].tolist()
    if missing:
        raise StabilityError(
            f"{len(missing)} scored patient(s) have no collection: {missing[:5]}. "
            f"The atlas writes '{STUDY}.HTA11_866' and the decomposition writes "
            f"'HTA11_866'; if neither form matches this is an identifier-space "
            f"mismatch rather than a patient with no collection."
        )
    return frame


def contrasts_within(
    values: pd.DataFrame, *, statistic: str, seed: int, label: str
) -> pd.DataFrame:
    """GUCA2A against every non-target gene, paired within patient."""
    rng = np.random.default_rng(seed)
    wide = values.pivot_table(index="patient_id", columns="gene", values=statistic)
    if PRIMARY_GENE not in wide.columns:
        return pd.DataFrame()
    rows = []
    for other in sorted(set(wide.columns) - TARGET_BLOCK):
        paired = (wide[PRIMARY_GENE] - wide[other]).dropna().to_numpy()
        n = int(paired.size)
        # `ratio` is summarised by the MEDIAN with a rank-based interval, the
        # same treatment adenoma_decomposition_scales gives it and for the same
        # reason: it is heavy-tailed, so a patient whose compositional term is
        # near zero sends its mean anywhere. Using a t-interval on it here was
        # an inconsistency with the pre-registered handling, and it produced
        # the only two sign flips in the first run — an artefact of the
        # summary, not a property of the data.
        if statistic in MEDIAN_STATISTICS:
            lo, hi = _median_interval(paired)
            centre = float(np.median(paired)) if n else float("nan")
            summary = "median"
        elif n >= 2:
            lo, hi = student_t_interval(paired, rng=rng)
            centre = float(paired.mean())
            summary = "mean"
        else:
            lo = hi = centre = float("nan")
            summary = "mean"
        rows.append({
            "half": label, "statistic": statistic, "summary": summary,
            "contrast": f"{PRIMARY_GENE} - {other}", "other": other,
            "n_patients": n,
            "mean": centre,
            "ci_low": lo, "ci_high": hi,
            "excludes_zero": bool(np.isfinite(lo) and (lo > 0 or hi < 0)),
        })
    return pd.DataFrame(rows)


def verdict(table: pd.DataFrame) -> dict:
    """§5's branches, taken against the table rather than by a reader."""
    load = table[(table["statistic"] == LOAD_BEARING)
                 & (table["half"].isin(PRIMARY_HALVES))]
    if load.empty:
        return {"verdict": "NOT COMPUTED", "detail": "no load-bearing rows"}

    per_half = load.groupby("half")["excludes_zero"].sum().to_dict()
    n_contrasts = load["contrast"].nunique()
    pooled_sign = np.sign(load.groupby("contrast")["mean"].mean())
    reversed_in = []
    for half in PRIMARY_HALVES:
        block = load[load["half"] == half].set_index("contrast")
        for contrast, mean in block["mean"].items():
            resolved = bool(block.loc[contrast, "excludes_zero"])
            if resolved and np.sign(mean) != pooled_sign.get(contrast, np.sign(mean)):
                reversed_in.append((half, contrast))

    if reversed_in:
        return {
            "verdict": "SIGN REVERSAL",
            "detail": (
                f"{reversed_in} reverse sign with an interval excluding zero. "
                f"**Width does not flip signs** (§6), so this is not a power "
                f"artefact: the pooled result would be a mixture of opposing "
                f"collections and the adenoma reading withdraws to the "
                f"m_T/m_N ratio table."
            ),
        }
    if all(per_half.get(h, 0) == n_contrasts for h in PRIMARY_HALVES):
        return {
            "verdict": "BATCH-DRIVENNESS EXCLUDED",
            "detail": (
                f"all {n_contrasts} contrasts exclude zero in both halves. The "
                f"result is not an artefact of one specimen collection. It "
                f"remains single-cohort, single-lab, single-platform — §1, and "
                f"this does not license removing that qualifier."
            ),
        }
    return {
        "verdict": "AMBIGUOUS AT THIS N",
        "detail": (
            f"holds in {per_half} of {n_contrasts}. §5 fixed this as ambiguous "
            f"BY DESIGN and it must not be read as a failure: at n=15 the "
            f"interval is 1.80x its width at n=43, so a real effect can miss. "
            f"Report both halves and the pooled result; claim nothing new."
        ),
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--split", type=Path, default=None)
    parser.add_argument("--obs-cache", type=Path,
                        default=INTERIM_DIR / "icbi_obs.parquet")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    path = args.split or newest("adenoma_decomposition")
    if path is None:
        raise SystemExit("no results/*/adenoma_decomposition.parquet")
    decomposition = pd.read_parquet(path)
    block = decomposition[(decomposition["granularity_rung"] == RUNG)
                          & (decomposition["weighting"] == WEIGHTING)]

    log.info("%s\nDIS/VAL STABILITY — a split, NOT a replication\n%s",
             "=" * 72, "=" * 72)
    log.info("  pre-registered in docs/prereg_disval_stability.md — read §1")
    log.info("  the same cells produced the pooled answer, so this CANNOT")
    log.info("  confirm it. It tests one failure mode: batch-drivenness.")
    log.info("  source: %s | %s, %s", path.parent.name, RUNG, WEIGHTING)

    assignment = assign_halves(collection_of(args.obs_cache),
                               block["patient_id"].unique())
    shared = assignment.loc[assignment["shared"], "patient_id"].tolist()
    usable = assignment[~assignment["shared"]]

    log.info("\n%s\nTHE SPLIT\n%s", "=" * 72, "=" * 72)
    log.info("%s", usable["collection"].value_counts().to_string())
    log.info("  excluded as cross-collection (§3): %s", shared or "none")

    values = scale_free(block).merge(
        usable[["patient_id", "collection"]], on="patient_id", how="inner")

    frames = []
    for label in (*PRIMARY_HALVES, "VUMC_HTAN_cohort3", "POOLED"):
        subset = values if label == "POOLED" else values[values["collection"] == label]
        if subset.empty:
            continue
        for statistic in STATISTICS:
            frames.append(contrasts_within(subset, statistic=statistic,
                                           seed=args.seed, label=label))
    table = pd.concat([f for f in frames if not f.empty], ignore_index=True)

    log.info("\n%s\nGUCA2A vs EACH CONTROL — %s (load-bearing, ac7eca1)\n%s",
             "=" * 72, LOAD_BEARING, "=" * 72)
    load = table[table["statistic"] == LOAD_BEARING]
    for label in (*PRIMARY_HALVES, "VUMC_HTAN_cohort3", "POOLED"):
        rows = load[load["half"] == label]
        if rows.empty:
            continue
        n = int(rows["n_patients"].max())
        log.info("\n-- %s (n=%d) --", label, n)
        log.info("%s", rows[["contrast", "mean", "ci_low", "ci_high",
                             "excludes_zero"]].to_string(index=False))

    outcome = verdict(table)
    log.info("\n%s\nVERDICT\n%s", "=" * 72, "=" * 72)
    log.info("  %s", outcome["verdict"])
    log.info("  %s", outcome["detail"])
    log.info(
        "\n  READ §6 BEFORE QUOTING A MISS. At n=15 the interval is 1.80x its "
        "width at\n  n=43, so a half-miss is WEAK evidence and a sign reversal "
        "is STRONG. That\n  asymmetry was fixed before these numbers existed.")

    meta = {
        "prereg": "docs/prereg_disval_stability.md",
        "source": f"{path.parent.name}/{path.name}",
        "what_this_is_not": (
            "an independent replication. Same lab, same platform, same "
            "population, and the same cells whose pooled answer is already "
            "known — a split of analysed data cannot confirm it. The "
            "single-cohort qualifier on the adenoma result STAYS."
        ),
        "what_this_is": (
            "a test of one failure mode: that the finding is an artefact of a "
            "single collection batch"
        ),
        "load_bearing_statistic": LOAD_BEARING,
        "statistic_fixed_in": "ac7eca1, before any of this",
        "primary_halves": list(PRIMARY_HALVES),
        "excluded_cross_collection": shared,
        "half_sizes": usable["collection"].value_counts().to_dict(),
        "width_penalty_vs_n43": {"n=15": 1.80, "n=13": 1.96},
        "reading_asymmetry": (
            "a half-miss is weak evidence (width), a sign reversal is strong "
            "(width does not flip signs). Fixed in §5/§6 before the numbers."
        ),
        "verdict": outcome,
        "exploratory": False,
        "pre_registered": True,
    }
    for frame, name in ((table, "disval_stability"),
                        (assignment, "disval_stability_assignment")):
        log.info("wrote %s", write_versioned_table(
            frame, name, seed=args.seed, results_dir=args.results_dir,
            allow_dirty=args.allow_dirty, extra_meta=meta,
        ))
    return 0


if __name__ == "__main__":
    sys.exit(main())
