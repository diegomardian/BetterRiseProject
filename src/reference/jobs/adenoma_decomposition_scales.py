"""Does the two-block reading depend on a scale nobody pre-specified?

    python -m src.reference.jobs.adenoma_decomposition_scales

WHY THIS EXISTS, STATED PLAINLY. `docs/prereg_adenoma_decomposition.md` §5 fixed
the comparator SET — score every pair of the six-gene panel, not the target
against each — and it did **not** fix the cross-gene STATISTIC. That is a real
gap in that document and its RESULT section records it as one.

It matters because the decomposition's raw terms are in each gene's own CP10K
units. ACTB sits near 25 CP10K and MS4A12 near 3, so `compositional` and
`intrinsic` are not comparable across genes as magnitudes — the same error the
detection-scale correction was written to stop, one estimand over. A scale-free
statistic is needed, and the one the RESULT used (the intrinsic **share**) was
chosen **after** seeing the output.

THIS JOB DOES NOT REPAIR THAT BY PICKING A WINNER. Choosing a statistic now,
with the answer visible, would be the same act with a longer paper trail. What
it does instead is what `adenoma_specificity_disagreements.parquet` already does
for the choice of detection statistic: **run the contrasts on every defensible
scale-free construction and report where they disagree.** If the two blocks
survive all of them, the conclusion demonstrably does not rest on the choice,
and that is reportable without any claim to have pre-specified it. Where they
part, the parting is the finding and no unqualified claim follows.

THE THREE STATISTICS, and why each is defensible rather than convenient:

``share_abs``      |i| / (|i| + |c|). Bounded [0, 1], finite whenever either
                   term is non-zero, and the one the RESULT quoted. Discards
                   sign, so a gene whose terms oppose each other reads the same
                   as one whose terms agree.
``share_signed``   i / (|i| + |c|). Bounded [-1, 1] and keeps the direction
                   ``share_abs`` throws away. Differs from it exactly where the
                   intrinsic term is positive.
``ratio``          i / c, which is the form `docs/NEXT_AVENUES.md` §1a states
                   the collapse in and therefore the one the identifiability
                   claim is made on. Unbounded and heavy-tailed — a patient
                   whose compositional term is near zero dominates any mean —
                   so it is summarised by the MEDIAN with a rank-based interval,
                   never by a t-interval on the raw values.

A fourth was considered and rejected: ``intrinsic / total``. The total passes
through zero whenever the two terms cancel, so it is not merely heavy-tailed
but undefined in the middle of its own range.
"""

from __future__ import annotations

import argparse
import logging
import sys
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from src.common.io import write_versioned_table
from src.common.paths import RESULTS_DIR
from src.common.provenance import DEFAULT_SEED
from src.reference.interval_calibration import NOMINAL_ALPHA, student_t_interval

log = logging.getLogger(__name__)

#: Panel order: controls, epithelial, identity, then the two targets. The blocks
#: the reading claims are {KRT8, ACTB, EPCAM, CDX2} and {MS4A12, GUCA2A}.
PANEL: tuple[str, ...] = ("KRT8", "ACTB", "EPCAM", "CDX2", "MS4A12", "GUCA2A")
TARGET_BLOCK: frozenset[str] = frozenset({"MS4A12", "GUCA2A"})

#: LOAD-BEARING, pre-registered in docs/prereg_decomposition_statistic.md and
#: committed BEFORE it was computed on anything. Chosen because the project
#: already made this decision one estimand over: detection_scale.py moved to a
#: LOG FOLD CHANGE because a difference of two proportions is not comparable
#: across genes at different baselines, and the decomposition presents the same
#: problem. It is also the only candidate unbounded in both directions, where
#: the shares compress at exactly the end the two targets occupy.
LOAD_BEARING = "log_ratio"

STATISTICS: dict[str, str] = {
    "log_ratio": "log(|i| / |c|) — LOAD-BEARING, pre-registered before computing",
    "share_abs": "|i| / (|i| + |c|), bounded [0,1] — what the RESULT quoted",
    "share_signed": "i / (|i| + |c|), bounded [-1,1] — keeps the sign",
    "ratio": "i / c — NEXT_AVENUES §1a's form; median, rank-based interval",
}

#: `ratio` is summarised by the median because its mean is not a statistic: a
#: patient whose compositional term is near zero sends it anywhere. The
#: interval is the order-statistic (sign-test) interval for a median, which
#: assumes only that the per-patient values are exchangeable.
MEDIAN_STATISTICS: frozenset[str] = frozenset({"ratio"})


class ScaleError(ValueError):
    """A decomposition frame that cannot carry a scale-free comparison."""


def newest(name: str) -> Path | None:
    matches = sorted(RESULTS_DIR.glob(f"*/{name}.parquet"))
    return matches[-1] if matches else None


def scale_free(split: pd.DataFrame) -> pd.DataFrame:
    """Per (patient, gene, rung) values on each scale-free statistic."""
    needed = {"patient_id", "gene", "granularity_rung", "weighting",
              "intrinsic", "compositional"}
    missing = sorted(needed - set(split.columns))
    if missing:
        raise ScaleError(f"the decomposition frame is missing {missing}")

    out = split.copy()
    for column in ("intrinsic", "compositional"):
        out[column] = pd.to_numeric(out[column], errors="coerce")
    # `intrinsic` is None wherever estimability said not_estimable (invariant 1).
    # Those rows carry no intrinsic term and must drop out here rather than
    # become a zero in a numerator.
    out = out.dropna(subset=["intrinsic", "compositional"])

    magnitude = out["intrinsic"].abs() + out["compositional"].abs()
    finite = magnitude.where(magnitude > 0)
    out["share_abs"] = out["intrinsic"].abs() / finite
    out["share_signed"] = out["intrinsic"] / finite
    out["ratio"] = (out["intrinsic"] / out["compositional"].where(
        out["compositional"] != 0)).replace([np.inf, -np.inf], np.nan)
    # The load-bearing one. Undefined where either term is exactly zero -- those
    # rows drop out and are counted, never floored to a small constant. At the
    # `epithelial` rung the compositional term is exactly zero for every gene by
    # construction, so log_ratio is undefined there for the whole panel. That is
    # stated in the prereg rather than discovered here.
    out["log_ratio"] = np.log(
        out["intrinsic"].abs().where(out["intrinsic"] != 0)
        / out["compositional"].abs().where(out["compositional"] != 0)
    ).replace([np.inf, -np.inf], np.nan)
    return out


def _median_interval(values: np.ndarray, alpha: float = NOMINAL_ALPHA):
    """Order-statistic (sign-test) interval for a median.

    Distribution-free: it assumes the per-patient values are exchangeable and
    nothing about their shape, which is the whole reason `ratio` is summarised
    this way. Returns ``(nan, nan)`` when n is too small for any pair of order
    statistics to reach the level — undefined, not widened until it fits.
    """
    values = np.sort(np.asarray(values, dtype=float))
    n = values.size
    if n < 6:
        return float("nan"), float("nan")
    k = int(stats.binom.ppf(alpha / 2, n, 0.5))
    if k < 1:
        return float("nan"), float("nan")
    return float(values[k - 1]), float(values[n - k])


def contrasts(frame: pd.DataFrame, *, statistic: str, seed: int) -> pd.DataFrame:
    """Every ordered pair on one statistic, paired within patient."""
    rng = np.random.default_rng(seed)
    rows = []
    for (rung, weighting), block in frame.groupby(
        ["granularity_rung", "weighting"], observed=True
    ):
        wide = block.pivot_table(index="patient_id", columns="gene",
                                 values=statistic)
        present = [g for g in PANEL if g in wide.columns]
        for a in present:
            for b in present:
                if a == b:
                    continue
                paired = (wide[a] - wide[b]).dropna().to_numpy()
                n = int(paired.size)
                record = {
                    "granularity_rung": rung, "weighting": weighting,
                    "statistic": statistic, "contrast": f"{a} - {b}",
                    "gene": a, "other": b, "n_patients": n,
                    "cross_block": (a in TARGET_BLOCK) != (b in TARGET_BLOCK),
                }
                if statistic in MEDIAN_STATISTICS:
                    lo, hi = _median_interval(paired)
                    record["centre"] = float(np.median(paired)) if n else float("nan")
                    record["summary"] = "median"
                elif n >= 2:
                    lo, hi = student_t_interval(paired, rng=rng)
                    record["centre"] = float(paired.mean())
                    record["summary"] = "mean"
                else:
                    lo = hi = float("nan")
                    record["centre"] = float("nan")
                    record["summary"] = "mean"
                record |= {
                    "ci_low": lo, "ci_high": hi,
                    "excludes_zero": bool(np.isfinite(lo) and np.isfinite(hi)
                                          and (lo > 0 or hi < 0)),
                }
                rows.append(record)
    return pd.DataFrame(rows)


def disagreements(all_contrasts: pd.DataFrame) -> pd.DataFrame:
    """Contrasts whose verdict is not the same on every statistic.

    THE POINT OF RUNNING THREE. If they always agreed the choice would not
    matter and the §5 gap would be harmless. Rows here are exactly the claims
    that depend on a statistic nobody pre-specified.
    """
    wide = all_contrasts.pivot_table(
        index=["granularity_rung", "weighting", "contrast", "cross_block"],
        columns="statistic", values="excludes_zero", aggfunc="first",
    )
    disputed = wide[wide.nunique(axis=1) > 1].reset_index()
    if disputed.empty:
        return disputed
    # A - B and B - A are one disagreement seen twice.
    pair = disputed["contrast"].str.split(" - ", expand=True)
    return disputed[pair[0] < pair[1]].reset_index(drop=True)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--split", type=Path, default=None)
    parser.add_argument("--weighting", default="doubly_robust")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    path = args.split or newest("adenoma_decomposition")
    if path is None:
        raise SystemExit("no results/*/adenoma_decomposition.parquet")
    split = pd.read_parquet(path)
    log.info("%s | %d rows", path.parent.name, len(split))
    log.info("\nTHE GAP THIS ADDRESSES: prereg §5 fixed the comparator SET and "
             "not the\nSTATISTIC, and the RESULT's share was chosen after "
             "seeing the output. This does\nNOT pick a winner — it reports "
             "whether the reading depends on the choice.\n")
    for name, why in STATISTICS.items():
        log.info("  %-13s %s", name, why)

    values = scale_free(split)
    everything = pd.concat(
        [contrasts(values, statistic=s, seed=args.seed) for s in STATISTICS],
        ignore_index=True,
    )

    shown = everything[everything["weighting"] == args.weighting]
    for rung in ("lineage", "best4"):
        block = shown[shown["granularity_rung"] == rung]
        if block.empty:
            continue
        log.info("\n%s\n%s — the two-block claim on each statistic\n%s",
                 "=" * 72, rung, "=" * 72)
        for statistic in STATISTICS:
            s = block[block["statistic"] == statistic]
            cross = s[s["cross_block"]]
            within = s[~s["cross_block"]]
            # Ordered pairs, so each unordered contrast appears twice.
            log.info("  %-13s cross-block %2d/%2d exclude zero | "
                     "within-block %2d/%2d",
                     statistic,
                     int(cross["excludes_zero"].sum()) // 2, len(cross) // 2,
                     int(within["excludes_zero"].sum()) // 2, len(within) // 2)
        gm = block[(block["contrast"] == "GUCA2A - MS4A12")]
        for _, r in gm.iterrows():
            log.info("    GUCA2A - MS4A12 on %-13s %s %+.3f [%+.3f, %+.3f] %s",
                     r["statistic"], r["summary"], r["centre"],
                     r["ci_low"], r["ci_high"],
                     "EXCLUDES ZERO" if r["excludes_zero"] else "contains zero")

    defined = values.groupby("granularity_rung")["log_ratio"].apply(
        lambda v: f"{v.notna().sum()}/{len(v)}")
    log.info("\n  log_ratio defined on (gene x patient) rows per rung: %s",
             defined.to_dict())

    disputed = disagreements(everything)
    log.info("\n%s\nWHERE THE STATISTICS DISAGREE\n%s", "=" * 72, "=" * 72)
    if disputed.empty:
        log.info("  Nowhere. Every contrast returns the same verdict on all "
                 "three, so the\n  reading does not depend on the statistic the "
                 "prereg failed to fix.")
    else:
        cross = disputed[disputed["cross_block"]]
        log.info("  %d contrast(s) disagree, %d of them CROSS-BLOCK:",
                 len(disputed), len(cross))
        log.info("%s", disputed.to_string(index=False))
        log.info("\n  A cross-block disagreement would mean the two-block "
                 "reading itself rests\n  on the choice of statistic. A "
                 "within-block one does not.")

    meta = {
        "purpose": (
            "prereg §5 fixed the comparator SET and not the STATISTIC, and the "
            "RESULT's intrinsic share was chosen after seeing the output. This "
            "reports whether the two-block reading depends on that choice. It "
            "deliberately does NOT designate a load-bearing statistic: doing so "
            "now, with the answer visible, would be the same act with a longer "
            "paper trail."
        ),
        "source": f"{path.parent.name}/{path.name}",
        "statistics": STATISTICS,
        "load_bearing": LOAD_BEARING,
        "load_bearing_prereg": "docs/prereg_decomposition_statistic.md",
        "median_statistics": sorted(MEDIAN_STATISTICS),
        "rejected_statistic": (
            "intrinsic / total — the total passes through zero whenever the two "
            "terms cancel, so it is undefined in the middle of its own range"
        ),
        "target_block": sorted(TARGET_BLOCK),
        "cross_block_disagreements": (
            disputed.loc[disputed["cross_block"], "contrast"].tolist()
            if not disputed.empty else []
        ),
        "exploratory": True,
        "pre_registered": False,
        "standing": (
            "The comparator SET is confirmatory (prereg §5). The SCALE is not, "
            "and this table is what stands in for that: agreement across every "
            "defensible construction, rather than a pre-specification that does "
            "not exist."
        ),
    }
    for frame, name in ((everything, "adenoma_decomposition_scales"),
                        (disputed, "adenoma_decomposition_scale_disagreements")):
        if frame.empty:
            continue
        log.info("wrote %s", write_versioned_table(
            frame, name, seed=args.seed, results_dir=args.results_dir,
            allow_dirty=args.allow_dirty, extra_meta=meta,
        ))
    return 0


if __name__ == "__main__":
    sys.exit(main())
