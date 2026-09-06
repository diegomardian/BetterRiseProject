"""The decomposition on adenoma — the project's deliverable, finally computable.

    python -m src.reference.jobs.adenoma_decomposition
    python -m src.reference.jobs.adenoma_decomposition --deltas <path>

A THIN LOCAL READ, in the same shape as ``adenoma_specificity``. No cluster and
no atlas: the cells were scored once by ``icbi_coexpression --arms adenoma
--collect-fractions``, and every quantity ``decompose_cohort`` needs is on the
committed rows. Pre-registered in ``docs/prereg_adenoma_decomposition.md``;
read §5 and §6 before reading any number this produces.

WHY THIS CAN RUN AT ALL, WHEN THE SAME ESTIMATOR COULD NOT ON CARCINOMA. The
Kitagawa ratio collapses when a gene's surviving per-cell mean goes to zero::

    intrinsic / compositional = (f_N / Δf) × (m_T/m_N − 1)

The bracket → −1 and the ratio becomes a property of the cell fractions,
identical for every gene on the same labels. On GSE178341 that constant is
−5.85. **On adenoma the bracket runs −0.05 to −0.63 and nothing is near −1**, so
the genes are distinguishable by it. That is measured, committed, and re-checked
(``docs/NEXT_AVENUES.md`` §1a).

TWO DENOMINATORS, ON PURPOSE — prereg §3.1 as amended. The mature fraction's
denominator is the RESOLVED epithelium (open decision #14: a cell that could not
be labelled is not a cell measured to be immature). The unresolved cells are
excluded by a depth cut that removes exactly ``DEPTH_QUANTILE`` of the
epithelium **and is endogenous** — a tumour arm that has lost expression carries
fewer counts, so more of its cells fall below the target. Measured on a
synthetic cohort with identical arms: shares of 0.229 normal against 0.270
tumour with the effect, ~0.25/0.25 without it.

So the whole decomposition is run **twice**, once against each denominator, and
a contrast whose sign or zero-exclusion flips between them is reported as
**denominator-dependent** rather than claimed. The second denominator is the one
decision #14 rejected, and it costs nothing to carry.

TWO INTERVALS, ALSO ON PURPOSE — prereg §4 as amended. ``ci_low``/``ci_high``
carry ``bootstrap_over_patients``'s percentile band because
``docs/open_decisions.md`` #10 settled that slot. That band is **miscalibrated
at these n** — 5.9% at n=44 and 7.1% at n=20 against a nominal 5%, by a closed
form with no data in it (``docs/HANDOFF.md`` §3a). So a Student-t interval is
carried in a **companion table** keyed by ``KEY_COLUMNS``, and every cross-rung
or cross-gene claim is made on that one. It is a companion rather than two extra
columns because ``coerce_results`` refuses any column outside the frozen schema.
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
from src.common.paths import RESULTS_DIR
from src.common.provenance import DEFAULT_SEED
from src.estimator.kitagawa import (
    CUTPOINTS,
    attach_intrinsic_ci,
    bootstrap_over_patients,
    decompose_cohort,
)
from src.harness.positivity import (
    COMPOSITIONAL_CUTPOINTS,
    classify_compositional_estimability,
)
from src.reference.interval_calibration import (
    NOMINAL_ALPHA,
    expected_false_positive_rate,
    student_t_interval,
)
from src.reference.jobs.coexpression_silencing import AXIS
from src.schema import coerce_results, write_results

log = logging.getLogger(__name__)

#: The two denominators. `resolved` is open decision #14 and is what the schema
#: table carries; `all_epithelial` is the denominator #14 rejected, run as the
#: sensitivity arm because the exclusion driving the difference is endogenous.
DENOMINATORS: dict[str, str] = {
    "resolved": "",
    "all_epithelial": "_all_epithelial",
}
PRIMARY_DENOMINATOR = "resolved"


class DecompositionError(ValueError):
    """The per-patient table cannot be read as a decomposition input."""


def newest(name: str) -> Path | None:
    matches = sorted(RESULTS_DIR.glob(f"*/{name}.parquet"))
    return matches[-1] if matches else None


def summary_frame(deltas: pd.DataFrame, *, denominator: str) -> pd.DataFrame:
    """The ten columns ``decompose_cohort`` requires, for one denominator.

    Renames rather than derives. Every quantity is already on the committed
    rows; if one is missing that is a scoring run without ``--collect-fractions``
    and the fix is upstream, not a default here.
    """
    suffix = DENOMINATORS[denominator]
    needed = {
        "patient_id", "study_id", "gene", "granularity_rung",
        "cp10k_normal", "cp10k_tumour", "n_tumour",
        f"frac_mature_normal{suffix}", f"frac_mature_tumour{suffix}",
    }
    missing = sorted(needed - set(deltas.columns))
    if missing:
        raise DecompositionError(
            f"the per-patient table is missing {missing}. The compositional arm "
            f"needs a mature FRACTION per arm, which only exists if the scoring "
            f"run passed --collect-fractions. Re-run the adenoma job; do not "
            f"default a fraction here."
        )
    out = deltas.loc[:, sorted(needed)].copy()
    return out.rename(columns={
        "cp10k_normal": "mean_normal",
        "cp10k_tumour": "mean_tumour",
        "n_tumour": "n_cells_mature",
        f"frac_mature_normal{suffix}": "frac_mature_normal",
        f"frac_mature_tumour{suffix}": "frac_mature_tumour",
    }).assign(labeling_axis=AXIS)


def compositional_estimability(deltas: pd.DataFrame) -> pd.DataFrame:
    """The compositional arm's own estimability call. Prereg §3.3.

    A DIFFERENT QUESTION FROM THE INTRINSIC ONE and it must not be folded into
    it. The intrinsic term gates on ``n_cells_mature`` under `CUTPOINTS`, which
    is PROVISIONAL; the compositional term gates on ``n_cells_resolved`` under
    `COMPOSITIONAL_CUTPOINTS`, which is decision #22 and is not. A row can have
    an estimable compositional term and an unestimable intrinsic one, and at a
    starved rung that is the ordinary case.
    """
    if "n_cells_resolved_normal" not in deltas.columns:
        raise DecompositionError(
            "no n_cells_resolved_* on the rows, so the compositional term's "
            "estimability cannot be called. It is not the same question as the "
            "intrinsic term's and it may not be inherited from it."
        )
    resolved = deltas[["n_cells_resolved_normal", "n_cells_resolved_tumour"]].min(axis=1)
    return pd.DataFrame({
        "patient_id": deltas["patient_id"].astype(str),
        "gene": deltas["gene"].astype(str),
        "granularity_rung": deltas["granularity_rung"].astype(str),
        "n_cells_resolved": resolved.astype(int),
        "compositional_estimability": [
            classify_compositional_estimability(int(n)) for n in resolved
        ],
    })


def student_t_companion(
    summary: pd.DataFrame, split: pd.DataFrame, *, seed: int,
) -> pd.DataFrame:
    """Per-(study, gene, rung, axis, weighting, term) Student-t interval.

    THE CALIBRATED INTERVAL, in its own table because the schema is frozen.
    Over PATIENTS, not cells — invariant 5 — so it is the same estimand as
    ``bootstrap_over_patients``'s band, reached by an interval that is actually
    5% at these n. The percentile band's real rate is carried on every row so
    the two can be compared without leaving the table.
    """
    rows = []
    keys = ["study_id", "gene", "granularity_rung", "labeling_axis", "weighting"]
    for key, block in split.groupby(keys, observed=True, dropna=False):
        for term in ("compositional", "intrinsic", "interaction"):
            values = pd.to_numeric(block[term], errors="coerce").dropna().to_numpy()
            n = int(values.size)
            record = dict(zip(keys, key, strict=True)) | {
                "term": term, "n_patients": n,
                "mean": float(values.mean()) if n else float("nan"),
            }
            if n >= 2:
                lo, hi = student_t_interval(
                    values, rng=np.random.default_rng(seed), alpha=NOMINAL_ALPHA
                )
                record |= {"t_ci_low": lo, "t_ci_high": hi,
                           "excludes_zero": bool(lo > 0 or hi < 0)}
            else:
                # One patient has no standard error. Undefined, not zero.
                record |= {"t_ci_low": float("nan"), "t_ci_high": float("nan"),
                           "excludes_zero": False}
            record["percentile_band_false_positive_rate"] = (
                expected_false_positive_rate(n) if n >= 2 else float("nan")
            )
            rows.append(record)
    return pd.DataFrame(rows)


def denominator_disagreements(companions: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Where the two denominators do not return the same reading. Prereg §3.1.

    The gate this replaces was a threshold on how unevenly the depth cut fell
    between arms. That gated on a diagnostic rather than on the harm, and at a
    tolerance inside the null's own noise band. This compares the ANSWERS.
    """
    keys = ["study_id", "gene", "granularity_rung", "labeling_axis",
            "weighting", "term"]
    primary = companions[PRIMARY_DENOMINATOR]
    other = next(k for k in companions if k != PRIMARY_DENOMINATOR)
    merged = primary.merge(
        companions[other], on=keys, suffixes=("", f"_{other}"), how="inner"
    )
    if merged.empty:
        return merged
    flipped_call = merged["excludes_zero"] != merged[f"excludes_zero_{other}"]
    # A SIGN FLIP ONLY COUNTS WHERE BOTH SIGNS ARE ESTABLISHED. The first
    # version compared np.sign directly and flagged 46 contrasts on real data,
    # nearly all of them at `epithelial`, where the compositional term is
    # EXACTLY 0.0 under the resolved denominator (delta-f is identically zero at
    # that rung) and about -0.001 under the other. np.sign(0.0) is 0 and
    # np.sign(-0.001) is -1, so the comparison fired on a difference of one
    # thousandth between two numbers that are both indistinguishable from zero.
    #
    # That is a check firing on noise -- the same defect as the threshold gate
    # this function was written to REPLACE, reappearing in its replacement. If
    # a term's interval contains zero its sign is not established, so a flip is
    # not a finding.
    both_resolved = merged["excludes_zero"] & merged[f"excludes_zero_{other}"]
    flipped_sign = both_resolved & (
        np.sign(merged["mean"]) != np.sign(merged[f"mean_{other}"])
    )
    keep = flipped_sign | flipped_call
    disputed = merged.loc[keep].copy()
    disputed["denominator_dependent_on"] = np.where(
        flipped_call[keep], "zero-exclusion", "sign"
    )
    return disputed.reset_index(drop=True)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--deltas", type=Path, default=None,
                        help="per-patient rows. Defaults to newest icbi_adenoma.")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--n-boot", type=int, default=2000)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    path = args.deltas or newest("icbi_adenoma")
    if path is None:
        raise SystemExit("no results/*/icbi_adenoma.parquet")
    deltas = pd.read_parquet(path)
    log.info("%s | %d rows, %d patients, rungs %s", path.parent.name,
             len(deltas), deltas["patient_id"].nunique(),
             sorted(deltas["granularity_rung"].unique()))

    log.info("\n%s\nWHAT IS FIXED BEFORE THE NUMBERS — prereg §3, §4\n%s",
             "=" * 72, "=" * 72)
    log.info("  denominators   %s (primary: %s, open decision #14)",
             sorted(DENOMINATORS), PRIMARY_DENOMINATOR)
    log.info("  intrinsic gate %s  <- PROVISIONAL", CUTPOINTS)
    log.info("  composit. gate %s", COMPOSITIONAL_CUTPOINTS)
    log.info("  schema band    percentile (open_decisions #10); calibrated "
             "interval in the companion")

    splits, companions = {}, {}
    for denominator in DENOMINATORS:
        summary = summary_frame(deltas, denominator=denominator)
        split = decompose_cohort(summary)
        splits[denominator] = split
        companions[denominator] = student_t_companion(
            summary, split, seed=args.seed
        )
        log.info("\n%s\nDENOMINATOR: %s\n%s", "=" * 72, denominator, "=" * 72)
        for rung in sorted(split["granularity_rung"].unique()):
            block = split[(split["granularity_rung"] == rung)
                          & (split["weighting"] == "doubly_robust")]
            mix = block.drop_duplicates("patient_id")["estimability"].value_counts()
            log.info("  %-15s n=%2d patients | intrinsic estimability %s",
                     rung, block["patient_id"].nunique(), dict(mix))

    primary_split = splits[PRIMARY_DENOMINATOR]
    bands = bootstrap_over_patients(
        summary_frame(deltas, denominator=PRIMARY_DENOMINATOR),
        n_boot=args.n_boot, seed=args.seed,
    )

    log.info("\n%s\nTHE SPLIT — %s denominator, doubly robust, Student-t\n%s",
             "=" * 72, PRIMARY_DENOMINATOR, "=" * 72)
    show = companions[PRIMARY_DENOMINATOR]
    show = show[show["weighting"] == "doubly_robust"]
    for rung in sorted(show["granularity_rung"].unique()):
        log.info("\n-- %s --", rung)
        log.info("%s", show[show["granularity_rung"] == rung][
            ["gene", "term", "n_patients", "mean", "t_ci_low", "t_ci_high",
             "excludes_zero"]].to_string(index=False))

    disputed = denominator_disagreements(companions)
    log.info("\n%s\nWHERE THE TWO DENOMINATORS DISAGREE\n%s", "=" * 72, "=" * 72)
    if disputed.empty:
        log.info("  nowhere — the reading does not depend on open decision #14.")
    else:
        log.info("  %d contrast(s) are DENOMINATOR-DEPENDENT and carry no "
                 "unqualified claim:", len(disputed))
        log.info("%s", disputed[["gene", "granularity_rung", "term", "weighting",
                                 "denominator_dependent_on"]].to_string(index=False))

    comp = compositional_estimability(deltas)

    meta = {
        "prereg": "docs/prereg_adenoma_decomposition.md",
        "source": f"{path.parent.name}/{path.name}",
        "denominators": {
            "resolved": "n_cells_epithelial - n_cells_unresolved (open decision #14, PRIMARY)",
            "all_epithelial": "n_cells_epithelial (the denominator #14 rejected; sensitivity arm)",
        },
        "why_two_denominators": (
            "the unresolved_depth exclusion is ENDOGENOUS -- an arm that has "
            "lost expression carries fewer counts, so more of its cells fall "
            "below the depth target. Measured on a synthetic cohort with "
            "identical arms: 0.229 normal / 0.270 tumour with the effect, "
            "~0.25/0.25 without. It did not propagate to the fraction there "
            "(0.500/0.502 against a true 0.500), and this run is where that is "
            "checked on real data instead of assumed from one fixture."
        ),
        "intrinsic_cutpoints": {"ok": CUTPOINTS.ok, "wide": CUTPOINTS.wide,
                                "source": CUTPOINTS.source},
        "compositional_cutpoints": {"ok": COMPOSITIONAL_CUTPOINTS.ok,
                                    "wide": COMPOSITIONAL_CUTPOINTS.wide,
                                    "source": COMPOSITIONAL_CUTPOINTS.source},
        "estimability_rule": (
            "two rules, reported separately, never folded: intrinsic gates on "
            "n_cells_mature (PROVISIONAL), compositional on n_cells_resolved "
            "(decision #22, not provisional)"
        ),
        "schema_interval": "percentile band, open_decisions #10",
        "companion_interval": (
            "Student-t, in adenoma_decomposition_t_intervals.parquet keyed by "
            "KEY_COLUMNS. The schema band's real false-positive rate is on "
            "every companion row -- 5.9% at n=44, 7.1% at n=20 (HANDOFF §3a). "
            "Every cross-rung or cross-gene claim is made on the companion."
        ),
        "denominator_dependent": (
            disputed[["gene", "granularity_rung", "term"]].to_dict("records")
            if not disputed.empty else []
        ),
        "exploratory": False,
        "pre_registered": True,
    }

    # attach_intrinsic_ci, NOT a hand-rolled merge. `bootstrap_over_patients`
    # is long-form by term, so merging it directly fans every patient row into
    # three and carries `term` and `n_boot` into a frozen schema that refuses
    # them. That function also encodes open_decisions #10's actual choice --
    # the schema's single ci slot carries the INTRINSIC term's band, because
    # estimability is defined for intrinsic and not for the other two.
    schema_frame = coerce_results(attach_intrinsic_ci(primary_split, bands))
    written = write_results(
        schema_frame, "adenoma_decomposition", seed=args.seed,
        results_dir=args.results_dir, allow_dirty=args.allow_dirty,
        notes=meta["companion_interval"],
    )
    log.info("\nwrote %s", written)

    for frame, name in (
        (companions[PRIMARY_DENOMINATOR], "adenoma_decomposition_t_intervals"),
        (companions["all_epithelial"], "adenoma_decomposition_t_intervals_all_epithelial"),
        (splits["all_epithelial"], "adenoma_decomposition_all_epithelial"),
        (disputed, "adenoma_decomposition_denominator_disagreements"),
        (comp, "adenoma_decomposition_compositional_estimability"),
    ):
        if frame.empty:
            continue
        log.info("wrote %s", write_versioned_table(
            frame, name, seed=args.seed, results_dir=args.results_dir,
            allow_dirty=args.allow_dirty, extra_meta=meta,
        ))
    return 0


if __name__ == "__main__":
    sys.exit(main())
