"""Re-read the adenoma contrasts on a scale that is comparable across genes.

    python -m src.reference.jobs.adenoma_specificity
    python -m src.reference.jobs.adenoma_specificity --deltas <path>

A THIN LOCAL READ over a committed per-patient table, in the same shape as
``coexpression_meta``. No cluster, no atlas: every quantity it needs -- both
arms' detection rates, both arms' cell counts, the CP10K means -- is already on
the committed rows, so the correction runs on a laptop from artifacts in git.

WHAT WAS WRONG WITH THE TABLE THIS REPLACES. ``icbi_adenoma_specificity.parquet``
compared six genes' DETECTION deltas to each other. Those genes sit at baseline
detection rates from 0.36 to 0.98, and the sensitivity of a proportion depends
on its baseline, so that table ranked genes substantially by how detectable they
were. Full argument in ``src/reference/detection_scale.py``.

It also only ever computed ``GUCA2A - X``. The conclusion drawn from it --
terminal differentiation down, intestinal identity retained -- is a claim about
where CDX2 sits relative to the CONTROLS, and no row said. It was read off the
fact that CDX2's delta (-0.075) looked small beside GUCA2A's (-0.174), which is
exactly the cross-gene magnitude comparison the scale does not support. Computed
properly on the detection statistic, ``CDX2 - KRT8`` EXCLUDES zero: it says the
opposite of the published claim.

WHAT SURVIVES. On the load-bearing scale the panel separates into two blocks --
{KRT8, ACTB, EPCAM, CDX2} and {MS4A12, GUCA2A} -- mutually indistinguishable
within each and jointly separated across, with CDX2 no longer distinguishable
from housekeeping. That is the same conclusion the prose reached, reached by a
comparison that supports it, and it is a sharper shape than the gradient: a
gradient is what a uniform thinning also produces, and two blocks is not.
"""

from __future__ import annotations

import argparse
import logging
import sys
from collections.abc import Sequence
from pathlib import Path

import pandas as pd

from src.common.io import write_versioned_table
from src.common.paths import RESULTS_DIR
from src.common.provenance import DEFAULT_SEED
from src.reference.detection_scale import (
    UNIFORM_THINNING_R2_CEILING,
    delta_cloglog,
    uniform_thinning_null,
)
from src.reference.jobs.icbi_coexpression import (
    LOAD_BEARING_STATISTIC,
    SPECIFICITY_PANEL,
    SPECIFICITY_STATISTICS,
    contrast_matrix,
    specificity,
)

log = logging.getLogger(__name__)


def newest(name: str) -> Path | None:
    matches = sorted(RESULTS_DIR.glob(f"*/{name}.parquet"))
    return matches[-1] if matches else None


def disagreements(contrasts: pd.DataFrame) -> pd.DataFrame:
    """Contrasts where the statistics do not return the same verdict.

    THE POINT OF REPORTING THREE. If they always agreed the choice would not
    matter and this table would be empty; it is not empty, and the rows in it
    are exactly the claims that depend on which statistic was used. A reader who
    quotes one number from a table like this should be able to see, in the same
    artifact, whether it would have survived a different scale.
    """
    if contrasts.empty:
        return pd.DataFrame()
    wide = contrasts.pivot_table(
        index=["granularity_rung", "contrast", "role_of_gene", "role_of_other"],
        columns="statistic", values="excludes_zero", aggfunc="first",
    )
    disputed = wide[wide.nunique(axis=1) > 1].reset_index()
    if disputed.empty:
        return disputed
    # A - B and B - A are the same disagreement seen twice; the matrix needs
    # both directions, a list of claims does not.
    pair = disputed["contrast"].str.split(" - ", expand=True)
    disputed = disputed[pair[0] < pair[1]].reset_index(drop=True)
    means = contrasts.pivot_table(
        index=["granularity_rung", "contrast"],
        columns="statistic", values="mean_difference", aggfunc="first",
    ).add_suffix("_mean").reset_index()
    return disputed.merge(means, on=["granularity_rung", "contrast"], how="left")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--deltas", type=Path, default=None,
        help="per-patient deltas. Defaults to the newest icbi_adenoma.",
    )
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--name", default="adenoma_specificity")
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    path = args.deltas or newest("icbi_adenoma")
    if path is None:
        raise SystemExit(
            "no results/*/icbi_adenoma.parquet. Run the adenoma job first, or "
            "pass --deltas."
        )
    deltas = pd.read_parquet(path)
    log.info("%s | %d rows, %d patients, rungs %s", path.parent.name + "/" + path.name,
             len(deltas), deltas["patient_id"].nunique(),
             sorted(deltas["granularity_rung"].unique()))

    deltas["delta_cloglog"] = delta_cloglog(deltas)

    log.info("\n%s\nTHE SCALE THIS RESTS ON, fixed before the numbers\n%s",
             "=" * 72, "=" * 72)
    for name, spec in SPECIFICITY_STATISTICS.items():
        log.info("  %-11s %-22s %s", name, spec["standing"], spec["scale"])
    log.info("  panel: %s", ", ".join(SPECIFICITY_PANEL))

    contrasts = specificity(deltas, seed=args.seed)
    for rung in sorted(contrasts["granularity_rung"].unique()):
        n = int(contrasts.loc[contrasts["granularity_rung"] == rung,
                              "n_patients"].max())
        log.info("\n%s\nPAIRWISE, %s, n=%d -- %s (load-bearing)\n%s",
                 "=" * 72, rung, n, LOAD_BEARING_STATISTIC, "=" * 72)
        log.info("rows minus columns; negative = row fell more; * excludes zero")
        log.info("%s", contrast_matrix(contrasts, rung).to_string())
        log.info("\n-- the same contrasts on detection (diagnostic only) --")
        log.info("%s", contrast_matrix(contrasts, rung, "detection").to_string())

    null_frames, null_verdicts = [], {}
    for rung, block in deltas.groupby("granularity_rung", observed=True):
        table, verdict = uniform_thinning_null(block)
        null_verdicts[str(rung)] = verdict
        if not table.empty:
            null_frames.append(table.assign(granularity_rung=rung))
        log.info("\n%s\nCOULD ONE UNIFORM THINNING PRODUCE THE GRADIENT? [%s]\n%s",
                 "=" * 72, rung, "=" * 72)
        log.info("  %s", verdict["verdict"])
        log.info("  %s", verdict["detail"])
        if not table.empty:
            log.info("%s", table.to_string(index=False))
    thinning = (pd.concat(null_frames, ignore_index=True)
                if null_frames else pd.DataFrame())

    disputed = disagreements(contrasts)
    log.info("\n%s\nWHERE THE STATISTICS DISAGREE\n%s", "=" * 72, "=" * 72)
    if disputed.empty:
        log.info("  nowhere -- every contrast returns the same verdict on all three.")
    else:
        log.info("%s", disputed.to_string(index=False))
        log.info(
            "\nThese are the claims that depend on which statistic was used. "
            "The reading\nrests on %r; the others are shown so that dependence "
            "is visible rather than\navailable to whoever quotes a number.",
            LOAD_BEARING_STATISTIC,
        )

    meta = {
        "source": str(path.parent.name + "/" + path.name),
        "specificity_statistics": SPECIFICITY_STATISTICS,
        "load_bearing_statistic": LOAD_BEARING_STATISTIC,
        "specificity_panel": list(SPECIFICITY_PANEL),
        "supersedes": (
            "results/2026-09-05_d869bdd/icbi_adenoma_specificity.parquet, which "
            "compared six genes' detection deltas to each other across baseline "
            "rates of 0.36 to 0.98, and reported only GUCA2A - X"
        ),
        "cross_gene_rule": (
            "every cross-gene statement is made on the load-bearing log "
            "fold-change scale; delta_detect is reported for comparability with "
            "earlier tables and is NOT read across genes"
        ),
        "uniform_thinning_null": null_verdicts,
        "uniform_thinning_r2_ceiling": UNIFORM_THINNING_R2_CEILING,
        "statistics_disagree_on": (
            disputed["contrast"].tolist() if not disputed.empty else []
        ),
        "exploratory": True,
        "pre_registered": False,
        "what_this_is_not": (
            "A new measurement. The cells were scored once, on the cluster; "
            "this re-reads the committed per-patient table on a scale that "
            "supports the comparison being made."
        ),
    }
    for frame, name in ((contrasts, args.name),
                        (thinning, f"{args.name}_thinning_null"),
                        (disputed, f"{args.name}_disagreements")):
        if frame.empty:
            continue
        written = write_versioned_table(
            frame, name, seed=args.seed, results_dir=args.results_dir,
            allow_dirty=args.allow_dirty, extra_meta=meta,
        )
        log.info("wrote %s", written)
    return 0


if __name__ == "__main__":
    sys.exit(main())
