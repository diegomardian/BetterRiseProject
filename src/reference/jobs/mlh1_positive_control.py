"""The instrument's first positive control: can it see a known silencing event?

    python -m src.reference.jobs.mlh1_positive_control --atlas <path>

WHAT THIS IS FOR. Every null this project has produced -- UNRESOLVED on three
cohorts, UNRESOLVED at thirteen studies, NOT SPECIFIC on the bulk CIMP screen,
not gene-specific on adenoma -- rests on a detection statistic whose sensitivity
to real transcriptional silencing has never been demonstrated. A negative result
from an instrument of unknown sensitivity is not evidence of absence. It is not
evidence of anything.

MLH1 is the one gene in this cohort where silencing is known from **an assay
rather than from expression**. Pelka's atlas rows carry
``MLH1_promoter_methylation_status`` for all 62 patients, patient-level,
from methylation calling; promoter hypermethylation silencing MLH1
transcription is established biology, not a hypothesis of this project. So
asking whether the statistic sees MLH1 fall in the mature cells of methylated
patients tests **the instrument, not the biology** -- and the answer reframes
every null above, in one direction or the other.

WHAT THIS IS NOT, AND THE THING THAT CHANGED. It is NOT the difference-in-
differences that ``docs/prereg_g2_mlh1.md`` pre-registered. That design needed a
mechanistic negative control -- ``mlh1_intact_mmrd``, patients who reach the
same MSI-H phenotype through MSH2/MSH6/PMS2 with MLH1 transcription untouched --
and the feasibility check says it does not exist at usable size. Only 29 of
Pelka's 62 patients survive this pipeline's own filters, and the intact-MMRd arm
is **four** of them. That is the same number the original prereg reached after
depth matching, arrived at here through a completely different route: the atlas
reprocessing rather than the GSE178341 path. Two pipelines agreeing on four says
it is a property of the cohort, not of anybody's filters.

Stratifying does not fix the dilution, because the stratum you would stratify
into has four patients in it.

SO THE CLAIM IS NARROWER AND IT STILL DISCRIMINATES. Within the methylated
stratum alone: does MLH1 fall in the mature cells of patients whose promoters
are known to be methylated? Pre-registered in
``docs/prereg_g2_mlh1_within_stratum.md``, with the interval, the power and the
falsifiers fixed before this ran.

THREE ARMS, WITH DIFFERENT STANDING, NEVER MERGED.

``mlh1_methylated`` (n=10)
    The reading. Powered.

``mlh1_unmethylated`` (n=19)
    Secondary and CONFOUNDED: it is mostly MMR-proficient patients, so it mixes
    methylation status with MSI status. Reported as specificity, never as the
    mechanistic control the prereg wanted, and the confound is on the row.

``mlh1_intact_mmrd`` (n=4)
    The pre-registered control, reported as UNDERPOWERED. At n=4 the project's
    own interval excludes zero 16-22% of the time under a true null, so this
    arm cannot support a verdict in either direction and none is taken from it.

THE INTERVAL IS NOT THE ONE THE REST OF THE REPOSITORY USES, ON PURPOSE. The
percentile bootstrap over patients is miscalibrated at these n -- 9.6% at n=10
against a nominal 5%, by a closed form with no data in it. This job reports the
Student-t interval, which is calibrated at every n measured. The measurement
that settles that is ``results/*/interval_calibration.parquet`` and it was
committed BEFORE this job ran; see ``src/reference/interval_calibration.py``.
"""

from __future__ import annotations

import argparse
import glob
import logging
import sys
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd

from src.common.io import write_versioned_table
from src.common.paths import INTERIM_DIR, RESULTS_DIR
from src.common.provenance import DEFAULT_SEED
from src.reference.detection_scale import cloglog_rate
from src.reference.interval_calibration import (
    CALIBRATED_METHOD,
    INTERVAL_METHODS,
    excludes_zero,
    expected_false_positive_rate,
)
from src.reference.jobs.coexpression_silencing import RUNG, premise_holds
from src.reference.jobs.icbi_coexpression import (
    DEFAULT_ATLAS,
    load_obs,
    study_deltas,
)

log = logging.getLogger(__name__)

#: The gene under test. Tier B of the frozen panel, and the only panel gene
#: whose silencing this cohort records from an assay.
TARGET = "MLH1"

STUDY = "Pelka_2021_Cell"

#: The arm the reading is about, and the two reported beside it.
PRIMARY_STRATUM = "mlh1_methylated"
SECONDARY_STRATUM = "mlh1_unmethylated"
UNDERPOWERED_STRATUM = "mlh1_intact_mmrd"

#: Patients whose NORMAL arm must carry at least one MLH1+ cell before the
#: reading is attempted.
#:
#: WHY THIS EXISTS AND WHY IT IS PRE-COMMITTED. ``cloglog_rate`` applies a
#: Jeffreys-style boundary correction so that a detection rate of exactly zero
#: is still transformable. That correction is what makes the statistic total,
#: and it is also what would let this job return a confident-looking number
#: computed almost entirely FROM the correction if MLH1 were undetected in most
#: patients. A delta between two boundary corrections is arithmetic about the
#: pseudocount, not a measurement of a gene.
#:
#: Five of ten is deliberately a majority. Below it the reading returns NOT
#: ESTIMABLE, which under CLAUDE.md invariant 1 is not the same as "no
#: silencing" and must never be written as zero.
MIN_PATIENTS_WITH_SIGNAL = 5


class ArmError(ValueError):
    """An arm that cannot be read as one."""


def strata_for(obs: pd.DataFrame) -> pd.DataFrame:
    """Per-patient methylation status, from the atlas AND from the prereg table.

    TWO SOURCES ON PURPOSE. ``docs/prereg_g2_mlh1.md`` §3 fixed the strata in
    week 0 from ``assign_mlh1_strata()`` over the clinical annotation; the atlas
    carries ``MLH1_promoter_methylation_status`` per cell, from the
    reprocessing. They are independent derivations of the same fact and they are
    compared here rather than one being trusted: on the 62 Pelka patients they
    agree exactly -- 22 ``meth`` against 22 ``mlh1_methylated``, no crossings.

    A disagreement would mean the arm this reading is about is not the arm the
    pre-registration named, which is the kind of thing that is invisible unless
    it is checked.
    """
    column = "MLH1_promoter_methylation_status"
    if column not in obs.columns:
        raise ArmError(
            f"{column} is not in the cached obs. It is what makes this a "
            f"positive control rather than another expression comparison; "
            f"without it there is no assay-derived truth to test against."
        )
    rows = obs[obs["study_id"] == STUDY].copy()
    rows["short_id"] = rows["patient_id"].astype(str).str.split(".").str[-1]
    rows[column] = rows[column].astype(str)

    records = []
    for patient, block in rows.groupby("short_id"):
        values = sorted(set(block[column][block[column] != "None"]))
        if len(values) > 1:
            raise ArmError(
                f"patient {patient} carries {values} for {column}. Methylation "
                f"status is a patient-level fact; more than one value means the "
                f"cells of one patient disagree about it and the arm cannot be "
                f"formed."
            )
        records.append({
            "short_id": patient,
            "atlas_methylation": values[0] if values else "NOT_ANNOTATED",
        })
    atlas = pd.DataFrame(records)

    cohort_paths = sorted(glob.glob(str(RESULTS_DIR / "*" / "cohort_table.parquet")))
    if not cohort_paths:
        raise ArmError("no results/*/cohort_table.parquet; the prereg strata are "
                       "the record and this reading is defined against them")
    cohort = pd.read_parquet(cohort_paths[-1])
    merged = atlas.merge(cohort[["patient_id", "mlh1_stratum", "matched"]],
                         left_on="short_id", right_on="patient_id", how="left")

    crossings = merged[
        ((merged["atlas_methylation"] == "meth")
         & (merged["mlh1_stratum"] != PRIMARY_STRATUM))
        | ((merged["atlas_methylation"] == "no_meth")
           & (merged["mlh1_stratum"] == PRIMARY_STRATUM))
    ]
    merged.attrs["crossings"] = crossings["short_id"].tolist()
    if len(crossings):
        log.warning(
            "  %d patient(s) where the atlas annotation and the pre-registered "
            "stratum disagree: %s. The reading is defined on the "
            "PRE-REGISTERED stratum; the disagreement is recorded, not resolved "
            "here.", len(crossings), crossings["short_id"].tolist())
    else:
        log.info("  atlas annotation and pre-registered strata agree on all "
                 "%d patients", len(merged))
    return merged


def arm_of(stratum: str) -> str:
    """Which reported arm a pre-registered stratum belongs to."""
    if stratum == PRIMARY_STRATUM:
        return PRIMARY_STRATUM
    if stratum == UNDERPOWERED_STRATUM:
        return UNDERPOWERED_STRATUM
    return SECONDARY_STRATUM


def arm_reading(
    deltas: pd.DataFrame, *, gene: str, arm: str, seed: int,
    method: str = CALIBRATED_METHOD,
) -> dict:
    """One gene in one arm: mean delta, interval, and what it may be read as.

    Returns a row rather than a verdict string, so that an arm the design calls
    underpowered still produces its numbers -- suppressing them would make the
    table's shape depend on the result, and a reader could not check the claim
    that the arm is uninformative.
    """
    block = deltas[deltas["gene"] == gene].drop_duplicates("patient_id")
    n = len(block)
    row = {
        "gene": gene, "arm": arm, "n_patients": n,
        "median_cells_per_arm": (float(block["n_tumour"].median()) if n else float("nan")),
        "detect_normal": float(block["detect_normal"].mean()) if n else float("nan"),
        "detect_tumour": float(block["detect_tumour"].mean()) if n else float("nan"),
        "patients_with_signal": (
            int((block["detect_normal"] > 0).sum()) if n else 0
        ),
        "interval_method": method,
        "percentile_false_positive_rate": (
            expected_false_positive_rate(n) if n >= 2 else float("nan")
        ),
    }
    if n < 2:
        row |= {"mean_delta_cloglog": float("nan"), "ci_low": float("nan"),
                "ci_high": float("nan"), "excludes_zero": False,
                "estimability": "not_estimable",
                "detail": f"{n} patient(s): no interval is defined"}
        return row

    values = np.asarray(
        cloglog_rate(block["detect_tumour"], block["n_tumour"])
        - cloglog_rate(block["detect_normal"], block["n_normal"])
    )
    rng = np.random.default_rng(seed)
    lo, hi = INTERVAL_METHODS[method](values, rng=rng)
    row |= {
        "mean_delta_cloglog": float(values.mean()),
        "ci_low": float(lo), "ci_high": float(hi),
        "excludes_zero": excludes_zero((lo, hi)),
        "estimability": "estimable",
        "detail": "",
    }
    return row


def instrument_verdict(primary: dict, premise: tuple[bool, str]) -> dict:
    """What the reading says about the instrument. Three outcomes, all real.

    The pre-committed consequences, from
    ``docs/prereg_g2_mlh1_within_stratum.md`` §5. Written as a function so the
    branch is taken by code against a table rather than by a reader against a
    paragraph.
    """
    held, why = premise
    if not held:
        return {
            "verdict": "UNINTERPRETABLE",
            "detail": (
                f"the premise does not hold in the methylated arm ({why}). The "
                f"cells compared are not established to be the same population, "
                f"so MLH1 moving or not moving inside them says nothing about "
                f"silencing. This is not a negative result."
            ),
        }
    if primary["estimability"] != "estimable":
        return {"verdict": "NOT ESTIMABLE", "detail": primary["detail"]}
    if primary["patients_with_signal"] < MIN_PATIENTS_WITH_SIGNAL:
        return {
            "verdict": "NOT ESTIMABLE",
            "detail": (
                f"{primary['patients_with_signal']} of {primary['n_patients']} "
                f"patients carry any MLH1+ cell in the normal arm, below "
                f"{MIN_PATIENTS_WITH_SIGNAL}. The delta would be mostly the "
                f"boundary correction rather than the data. Under invariant 1 "
                f"this is NOT zero silencing and must not be written as zero."
            ),
        }
    if primary["excludes_zero"] and primary["mean_delta_cloglog"] < 0:
        return {
            "verdict": "INSTRUMENT SEES KNOWN SILENCING",
            "detail": (
                f"MLH1 falls {primary['mean_delta_cloglog']:+.3f} "
                f"[{primary['ci_low']:+.3f}, {primary['ci_high']:+.3f}] in the "
                f"mature cells of patients whose promoters are methylated. The "
                f"statistic has demonstrated sensitivity to a silencing event "
                f"known from an assay, so this project's nulls become evidence "
                f"of absence rather than absence of evidence -- WITHOUT a "
                f"mechanistic control arm behind them, which this cohort "
                f"cannot supply."
            ),
        }
    if primary["excludes_zero"]:
        return {
            "verdict": "WRONG DIRECTION",
            "detail": (
                f"MLH1 RISES {primary['mean_delta_cloglog']:+.3f} "
                f"[{primary['ci_low']:+.3f}, {primary['ci_high']:+.3f}] where "
                f"methylation predicts a fall. §5 names this a falsifier: an "
                f"instrument that fires in the wrong direction on a known event "
                f"is not a calibrated instrument, and the nulls are not "
                f"rehabilitated by it."
            ),
        }
    return {
        "verdict": "INSTRUMENT DOES NOT SEE IT",
        "detail": (
            f"MLH1 {primary['mean_delta_cloglog']:+.3f} "
            f"[{primary['ci_low']:+.3f}, {primary['ci_high']:+.3f}], containing "
            f"zero, in patients whose promoters are methylated. At the "
            f"pre-registered power this is informative against STRONG silencing "
            f"and NOT informative against moderate silencing -- read the power "
            f"table before reading this as a null. Its consequence is that the "
            f"project's other nulls stay uninformative: an instrument that "
            f"cannot see a known event cannot be cited for not seeing an "
            f"unknown one."
        ),
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", type=Path, default=DEFAULT_ATLAS)
    parser.add_argument("--obs-cache", type=Path,
                        default=INTERIM_DIR / "icbi_obs.parquet")
    parser.add_argument("--rung", default=RUNG)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--max-patients", type=int, default=None)
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    log.info("%s\nMLH1 POSITIVE CONTROL -- %s, rung %s\n%s",
             "=" * 72, STUDY, args.rung, "=" * 72)
    log.info("pre-registered in docs/prereg_g2_mlh1_within_stratum.md")

    obs = load_obs(args.obs_cache, args.atlas)
    strata = strata_for(obs)

    deltas, report = study_deltas(
        args.atlas, obs, STUDY, seed=args.seed, rung=args.rung,
        reading="carcinoma", max_patients=args.max_patients,
        extra_genes=(TARGET,),
    )
    if deltas.empty:
        log.error("no patient scored; nothing to read")
        return 3

    deltas["short_id"] = deltas["patient_id"].astype(str).str.split(".").str[-1]
    deltas = deltas.merge(strata, on="short_id", how="left")
    unassigned = sorted(deltas.loc[deltas["mlh1_stratum"].isna(), "short_id"].unique())
    if unassigned:
        # Named, not dropped. A patient with no stratum is a patient this
        # reading cannot place, and how many there are bounds the reading.
        log.warning("  %d scored patient(s) carry no pre-registered stratum: %s",
                    len(unassigned), unassigned)
    deltas["arm"] = deltas["mlh1_stratum"].map(arm_of)
    deltas["granularity_rung"] = args.rung

    log.info("\nARMS ACTUALLY SCORED")
    counts = (deltas.drop_duplicates("short_id")
              .groupby(["arm", "mlh1_stratum"], dropna=False).size())
    log.info("%s", counts.to_string())

    primary_block = deltas[deltas["arm"] == PRIMARY_STRATUM]
    log.info("\n%s\nTHE PREMISE, IN THE METHYLATED ARM ONLY\n%s", "=" * 72, "=" * 72)
    log.info("  The controls have to hold in the arm the reading is about. A "
             "premise that\n  holds over all 29 patients does not license a "
             "claim about these 10.")
    premise = premise_holds(primary_block, seed=args.seed)
    log.info("  %s", premise[1])
    log.info(
        "  NOTE: premise_holds uses the percentile bootstrap, and at n=%d that "
        "interval\n  is %.0f%% narrower than the calibrated one. For an "
        "EQUIVALENCE test a narrower\n  interval is more likely to fit inside "
        "the tolerance, so this check is\n  anti-conservative here: it errs "
        "toward HOLDS. Read a HOLDS at this n as\n  weaker than the same word "
        "at n=44.", primary_block["patient_id"].nunique(),
        100 * (1 - 0.822))

    rows = []
    for arm in (PRIMARY_STRATUM, SECONDARY_STRATUM, UNDERPOWERED_STRATUM):
        block = deltas[deltas["arm"] == arm] if arm != UNDERPOWERED_STRATUM else \
            deltas[deltas["mlh1_stratum"] == UNDERPOWERED_STRATUM]
        for gene in sorted(block["gene"].unique()):
            rows.append(arm_reading(block, gene=gene, arm=arm, seed=args.seed)
                        | {"granularity_rung": args.rung,
                           "study_id": STUDY,
                           "standing": {
                               PRIMARY_STRATUM: "primary, powered",
                               SECONDARY_STRATUM: "secondary, CONFOUNDED with MMR status",
                               UNDERPOWERED_STRATUM: "pre-registered control, UNDERPOWERED",
                           }[arm]})
    readings = pd.DataFrame(rows)

    for arm in (PRIMARY_STRATUM, SECONDARY_STRATUM, UNDERPOWERED_STRATUM):
        block = readings[readings["arm"] == arm]
        if block.empty:
            continue
        log.info("\n%s\n%s -- %s\n%s", "=" * 72, arm,
                 block["standing"].iloc[0], "=" * 72)
        show = block[["gene", "n_patients", "detect_normal", "detect_tumour",
                      "mean_delta_cloglog", "ci_low", "ci_high",
                      "excludes_zero", "patients_with_signal"]]
        log.info("%s", show.to_string(index=False))

    primary = readings[(readings["arm"] == PRIMARY_STRATUM)
                       & (readings["gene"] == TARGET)]
    if primary.empty:
        log.error("%s produced no row in the primary arm", TARGET)
        return 3
    verdict = instrument_verdict(primary.iloc[0].to_dict(), premise)
    log.info("\n%s\nVERDICT ON THE INSTRUMENT\n%s", "=" * 72, "=" * 72)
    log.info("  %s", verdict["verdict"])
    log.info("  %s", verdict["detail"])
    log.info(
        "\n  WHAT IS NOT CLAIMED EITHER WAY. This is a positive control without "
        "a\n  mechanistic negative control -- the n=4 intact-MMRd arm cannot "
        "supply one.\n  So a fall here shows the statistic responds to a gene "
        "whose promoter is\n  methylated; it does not show the response is "
        "SPECIFIC to methylation.")

    meta = {
        "prereg": "docs/prereg_g2_mlh1_within_stratum.md",
        "supersedes_design": (
            "docs/prereg_g2_mlh1.md's difference-in-differences, which needed "
            "mlh1_intact_mmrd as a mechanistic negative control. That arm is "
            "n=4 here and n=4 on the original GSE178341 path after depth "
            "matching -- two independent pipelines, so it is a property of the "
            "cohort. The DiD is not available on this data."
        ),
        "target": TARGET,
        "study": STUDY,
        "granularity_rung": args.rung,
        "interval_method": CALIBRATED_METHOD,
        "why_not_the_usual_interval": (
            "the percentile bootstrap over patients, used everywhere else in "
            "this repository, has a false-positive rate of "
            f"{100 * expected_false_positive_rate(10):.1f}% at n=10 against a "
            "nominal 5% -- z*sqrt((n-1)/n)/t(n-1), a function of n alone. "
            "Measured and committed BEFORE this ran: "
            "results/*/interval_calibration.parquet."
        ),
        "arms": {
            PRIMARY_STRATUM: "primary, powered",
            SECONDARY_STRATUM: (
                "secondary, CONFOUNDED: mostly MMR-proficient patients, so it "
                "mixes methylation status with MSI status. Not the mechanistic "
                "control the original prereg wanted."
            ),
            UNDERPOWERED_STRATUM: (
                "pre-registered control, UNDERPOWERED. At n=4 the percentile "
                f"interval fires {100 * expected_false_positive_rate(4):.0f}% "
                "of the time under a true null. No verdict is taken from it."
            ),
        },
        "min_patients_with_signal": MIN_PATIENTS_WITH_SIGNAL,
        "premise_in_primary_arm": {"holds": bool(premise[0]), "detail": premise[1]},
        "premise_caveat": (
            "premise_holds is an EQUIVALENCE test on the percentile bootstrap. "
            "A narrower interval fits inside the tolerance more easily, so at "
            "n=10 that check errs toward HOLDS. Anti-conservative here."
        ),
        "instrument_verdict": verdict,
        "strata_sources_agree": not strata.attrs.get("crossings"),
        "strata_crossings": strata.attrs.get("crossings", []),
        "patients_without_stratum": unassigned,
        "study_deltas_report": report,
        "what_this_cannot_say": (
            "that the response is SPECIFIC to promoter methylation. That needs "
            "the mechanistic negative control, which this cohort does not have "
            "at usable size."
        ),
        "exploratory": False,
        "pre_registered": True,
    }
    written = write_versioned_table(
        readings, "mlh1_positive_control", seed=args.seed,
        results_dir=args.results_dir, allow_dirty=args.allow_dirty,
        extra_meta=meta,
    )
    log.info("\nwrote %s", written)

    per_patient = deltas.drop(columns=["patient_id_cohort"], errors="ignore")
    written = write_versioned_table(
        per_patient, "mlh1_positive_control_per_patient", seed=args.seed,
        results_dir=args.results_dir, allow_dirty=args.allow_dirty,
        extra_meta=meta,
    )
    log.info("wrote %s", written)
    return 0


if __name__ == "__main__":
    sys.exit(main())
