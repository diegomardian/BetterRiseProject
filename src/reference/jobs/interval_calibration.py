"""What the project's interval actually does at the n it is about to run at.

    python -m src.reference.jobs.interval_calibration

A THIN LOCAL RUN. No cluster and no atlas: the cohorts it calibrates against are
read off committed per-patient tables, and everything else is simulation. It
emits three tables:

``interval_heterogeneity``
    Between-patient SD of the per-patient delta, per gene, net of binomial
    sampling noise, on Pelka. This is the number a power calculation needs and
    the number the first version of the MLH1 power statement assumed was zero.

``interval_calibration``
    False-positive rate of the percentile bootstrap, BCa and the Student-t
    interval under a true null, at each n this project reports at.

``mlh1_power``
    Power of the MLH1 positive control at several silencing depths, on the
    calibrated interval, with that interval's own false-positive rate on every
    row.

WHY IT IS RUN BEFORE THE MLH1 READING RATHER THAN BESIDE IT. The point of a
pre-registration is that the analysis is fixed before the numbers exist. An
interval's calibration is part of the analysis: "the 95% interval excluded zero"
means something different if the interval is really an 89% one. Measuring that
after seeing the result would make the choice of interval a free parameter, and
this project has spent four routes' worth of effort on not having free
parameters at the end.

WHAT IT FOUND, so a reader knows what to look for. The percentile bootstrap --
used by ``premise_holds``, ``summarise``, ``specificity`` and
``control_log2_interval`` -- excludes zero about 9-11% of the time under a true
null at n=10, against a nominal 5%. BCa is no better. The Student-t interval is
calibrated at every n measured. Full argument in
``src/reference/interval_calibration.py``.
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
from src.common.paths import RESULTS_DIR
from src.common.provenance import DEFAULT_SEED
from src.reference.interval_calibration import (
    CALIBRATED_METHOD,
    INTERVAL_METHODS,
    MISCALIBRATION_TOLERANCE,
    NOMINAL_ALPHA,
    calibration_table,
    check_power_carries_its_own_calibration,
    expected_false_positive_rate,
    heterogeneity_tau,
    power_curve,
    width_ratio,
)

log = logging.getLogger(__name__)

#: MLH1's per-cell mean in the mature cells of the normal arm, counts per 10k.
#: From ``docs/prereg_g2_mlh1.md`` §"the expression scale explains it", measured
#: on the depth-matched decomposition. About one count per 250,000 UMIs, which
#: is the entire reason this reading needs a power calculation at all.
MLH1_CP10K_NORMAL = 0.039

#: A common gene's abundance on the same cohort, as the contrast case. GUCA2A's
#: baseline detection on Pelka is ~0.51 against MLH1's ~0.03: if the percentile
#: bootstrap only misbehaved on rare transcripts the fix would be a floor on
#: abundance rather than a change of interval, so both are measured.
COMMON_CP10K_NORMAL = 3.0

#: Fold changes to report power at. 0.5 and 0.25 bracket what "silencing" could
#: plausibly mean for a promoter-methylated gene; 1.0 is carried so the null
#: sits in the same table as the alternatives it is compared against.
FOLD_CHANGES = (1.0, 0.5, 0.25, 0.1)

#: Between-patient heterogeneity to report power across. 0.0 is the homogeneous
#: model -- kept only to show what it overstates. 0.2 is what Pelka's control
#: genes measure. 0.4 is the value the trend across baseline detection suggests
#: for a gene as rare as MLH1, and is the pessimistic end.
TAUS = (0.0, 0.2, 0.4)


def newest(name: str) -> Path | None:
    matches = sorted(RESULTS_DIR.glob(f"*/{name}.parquet"))
    return matches[-1] if matches else None


def pelka_with_strata() -> pd.DataFrame:
    """Pelka's committed per-patient rows, joined to the pre-registered strata.

    The join is on the short patient id: the atlas writes
    ``Pelka_2021_Cell.C110`` and the cohort table writes ``C110``. That is an
    identifier-space difference, and this repository has been bitten by four of
    them -- every time the symptom was an empty intersection reported as a
    finding. So the overlap is asserted rather than assumed.
    """
    deltas_path = newest("icbi_coexpression")
    if deltas_path is None:
        raise SystemExit("no results/*/icbi_coexpression.parquet")
    cohort_paths = sorted(glob.glob(str(RESULTS_DIR / "*" / "cohort_table.parquet")))
    if not cohort_paths:
        raise SystemExit("no results/*/cohort_table.parquet")

    deltas = pd.read_parquet(deltas_path)
    pelka = deltas[deltas["study_id"] == "Pelka_2021_Cell"].copy()
    if pelka.empty:
        raise SystemExit(f"{deltas_path} carries no Pelka_2021_Cell rows")
    pelka["short_id"] = pelka["patient_id"].astype(str).str.split(".").str[-1]

    cohort = pd.read_parquet(cohort_paths[-1])
    merged = pelka.merge(
        cohort[["patient_id", "mlh1_stratum", "matched"]],
        left_on="short_id", right_on="patient_id", suffixes=("", "_cohort"),
    )
    if merged.empty:
        raise SystemExit(
            f"joining {deltas_path.name} to {Path(cohort_paths[-1]).name} on the "
            f"short patient id produced no rows. The atlas writes "
            f"'Pelka_2021_Cell.C110' and the cohort table writes 'C110'; if "
            f"either has changed shape this is an identifier-space mismatch, "
            f"not a cohort with no methylated patients in it."
        )
    log.info("%s: %d Pelka patients scored, %d joined to a stratum",
             deltas_path.parent.name, pelka["short_id"].nunique(),
             merged["short_id"].nunique())
    return merged


def cohort_vectors(frame: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """One row per patient, as ``(mature cells per arm, median depth)``.

    The arms are depth-matched upstream, so ``n_normal == n_tumour`` and the two
    depths agree; the tumour count and the normal depth are taken and that
    equality is asserted rather than relied on silently.
    """
    one = frame.drop_duplicates("short_id")
    n_cells = one["n_tumour"].to_numpy().astype(int)
    if not np.array_equal(n_cells, one["n_normal"].to_numpy().astype(int)):
        raise SystemExit(
            "the two arms do not carry the same cell counts, so simulating them "
            "with one vector would model a design that was not run. Depth "
            "matching is supposed to equalise them."
        )
    return n_cells, one["depth_normal"].to_numpy().astype(float)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--trials", type=int, default=1500)
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    pelka = pelka_with_strata()

    # -- 1. heterogeneity, measured on the six committed genes ---------------
    log.info("\n%s\nBETWEEN-PATIENT HETEROGENEITY, Pelka, net of sampling noise\n%s",
             "=" * 72, "=" * 72)
    tau_table = heterogeneity_tau(pelka.drop_duplicates(["short_id", "gene"]),
                                  seed=args.seed)
    log.info("%s", tau_table.to_string(index=False))
    log.info(
        "\nEvery gene's spread exceeds what binomial sampling explains, and tau "
        "rises as\nbaseline detection falls. A power calculation with tau = 0 "
        "is a calculation about\na cohort more homogeneous than this one.")

    # -- 2. the cohorts this project reports at ------------------------------
    meth = pelka[pelka["mlh1_stratum"] == "mlh1_methylated"]
    unmeth = pelka[pelka["mlh1_stratum"] != "mlh1_methylated"]
    intact = pelka[pelka["mlh1_stratum"] == "mlh1_intact_mmrd"]

    cohorts = {
        "mlh1_methylated": cohort_vectors(meth),
        "mlh1_unmethylated": cohort_vectors(unmeth),
        "mlh1_intact_mmrd": cohort_vectors(intact),
    }
    adenoma_path = newest("icbi_adenoma")
    if adenoma_path is not None:
        adenoma = pd.read_parquet(adenoma_path)
        adenoma["short_id"] = adenoma["patient_id"].astype(str)
        for rung, block in adenoma.groupby("granularity_rung", observed=True):
            cohorts[f"adenoma_{rung}"] = cohort_vectors(block)

    log.info("\ncohorts calibrated against:")
    for name, (n_cells, depth) in cohorts.items():
        log.info("  %-22s n=%2d patients, cells/arm median %5.0f (min %d), "
                 "depth median %6.0f", name, len(n_cells), np.median(n_cells),
                 n_cells.min(), np.median(depth))

    # -- 3. calibration under a true null ------------------------------------
    log.info("\n%s\nFALSE-POSITIVE RATE UNDER A TRUE NULL. Nominal %.1f%%.\n%s",
             "=" * 72, 100 * NOMINAL_ALPHA, "=" * 72)
    calibration = calibration_table(
        cohorts=cohorts,
        abundances={"rare (MLH1)": MLH1_CP10K_NORMAL,
                    "common": COMMON_CP10K_NORMAL},
        taus=(0.0, 0.2), seed=args.seed, n_trials=args.trials,
    )
    wide = calibration.pivot_table(
        index=["cohort", "n_patients", "abundance", "tau"],
        columns="method", values="false_positive_rate",
    )
    log.info("%s", (100 * wide).round(1).to_string())

    worst = calibration.sort_values("false_positive_rate", ascending=False)
    bad = worst[worst["verdict"] == "MISCALIBRATED"]
    log.info("\n  %d of %d (cohort, abundance, tau, method) cells MISCALIBRATED",
             len(bad), len(calibration))
    for method in INTERVAL_METHODS:
        block = calibration[calibration["method"] == method]
        log.info("  %-11s worst %.1f%%, %d/%d miscalibrated", method,
                 100 * block["false_positive_rate"].max(),
                 int((block["verdict"] == "MISCALIBRATED").sum()), len(block))
    log.info(
        "\n  The bar is nominal + %.0fpp. It is small-n on BOTH statistics and "
        "BOTH\n  abundances, so it is not a property of MLH1 being rare, and "
        "BCa -- whose bias\n  and acceleration are estimated from the same ten "
        "numbers -- does not repair it.", 100 * MISCALIBRATION_TOLERANCE)

    log.info("\n%s\nTHE SIMULATION AGAINST THE CLOSED FORM\n%s", "=" * 72, "=" * 72)
    log.info("  percentile bootstrap of a mean ~ mean +/- z*s*sqrt((n-1)/n)/sqrt(n)")
    log.info("  the calibrated interval       =  mean +/- t(n-1)*s/sqrt(n)")
    log.info("  so the rate is a function of n ALONE -- no gene, no scale, no data\n")
    pct = calibration[calibration["method"] == "percentile"]
    check = (pct.groupby(["cohort", "n_patients"], as_index=False)
             .agg(simulated=("false_positive_rate", "mean"),
                  closed_form=("closed_form_rate", "first"),
                  width_ratio=("width_ratio_vs_t", "first"))
             .sort_values("n_patients"))
    check["skew_excess"] = check["simulated"] - check["closed_form"]
    log.info("%s", check.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    log.info(
        "\n  The closed form is a FLOOR: it assumes the per-patient values are "
        "normal\n  enough for the bootstrap mean to be. The excess above it is "
        "the skew of a\n  rare transcript's delta. Two independent routes to "
        "one number -- where they\n  agree the mechanism is understood, and "
        "where they part the gap is readable.")
    log.info(
        "\n  THIS GENERALISES BEYOND MLH1. Every percentile-bootstrap interval "
        "in this\n  repository -- premise_holds, summarise, specificity, "
        "control_log2_interval --\n  carries this, at whatever n its table ran "
        "at: %s.",
        ", ".join(f"n={int(r.n_patients)} -> {100 * r.closed_form:.1f}%"
                  for r in check.itertuples()))

    # -- 4. power for the MLH1 reading ---------------------------------------
    log.info("\n%s\nPOWER, MLH1 within the methylated stratum, %s interval\n%s",
             "=" * 72, CALIBRATED_METHOD, "=" * 72)
    n_cells, depth = cohorts["mlh1_methylated"]
    power = power_curve(
        n_cells=n_cells, depth=depth, cp10k=MLH1_CP10K_NORMAL,
        fold_changes=FOLD_CHANGES, taus=TAUS, method=CALIBRATED_METHOD,
        seed=args.seed, n_trials=args.trials, cohort="mlh1_methylated",
    )
    check_power_carries_its_own_calibration(power)
    log.info("%s", power.pivot_table(index="silencing_pct", columns="tau",
                                     values="power").mul(100).round(1).to_string())
    log.info("\n  false-positive rate of this interval on this cohort: %s",
             ", ".join(f"tau={t}: {100 * r:.1f}%" for t, r in
                       power.drop_duplicates("tau")[["tau", "false_positive_rate"]]
                       .to_numpy()))
    log.info(
        "\n  READ THE tau COLUMN BEFORE QUOTING A NUMBER. tau=0 is the "
        "homogeneous model\n  and it overstates power; Pelka's control genes "
        "measure tau ~ 0.2 and the trend\n  across baseline detection puts a "
        "gene as rare as MLH1 nearer 0.4.")

    meta = {
        "purpose": (
            "calibrate the interval the MLH1 positive control will report, "
            "BEFORE that reading is run. The interval is part of the analysis; "
            "choosing it after seeing the result would make it a free parameter."
        ),
        "source_deltas": str(newest("icbi_coexpression")),
        "mlh1_cp10k_normal": MLH1_CP10K_NORMAL,
        "mlh1_cp10k_source": "docs/prereg_g2_mlh1.md, depth-matched decomposition",
        "nominal_alpha": NOMINAL_ALPHA,
        "miscalibration_tolerance": MISCALIBRATION_TOLERANCE,
        "calibrated_method": CALIBRATED_METHOD,
        "methods_measured": sorted(INTERVAL_METHODS),
        "taus": list(TAUS),
        "fold_changes": list(FOLD_CHANGES),
        "n_trials": int(args.trials),
        "closed_form": (
            "the percentile bootstrap of a mean is approximately "
            "mean +/- z*s*sqrt((n-1)/n)/sqrt(n) against the correct "
            "mean +/- t(n-1)*s/sqrt(n), so its false-positive rate is "
            "P(|t(n-1)| > z*sqrt((n-1)/n)) -- a function of n alone. It is a "
            "floor: skew in the per-patient values adds to it."
        ),
        "closed_form_rates": {
            str(n): expected_false_positive_rate(n)
            for n in sorted({int(len(c[0])) for c in cohorts.values()})
        },
        "width_ratios_vs_t": {
            str(n): width_ratio(n)
            for n in sorted({int(len(c[0])) for c in cohorts.values()})
        },
        "generative_model": (
            "Poisson thinning p = 1 - exp(-mu) with mu = cp10k/1e4 * depth, "
            "per-patient cell counts and depths taken from the real cohort, "
            "per-patient true log fold change ~ Normal(log fc, tau^2), scored "
            "through cloglog_rate including its boundary rule"
        ),
        "what_this_is_not": (
            "a measurement of MLH1. No target-gene expression is read anywhere "
            "in this job; the only MLH1 quantity used is its already-published "
            "normal-arm abundance, which sets the detection rate."
        ),
        "power_rule": (
            "every row of mlh1_power carries the false-positive rate of its OWN "
            "method, measured on the same generator. Power quoted from one "
            "interval beside another interval's coverage overstates the design."
        ),
        "exploratory": False,
        "pre_registered": True,
        "prereg": "docs/prereg_g2_mlh1_within_stratum.md",
    }
    for frame, name in ((tau_table, "interval_heterogeneity"),
                        (calibration, "interval_calibration"),
                        (power, "mlh1_power")):
        written = write_versioned_table(
            frame, name, seed=args.seed, results_dir=args.results_dir,
            allow_dirty=args.allow_dirty, extra_meta=meta,
        )
        log.info("wrote %s", written)
    return 0


if __name__ == "__main__":
    sys.exit(main())
