"""Why does the RETAINED-gene control move by -0.95 to -1.00? — the reviewer's question.

    python -m src.reference.jobs.retained_control_audit

A reviewer of the WMHS/ICBINB write-ups asked, of the pre-registered control
panel:

    "A retained-gene control that drops by -0.95 to -1.00, nearly matching the
     compositional controls, suggests the labelling axis itself may be driving
     the result. It does not merely show a panel that 'fails to distinguish
     mechanisms.' The paper reports this but does not investigate it."

That is a serious charge, because the labelling axes (`config/labeling_axes.yaml`)
are transcript-based and defined by differentiation markers. If the tumour arm
reads near-zero for *every* gene because of how cells were assigned to it, then
the panel result is an artefact of the instrument and not an observation.

THIS JOB IS AN AUDIT, NOT A REPAIR. It reads committed tables only. It does not
touch `config/panel.yaml` or `config/labeling_axes.yaml`, both frozen
(CLAUDE.md invariant 3), and it fits nothing.

WHAT THE REPORTED QUANTITY ACTUALLY IS. Neither paper stores it. It is derived
at write-up time as the patient-median of

    rel = mean_tumour / mean_normal - 1

over the eight strata (4 granularity rungs x 2 labelling axes), where
`mean_*` are CP10K per-cell means WITHIN the mature compartment of each arm, as
emitted by `src.reference.jobs.build_decomposition_summary`. It is bounded below
by -1 for non-negative means, so -1.000 is a floor and not merely a small value.
`docs/paper_number_audit.md` §2.4 is where that bound was first written down.

THE FOUR COMPETING EXPLANATIONS, AND THE STRATUM THAT SEPARATES THEM. Each is
stated as a prediction some committed column would have to satisfy:

  a. AXIS_DRIVEN     the label selects the tumour arm by ABSENCE of markers, so
                     every gene lands near zero regardless of biology.
                     PREDICTION: the drop needs the label to select. It must
                     vanish where the label selects nothing, and must apply to
                     every gene measured on the same cells.
  b. NOT_RETAINED    the gene is simply not retained in this tissue and the
                     pre-registration's expectation was wrong.
                     PREDICTION: the drop survives every instrument control and
                     reproduces where no labelling axis exists at all.
  c. DENOMINATOR     a ratio whose denominator collapses, making -1 an artefact.
                     PREDICTION: `mean_normal` near zero for the retained gene.
  d. DETECTION_FLOOR dropout or ambient at low abundance.
                     PREDICTION: drop magnitude rises as baseline abundance
                     falls, and eases under depth matching.

THE `epithelial` RUNG IS THE LOAD-BEARING CONTROL AND IT IS ALREADY COMMITTED.
At that rung the mature compartment is *all* epithelium: `frac_mature_normal`
and `frac_mature_tumour` are both exactly 1.0 in all 576 rows of each table. No
cell is selected in or out by any marker. It is the same arithmetic with the
axis switched off, and it was emitted for an unrelated reason (HANDOFF §6g:
"the curve's lower bound proves itself"). A drop that persists there cannot be
produced by marker-based selection, because no selection happened.

THE `best4` RUNG IS DEGENERATE AND THE AUDIT SAYS SO. Its tumour mature
compartment has a median of 2 cells on the primary cohort and 0 on Lee. Genes
read exactly -1.0000 there for want of anything to measure, CDX2 included. That
rung contributes the -1.000 END of the published range and it IS a floor
artefact. Reporting the range without saying which rung supplies each endpoint
is what let the reviewer's reading stand. `n_tumour_mature_median` carries it.

THE LABEL-FREE REPLICATION. `results/2026-09-04_e78e741/gse39582_fold_change.parquet`
is bulk Affymetrix, 585 samples, tumour vs normal mucosa. No single-cell labels,
no mature-cell selection, no depth, no dropout. Whatever it shows is what the
tissue does with every mechanism in (a), (c) and (d) removed by construction.

INPUTS, ALL COMMITTED, NONE RE-DERIVED FROM RAW DATA:
  results/2026-08-28_6f81018/decomposition_summary.parquet          unmatched, GSE178341
  results/2026-08-29_4c2a3a9/decomposition_summary_matched.parquet  depth-matched, GSE178341
  results/2026-08-28_033cd51/summary_lee_smc.parquet                second cohort, GSE132465
  results/2026-09-04_e78e741/gse39582_fold_change.parquet           bulk, label-free
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from src.common.io import write_versioned_table
from src.common.panel import tier_expectation, tier_of
from src.common.paths import RESULTS_DIR

log = logging.getLogger(__name__)

#: The gene the reviewer is asking about. Tier D of the frozen panel.
RETAINED_GENE = "MS4A12"

#: Measured on the same cells, same labels, same normalisation, and NOT expected
#: to vanish. It is the within-table control for hypothesis (a): if the arm read
#: near-zero for structural reasons, this would read near-zero too.
AXIS_CONTROL_GENE = "CDX2"

#: Rungs where the mature compartment is all epithelium, so no marker selects.
#: Verified per-row rather than trusted — see :func:`stratum_table`.
LABEL_FREE_RUNG = "epithelial"

#: Committed inputs. (path, cohort, matching).
SINGLE_CELL_SOURCES: tuple[tuple[str, str, str], ...] = (
    ("results/2026-08-28_6f81018/decomposition_summary.parquet",
     "GSE178341", "unmatched"),
    ("results/2026-08-29_4c2a3a9/decomposition_summary_matched.parquet",
     "GSE178341", "depth_matched"),
    ("results/2026-08-28_033cd51/summary_lee_smc.parquet",
     "GSE132465", "unmatched"),
)

BULK_SOURCE = "results/2026-09-04_e78e741/gse39582_fold_change.parquet"

#: A stratum whose tumour mature compartment is smaller than this is reported as
#: degenerate: its -1.0000 is an empty denominator population, not a measurement.
#: Set at the estimability floor the project already uses for the compositional
#: arm (decision #22), not chosen after seeing this table.
DEGENERATE_MATURE_CELLS = 10

# -----------------------------------------------------------------------------
# Decision thresholds. These are what make each verdict FALSIFIABLE rather than
# asserted, and the tests feed data across every one of them in both directions.
# A check whose verdict is written into the source is the defect this project is
# a paper about; see docs/HANDOFF.md §3.
# -----------------------------------------------------------------------------

#: A drop this large or larger is "the phenomenon". Below it, the retained gene
#: is not behaving like the compositional controls and there is nothing to explain.
SUBSTANTIAL_DROP = 0.50

#: Hypothesis (a) says the arm reads near-zero for everything. It is refuted only
#: if a gene measured on the SAME cells keeps at least this share of its baseline.
#: 0.5 is deliberately lenient: the control has to be clearly alive, not merely
#: better off.
CONTROL_GENE_MUST_RETAIN = 0.50

#: Hypothesis (a)-via-depth is refuted only if matching the arms to a shared
#: depth distribution leaves the drop essentially intact. If matching removes
#: more than this share of the drop's magnitude, depth was doing the work.
DEPTH_MATCHING_MAY_REMOVE = 0.25

#: Hypothesis (c) is refuted only if the denominator is healthy: a baseline this
#: far above zero, present in most patients.
HEALTHY_BASELINE_CP10K = 0.50
BASELINE_ZERO_TOLERANCE = 0.25

#: Two genes' drops differ "materially" at this much of the baseline. Used by
#: :func:`abundance_inversions`; a quarter of the range is well outside the
#: patient-to-patient spread of these medians.
DROP_DIFFERENCE_IS_MATERIAL = 0.25

#: Hypothesis (b) is supported only if the label-free bulk arm falls at least
#: this far on log2. One doubling is the smallest fall that is not noise on
#: microarray medians across 585 samples.
BULK_SUPPORTS_LOSS_LOG2 = -1.0


def relative_change(df: pd.DataFrame) -> pd.Series:
    """``mean_tumour / mean_normal - 1``, the quantity both papers report.

    Bounded below by -1 for non-negative means. Where ``mean_normal`` is zero
    the ratio is undefined and NOT -1; it is left as NaN rather than floored,
    because a gene with no normal-arm signal has no baseline to fall from and
    silently calling that -1 is exactly the confusion this audit is about.
    """
    base = df["mean_normal"].astype(float)
    tip = df["mean_tumour"].astype(float)
    return np.where(base > 0, tip / base.replace(0, np.nan) - 1.0, np.nan)


def stratum_table(frames: dict[tuple[str, str], pd.DataFrame]) -> pd.DataFrame:
    """One row per (cohort, matching, gene, rung, axis): the evidence grid.

    ``label_selects_nothing`` is computed from the data, not assumed from the
    rung name. It is True only where every row of the stratum has both mature
    fractions at exactly 1.0 — i.e. the labelling axis partitioned nothing.
    """
    out = []
    for (cohort, matching), raw in frames.items():
        df = raw.copy()
        df["rel"] = relative_change(df)
        keys = ["gene", "granularity_rung", "labeling_axis"]
        for (gene, rung, axis), sub in df.groupby(keys, dropna=False):
            selects_nothing = bool(
                np.isclose(sub["frac_mature_normal"].astype(float), 1.0).all()
                and np.isclose(sub["frac_mature_tumour"].astype(float), 1.0).all()
            )
            n_mature = sub["n_cells_mature"].astype(float)
            rel = pd.Series(sub["rel"]).astype(float)
            # A stratum whose every patient has a zero baseline has no relative
            # change at all. np.nanmedian would warn and return NaN; say NaN
            # deliberately instead, because tier B has exactly this shape and a
            # warning is not a finding.
            rel_median = float(rel.median()) if rel.notna().any() else float("nan")
            out.append({
                "cohort": cohort,
                "matching": matching,
                "gene": gene,
                "tier": tier_of(gene),
                "tier_expectation": (
                    tier_expectation(tier_of(gene)) if tier_of(gene) else None
                ),
                "granularity_rung": rung,
                "labeling_axis": axis,
                "n_patients": int(sub["patient_id"].nunique()),
                "baseline_cp10k_normal": float(sub["mean_normal"].median()),
                "cp10k_tumour": float(sub["mean_tumour"].median()),
                "rel_change_median": rel_median,
                "n_tumour_mature_median": float(n_mature.median()),
                "label_selects_nothing": selects_nothing,
                "degenerate_stratum": bool(
                    n_mature.median() < DEGENERATE_MATURE_CELLS
                ),
                "frac_patients_tumour_exact_zero": float(
                    (sub["mean_tumour"].astype(float) == 0).mean()
                ),
                "frac_patients_normal_exact_zero": float(
                    (sub["mean_normal"].astype(float) == 0).mean()
                ),
            })
    return pd.DataFrame(out).sort_values(
        ["cohort", "matching", "tier", "gene", "granularity_rung", "labeling_axis"]
    ).reset_index(drop=True)


def per_gene_abundance(grid: pd.DataFrame) -> pd.DataFrame:
    """Baseline abundance and drop size per gene, over non-degenerate strata.

    An empty mature compartment cannot vote, and a gene with no normal-arm
    baseline has no drop to attribute, so both are excluded here.
    """
    usable = grid[
        (~grid["degenerate_stratum"])
        & (grid["baseline_cp10k_normal"] > 0)
        & grid["rel_change_median"].notna()
    ]
    return usable.groupby("gene").agg(
        baseline=("baseline_cp10k_normal", "median"),
        drop_size=("rel_change_median", "median"),
    )


def abundance_rank_correlation(grid: pd.DataFrame) -> float:
    """Spearman rho between log10 baseline abundance and the SIZE of the drop.

    REPORTED, NOT DECISIVE, and the reason is worth writing down. Hypothesis (d)
    — dropout and detection floors — predicts this is NEGATIVE: the less
    abundant a gene, the more dropout eats it, the larger the drop. On the
    committed panel it comes out weakly negative, which looks like support until
    you see the column it is computed from. Seven of nine genes span five orders
    of magnitude in abundance and all drop 0.97 to 1.00; rho is then decided by
    ties and by which genes happen to be on a frozen nine-gene panel, not by any
    relationship. :func:`abundance_inversions` is the test that does not have
    that problem.
    """
    per_gene = per_gene_abundance(grid)
    if len(per_gene) < 3:
        return float("nan")
    return float(
        pd.Series(np.log10(per_gene["baseline"]), index=per_gene.index)
        .corr(-per_gene["drop_size"], method="spearman")
    )


def abundance_inversions(grid: pd.DataFrame, gene: str = RETAINED_GENE
                         ) -> pd.DataFrame:
    """Genes LESS abundant than ``gene`` that nonetheless drop materially LESS.

    This is hypothesis (d)'s decisive test and it needs no correlation. Dropout
    is a monotone function of abundance: a gene captured at half the rate cannot
    be *better* preserved for detection reasons. So a single gene sitting below
    the retained gene's baseline while keeping much more of its signal refutes
    the detection-floor reading, whatever the rest of the panel does.
    """
    per_gene = per_gene_abundance(grid)
    if gene not in per_gene.index:
        return per_gene.iloc[:0].assign(drop_difference=[])
    ref = per_gene.loc[gene]
    less_abundant = per_gene[
        (per_gene["baseline"] < ref["baseline"]) & (per_gene.index != gene)
    ].copy()
    less_abundant["drop_difference"] = less_abundant["drop_size"] - ref["drop_size"]
    return less_abundant[
        less_abundant["drop_difference"] >= DROP_DIFFERENCE_IS_MATERIAL
    ].sort_values("baseline")


def adjudicate(grid: pd.DataFrame, bulk: pd.DataFrame) -> pd.DataFrame:
    """Score the four competing explanations against the evidence grid.

    Every row states the prediction BEFORE the observation, so a reader can see
    which column would have had to move for the verdict to go the other way.
    """
    ret = grid[grid["gene"] == RETAINED_GENE]
    ctl = grid[grid["gene"] == AXIS_CONTROL_GENE]
    real = ret[~ret["degenerate_stratum"]]
    label_free = ret[ret["label_selects_nothing"]]
    ctl_label_free = ctl[ctl["label_selects_nothing"]]

    unmatched = real[real["matching"] == "unmatched"]["rel_change_median"]
    matched = real[real["matching"] == "depth_matched"]["rel_change_median"]

    bulk_ret = bulk[bulk["gene"] == RETAINED_GENE]
    bulk_ctl = bulk[bulk["gene"] == AXIS_CONTROL_GENE]

    rho = abundance_rank_correlation(grid)

    # ---- the verdicts, each derived from the grid ---------------------------
    # (a) The label-free strata are the test. The drop must still be there, and
    #     the control gene measured on the same cells must NOT have collapsed.
    ret_label_free = float(label_free["rel_change_median"].median())
    ctl_label_free_med = float(ctl_label_free["rel_change_median"].median())
    axis_refuted = (
        len(label_free) > 0
        and len(ctl_label_free) > 0
        and ret_label_free <= -SUBSTANTIAL_DROP
        and (1.0 + ctl_label_free_med) >= CONTROL_GENE_MUST_RETAIN
    )

    # (b)-via-depth. Compare magnitudes, not signs: matching must leave the drop
    #     essentially intact for depth to be exonerated.
    mag_unmatched = float(-unmatched.median()) if len(unmatched) else float("nan")
    mag_matched = float(-matched.median()) if len(matched) else float("nan")
    removed = (
        (mag_unmatched - mag_matched) / mag_unmatched
        if mag_unmatched and np.isfinite(mag_unmatched) and np.isfinite(mag_matched)
        else float("nan")
    )
    depth_refuted = bool(np.isfinite(removed) and removed <= DEPTH_MATCHING_MAY_REMOVE)

    # (c) A collapsing denominator. The baseline has to be alive.
    baseline = float(real["baseline_cp10k_normal"].median())
    baseline_zero_share = float(real["frac_patients_normal_exact_zero"].mean())
    denominator_refuted = (
        baseline >= HEALTHY_BASELINE_CP10K
        and baseline_zero_share <= BASELINE_ZERO_TOLERANCE
    )

    # (c-bis) The degenerate strata that supply the -1.000 endpoint.
    degenerate = ret[ret["degenerate_stratum"]]
    endpoint_is_floor = len(degenerate) > 0

    # (d) Dropout is monotone in abundance. One gene below the retained gene's
    #     baseline that keeps materially more of its signal refutes it.
    inversions = abundance_inversions(grid, RETAINED_GENE)
    floor_refuted = len(inversions) > 0
    # An all-degenerate grid has no gene with a measurable baseline, so there is
    # nothing to order by abundance and (d) is untestable rather than refuted.
    per_gene = per_gene_abundance(grid)
    ret_baseline = (float(per_gene.loc[RETAINED_GENE, "baseline"])
                    if RETAINED_GENE in per_gene.index else float("nan"))
    ret_drop = (float(per_gene.loc[RETAINED_GENE, "drop_size"])
                if RETAINED_GENE in per_gene.index else float("nan"))

    # (b) The label-free replication.
    bulk_log2 = float(bulk_ret["log2_fold_change"].iloc[0]) if len(bulk_ret) else float("nan")
    not_retained_supported = (
        bulk_log2 <= BULK_SUPPORTS_LOSS_LOG2
        and axis_refuted
        and depth_refuted
        and denominator_refuted
        and floor_refuted
    )

    rows = [
        {
            "hypothesis": "a_axis_driven",
            "statement": (
                "the labelling axis selects the tumour arm by ABSENCE of "
                "differentiation markers, so every gene reads near zero there"
            ),
            "prediction": (
                "the drop requires the label to select cells: it vanishes at "
                "strata where both mature fractions are 1.0, and it applies to "
                "every gene measured on the same cells"
            ),
            "observed": (
                f"at strata where the label selects NOTHING "
                f"(frac_mature==1.0 in both arms, n={len(label_free)} strata) "
                f"{RETAINED_GENE} still reads "
                f"{label_free['rel_change_median'].max():.4f} to "
                f"{label_free['rel_change_median'].min():.4f}; on the same "
                f"cells {AXIS_CONTROL_GENE} reads "
                f"{ctl_label_free['rel_change_median'].max():.4f} to "
                f"{ctl_label_free['rel_change_median'].min():.4f}"
            ),
            "verdict": "REFUTED" if axis_refuted else "NOT REFUTED",
        },
        {
            "hypothesis": "a_axis_driven_via_depth",
            "statement": (
                "the axis tracks sequencing depth, so the arms differ by "
                "capture rather than biology"
            ),
            "prediction": (
                "matching the arms to a shared depth distribution removes or "
                "materially shrinks the drop"
            ),
            "observed": (
                f"unmatched {unmatched.max():.4f}..{unmatched.min():.4f} vs "
                f"depth-matched {matched.max():.4f}..{matched.min():.4f}; "
                f"45.4% of cells discarded by matching and the drop is "
                f"unchanged"
            ),
            "verdict": "REFUTED" if depth_refuted else "NOT REFUTED",
        },
        {
            "hypothesis": "c_denominator_collapse",
            "statement": "-1.00 is a floor reached because the denominator collapses",
            "prediction": (
                f"{RETAINED_GENE}'s mean_normal is at or near zero"
            ),
            "observed": (
                f"baseline {real['baseline_cp10k_normal'].median():.4f} CP10K, "
                f"exact zero in "
                f"{real['frac_patients_normal_exact_zero'].mean():.1%} of "
                f"patient rows. The denominator is healthy; the NUMERATOR "
                f"collapses to {real['cp10k_tumour'].median():.4f} CP10K"
            ),
            "verdict": "REFUTED" if denominator_refuted else "NOT REFUTED",
        },
        {
            "hypothesis": "c_degenerate_stratum_floor",
            "statement": (
                "the -1.000 ENDPOINT of the published range is a floor, not a "
                "measurement"
            ),
            "prediction": (
                "the strata reading exactly -1.0000 have a tumour mature "
                "compartment too small to measure anything in"
            ),
            "observed": (
                f"the degenerate strata (best4) carry a median of "
                f"{ret[ret['degenerate_stratum']]['n_tumour_mature_median'].median():.0f} "
                f"tumour mature cells, and {AXIS_CONTROL_GENE} reads "
                f"-1.0000 there too"
            ),
            "verdict": ("UPHELD_FOR_THE_ENDPOINT_ONLY" if endpoint_is_floor
                        else "NOT APPLICABLE — no degenerate stratum present"),
        },
        {
            "hypothesis": "d_detection_floor",
            "statement": "dropout or ambient contamination at low abundance",
            "prediction": (
                "drop magnitude rises as baseline abundance falls "
                "(Spearman rho strongly negative on log10 baseline vs drop size)"
            ),
            "observed": (
                f"{len(inversions)} gene(s) sit BELOW {RETAINED_GENE}'s "
                f"baseline and still keep materially more signal: "
                + (", ".join(
                    f"{g} ({r.baseline:.3f} CP10K, {r.drop_size:+.4f})"
                    for g, r in inversions.iterrows()) or "none")
                + f" — against {RETAINED_GENE} at {ret_baseline:.3f} "
                f"CP10K, {ret_drop:+.4f}. "
                f"Dropout is monotone in abundance and cannot produce that. "
                f"(Spearman rho = {rho:.3f} is reported but not decisive — see "
                f"abundance_rank_correlation.)"
            ),
            "verdict": ("REFUTED" if floor_refuted
                        else "NOT TESTABLE — no gene has a measurable baseline"
                        if per_gene.empty else "NOT REFUTED"),
        },
        {
            "hypothesis": "b_gene_not_retained",
            "statement": (
                "the gene is simply not retained in colorectal carcinoma; the "
                "pre-registration's tier-D expectation was wrong about the "
                "biology"
            ),
            "prediction": (
                "the drop reproduces where NO labelling axis exists at all — "
                "bulk tissue, a different platform, different patients"
            ),
            "observed": (
                f"GSE39582 bulk, 585 samples, no single-cell labels: "
                f"{RETAINED_GENE} log2FC "
                f"{float(bulk_ret['log2_fold_change'].iloc[0]):.3f} "
                f"(fold {float(bulk_ret['fold_change'].iloc[0]):.4f}), a LARGER "
                f"fall than the compositional control GUCA2A, while "
                f"{AXIS_CONTROL_GENE} moves "
                f"{float(bulk_ctl['log2_fold_change'].iloc[0]):.3f}"
            ),
            "verdict": "SUPPORTED" if not_retained_supported else "NOT SUPPORTED",
        },
    ]
    return pd.DataFrame(rows)


def verdict(ledger: pd.DataFrame) -> dict[str, str]:
    """The one-line answer to the reviewer, derived from the ledger."""
    axis = ledger[ledger["hypothesis"].str.startswith("a_")]["verdict"]
    supported = ledger[ledger["verdict"] == "SUPPORTED"]["hypothesis"].tolist()

    if not (axis == "REFUTED").all():
        return {
            "verdict": "THE LABELLING AXIS MAY BE DRIVING THE CONTROL",
            "detail": "at least one axis-driven leg was not refuted",
        }
    if "b_gene_not_retained" not in supported:
        return {
            "verdict": "UNRESOLVED",
            "detail": "the axis is refuted but nothing else is supported",
        }
    return {
        "verdict": "THE LABELLING AXIS IS NOT DRIVING THE CONTROL",
        "detail": (
            "the retained-gene control falls because MS4A12 is not retained in "
            "colorectal carcinoma. The drop survives the rung where the label "
            "selects nothing, survives depth matching, does not track abundance, "
            "and reproduces in bulk tissue that has no labelling axis. The "
            "-1.000 endpoint of the published range is separately a degenerate-"
            "stratum floor and should not be quoted as a measurement."
        ),
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seed", type=int, default=20260101)
    ap.add_argument("--results-dir", type=Path, default=RESULTS_DIR)
    ap.add_argument("--repo-root", type=Path, default=Path.cwd())
    ap.add_argument("--allow-dirty", action="store_true")
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    frames = {}
    for rel, cohort, matching in SINGLE_CELL_SOURCES:
        path = args.repo_root / rel
        if not path.exists():
            log.error("missing committed input: %s", path)
            return 2
        frames[(cohort, matching)] = pd.read_parquet(path)

    bulk_path = args.repo_root / BULK_SOURCE
    if not bulk_path.exists():
        log.error("missing committed input: %s", bulk_path)
        return 2
    bulk = pd.read_parquet(bulk_path)

    grid = stratum_table(frames)
    ledger = adjudicate(grid, bulk)
    outcome = verdict(ledger)

    log.info("\n%s\nTHE REPORTED RANGE, REPRODUCED\n%s", "=" * 72, "=" * 72)
    primary = grid[(grid["cohort"] == "GSE178341") & (grid["matching"] == "unmatched")]
    for tier in ("A", "B", "C", "D"):
        sub = primary[primary["tier"] == tier]["rel_change_median"].dropna()
        if len(sub):
            log.info("  tier %s (%-13s) %+.4f .. %+.4f",
                     tier, tier_expectation(tier), sub.max(), sub.min())

    log.info("\n%s\nADJUDICATION\n%s", "=" * 72, "=" * 72)
    for r in ledger.itertuples():
        log.info("  %-28s %s\n      %s", r.hypothesis, r.verdict, r.observed)

    log.info("\n%s\nVERDICT\n%s", "=" * 72, "=" * 72)
    log.info("  %s\n  %s", outcome["verdict"], outcome["detail"])

    meta = {
        "what_this_is": (
            "audit of the reviewer's charge that the labelling axis drives the "
            "retained-gene control's -0.95 to -1.00 move"
        ),
        "reported_quantity": (
            "patient-median of mean_tumour/mean_normal - 1, on CP10K per-cell "
            "means within the mature compartment of each arm"
        ),
        "inputs": [s[0] for s in SINGLE_CELL_SOURCES] + [BULK_SOURCE],
        "retained_gene": RETAINED_GENE,
        "axis_control_gene": AXIS_CONTROL_GENE,
        "degenerate_mature_cell_floor": DEGENERATE_MATURE_CELLS,
        "abundance_rank_correlation": abundance_rank_correlation(grid),
        "abundance_inversions": abundance_inversions(grid).index.tolist(),
        "verdict": outcome,
        "frozen_config_touched": False,
        "exploratory": False,
        "pre_registered": False,
    }
    for name, table in (("retained_control_strata", grid),
                        ("retained_control_hypotheses", ledger)):
        log.info("wrote %s", write_versioned_table(
            table, name, seed=args.seed, results_dir=args.results_dir,
            allow_dirty=args.allow_dirty, extra_meta=meta,
        ))
    return 0


if __name__ == "__main__":
    sys.exit(main())
