"""Does the differentiation fall track Wnt, inside surviving mature cells?

    python -m src.reference.jobs.wnt_mechanism --atlas <path>

Pre-registered in ``docs/prereg_wnt_mechanism.md``. Read §5 before reading any
number this produces: the branches and their consequences are fixed there.

WHAT THIS FOLLOWS. Avenue A established *that* GUCA2A's loss in adenoma is
cell-intrinsic — the cells are still there and have turned their output down
(``docs/prereg_adenoma_decomposition.md`` RESULT). It says nothing about why.
This is Stage 3's per-cell test, never run on human tissue, and it is the
cheapest mechanism question available: the cells are already scored, the atlas
is already on disk, and invariant 8 fixed the signature in week 0.

THE ESTIMAND, and both conditioners are load-bearing. Per patient, within the
mature cells of one arm: the **partial Spearman** correlation between the
per-cell Wnt-target score and each panel gene's per-cell CP10K, conditioning on
**the maturity score** and **log library depth**. Then aggregate across patients
with a Student-t interval — the patient is the unit of inference (invariant 5),
and the interval is not the percentile bootstrap (``docs/HANDOFF.md`` §3a).

Without the maturity conditioner the correlation recovers residual
differentiation state inside the mature bin and reports it as Wnt. Without the
depth conditioner it recovers library size. The housekeeping genes are scored
through the identical path and are the empirical floor: whatever ACTB and KRT8
return is how much correlation survives conditioning for technical reasons, and
a target gene's value means only the amount by which it exceeds that.

WHAT IT NEEDS. A cluster pass over the ICBI atlas — **no new data and no
download.** No committed table carries per-cell values; they are all per
(patient, gene) aggregates, which is why this cannot run on a laptop.
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
from src.reference.interval_calibration import student_t_interval
from src.reference.jobs.coexpression_silencing import (
    AXIS,
    DEPTH_QUANTILE,
    GENE_ROLES,
    RUNG,
)
from src.reference.jobs.icbi_coexpression import (
    ARM_MAPS,
    BATCH_KEY,
    DEFAULT_ATLAS,
    eligible_patients,
    load_obs,
    mature_bin,
)
from src.reference.wnt_score import (
    MIN_CELLS_FOR_CORRELATION,
    SIGNATURE,
    assert_no_signature_leakage,
    partial_spearman,
    wnt_score,
    wnt_stem_verdict,
)

log = logging.getLogger(__name__)

STUDY = "Chen_2021_Cell"

#: `lineage` only. Not `best4`: its intrinsic arm carries no claim
#: (prereg_adenoma_decomposition.md Amendment 3 and RESULT), and running a
#: mechanism test against a rung whose result is retracted would be asking why
#: something is true that has not been established.
RUNGS: tuple[str, ...] = (RUNG,)


def per_patient(
    atlas: Path, obs: pd.DataFrame, *, seed: int, rung: str,
    max_patients: int | None = None,
) -> tuple[pd.DataFrame, dict]:
    """One row per (patient, arm, gene). Returns (rows, report)."""
    from src.common.panel import load_axes, panel_genes, tier_genes
    from src.reference.icbi_slice import assert_raw_counts, read_cells, read_var
    from src.reference.labels import (
        assign_labels,
        label_column,
        maturity_score,
    )
    from src.reference.qc import apply_qc, cell_qc_metrics, qc_thresholds

    axis_markers = [g for a in load_axes()["axes"].values()
                    for g in (a.get("genes") or [])]
    # THE CIRCULARITY GUARD, before any cell is read. A signature containing a
    # label marker makes the correlation with maturity true by construction.
    assert_no_signature_leakage(panel_genes(), axis_markers)

    var = read_var(atlas)
    symbols = var["gene_symbol"]
    rows_obs, patients = eligible_patients(obs, STUDY, reading="adenoma")
    if max_patients:
        patients = patients[:max_patients]

    report = {
        "study_id": STUDY, "granularity_rung": rung, "labeling_axis": AXIS,
        "signature": list(SIGNATURE),
        "n_patients_eligible": len(patients),
        "conditioners": ["maturity_score", "log_depth"],
        "leakage_checked_against": {"panel": sorted(panel_genes()),
                                    "axis_markers": sorted(set(axis_markers))},
    }
    out: list[dict] = []
    skipped: list[dict] = []

    for i, patient in enumerate(patients, 1):
        block_obs = rows_obs[rows_obs["patient_id"].astype(str) == patient]
        block = read_cells(atlas, block_obs.index.to_numpy(), n_genes=len(symbols))
        assert_raw_counts(block, context=f"{STUDY}/{patient}")

        metrics = cell_qc_metrics(block, symbols, batch=block_obs[BATCH_KEY])
        keep = apply_qc(metrics, qc_thresholds(metrics)).to_numpy()
        if keep.sum() < 50:
            skipped.append({"patient_id": patient, "reason": "QC"})
            continue

        kept = block[keep]
        tissue = block_obs["tissue"].to_numpy()[keep]
        comp = block_obs["compartment"].to_numpy()[keep]
        try:
            labels = assign_labels(
                kept, symbols, compartment=comp,
                sample_id=block_obs[BATCH_KEY].to_numpy()[keep],
                target_genes=sorted(tier_genes("A")), tissue=tissue,
                patient_id=block_obs["patient_id"].to_numpy()[keep],
                depth_quantile=DEPTH_QUANTILE, seed=seed,
                index=pd.Index(block_obs.index.to_numpy()[keep]),
            )
        except Exception as exc:  # noqa: BLE001 — same skip-and-name as study_deltas
            skipped.append({"patient_id": patient, "reason": f"label: {exc}"})
            continue

        # The maturity SCORE, not the label. The label is a bin; the score is
        # the continuous quantity the conditioning needs, and it comes from the
        # same axis the labeller used so the two cannot drift apart.
        maturity = maturity_score(kept, symbols, AXIS,
                                  target_genes=sorted(tier_genes("A")), seed=seed)
        call = labels[label_column(AXIS, rung)].astype(str).to_numpy()
        depth = np.asarray(kept.sum(axis=1), dtype=float).ravel()
        score, score_report = wnt_score(kept, symbols, depth=depth)

        for arm in ("normal", "tumour"):
            mask = (call == mature_bin(rung)) & (tissue == arm)
            n = int(mask.sum())
            if n < MIN_CELLS_FOR_CORRELATION:
                skipped.append({"patient_id": patient, "arm": arm,
                                "reason": f"{n} mature cells"})
                continue
            conditioners = np.column_stack([
                maturity[mask], np.log(np.where(depth[mask] > 0, depth[mask], np.nan)),
            ])
            wnt_here = score[mask]
            base = {
                "study_id": STUDY, "patient_id": patient, "arm": arm,
                "granularity_rung": rung, "labeling_axis": AXIS,
                "n_cells": n, "n_signature_genes": score_report["n_genes_used"],
                # Invariant 8's REQUIRED report: is the Wnt score just maturity?
                "wnt_stem_correlation": partial_spearman(
                    wnt_here, maturity[mask],
                    np.log(np.where(depth[mask] > 0, depth[mask], np.nan))[:, None],
                ),
            }
            for gene, role in GENE_ROLES.items():
                hit = np.flatnonzero(symbols.to_numpy() == gene)
                if not hit.size:
                    continue
                column = kept[mask][:, hit[0]]
                values = np.asarray(column.todense()).ravel().astype(float)
                cp10k = values / np.where(depth[mask] > 0, depth[mask], np.nan) * 1e4
                out.append(base | {
                    "gene": gene, "role": role,
                    "partial_rho": partial_spearman(wnt_here, cp10k, conditioners),
                    "unconditioned_rho": partial_spearman(
                        wnt_here, cp10k, np.zeros((n, 1))),
                })
        log.info("  [%d/%d] %s scored", i, len(patients), patient)

    frame = pd.DataFrame(out)
    report["n_patients_scored"] = (
        int(frame["patient_id"].nunique()) if not frame.empty else 0)
    report["skipped"] = skipped
    return frame, report


def summarise(rows: pd.DataFrame, *, seed: int) -> pd.DataFrame:
    """Across patients, per (arm, gene). Student-t, patients as the unit."""
    rng = np.random.default_rng(seed)
    out = []
    for (arm, gene, role), block in rows.groupby(["arm", "gene", "role"],
                                                 observed=True):
        for column in ("partial_rho", "unconditioned_rho"):
            values = block[column].dropna().to_numpy()
            n = int(values.size)
            lo, hi = (student_t_interval(values, rng=rng) if n >= 2
                      else (float("nan"), float("nan")))
            out.append({
                "arm": arm, "gene": gene, "role": role, "statistic": column,
                "n_patients": n,
                "mean_rho": float(values.mean()) if n else float("nan"),
                "ci_low": lo, "ci_high": hi,
                "excludes_zero": bool(np.isfinite(lo) and (lo > 0 or hi < 0)),
            })
    return pd.DataFrame(out)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", type=Path, default=DEFAULT_ATLAS)
    parser.add_argument("--obs-cache", type=Path,
                        default=INTERIM_DIR / "icbi_obs.parquet")
    parser.add_argument("--rung", default=RUNG, choices=list(RUNGS))
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--max-patients", type=int, default=None)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    log.info("%s\nWNT MECHANISM TEST — %s, rung %s\n%s",
             "=" * 72, STUDY, args.rung, "=" * 72)
    log.info("pre-registered in docs/prereg_wnt_mechanism.md — read §5 first")
    log.info("  signature (invariant 8): %s", ", ".join(SIGNATURE))
    log.info("  conditioners: maturity score AND log depth. Neither optional —")
    log.info("  without them the correlation recovers maturity, or library size.")
    log.info("  ACTB and KRT8 are the technical floor, not a sanity check.")

    obs = load_obs(args.obs_cache, args.atlas)
    obs["tissue"] = obs["sample_type"].map(ARM_MAPS["adenoma"])
    rows, report = per_patient(args.atlas, obs, seed=args.seed, rung=args.rung,
                               max_patients=args.max_patients)
    if rows.empty:
        log.error("no patient produced a correlation")
        return 3

    summary = summarise(rows, seed=args.seed)

    log.info("\n%s\nIS THE WNT SCORE JUST MATURITY? (invariant 8's required report)\n%s",
             "=" * 72, "=" * 72)
    stem = rows.drop_duplicates(["patient_id", "arm"])["wnt_stem_correlation"]
    verdict = wnt_stem_verdict(float(stem.mean()))
    log.info("  mean partial r(Wnt, maturity) over %d patient-arms = %+.3f",
             len(stem), stem.mean())
    log.info("  %s — %s", verdict["verdict"], verdict["detail"])

    for arm in ("tumour", "normal"):
        block = summary[(summary["arm"] == arm)
                        & (summary["statistic"] == "partial_rho")]
        if block.empty:
            continue
        log.info("\n%s\n%s ARM — partial rho, conditioned on maturity and depth\n%s",
                 "=" * 72, "polyp" if arm == "tumour" else "normal", "=" * 72)
        order = ["ACTB", "KRT8", "EPCAM", "CDX2", "MS4A12", "GUCA2A"]
        block = block.set_index("gene").reindex([g for g in order
                                                 if g in set(block["gene"])])
        log.info("%s", block[["role", "n_patients", "mean_rho", "ci_low",
                              "ci_high", "excludes_zero"]].to_string())
        floor = block.loc[block["role"] == "control", "mean_rho"]
        if len(floor):
            log.info("\n  technical floor (controls): %+.3f to %+.3f. A target's "
                     "value means\n  only the amount by which it exceeds this.",
                     floor.min(), floor.max())

    meta = {
        "prereg": "docs/prereg_wnt_mechanism.md",
        "signature": list(SIGNATURE),
        "invariant_8": (
            "target signature, not CTNNB1/TCF7L2; ASCL2 and LGR5 excluded "
            "because the stem axis IS in play and they are its markers"
        ),
        "conditioners": ["maturity_score (stem_pole)", "log library depth"],
        "why_conditioned": (
            "the mature bin is a bin and not a point, so an unconditioned "
            "correlation recovers residual maturity; and two per-cell scores in "
            "one cell correlate through library size even after CP10K"
        ),
        "technical_floor": "ACTB and KRT8, scored through the identical path",
        "wnt_stem_verdict": verdict,
        "interval": "Student-t over patients (HANDOFF §3a), not the bootstrap",
        "rung": args.rung,
        "why_not_best4": (
            "best4's intrinsic arm carries no claim after Amendment 3; a "
            "mechanism test there would ask why something is true that has not "
            "been established"
        ),
        "report": report,
        "exploratory": True,
        "pre_registered": True,
        "what_this_cannot_show": (
            "direction. A within-cell association is equally consistent with "
            "Wnt suppressing the programme and with less-differentiated cells "
            "carrying more Wnt tone for other reasons."
        ),
    }
    for frame, name in ((rows, "wnt_mechanism_per_patient"),
                        (summary, "wnt_mechanism_summary")):
        log.info("wrote %s", write_versioned_table(
            frame, name, seed=args.seed, results_dir=args.results_dir,
            allow_dirty=args.allow_dirty, extra_meta=meta,
        ))
    return 0


if __name__ == "__main__":
    sys.exit(main())
