#!/usr/bin/env python
"""The split as a function of the maturity cut, not at four fixed rungs. §6.2.

**What §6.2 actually asks for:**

    Report the split as a curve across four resolutions plus a continuous
    maturation-score version. **A single point estimate would present a
    modelling choice as a measurement.**

The four rungs are four choices of where to cut the maturity score. Three of them
are disqualified, each for a reason tied to *that particular cut*:

    epithelial       the cut is at 0% — every scored cell is mature, so
                     Δfraction is identically zero
    crypt_position   the tertile boundaries collide on ~90% of patients (#42)
    best4            the cut is at 95%, so no patient clears the estimability
                     floor (G4, #48)

**So the disqualifications are properties of three specific thresholds, not of
the estimand.** This sweeps the cut continuously and reports the decomposition at
each, which shows where the estimate is stable and where it falls apart.

WHAT THIS IS NOT
----------------
**Not a new estimator.** `decompose()` is W4's and is called per row, unmodified.
Nothing here re-derives a term or a cutpoint. A genuinely continuous
decomposition — regressing expression on a maturation score — would be a
different estimand, and designing one now, having seen which rungs failed, is the
move this project refuses everywhere.

**Not a fifth rung.** The frozen schema pins `granularity_rung` to the four, and
this writes its own table with `threshold` as a column rather than pretending to
be one. The pre-committed thresholds are the rung definitions; this is context
around them.

**Not a menu.** A curve invites reading the most favourable point off it. The
threshold that counts was fixed at the rung definitions before any of this ran.
If the estimate only holds near one cut, that is the finding.

    qsub src/reference/jobs/run_threshold_sweep.sh
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import sparse  # noqa: E402

from src.common.io import write_versioned_table  # noqa: E402
from src.common.panel import tier_genes  # noqa: E402
from src.common.provenance import DEFAULT_SEED, set_global_seeds  # noqa: E402
from src.estimator.kitagawa import decompose  # noqa: E402
from src.harness.depth_confound import match_arm_depth  # noqa: E402
from src.harness.positivity import classify_estimability  # noqa: E402
from src.reference.ingest import (  # noqa: E402
    assign_compartments,
    read_gse178341,
    read_gse178341_clusters,
    read_gse178341_index,
    read_gse178341_metadata,
)
from src.reference.labels import (  # noqa: E402
    TRANSCRIPT_AXES,
    _positions,
    maturity_score,
)
from src.reference.qc import apply_qc, cell_qc_metrics, qc_thresholds  # noqa: E402

STUDY_ID = "GSE178341"
UNSORTED = "unsorted"
DEPTH_QUANTILE = 0.25
CP10K = 1e4
WEIGHTINGS = ("normal", "tumour", "doubly_robust")

#: Where the four rungs sit, so the curve can be read against them rather than
#: instead of them. `epithelial` is the cut at 0 and `best4` at 0.95.
RUNG_THRESHOLDS = {"epithelial": 0.0, "lineage": 0.50,
                   "crypt_position": 0.667, "best4": 0.95}

#: 5% to 95%. Below 5% the mature set is everything and Δfraction cannot move;
#: above 95% nothing clears the estimability floor. Both ends are included
#: deliberately — the curve should show its own failure modes.
THRESHOLDS = np.round(np.linspace(0.05, 0.95, 19), 3)

TARGET_TIERS = ("A", "B", "C", "D")


def main() -> int:
    set_global_seeds(DEFAULT_SEED)
    parser = argparse.ArgumentParser()
    parser.add_argument("--patients", nargs="*", default=None)
    parser.add_argument("--no-match-depth", action="store_true")
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args()
    match_depth = not args.no_match_depth

    data = Path(os.environ.get("BRP_DATA_DIR", "data")) / "raw" / "GSE178341"
    h5 = data / "GSE178341_crc10x_full_c295v4_submit.h5"
    clusters = read_gse178341_clusters(
        data / "GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz")
    metadata = read_gse178341_metadata(
        data / "GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz")
    obs, _var = read_gse178341_index(h5)
    patients = args.patients or sorted(obs["patient_id"].unique())
    genes = sorted({g for t in TARGET_TIERS for g in tier_genes(t)})
    print(f"{len(patients)} patients · {len(THRESHOLDS)} thresholds · "
          f"{len(genes)} genes · depth-matched: {match_depth}")

    rows = []
    for i, patient in enumerate(patients, 1):
        adata = read_gse178341(h5, patients=[patient])
        compartment = assign_compartments(clusters).reindex(adata.obs.index)
        metrics = cell_qc_metrics(adata.X, adata.var["gene_symbol"],
                                  batch=adata.obs["sample_id"])
        keep = apply_qc(metrics, qc_thresholds(metrics)).to_numpy()
        joined = adata.obs.join(metadata, how="left")
        keep &= (joined["PROCESSING_TYPE"] == UNSORTED).to_numpy()
        if keep.sum() < 50:
            continue

        tissue = adata.obs["tissue"].to_numpy()[keep]
        comp = compartment.to_numpy()[keep]
        epi = comp == "epithelial"
        is_n, is_t = tissue == "normal", tissue == "tumour"
        if not (epi & is_n).any() or not (epi & is_t).any():
            print(f"[{i}/{len(patients)}] {patient} — no paired epithelium")
            continue

        block = adata.X[keep]
        names = adata.var["gene_symbol"]
        totals = np.asarray(block.sum(axis=1), dtype=float).ravel()

        # The same depth floor the labels use: cells below it are not scored.
        target = float(np.quantile(totals[epi], DEPTH_QUANTILE))
        scorable = epi & (totals >= target)

        # Matched on the SCORED EPITHELIUM, per #24.1 and the correction in #57.
        matched = np.ones(len(tissue), bool)
        if match_depth:
            idx = np.flatnonzero(scorable)
            if idx.size and len(set(tissue[scorable].tolist())) == 2:
                km = match_arm_depth(totals[scorable], tissue[scorable],
                                     seed=DEFAULT_SEED)
                matched = np.zeros(len(tissue), bool)
                matched[idx[km]] = True
            else:
                matched = np.zeros(len(tissue), bool)
        usable = scorable & matched
        if usable.sum() < 50 or not (usable & is_n).any() or not (usable & is_t).any():
            print(f"[{i}/{len(patients)}] {patient} — too few after matching")
            continue

        positions, found = _positions(names, genes)
        if len(found) < len(genes):
            continue
        sub = block[:, positions]
        sub = sub.toarray() if sparse.issparse(sub) else np.asarray(sub)
        tot = totals.copy()
        tot[tot == 0] = np.nan
        expr = sub / tot[:, None] * CP10K

        for axis in TRANSCRIPT_AXES:
            score = maturity_score(
                block, names, axis, target_genes=sorted(tier_genes("A")),
                normalise=True, depth_target=target, seed=DEFAULT_SEED,
                totals=totals,
            )
            ref = usable & is_n
            for q in THRESHOLDS:
                # THE CUT COMES FROM THE NORMAL ARM and is applied as an
                # absolute value to the tumour — decision #13's reference-arm
                # rule. A within-arm quantile would pin Δfraction at q in both
                # arms by construction and the compositional term could not move.
                cut = float(np.quantile(score[ref], q))
                mature = usable & (score >= cut)
                den_n, den_t = (usable & is_n).sum(), (usable & is_t).sum()
                mat_n, mat_t = mature & is_n, mature & is_t
                frac_n, frac_t = mat_n.sum() / den_n, mat_t.sum() / den_t

                for j, gene in enumerate(found):
                    if not mat_n.any() or not mat_t.any():
                        continue
                    mn = float(np.nanmean(expr[mat_n, j]))
                    mt = float(np.nanmean(expr[mat_t, j]))
                    n_mat = int(mat_t.sum())
                    est = classify_estimability(n_mat)
                    for w in WEIGHTINGS:
                        d = decompose(frac_n, frac_t, mn, mt,
                                      n_cells_mature=n_mat, weighting=w)
                        rows.append({
                            "patient_id": patient, "study_id": STUDY_ID,
                            "gene": gene, "labeling_axis": axis,
                            "threshold": float(q), "weighting": w,
                            "compositional": d.compositional,
                            # Invariant 1 at the source: not_estimable means
                            # None here, never 0.0, before anything downstream
                            # gets a chance to coerce it.
                            "intrinsic": d.intrinsic if est != "not_estimable" else None,
                            "interaction": d.interaction,
                            "estimability": est,
                            "frac_mature_normal": float(frac_n),
                            "frac_mature_tumour": float(frac_t),
                            "n_cells_mature": n_mat,
                            "depth_matched": bool(match_depth),
                        })
        print(f"[{i}/{len(patients)}] {patient} — {len(THRESHOLDS)} thresholds")

    if not rows:
        raise SystemExit("no patient produced a sweep")
    out = pd.DataFrame(rows)

    print("\n" + "=" * 68)
    print("THE SPLIT AS A FUNCTION OF THE CUT — GUCA2A, doubly robust")
    print("=" * 68)
    g = out[(out.gene == "GUCA2A") & (out.weighting == "doubly_robust")
            & (out.estimability != "not_estimable")]
    if len(g):
        curve = g.groupby(["labeling_axis", "threshold"]).agg(
            n_patients=("patient_id", "nunique"),
            compositional=("compositional", "median"),
            intrinsic=("intrinsic", "median"),
        ).round(3)
        print(curve.to_string())
    print("\n  rung thresholds for reference: " +
          ", ".join(f"{k} {v}" for k, v in RUNG_THRESHOLDS.items()))
    print("\n  A CURVE IS NOT A MENU. The thresholds that count were fixed at the")
    print("  rung definitions before any of this ran; this is context around")
    print("  them, and reading the most favourable point off it would be the")
    print("  move the pre-registration exists to prevent.")

    path = write_versioned_table(
        out, "threshold_sweep", seed=DEFAULT_SEED, allow_dirty=args.allow_dirty,
        notes=(
            "§6.2's curve across resolutions: the Kitagawa split recomputed at 19 "
            "maturity cut points from 0.05 to 0.95, rather than at the four "
            "frozen rungs. decompose() is W4's and is called per row, unmodified "
            "— this is not a new estimator and not a fifth rung. Cuts come from "
            "each patient's normal arm and are applied absolutely to the tumour "
            "(#13). Depth-matched on the scored epithelium (#24.1). Written as "
            "its own table because the frozen schema pins granularity_rung to "
            "the four rungs, and a threshold is not one. The pre-committed "
            "thresholds are the rung definitions; this is context around them, "
            "not a menu to choose from."
        ),
        extra_meta={
            "thresholds": [float(x) for x in THRESHOLDS],
            "rung_thresholds": RUNG_THRESHOLDS,
            "depth_matched": bool(match_depth),
            "n_patients": int(out.patient_id.nunique()),
        },
    )
    print(f"\nwrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
