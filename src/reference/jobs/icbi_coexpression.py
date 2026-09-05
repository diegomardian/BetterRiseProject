"""The coexpression reading, per ICBI study. Piece 2 of the B path.

    python -m src.reference.jobs.icbi_coexpression --study Pelka_2021_Cell
    python -m src.reference.jobs.icbi_coexpression --all

Mirrors ``coexpression_silencing.gse178341_deltas`` exactly -- same QC, same
labeller, same ``rows_for_patient``, same ``summarise`` -- with a different
loader. Nothing about the statistic changes; only where the cells come from.

WHY PER STUDY, AND ONE AT A TIME. Invariant 4: estimate per study, never pool.
Beyond that it is a disk constraint. The 136 qualifying patients hold 1,398,508
cells across all compartments, which is ~9.85 GB gzipped against ~10 GB free.
So each study is extracted, analysed, written and DISCARDED before the next
begins; peak disk is the largest single study and nothing intermediate is kept.

THE QC POPULATION IS THE FULL PATIENT BLOCK, NOT THE EPITHELIUM.
``cell_qc_metrics`` computes per-batch MAD thresholds over every compartment,
and the epithelial subset is only the SCORING population. Extracting epithelium
alone would change the QC basis and make a Pelka comparison diverge from the
committed GSE178341 result for a reason that has nothing to do with the
adaptation being right or wrong.

TWO DECISIONS RECORDED RATHER THAN ASSUMED.

*The sorted-fraction filter.* GSE178341 filters ``PROCESSING_TYPE == UNSORTED``.
That column does not exist here; the analogue is
``enrichment_cell_types == "naive"``. It is not cosmetic -- Pelka is 210,667
naive against 130,019 CD45+, and a CD45-sorted fraction is immune-enriched by
construction, so leaving it in makes every compartment fraction a statement
about the sort rather than the tissue.

*The QC batch key.* ``sample_id``, matching GSE178341. The atlas does carry it
(median 3 samples per patient in Pelka), so no patient-by-tissue fallback is
needed. Recorded because if it ever were, the MAD thresholds would be computed
over a different unit and the numbers would move for that reason alone.
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
from src.reference.icbi_slice import (
    COUNTS_LAYER,
    SliceError,
    arms,
    assert_raw_counts,
    compartments,
    read_cells,
    read_var,
)
from src.reference.jobs.coexpression_silencing import (
    AXIS,
    CONTROL_LOG2_TOLERANCE,
    DEPTH_QUANTILE,
    GENE_ROLES,
    MATURE_BIN,
    MIN_CELLS_BOTH_ARMS,
    N_BOOTSTRAP,
    RUNG,
    premise_holds,
    rows_for_patient,
    summarise,
)
from src.reference.labels import label_column

log = logging.getLogger(__name__)

#: The unsorted fraction. GSE178341's UNSORTED, in this atlas's vocabulary.
NAIVE = "naive"

#: The QC batch unit. Same as GSE178341's.
BATCH_KEY = "sample_id"

#: Minimum epithelial cells per arm for a patient to be worth loading at all.
#: Matches the feasibility table that chose the 14 studies.
MIN_EPITHELIAL_PER_ARM = 100

DEFAULT_ATLAS = Path("/project/rise-batteries/bode/icbi/final_crc_atlas-adata.h5ad")


def load_obs(cache: Path, atlas: Path) -> pd.DataFrame:
    """The cached per-cell metadata, checked against the atlas it indexes.

    Row *i* of this frame must be row *i* of the atlas -- every slice is by
    positional index. A mismatched cache would hand every cell someone else's
    metadata, silently, so the lengths are compared rather than trusted.
    """
    import h5py

    if not cache.exists():
        raise SliceError(
            f"{cache} not found. Run\n"
            f"    python -m src.reference.jobs.pull_icbi_metadata\n"
            f"first -- it reads /obs over range requests and caches it."
        )
    obs = pd.read_parquet(cache)
    with h5py.File(str(atlas), "r") as h5:
        n_cells = int(len(h5[f"{COUNTS_LAYER}/indptr"][:]) - 1)
    if len(obs) != n_cells:
        raise SliceError(
            f"the cached obs has {len(obs):,} rows and the atlas has "
            f"{n_cells:,}. Every slice here is by POSITION, so a mismatch gives "
            f"each cell another cell's metadata. Re-pull the cache against this "
            f"atlas build."
        )
    return obs.reset_index(drop=True)


def eligible_patients(obs: pd.DataFrame, study_id: str) -> tuple[pd.DataFrame, list[str]]:
    """The naive, two-armed cells of one study, and the patients worth loading."""
    rows = obs[obs["study_id"] == study_id].copy()
    if rows.empty:
        raise SliceError(f"no cells for study {study_id!r}")

    before = len(rows)
    rows = rows[rows["enrichment_cell_types"].astype(str) == NAIVE]
    log.info("  %s fraction: %d of %d cells", NAIVE, len(rows), before)
    if rows.empty:
        return rows, []

    rows["tissue"] = arms(rows["sample_type"])
    rows["compartment"] = compartments(rows["atlas_cell_type_coarse"])
    rows = rows[rows["tissue"].notna()]

    epithelial = rows[rows["compartment"] == "epithelial"]
    per_arm = epithelial.groupby(["patient_id", "tissue"]).size().unstack(fill_value=0)
    for arm in ("normal", "tumour"):
        if arm not in per_arm.columns:
            per_arm[arm] = 0
    keep = per_arm[
        (per_arm["normal"] >= MIN_EPITHELIAL_PER_ARM)
        & (per_arm["tumour"] >= MIN_EPITHELIAL_PER_ARM)
    ]
    return rows, sorted(keep.index.astype(str))


def study_deltas(
    atlas: Path,
    obs: pd.DataFrame,
    study_id: str,
    *,
    seed: int = DEFAULT_SEED,
    max_patients: int | None = None,
) -> tuple[pd.DataFrame, dict]:
    """Per-patient, per-gene deltas for one study. Returns (deltas, report)."""
    from src.common.panel import tier_genes
    from src.reference.labels import assign_labels
    from src.reference.qc import apply_qc, cell_qc_metrics, qc_thresholds

    var = read_var(atlas)
    symbols = var["gene_symbol"]
    missing = sorted(set(GENE_ROLES) - set(symbols))
    if missing:
        # Skip-and-report: a study whose slice lacks a gene is a smaller
        # reading, not a failed run.
        log.warning("  genes absent from /var: %s", missing)

    rows, patients = eligible_patients(obs, study_id)
    if max_patients:
        patients = patients[:max_patients]
    report = {
        "study_id": study_id, "n_patients_eligible": len(patients),
        "genes_absent_from_var": missing, "batch_key": BATCH_KEY,
        "sorted_fraction_filter": f"enrichment_cell_types == {NAIVE!r}",
    }
    if not patients:
        log.warning("  no patient clears %d epithelial cells per arm",
                    MIN_EPITHELIAL_PER_ARM)
        return pd.DataFrame(), report

    log.info("  %d patient(s) clear %d epithelial cells per arm",
             len(patients), MIN_EPITHELIAL_PER_ARM)

    out: list[dict] = []
    for i, patient in enumerate(patients, 1):
        block_obs = rows[rows["patient_id"].astype(str) == patient]
        # ALL compartments: the QC population is the full patient block.
        row_index = block_obs.index.to_numpy()
        block = read_cells(atlas, row_index, n_genes=len(symbols))
        assert_raw_counts(block, context=f"{study_id}/{patient}")

        metrics = cell_qc_metrics(block, symbols, batch=block_obs[BATCH_KEY])
        keep = apply_qc(metrics, qc_thresholds(metrics)).to_numpy()
        if keep.sum() < 50:
            log.info("  [%d/%d] %s — %d cells survive QC, skipped",
                     i, len(patients), patient, int(keep.sum()))
            continue

        tissue = block_obs["tissue"].to_numpy()[keep]
        comp = block_obs["compartment"].to_numpy()[keep]
        kept = block[keep]
        labels = assign_labels(
            kept, symbols, compartment=comp,
            sample_id=block_obs[BATCH_KEY].to_numpy()[keep],
            target_genes=sorted(tier_genes("A")), tissue=tissue,
            patient_id=block_obs["patient_id"].to_numpy()[keep],
            depth_quantile=DEPTH_QUANTILE, seed=seed,
            index=pd.Index(block_obs.index.to_numpy()[keep]),
        )
        call = labels[label_column(AXIS, RUNG)].astype(str).to_numpy()
        mature = (call == MATURE_BIN) & np.isin(tissue, ["normal", "tumour"])
        if mature.sum() < MIN_CELLS_BOTH_ARMS:
            log.info("  [%d/%d] %s — %d mature cells, below the floor",
                     i, len(patients), patient, int(mature.sum()))
            continue

        mature_block = kept[mature]
        depth = np.asarray(mature_block.sum(axis=1), dtype=float).ravel()
        counts = {}
        for gene in GENE_ROLES:
            hit = np.where(symbols.to_numpy() == gene)[0]
            if len(hit) == 0:
                continue
            counts[gene] = np.asarray(
                mature_block[:, hit[0]].todense()
            ).ravel().astype(float)

        out.extend(rows_for_patient(
            study_id=study_id, patient=patient, counts=counts,
            depth=depth, tissue=tissue[mature], seed=seed,
        ))
        log.info("  [%d/%d] %s — %d mature cells scored",
                 i, len(patients), patient, int(mature.sum()))

    frame = pd.DataFrame(out)
    report["n_patients_scored"] = int(frame["patient_id"].nunique()) if not frame.empty else 0
    return frame, report


# ---------------------------------------------------------------------------
# The validation bar, pre-committed and checked in code


#: The committed GSE178341 result this run must reproduce. Pelka_2021_Cell IS
#: GSE178341, so the ICBI path has ground truth to be checked against -- the
#: only study of the fourteen that does.
#:
#: NOT bit-identical, and the bar says so: the atlas reprocessed every study
#: through its own QC, normalisation and doublet calling. What must survive is
#: the verdict and an overlapping interval.
VALIDATION = {
    "study_id": "Pelka_2021_Cell",
    "against": "results/2026-09-04_975cf5c/coexpression_silencing.parquet",
    "committed_study_id": "GSE178341",
    "requirements": [
        "the premise verdict word is the same (UNRESOLVED)",
        "the ACTB log2 control interval overlaps the committed one",
        "the GUCA2A detection delta is within 0.15 of the committed one",
    ],
}
GUCA2A_DELTA_TOLERANCE = 0.15


def control_log2_interval(
    deltas: pd.DataFrame, gene: str, *, seed: int = DEFAULT_SEED
) -> tuple[float, float, float]:
    """Bootstrap interval on a control's log2 ratio. Returns (mean, lo, hi).

    THE STATISTIC THE PREMISE ACTUALLY USES for a saturated control. ACTB and
    KRT8 sit at 0.99-1.00 detection in both arms, so `premise_holds` switches
    them to log2 expression against a 0.5 tolerance -- and their DETECTION
    interval is about +/-0.01 whatever happened, which makes an overlap test on
    it very nearly vacuous.

    Deliberately mirrors `premise_holds`' own bootstrap rather than parsing its
    reading string, and `test_the_log2_interval_matches_premise_holds` pins the
    two to agree. Two code paths that agree by review is how a validation bar
    ends up measuring something other than what it claims.
    """
    values = deltas.loc[deltas["gene"] == gene, "log2_cp10k_ratio"]
    values = values.dropna().to_numpy(dtype=float)
    if values.size < 3:
        return float("nan"), float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    draws = rng.choice(
        values, size=(N_BOOTSTRAP, values.size), replace=True
    ).mean(axis=1)
    lo, hi = (float(x) for x in np.percentile(draws, [2.5, 97.5]))
    return float(values.mean()), lo, hi


def verdict_word(reading: str) -> str:
    """The state at the front of a premise reading: UNRESOLVED / REFUSED / ...

    The rest of the string carries numbers that will not survive the atlas's
    reprocessing. The word is what has to match.
    """
    return str(reading).split(":", 1)[0].strip().upper()


def check_against_committed(deltas: pd.DataFrame, seed: int = DEFAULT_SEED) -> dict:
    """Compare a Pelka run against the committed GSE178341 result. Pass/fail.

    Runs IN the job rather than in someone's head afterwards -- the bar is in
    ``VALIDATION`` above, written before the number exists.
    """
    committed_path = RESULTS_DIR.parent / VALIDATION["against"]
    if not committed_path.exists():
        return {"verdict": "SKIPPED", "detail": f"{committed_path} not committed"}

    committed = pd.read_parquet(committed_path)
    committed = committed[committed["study_id"] == VALIDATION["committed_study_id"]]
    if committed.empty:
        return {"verdict": "SKIPPED",
                "detail": f"no {VALIDATION['committed_study_id']} rows to compare"}

    got = summarise(deltas, seed=seed)
    old = summarise(committed, seed=seed)
    checks, failures = {}, []

    for frame, label in ((old, "committed"), (got, "icbi")):
        if frame.empty:
            return {"verdict": "FAIL", "detail": f"{label} summary is empty"}

    def row(frame, gene, statistic):
        hit = frame[(frame["gene"] == gene) & (frame["statistic"] == statistic)]
        return hit.iloc[0] if len(hit) else None

    # (1) The premise verdict. THE ONE THAT MATTERS: it decides whether the
    # whole reading is interpretable, and it is the thing a detection-only
    # check cannot see moving.
    holds_new, reading_new = premise_holds(deltas, seed=seed)
    holds_old, reading_old = premise_holds(committed, seed=seed)
    checks["premise_icbi"] = reading_new
    checks["premise_committed"] = reading_old
    checks["premise_verdict_icbi"] = verdict_word(reading_new)
    checks["premise_verdict_committed"] = verdict_word(reading_old)
    if verdict_word(reading_new) != verdict_word(reading_old):
        failures.append(
            f"the premise verdict CHANGED: committed "
            f"{verdict_word(reading_old)}, icbi {verdict_word(reading_new)}. "
            f"That gate decides whether the reading is interpretable at all."
        )
    if holds_new != holds_old:
        failures.append(f"premise_holds returned {holds_new} against {holds_old}")

    # (2) ACTB on LOG2, which is the statistic the premise uses for it. Its
    # detection interval is ~+/-0.01 either way, so an overlap test there
    # cannot detect the shift that flips the verdict.
    mean_new, lo_new, hi_new = control_log2_interval(deltas, "ACTB", seed=seed)
    mean_old, lo_old, hi_old = control_log2_interval(committed, "ACTB", seed=seed)
    if not np.isfinite(lo_new) or not np.isfinite(lo_old):
        failures.append("ACTB log2 interval undefined on one side")
    else:
        overlap = (lo_new <= hi_old) and (lo_old <= hi_new)
        checks["actb_log2_overlaps"] = bool(overlap)
        checks["actb_log2_icbi"] = [round(mean_new, 3), round(lo_new, 3), round(hi_new, 3)]
        checks["actb_log2_committed"] = [round(mean_old, 3), round(lo_old, 3), round(hi_old, 3)]
        checks["control_log2_tolerance"] = CONTROL_LOG2_TOLERANCE
        if not overlap:
            failures.append(
                f"ACTB log2 intervals disjoint: icbi "
                f"{checks['actb_log2_icbi']} vs committed "
                f"{checks['actb_log2_committed']}"
            )

    g_new, g_old = row(got, "GUCA2A", "detection"), row(old, "GUCA2A", "detection")
    if g_new is None or g_old is None:
        failures.append("GUCA2A detection row missing from one side")
    else:
        drift = abs(float(g_new["mean_delta"]) - float(g_old["mean_delta"]))
        checks["guca2a_drift"] = drift
        checks["guca2a_icbi"] = float(g_new["mean_delta"])
        checks["guca2a_committed"] = float(g_old["mean_delta"])
        if drift > GUCA2A_DELTA_TOLERANCE:
            failures.append(
                f"GUCA2A delta moved {drift:.3f}, beyond the "
                f"{GUCA2A_DELTA_TOLERANCE} bar"
            )

    return {
        "verdict": "FAIL" if failures else "PASS",
        "detail": "; ".join(failures) if failures else "every requirement met",
        **checks,
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", type=Path, default=DEFAULT_ATLAS)
    parser.add_argument("--obs-cache", type=Path,
                        default=INTERIM_DIR / "icbi_obs.parquet")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--study")
    group.add_argument("--all", action="store_true")
    parser.add_argument("--max-patients", type=int, default=None)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    if not args.atlas.exists():
        raise SystemExit(
            f"{args.atlas} not found. Fetch it with\n"
            f"    qsub src/reference/jobs/fetch_icbi_atlas.sh\n"
            f"(30.44 GiB; it lives on /project, not under BRP_DATA_DIR)."
        )
    obs = load_obs(args.obs_cache, args.atlas)
    log.info("atlas %s | obs %s rows", args.atlas.name, f"{len(obs):,}")

    if args.all:
        candidates = sorted(RESULTS_DIR.glob("*/icbi_premise_candidate_studies.parquet"))
        if not candidates:
            raise SystemExit(
                "no icbi_premise_candidate_studies.parquet. Run\n"
                "    python -m src.reference.jobs.icbi_premise_feasibility"
            )
        studies = pd.read_parquet(candidates[-1])["study_id"].tolist()
    else:
        studies = [args.study]
    log.info("studies: %s", ", ".join(studies))

    frames, reports = [], []
    for study in studies:
        log.info("\n=== %s ===", study)
        try:
            deltas, report = study_deltas(
                args.atlas, obs, study, seed=args.seed,
                max_patients=args.max_patients,
            )
        except SliceError as exc:
            log.error("  REFUSED: %s", exc)
            reports.append({"study_id": study, "error": str(exc)})
            continue
        if not deltas.empty:
            holds, reading = premise_holds(deltas, seed=args.seed)
            report["premise_holds"] = bool(holds)
            report["premise_reading"] = reading
            log.info("  premise: %s", reading)
            frames.append(deltas)
        reports.append(report)

    if not frames:
        log.error("\nno study produced any scored patient.")
        return 3

    deltas = pd.concat(frames, ignore_index=True)
    summary = summarise(deltas, seed=args.seed)

    validation = {}
    if VALIDATION["study_id"] in set(deltas["study_id"]):
        subset = deltas[deltas["study_id"] == VALIDATION["study_id"]]
        validation = check_against_committed(subset, seed=args.seed)
        log.info("\n%s", "=" * 68)
        log.info("VALIDATION against the committed GSE178341 result: %s",
                 validation["verdict"])
        log.info("  %s", validation["detail"])
        for key in ("actb_icbi", "actb_committed", "guca2a_icbi",
                    "guca2a_committed", "guca2a_drift"):
            if key in validation:
                log.info("  %-18s %s", key, validation[key])
        log.info("%s", "=" * 68)

    for frame, name in ((deltas, "icbi_coexpression"),
                        (summary, "icbi_coexpression_summary")):
        path = write_versioned_table(
            frame, name, seed=args.seed,
            results_dir=args.results_dir, allow_dirty=args.allow_dirty,
            extra_meta={
                "atlas": args.atlas.name,
                "counts_layer": COUNTS_LAYER,
                "counts_layer_note": (
                    "/X is log1p-normalised on this atlas; reading it would "
                    "give log values to a detection statistic that reports "
                    "them without complaint."
                ),
                "studies": studies,
                "per_study_report": reports,
                "batch_key": BATCH_KEY,
                "sorted_fraction_filter": f"enrichment_cell_types == {NAIVE!r}",
                "qc_population": (
                    "the full patient block, ALL compartments -- the epithelial "
                    "subset is only the scoring population"
                ),
                "min_epithelial_per_arm": MIN_EPITHELIAL_PER_ARM,
                "validation": validation or VALIDATION,
                "premise_reading": {
                    r["study_id"]: r.get("premise_reading", "not scored")
                    for r in reports if "study_id" in r
                },
                "exploratory": True,
                "pre_registered": False,
                "what_this_is_not": (
                    "Confirmatory. Built after the week-0 falsification rule had "
                    "already forbidden a biological claim from the decomposition."
                ),
            },
        )
        log.info("wrote %s", path)

    if validation.get("verdict") == "FAIL":
        log.error("\nThe Pelka validation FAILED. Do not run the other thirteen "
                  "until the divergence is understood.")
        return 4
    return 0


if __name__ == "__main__":
    sys.exit(main())
