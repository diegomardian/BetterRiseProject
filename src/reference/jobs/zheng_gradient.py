"""Descriptive normal-to-polyp-to-carcinoma trajectories in Zheng_2022.

Run on the cluster holding the ICBI atlas:

    python -m src.reference.jobs.zheng_gradient \
      --atlas /project/rise-batteries/bode/icbi/final_crc_atlas-adata.h5ad

This intentionally writes three individual patient trajectories.  It performs
no pooled inference, interval, stage contrast, premise verdict, or model.
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
from src.common.label_provenance import Measurement, check_no_circular_claim, provenance_meta
from src.common.paths import INTERIM_DIR
from src.common.provenance import DEFAULT_SEED
from src.reference.icbi_slice import (
    assert_raw_counts,
    compartments,
    read_cells,
    read_var,
)
from src.reference.jobs.coexpression_silencing import DETECTION_MIN_UMI, GENE_ROLES
from src.reference.jobs.icbi_coexpression import BATCH_KEY, DEFAULT_ATLAS, NAIVE, load_obs
from src.reference.qc import apply_qc, cell_qc_metrics, qc_thresholds

log = logging.getLogger(__name__)

STUDY_ID = "Zheng_2022_Signal_Transduct_Target_Ther"
STAGES = ("normal", "polyp", "carcinoma")
STAGE_MAP = {
    "adjacent normal": "normal",
    "healthy normal": "normal",
    "polyp": "polyp",
    "primary tumor": "carcinoma",
}
MIN_EPITHELIAL_PER_STAGE = 100
EXPECTED_PATIENTS = 3
PANEL = tuple(GENE_ROLES)

LABEL_PROVENANCE = Measurement(
    modality="sample_annotation",
    assay="ICBI coarse cell-type and sample-type annotations",
    genes=(),
)
CLAIM_PROVENANCE = Measurement(
    modality="transcript",
    assay="ICBI raw UMI counts layer",
    genes=PANEL,
)


class ZhengGradientError(ValueError):
    """The fixed Zheng gradient input no longer has its verified shape."""


def validate_specification() -> tuple[str, ...]:
    """Refuse label leakage before opening the atlas counts layer."""
    return check_no_circular_claim(labels=LABEL_PROVENANCE, claim=CLAIM_PROVENANCE)


def panel_indices(symbols: Sequence[str]) -> dict[str, int]:
    """Locate the frozen panel exactly once; a missing symbol is not a zero."""
    values = np.asarray([str(value) for value in symbols])
    index: dict[str, int] = {}
    missing: list[str] = []
    duplicated: list[str] = []
    for gene in PANEL:
        hits = np.flatnonzero(values == gene)
        if hits.size == 0:
            missing.append(gene)
        elif hits.size > 1:
            duplicated.append(gene)
        else:
            index[gene] = int(hits[0])
    if missing or duplicated:
        raise ZhengGradientError(
            "frozen panel lookup refused: "
            + (f"missing symbols {missing}. " if missing else "")
            + (f"duplicated symbols {duplicated}." if duplicated else "")
        )
    return index


def eligible_gradient_rows(obs: pd.DataFrame) -> tuple[pd.DataFrame, tuple[str, ...]]:
    """Return exactly the pre-specified three-stage patient cohort."""
    rows = obs[
        (obs["study_id"].astype(str) == STUDY_ID)
        & (obs["enrichment_cell_types"].astype(str) == NAIVE)
    ].copy()
    if rows.empty:
        raise ZhengGradientError(f"no naive cells for {STUDY_ID}")
    rows["stage"] = rows["sample_type"].astype(str).map(STAGE_MAP)
    rows["compartment"] = compartments(rows["atlas_cell_type_coarse"])
    rows = rows[rows["stage"].notna()]
    epithelial = rows[rows["compartment"] == "epithelial"]
    counts = epithelial.groupby(["patient_id", "stage"], observed=True).size().unstack(fill_value=0)
    for stage in STAGES:
        if stage not in counts:
            counts[stage] = 0
    patients = tuple(sorted(
        counts.index[
            (counts[list(STAGES)] >= MIN_EPITHELIAL_PER_STAGE).all(axis=1)
        ].astype(str)
    ))
    if len(patients) != EXPECTED_PATIENTS:
        raise ZhengGradientError(
            f"{len(patients)} patients with all three stages at >= "
            f"{MIN_EPITHELIAL_PER_STAGE} epithelial cells, expected "
            f"{EXPECTED_PATIENTS}. Do not change the cohort after inspection."
        )
    return rows[rows["patient_id"].astype(str).isin(patients)], patients


def _values(matrix, column: int) -> np.ndarray:
    selected = matrix[:, column]
    return (
        np.asarray(selected.toarray()).ravel()
        if hasattr(selected, "toarray") else np.asarray(selected).ravel()
    )


def summarise_patient(
    *, patient_id: str, stage: np.ndarray, compartment: np.ndarray,
    matrix, metrics: pd.DataFrame, keep: np.ndarray, indices: dict[str, int],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Emit all stage rows for one patient, including post-QC zeroes."""
    output: list[dict[str, object]] = []
    availability: list[dict[str, object]] = []
    epithelial = compartment == "epithelial"
    for current in STAGES:
        before = epithelial & (stage == current)
        after = before & keep
        depth = metrics.loc[after, "n_counts"].to_numpy(dtype=float)
        n_genes = metrics.loc[after, "n_genes"].to_numpy(dtype=float)
        availability.append({
            "study_id": STUDY_ID,
            "patient_id": patient_id,
            "stage": current,
            "n_epithelial_before_qc": int(before.sum()),
            "n_epithelial_after_qc": int(after.sum()),
            "n_samples_before_qc": int(metrics.loc[before, "batch"].nunique()),
            "n_samples_after_qc": int(metrics.loc[after, "batch"].nunique()),
            "median_counts_after_qc": float(np.median(depth)) if depth.size else float("nan"),
            "median_genes_after_qc": float(np.median(n_genes)) if n_genes.size else float("nan"),
            "stage_measurable_after_qc": bool(depth.size > 0),
        })
        for gene, column in indices.items():
            values = _values(matrix[after], column)
            output.append({
                "study_id": STUDY_ID,
                "patient_id": patient_id,
                "stage": current,
                "gene": gene,
                "role": GENE_ROLES[gene],
                "n_epithelial_before_qc": int(before.sum()),
                "n_epithelial_after_qc": int(after.sum()),
                "n_detected": int((values >= DETECTION_MIN_UMI).sum()),
                "detection": (
                    float((values >= DETECTION_MIN_UMI).mean())
                    if values.size else float("nan")
                ),
                "cp10k_mean": (
                    float(np.mean(values / depth * 1e4)) if depth.size else float("nan")
                ),
                "median_counts_after_qc": (
                    float(np.median(depth)) if depth.size else float("nan")
                ),
                "median_genes_after_qc": (
                    float(np.median(n_genes)) if n_genes.size else float("nan")
                ),
                "stage_measurable_after_qc": bool(depth.size > 0),
            })
    return output, availability


def run_gradient(atlas: Path, obs: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Read the fixed three-patient cohort one patient at a time."""
    rows, patients = eligible_gradient_rows(obs)
    var = read_var(atlas)
    indices = panel_indices(var["gene_symbol"])
    trajectories: list[dict[str, object]] = []
    availability: list[dict[str, object]] = []
    for number, patient in enumerate(patients, start=1):
        block_obs = rows[rows["patient_id"].astype(str) == patient]
        matrix = read_cells(atlas, block_obs.index.to_numpy(), n_genes=len(var))
        assert_raw_counts(matrix, context=f"{STUDY_ID}/{patient}")
        metrics = cell_qc_metrics(matrix, var["gene_symbol"], batch=block_obs[BATCH_KEY])
        keep = apply_qc(metrics, qc_thresholds(metrics)).to_numpy()
        rows_out, availability_out = summarise_patient(
            patient_id=patient,
            stage=block_obs["stage"].to_numpy(),
            compartment=block_obs["compartment"].to_numpy(),
            matrix=matrix,
            metrics=metrics,
            keep=keep,
            indices=indices,
        )
        trajectories.extend(rows_out)
        availability.extend(availability_out)
        log.info("[%d/%d] %s — %d/%d cells survive QC", number, len(patients), patient,
                 int(keep.sum()), len(keep))
    report = {
        "study_id": STUDY_ID,
        "patients": list(patients),
        "n_patients": len(patients),
        "stages": list(STAGES),
        "min_epithelial_per_stage_before_qc": MIN_EPITHELIAL_PER_STAGE,
    }
    return pd.DataFrame(trajectories), pd.DataFrame(availability), report


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", type=Path, default=DEFAULT_ATLAS)
    parser.add_argument("--obs-cache", type=Path, default=INTERIM_DIR / "icbi_obs.parquet")
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    validate_specification()
    if not args.atlas.exists():
        raise SystemExit(f"{args.atlas} not found")

    obs = load_obs(args.obs_cache, args.atlas)
    trajectories, availability, report = run_gradient(args.atlas, obs)
    meta = {
        "prereg": "docs/prereg_zheng_gradient.md",
        "report": report,
        "counts_layer": "layers/counts",
        "panel": list(PANEL),
        "detection_min_umi": DETECTION_MIN_UMI,
        "population": (
            "naive, coarse-annotation epithelial cells; no transcript-derived mature label"
        ),
        "inference": (
            "none: individual trajectories only; no pooling, contrast, interval, test, or model"
        ),
        "pre_registered": True,
        "exploratory": False,
        **provenance_meta(LABEL_PROVENANCE, CLAIM_PROVENANCE),
    }
    for frame, name in (
        (trajectories, "zheng_three_stage_trajectories"),
        (availability, "zheng_three_stage_availability"),
    ):
        path = write_versioned_table(
            frame, name, seed=args.seed, results_dir=args.results_dir,
            allow_dirty=args.allow_dirty, extra_meta=meta,
        )
        log.info("wrote %s", path)
    log.info("DESCRIPTIVE ONLY — three patient trajectories; no pooled inference.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
