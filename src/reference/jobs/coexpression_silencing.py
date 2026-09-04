"""Is the marker gone because the cells are gone, or because the cells went quiet?

WHY THIS EXISTS. The decomposition in ``src/estimator/kitagawa.py`` cannot
answer that question on this panel, and the reason is algebraic rather than
statistical. Writing ``i/c`` for the ratio of the intrinsic to the compositional
term,

    i / c = (f_N / delta_f) * (m_T/m_N - 1)

so as a gene's surviving per-cell mean goes to zero the bracket goes to -1 and
the ratio collapses onto ``-(f_N / delta_f)`` -- a property of the *cell
fractions*, identical for every gene scored on the same labels. On the primary
cohort that constant is -5.85, and GUCA2A (5.67), GUCA2B (5.80), OTOP2 (5.85),
CA7 (5.83) and MS4A12 (5.58) are not distinguishable by it. Tier A was
pre-registered as compositional and tier D as neither; both come back ~99%
intrinsic. That is the week-0 falsification rule firing, and this project's
pre-committed response to it stands: **the decomposition supports no
gene-specific mechanism claim.**

WHAT THIS MEASURES INSTEAD. Not a variance split. A per-cell detection rate,
inside a population fixed *before* either gene is looked at, compared between
arms of the same patient at matched sequencing depth. If a marker's cells were
destroyed, the cells that remain and still read as epithelial should carry the
marker at the rate they always did. If the cells are present and quiet, the
marker falls while everything else about them does not. None of the algebra
above applies, because nothing here is divided by ``delta_f``.

THE PREMISE IS A CONTROL, NOT AN ASSUMPTION. The reading only means anything if
the cells scored in the diseased arm really are still epithelial cells. So
housekeeping and structural genes are scored the same way in the same cells, and
:func:`premise_holds` **refuses the comparison** if they move: a population whose
ACTB has shifted is not the same population, and a marker falling inside it says
nothing. A check that cannot fail is worse than no check, and the failure mode
of this one is a silently different population.

WHAT IT IS NOT. **Exploratory, and post-hoc.** It was not pre-registered, it was
built after the falsification rule had already fired, and reaching a mechanism
claim by a second route after a pre-committed rule forbade the first is exactly
what pre-registration exists to catch. Nothing here is confirmatory. It is a
hypothesis with a measurement attached, and it should be described that way
wherever it is cited. It also does not show silencing is the *only* mechanism:
identity markers fall too, so real de-differentiation is happening alongside.

    python -m src.reference.jobs.coexpression_silencing
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
from src.harness.depth_confound import match_arm_depth
from src.reference.labels import label_column

log = logging.getLogger(__name__)

#: Set from --allow-dirty. Default False.
ALLOW_DIRTY = False

DEFAULT_SEED = 20260904

#: A cell counts as expressing a gene at one UMI or more. Detection rather than
#: mean expression, because detection is what survives a depth difference once
#: the arms are matched, and because a mean over mostly-zero counts is dominated
#: by the few cells that fired.
DETECTION_MIN_UMI = 1

#: Cells a patient needs across both arms before depth matching, and per arm
#: after it. Matching discards cells, so both floors are checked.
MIN_CELLS_BOTH_ARMS = 40
MIN_CELLS_PER_ARM = 15

#: CDX2-positive cells an arm needs before the conditional reading is reported.
MIN_CONDITIONING_CELLS = 10

N_BOOTSTRAP = 10_000

#: The gene tiers, and what each one is doing here. ``control`` is load-bearing:
#: if those move, the population is not the same and nothing else on this table
#: can be read.
GENE_ROLES: dict[str, str] = {
    "ACTB": "control",
    "KRT8": "control",
    "EPCAM": "epithelial",
    "CDX2": "identity",
    "MS4A12": "identity",
    "GUCA2A": "target",
}

#: How far a control gene may move before the premise is refused. Detection
#: rates are proportions, so this is in percentage points.
CONTROL_TOLERANCE = 0.10

#: Patients the premise needs before "the controls did not move" means anything.
#: On one patient a control's mean IS that patient, and a tolerance compared
#: against a single noisy number reports clean whatever the data does -- which
#: is the failure this project is named after. Below this the premise is
#: UNDEFINED, not satisfied.
MIN_PREMISE_PATIENTS = 3

AXIS, RUNG, MATURE_BIN = "stem_pole", "lineage", "differentiated"

#: GSE178341 QC settings, identical to build_decomposition_summary's so the
#: population scored here is the one the decomposition is computed on.
UNSORTED = "unsorted"
DEPTH_QUANTILE = 0.25
CONDITION_ON = "CDX2"
COHORTS = ("smc", "kul3")

#: The primary cohort. Cluster only -- 371k cells, and the loader streams one
#: patient at a time because the whole deposit does not belong in memory.
PRIMARY = "gse178341"


def premise_holds(
    deltas: pd.DataFrame, tolerance: float = CONTROL_TOLERANCE
) -> tuple[bool, str]:
    """Whether the diseased cells are still the same kind of cell.

    The whole reading rests on this. A marker falling inside a population that
    has itself changed is not silencing; it is a different population. So the
    control genes are scored in the same cells and must not have moved.

    Takes the per-patient rows, not a pre-aggregated mean, because how many
    patients stand behind a control is part of whether it says anything: on one
    patient the mean is that patient, and comparing a single noisy number
    against a tolerance reports clean whatever happened.

    Returns ``(holds, reading)``. When it does not hold the caller must not
    report a mechanism, and the reading says why.
    """
    controls = [g for g, role in GENE_ROLES.items() if role == "control"]
    present = deltas.loc[deltas["gene"].isin(controls)]
    if present.empty:
        return False, "UNDEFINED: no control gene was scored, so the premise is untested"

    if "patient_id" in present.columns:
        n_patients = int(present["patient_id"].nunique())
    else:
        n_patients = int(present.groupby("gene").size().max())
    if n_patients < MIN_PREMISE_PATIENTS:
        return False, (
            f"UNDEFINED: {n_patients} patient(s) behind the controls, below the "
            f"{MIN_PREMISE_PATIENTS} needed to tell a shift from noise. The "
            f"premise is untested here, not satisfied."
        )

    means = present.groupby("gene", as_index=False)["delta_detect"].mean()
    worst = means.reindex(means["delta_detect"].abs().sort_values().index).iloc[-1]
    if abs(worst["delta_detect"]) > tolerance:
        return False, (
            f"REFUSED: control gene {worst['gene']} moved by "
            f"{worst['delta_detect']:+.3f}, beyond the {tolerance} tolerance. The "
            f"cells scored in the two arms are not the same population, so a "
            f"marker falling inside them says nothing about silencing."
        )
    return True, (
        f"holds over {n_patients} patients: worst control ({worst['gene']}) moved "
        f"{worst['delta_detect']:+.3f}, within {tolerance}"
    )


def _detection(counts: np.ndarray) -> float:
    return float((counts >= DETECTION_MIN_UMI).mean())


def rows_for_patient(
    *, study_id: str, patient: str, counts: dict[str, np.ndarray],
    depth: np.ndarray, tissue: np.ndarray, seed: int,
) -> list[dict]:
    """One row per gene for one patient, arms depth-matched.

    ``counts`` holds raw integer counts per gene, already restricted to the
    mature cells of this patient in both arms. Matching happens on exactly that
    set: matching a wider population and intersecting afterwards leaves the
    analysed subset unmatched, which `build_decomposition_summary` learned the
    hard way and records in its own comments.
    """
    if len(depth) < MIN_CELLS_BOTH_ARMS or len(set(tissue.tolist())) < 2:
        return []
    keep = match_arm_depth(depth, tissue, seed=seed)
    arm, depth = tissue[keep], depth[keep]
    if min((arm == "normal").sum(), (arm == "tumour").sum()) < MIN_CELLS_PER_ARM:
        return []

    conditioner = counts[CONDITION_ON][keep] >= DETECTION_MIN_UMI
    out = []
    for gene, values in counts.items():
        v = values[keep]
        row = {"study_id": study_id, "patient_id": str(patient),
               "gene": gene, "role": GENE_ROLES[gene]}
        for name in ("normal", "tumour"):
            m = arm == name
            row[f"n_{name}"] = int(m.sum())
            row[f"depth_{name}"] = float(np.median(depth[m]))
            row[f"detect_{name}"] = _detection(v[m])
            both = m & conditioner
            row[f"n_conditioned_{name}"] = int(both.sum())
            row[f"detect_given_{CONDITION_ON.lower()}_{name}"] = (
                _detection(v[both]) if both.sum() >= MIN_CONDITIONING_CELLS else np.nan
            )
        row["depth_ratio"] = row["depth_tumour"] / row["depth_normal"]
        row["delta_detect"] = row["detect_tumour"] - row["detect_normal"]
        row["delta_given_conditioner"] = (
            row[f"detect_given_{CONDITION_ON.lower()}_tumour"]
            - row[f"detect_given_{CONDITION_ON.lower()}_normal"]
        )
        out.append(row)
    return out


def per_patient_deltas(cohort, *, seed: int = DEFAULT_SEED) -> pd.DataFrame:
    """The Lee cohorts, which carry cells, labels and raw counts in one object."""
    cells = cohort.cells
    labels = cohort.labels[label_column(AXIS, RUNG)].reindex(cells.index).astype(str)
    raw = cohort.raw_counts.reindex(cells.index)
    depth = cohort.n_counts.reindex(cells.index)
    genes = [g for g in GENE_ROLES if g in raw.columns]

    rows: list[dict] = []
    for patient in sorted(cells["patient_id"].unique()):
        sel = (
            (labels == MATURE_BIN)
            & (cells["patient_id"] == patient)
            & cells["tissue"].isin(["normal", "tumour"])
        ).to_numpy()
        rows.extend(rows_for_patient(
            study_id=cohort.study_id, patient=patient,
            counts={g: raw[g].to_numpy(dtype=float)[sel] for g in genes},
            depth=depth.to_numpy(dtype=float)[sel],
            tissue=cells["tissue"].to_numpy()[sel], seed=seed,
        ))
    return pd.DataFrame(rows)


def gse178341_deltas(*, seed: int = DEFAULT_SEED, patients=None, data_dir=None) -> pd.DataFrame:
    """The primary cohort, streamed one patient at a time.

    Thirty-two patients against ten and six, so this is the read that decides
    whether the two Lee cohorts disagree because they are small or because they
    disagree. It is also the only one of the three that needs the cluster: the
    deposit is 371k cells and the loader materialises one patient at a time for
    exactly that reason.

    The QC and labelling path is `build_decomposition_summary`'s, unchanged, so
    the population scored here is the one the decomposition is computed on.
    """
    import os

    from scipy import sparse

    from src.common.panel import tier_genes
    from src.reference.ingest import (
        assign_compartments,
        read_gse178341,
        read_gse178341_clusters,
        read_gse178341_index,
        read_gse178341_metadata,
    )
    from src.reference.labels import assign_labels
    from src.reference.qc import apply_qc, cell_qc_metrics, qc_thresholds

    base = Path(data_dir or os.environ.get("BRP_DATA_DIR", "data")) / "raw" / "GSE178341"
    h5 = base / "GSE178341_crc10x_full_c295v4_submit.h5"
    clusters = read_gse178341_clusters(
        base / "GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz"
    )
    metadata = read_gse178341_metadata(
        base / "GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz"
    )
    obs, _ = read_gse178341_index(h5)
    todo = patients or sorted(obs["patient_id"].unique())

    rows: list[dict] = []
    for i, patient in enumerate(todo, 1):
        adata = read_gse178341(h5, patients=[patient])
        compartment = assign_compartments(clusters).reindex(adata.obs.index)
        metrics = cell_qc_metrics(adata.X, adata.var["gene_symbol"],
                                  batch=adata.obs["sample_id"])
        keep = apply_qc(metrics, qc_thresholds(metrics)).to_numpy()
        keep &= (adata.obs.join(metadata, how="left")["PROCESSING_TYPE"] == UNSORTED).to_numpy()
        tissue = adata.obs["tissue"].to_numpy()[keep]
        comp = compartment.to_numpy()[keep]
        if keep.sum() < 50 or not ((comp == "epithelial") & (tissue == "normal")).any():
            log.info("[%d/%d] %s — skipped (no reference arm or too few cells)",
                     i, len(todo), patient)
            continue

        names = adata.var["gene_symbol"]
        block = adata.X[keep]
        labels = assign_labels(
            block, names, compartment=comp,
            sample_id=adata.obs["sample_id"].to_numpy()[keep],
            target_genes=sorted(tier_genes("A")), tissue=tissue,
            patient_id=adata.obs["patient_id"].to_numpy()[keep],
            depth_quantile=DEPTH_QUANTILE, seed=seed, index=adata.obs.index[keep],
        )
        call = labels[label_column(AXIS, RUNG)].astype(str).to_numpy()
        mature = (call == MATURE_BIN) & np.isin(tissue, ["normal", "tumour"])
        if mature.sum() < MIN_CELLS_BOTH_ARMS:
            log.info("[%d/%d] %s — %d mature cells, below the floor",
                     i, len(todo), patient, int(mature.sum()))
            continue

        depth = np.asarray(block.sum(axis=1), dtype=float).ravel()
        counts = {}
        for gene in GENE_ROLES:
            hit = np.where(names.to_numpy() == gene)[0]
            if len(hit) == 0:
                continue
            column = block[:, hit[0]]
            column = (column.toarray().ravel() if sparse.issparse(column)
                      else np.asarray(column).ravel())
            counts[gene] = column[mature]
        if CONDITION_ON not in counts:
            log.info("[%d/%d] %s — %s absent, cannot condition",
                     i, len(todo), patient, CONDITION_ON)
            continue

        found = rows_for_patient(
            study_id="GSE178341", patient=patient, counts=counts,
            depth=depth[mature], tissue=tissue[mature], seed=seed,
        )
        log.info("[%d/%d] %s — %d mature cells, %d gene rows", i, len(todo), patient,
                 int(mature.sum()), len(found))
        rows.extend(found)
    return pd.DataFrame(rows)


def summarise(deltas: pd.DataFrame, *, seed: int = DEFAULT_SEED) -> pd.DataFrame:
    """Per (study, gene): the paired change, bootstrapped over PATIENTS.

    Per study and never pooled (CLAUDE.md invariant 4), and the resample is over
    patients rather than cells (invariant 5) because the patient is the unit of
    inference. Two studies is not a meta-analysis and none is attempted here.
    """
    rng = np.random.default_rng(seed)
    out: list[dict] = []
    for (study, gene), block in deltas.groupby(["study_id", "gene"], sort=True):
        for column, label in (("delta_detect", "detection"),
                              ("delta_given_conditioner", f"detection|{CONDITION_ON}+")):
            values = block[column].dropna().to_numpy(dtype=float)
            row = {
                "study_id": study, "gene": gene, "role": GENE_ROLES.get(gene, "?"),
                "statistic": label, "n_patients": int(len(values)),
                "mean_delta": float(values.mean()) if len(values) else np.nan,
            }
            if len(values) >= 3:
                draws = rng.choice(
                    values, size=(N_BOOTSTRAP, len(values)), replace=True
                ).mean(axis=1)
                lo, hi = np.percentile(draws, [2.5, 97.5])
                row["ci_low"], row["ci_high"] = float(lo), float(hi)
                row["excludes_zero"] = bool(row["ci_low"] * row["ci_high"] > 0)
            else:
                row["ci_low"] = row["ci_high"] = np.nan
                row["excludes_zero"] = False
            out.append(row)
    return pd.DataFrame(out)


def main(argv: Sequence[str] | None = None) -> int:
    global ALLOW_DIRTY
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohorts", nargs="+", default=list(COHORTS),
                        choices=[*COHORTS, "gse178341"])
    parser.add_argument("--patients", nargs="*", default=None)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--raw-dir", type=Path, default=None)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument(
        "--allow-dirty", action="store_true",
        help="write from a dirty tree; the recorded sha will not reproduce it",
    )
    args = parser.parse_args(argv)
    ALLOW_DIRTY = args.allow_dirty
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    from src.estimator.lee_io import load_lee_cohort

    frames = []
    for which in args.cohorts:
        if which == "gse178341":
            log.info("streaming GSE178341, one patient at a time …")
            frames.append(gse178341_deltas(seed=args.seed, patients=args.patients,
                                           data_dir=args.raw_dir))
            continue
        log.info("loading Lee/%s …", which.upper())
        cohort = load_lee_cohort(
            which,
            target_genes=["GUCA2A"],
            axes=(AXIS,),
            rungs=(RUNG,),
            raw_dir=args.raw_dir,
            extra_genes=[g for g in GENE_ROLES if g != "GUCA2A"],
            keep_raw_counts=True,
        )
        frames.append(per_patient_deltas(cohort, seed=args.seed))
    deltas = pd.concat(frames, ignore_index=True)
    if deltas.empty:
        log.error("no patient met the cell floors — nothing to report")
        return 1

    summary = summarise(deltas, seed=args.seed)

    readings: dict[str, str] = {}
    for study, block in deltas.groupby("study_id"):
        holds, reading = premise_holds(block)
        readings[str(study)] = reading
        log.info("%s premise %s", study, reading)
        detect = summary[(summary.study_id == study) & (summary.statistic == "detection")]
        for _, r in detect.sort_values("mean_delta").iterrows():
            log.info("  %-8s %-11s %+0.3f  [%+0.3f, %+0.3f]%s",
                     r.gene, r.role, r.mean_delta, r.ci_low, r.ci_high,
                     "" if holds else "   (premise refused — not a mechanism claim)")

    meta = {
        "exploratory": True,
        "pre_registered": False,
        "what_this_is_not": (
            "Post-hoc. Built after the week-0 falsification rule had already "
            "fired and forbidden a biological claim from the decomposition. "
            "Reaching a mechanism claim by a second route afterwards is what "
            "pre-registration exists to catch, so nothing here is confirmatory."
        ),
        "population": f"label_{AXIS}_{RUNG} == {MATURE_BIN!r}, depth-matched within patient",
        "detection_min_umi": DETECTION_MIN_UMI,
        "conditioned_on": CONDITION_ON,
        "control_tolerance": CONTROL_TOLERANCE,
        "gene_roles": GENE_ROLES,
        "premise_reading": readings,
        "n_bootstrap": N_BOOTSTRAP,
    }
    for frame, name in ((deltas, "coexpression_silencing"),
                        (summary, "coexpression_silencing_summary")):
        path = write_versioned_table(
            frame, name, seed=args.seed, results_dir=args.results_dir,
            allow_dirty=ALLOW_DIRTY, extra_meta=meta,
        )
        log.info("wrote %s (%d rows)", path, len(frame))
    return 0


if __name__ == "__main__":
    sys.exit(main())
