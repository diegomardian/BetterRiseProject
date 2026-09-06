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
from src.reference.detection_scale import (
    UNIFORM_THINNING_R2_CEILING,
    delta_cloglog,
    uniform_thinning_null,
)
from src.reference.icbi_slice import (
    ADENOMA_TISSUE_MAP,
    COUNTS_LAYER,
    TISSUE_MAP,
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
    MIN_CELLS_BOTH_ARMS,
    N_BOOTSTRAP,
    RUNG,
    premise_holds,
    rows_for_patient,
    summarise,
)
from src.reference.labels import LabelError, label_column

log = logging.getLogger(__name__)

#: The unsorted fraction. GSE178341's UNSORTED, in this atlas's vocabulary.
NAIVE = "naive"

#: The QC batch unit. Same as GSE178341's.
BATCH_KEY = "sample_id"

#: Minimum epithelial cells per arm for a patient to be worth loading at all.
#: Matches the feasibility table that chose the 14 studies.
MIN_EPITHELIAL_PER_ARM = 100

DEFAULT_ATLAS = Path("/project/rise-batteries/bode/icbi/final_crc_atlas-adata.h5ad")

#: The two readings. `carcinoma` is what produced the committed 13-study result
#: and its mapping is unchanged; `adenoma` compares a polyp against the same
#: patient's own normal.
ARM_MAPS: dict[str, dict[str, str]] = {
    "carcinoma": TISSUE_MAP,
    "adenoma": ADENOMA_TISSUE_MAP,
}


def mature_bin(rung: str) -> str:
    """The bin this rung calls mature, read off RUNG_SPECS rather than fixed.

    `lineage` -> differentiated, `best4` -> best4. Hard-coding one rung's bin is
    how a run at another rung silently scores zero mature cells, which reads
    exactly like the finding this project keeps reporting.
    """
    from src.reference.labels import RUNG_SPECS

    return RUNG_SPECS[rung].mature


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


def eligible_patients(
    obs: pd.DataFrame, study_id: str, *, reading: str = "carcinoma"
) -> tuple[pd.DataFrame, list[str]]:
    """The naive, two-armed cells of one study, and the patients worth loading.

    ``reading`` picks the arm map. Under `adenoma` the diseased arm is the
    polyp and either normal label serves as the reference -- what makes it the
    patient's own is this function grouping by `patient_id`, not the label.
    """
    rows = obs[obs["study_id"] == study_id].copy()
    if rows.empty:
        raise SliceError(f"no cells for study {study_id!r}")

    before = len(rows)
    rows = rows[rows["enrichment_cell_types"].astype(str) == NAIVE]
    log.info("  %s fraction: %d of %d cells", NAIVE, len(rows), before)
    if rows.empty:
        return rows, []

    rows["tissue"] = arms(rows["sample_type"], ARM_MAPS[reading])
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


def arm_fractions(
    labels: pd.DataFrame, *, tissue, patient_id, rung: str
) -> dict | None:
    """The compositional arm's inputs for one patient, pivoted normal/tumour.

    Returns the four numbers `decompose_cohort` and the depth gate need, or
    ``None`` when either arm has no resolved epithelium -- in which case the
    fraction is not a fraction and the caller must count the patient rather
    than default it.

    EVERY NUMBER HERE COMES FROM `mature_cell_counts`, WHICH ALREADY DERIVES
    THEM. The denominator is `n_cells_resolved = n_cells_epithelial -
    n_cells_unresolved`, settled as open decision #14 on the ground that a cell
    which could not be labelled is not a cell measured to be immature. This
    function pivots, it does not decide.

    `unresolved_share_*` is carried because nothing has ever compared it
    BETWEEN ARMS. The `unresolved_depth` cut is applied per patient over both
    arms pooled and before depth matching, and it removes exactly
    `DEPTH_QUANTILE` of the epithelium by construction -- every patient in the
    MLH1 run reported 25.0%. If that quarter is drawn unevenly from the two
    arms, the two denominators differ for a purely technical reason and the
    compositional term is measuring sequencing depth. The gate lives in the
    reading job; the measurement it needs is emitted here.
    """
    from src.reference.labels import mature_cell_counts

    counts = mature_cell_counts(
        labels, patient_id=patient_id, tissue=tissue, axes=[AXIS], rungs=[rung],
    )
    by_arm = counts.set_index("tissue")
    if not {"normal", "tumour"}.issubset(by_arm.index):
        return None
    out: dict = {}
    for arm in ("normal", "tumour"):
        row = by_arm.loc[arm]
        resolved = int(row["n_cells_resolved"])
        if resolved <= 0:
            return None
        epithelial = int(row["n_cells_epithelial"])
        out[f"frac_mature_{arm}"] = float(row["mature_fraction"])
        out[f"n_cells_resolved_{arm}"] = resolved
        out[f"n_cells_epithelial_{arm}"] = epithelial
        out[f"unresolved_share_{arm}"] = float(row["unresolved_fraction"])
        # THE SENSITIVITY ARM. Open decision #14 put the unresolved cells
        # OUTSIDE the denominator -- a cell that could not be labelled is not a
        # cell measured to be immature. This is the same fraction computed as
        # if #14 had gone the other way, and it is free: it needs no extra
        # counting, only the denominator the decision rejected.
        #
        # It exists because the unresolved share is ENDOGENOUS. Measured on a
        # synthetic cohort whose two arms have identical composition by
        # construction: when the targets collapse in the tumour arm its cells
        # carry fewer counts, more of them fall below the depth target, and the
        # tumour arm's unresolved share rises to 0.270 against the normal arm's
        # 0.229. With the collapse removed both sit at ~0.25. So the exclusion
        # is driven by the effect under study, not independent of it.
        #
        # What that fixture ALSO showed is that it does not propagate: the
        # mature fraction came out at 0.500/0.502 against a true 0.500. One
        # fixture is not a proof, so the two denominators are both carried and
        # the comparison is made on the real data rather than assumed from this
        # one.
        out[f"frac_mature_{arm}_all_epithelial"] = (
            float(row["n_cells_mature"]) / epithelial if epithelial > 0 else float("nan")
        )
    out["unresolved_arm_gap"] = abs(
        out["unresolved_share_normal"] - out["unresolved_share_tumour"]
    )
    return out


def study_deltas(
    atlas: Path,
    obs: pd.DataFrame,
    study_id: str,
    *,
    seed: int = DEFAULT_SEED,
    max_patients: int | None = None,
    rung: str = RUNG,
    reading: str = "carcinoma",
    extra_genes: Sequence[str] = (),
    collect_fractions: bool = False,
) -> tuple[pd.DataFrame, dict]:
    """Per-patient, per-gene deltas for one study. Returns (deltas, report).

    ``extra_genes`` scores genes outside the standing six-gene panel on the
    SAME cells, QC, labels and depth matching. The MLH1 positive control is
    what it exists for: a gene whose silencing is known from an assay rather
    than from expression, scored by the instrument whose sensitivity is in
    question. They are additional to `GENE_ROLES`, never a replacement -- the
    controls have to be on the same rows or the premise cannot be read.

    An extra gene is added to ``target_genes`` as well as to the count matrix,
    so the leakage guard (invariant 2) is ASKED about it rather than merely not
    tripped by it. MLH1 is in no labelling axis today; a future axis that
    picked it up would be caught here instead of silently making a silenced
    cell unreadable as a silenced cell.
    """
    from src.common.panel import tier_genes
    from src.reference.labels import assign_labels
    from src.reference.qc import apply_qc, cell_qc_metrics, qc_thresholds

    extra_genes = tuple(dict.fromkeys(extra_genes))
    overlap = sorted(set(extra_genes) & set(GENE_ROLES))
    if overlap:
        raise SliceError(
            f"{overlap} are already on the standing panel. Scoring one twice "
            f"would put two rows for the same gene on one patient, and the "
            f"second would silently win every groupby downstream."
        )
    roles = {**GENE_ROLES, **{g: "positive_control" for g in extra_genes}}

    var = read_var(atlas)
    symbols = var["gene_symbol"]
    missing = sorted(set(roles) - set(symbols))
    if missing:
        # Skip-and-report: a study whose slice lacks a gene is a smaller
        # reading, not a failed run.
        log.warning("  genes absent from /var: %s", missing)

    rows, patients = eligible_patients(obs, study_id, reading=reading)
    if max_patients:
        patients = patients[:max_patients]
    report = {
        "study_id": study_id, "reading": reading, "granularity_rung": rung,
        "mature_bin": mature_bin(rung),
        "n_patients_eligible": len(patients),
        "genes_absent_from_var": missing, "batch_key": BATCH_KEY,
        "sorted_fraction_filter": f"enrichment_cell_types == {NAIVE!r}",
        "extra_genes": list(extra_genes),
        "gene_roles": roles,
    }
    absent_extras = sorted(set(extra_genes) - set(symbols))
    if absent_extras:
        # A CONTROL THAT IS NOT IN THE MATRIX IS NOT A CONTROL THAT DID NOT
        # MOVE. The standing panel tolerates a missing gene because the reading
        # is still a reading without it; an extra gene IS the reading, so its
        # absence has to stop the run rather than produce a table with one
        # fewer column and the same name.
        raise SliceError(
            f"{absent_extras} absent from /var, and they are what this run is "
            f"for. The atlas indexes /var/_index by Ensembl id and carries "
            f"symbols in a separate column -- if the symbol lookup is empty, "
            f"suspect the identifier space before concluding the gene is not "
            f"measured."
        )
    if not patients:
        log.warning("  no patient clears %d epithelial cells per arm",
                    MIN_EPITHELIAL_PER_ARM)
        return pd.DataFrame(), report

    log.info("  %d patient(s) clear %d epithelial cells per arm",
             len(patients), MIN_EPITHELIAL_PER_ARM)

    out: list[dict] = []
    unlabellable: list[dict] = []
    unfractionable: list[dict] = []
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
        try:
            labels = assign_labels(
                kept, symbols, compartment=comp,
                sample_id=block_obs[BATCH_KEY].to_numpy()[keep],
                target_genes=sorted(set(tier_genes("A")) | set(extra_genes)),
                tissue=tissue,
                patient_id=block_obs["patient_id"].to_numpy()[keep],
                depth_quantile=DEPTH_QUANTILE, seed=seed,
                index=pd.Index(block_obs.index.to_numpy()[keep]),
            )
        except LabelError as exc:
            # A PATIENT THIS COHORT CANNOT LABEL IS NOT A FAILED STUDY.
            # The commonest case is a tumour arm so much deeper than its own
            # normal arm that no reference cell clears the matching threshold
            # -- Pelka's C166 has a depth target of 27,146 against 949
            # epithelial cells. That is a property of one patient's library
            # preparation, and the two conditions already handled here by
            # skipping (too few QC survivors, too few mature cells) are the
            # same kind of thing.
            #
            # Counted and named, never silently dropped: a study that loses
            # patients this way is a smaller reading, and how much smaller has
            # to be visible in the sidecar.
            unlabellable.append({"patient_id": patient, "reason": str(exc)})
            log.warning("  [%d/%d] %s — cannot label: %s",
                        i, len(patients), patient, exc)
            continue
        call = labels[label_column(AXIS, rung)].astype(str).to_numpy()
        mature = (call == mature_bin(rung)) & np.isin(tissue, ["normal", "tumour"])
        if mature.sum() < MIN_CELLS_BOTH_ARMS:
            log.info("  [%d/%d] %s — %d mature cells, below the floor",
                     i, len(patients), patient, int(mature.sum()))
            continue

        mature_block = kept[mature]
        depth = np.asarray(mature_block.sum(axis=1), dtype=float).ravel()
        counts = {}
        for gene in roles:
            hit = np.where(symbols.to_numpy() == gene)[0]
            if len(hit) == 0:
                continue
            counts[gene] = np.asarray(
                mature_block[:, hit[0]].todense()
            ).ravel().astype(float)

        patient_rows = rows_for_patient(
            study_id=study_id, patient=patient, counts=counts,
            depth=depth, tissue=tissue[mature], seed=seed, roles=roles,
        )
        if collect_fractions and patient_rows:
            # THE COMPOSITIONAL ARM. `decompose_cohort` needs a mature FRACTION
            # per arm, and nothing in this path has ever emitted a denominator
            # -- the mask is computed and the population it came from is thrown
            # away. `mature_cell_counts` already derives all of it (open
            # decision #14: the denominator is the RESOLVED epithelium, because
            # a cell that could not be labelled is not a cell measured to be
            # immature), so this CALLS it rather than re-deriving anything.
            #
            # Attached to every gene row of this patient because that is the
            # shape decompose_cohort consumes: one row per (patient, gene,
            # rung, axis) carrying both fractions.
            fractions = arm_fractions(
                labels, tissue=tissue,
                patient_id=block_obs["patient_id"].to_numpy()[keep],
                rung=rung,
            )
            if fractions is None:
                # Both arms must have a resolved epithelium or the fraction is
                # not a fraction. Skipped and counted, never defaulted to zero.
                unfractionable.append({"patient_id": patient,
                                       "reason": "an arm has no resolved epithelium"})
                log.warning("  [%d/%d] %s — no resolved epithelium in one arm; "
                            "compositional term not formed", i, len(patients), patient)
            else:
                for row in patient_rows:
                    row.update(fractions)
        out.extend(patient_rows)
        log.info("  [%d/%d] %s — %d mature cells scored",
                 i, len(patients), patient, int(mature.sum()))

    frame = pd.DataFrame(out)
    report["n_patients_scored"] = int(frame["patient_id"].nunique()) if not frame.empty else 0
    report["n_patients_unlabellable"] = len(unlabellable)
    report["unlabellable"] = unlabellable
    report["collect_fractions"] = bool(collect_fractions)
    report["n_patients_unfractionable"] = len(unfractionable)
    report["unfractionable"] = unfractionable
    if unlabellable:
        log.warning("  %d of %d eligible patient(s) could not be labelled",
                    len(unlabellable), len(patients))
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

#: THE ADENOMA READ BAR, written before the numbers exist.
#:
#: There is no committed Chen_2021 result, so the Pelka-style validation --
#: reproduce a known answer -- does not apply. What replaces it is a statement
#: of how the output will be read, recorded in the run's own sidecar so it
#: cannot be composed after the fact. That is the same discipline the Pelka bar
#: followed, minus the baseline.
ADENOMA_READ_BAR = {
    "primary_study": "Chen_2021_Cell",
    "not_a_meta_analysis": (
        "Two studies carry a polyp arm and MIN_STUDIES is 3, so this is NOT "
        "meta-analysed. Chen_2021_Cell (44 paired patients) is the primary "
        "per-study reading -- a single 44-patient study IS the power argument. "
        "Zheng_2022 (3 patients) is a gradient companion, reported beside it "
        "and never pooled with it."
    ),
    "rungs": ["lineage", "best4"],
    "why_two_rungs": (
        "lineage for comparability with the 13 carcinoma studies; best4 because "
        "it is the resolution the question is actually posed at and the one "
        "carcinoma structurally could not support -- a median of 3 mature cells "
        "there. Adenomas retain differentiated epithelium, so best4 may be "
        "estimable for the first time. That measurement is the whole reason "
        "path C exists."
    ),
    "how_this_will_be_read": [
        "1. The premise verdict comes FIRST and gates everything. HOLDS licenses "
        "the detection reading; REFUSED or UNRESOLVED does not, at either rung.",
        "2. Only if the premise holds: the GUCA2A detection delta within the "
        "mature population. A fall there, in cells still called mature, is "
        "consistent with silencing.",
        "3. The best4 mature-cell counts are reported whatever the verdict, "
        "because 'best4 is estimable in adenoma' is itself a result and is the "
        "measurement carcinoma could not take.",
    ],
    "what_a_positive_reading_would_and_would_not_mean": (
        "It would be consistency with silencing, not proof of it. Survivorship "
        "-- GUCA2A-high cells preferentially destroyed -- is not "
        "transcript-detectable, and no amount of scale or resolution changes "
        "that. A premise that does not hold is an honest negative at a "
        "resolution nobody has measured, which is a result either way."
    ),
    "known_limits_before_the_run": [
        "TruDrop: median 1,881 epithelial genes/cell against Pelka's 2,630. "
        "Shallower cuts both ways -- less control saturation, less GUCA2A "
        "sensitivity.",
        "No within-patient polyp->carcinoma gradient exists in Chen_2021 "
        "(0 patients hold both), so timing evidence is polyp-versus-normal only.",
        "The reference arm is labelled `healthy normal` and is the patient's "
        "own mucosa; pairing is established by patient_id, not by the label.",
    ],
}


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


#: WHICH STATISTIC THE READING RESTS ON, fixed before the re-run.
#:
#: The three below do not agree, so leaving the choice open until the numbers
#: are in front of you is choosing on the answer. Recorded here and copied into
#: every sidecar this job writes.
#:
#: `cloglog` is LOAD-BEARING. It is a monotone transform of the detection rate
#: the project already chose, so it changes the scale without changing the
#: measurement, and on that scale a difference between two genes is a ratio of
#: fold changes rather than a comparison of two proportions with different room
#: to move.
#:
#: `log2_cp10k` CORROBORATES. It reaches the same quantity by a different route
#: -- a mean, not a rate -- so agreement between the two is evidence the
#: conclusion does not rest on the Poisson model being right. Where they
#: disagree, say so; do not pick.
#:
#: `detection` is DIAGNOSTIC ONLY. Kept because it is what every previously
#: committed table reports and dropping it would make this run incomparable to
#: them, and because the gap between it and the other two IS the finding here.
#: It is not read across genes.
#:
#: HONESTY ABOUT WHAT "PRE-COMMITTED" MEANS HERE. This rule is pre-committed
#: relative to any future re-run and to the MLH1 reading. It is NOT pre-committed
#: relative to the adenoma numbers: all three statistics were computed
#: exploratorily on 2026-09-05, before this constant was written, and they were
#: seen to agree on the load-bearing contrast (GUCA2A - MS4A12 contains zero on
#: every one of them, at both rungs). Claiming a blind pre-commitment that did
#: not happen would be the same defect this repository is about, so the sequence
#: is recorded instead of dressed up.
SPECIFICITY_STATISTICS: dict[str, dict[str, str]] = {
    "cloglog": {
        "column": "delta_cloglog",
        "standing": "load-bearing",
        "scale": "log fold change (natural log), from detection under a Poisson model",
        "why": (
            "comparable across genes of any baseline abundance, which the raw "
            "detection delta is not"
        ),
    },
    "log2_cp10k": {
        "column": "log2_cp10k_ratio",
        "standing": "corroborating",
        "scale": "log2 fold change of the per-cell CP10K mean",
        "why": (
            "independent of the Poisson model, so agreement with cloglog is "
            "evidence the reading does not rest on that model"
        ),
    },
    "detection": {
        "column": "delta_detect",
        "standing": "diagnostic only -- NOT read across genes",
        "scale": "difference of two proportions",
        "why": (
            "the statistic every earlier table reports, kept for comparability. "
            "Its sensitivity varies with baseline abundance, so ranking genes "
            "by it ranks them substantially by how detectable they were."
        ),
    },
}

LOAD_BEARING_STATISTIC = "cloglog"

#: Every gene against every other, not just the target against each. The tier
#: claim -- terminal differentiation down, intestinal identity retained --
#: is a statement about where CDX2 sits relative to the CONTROLS, and a table
#: of `GUCA2A - X` rows cannot express it. Reading it off the fact that CDX2's
#: delta looked small next to GUCA2A's is the cross-gene magnitude comparison
#: this module now exists to stop.
SPECIFICITY_PANEL: tuple[str, ...] = ("KRT8", "ACTB", "EPCAM", "CDX2",
                                      "MS4A12", "GUCA2A")


def specificity(deltas: pd.DataFrame, *, seed: int = DEFAULT_SEED) -> pd.DataFrame:
    """Every pairwise contrast in the panel, on all three statistics.

    The premise check compares housekeeping between arms, which asks whether
    these are still cells. It cannot ask whether they are still equally MATURE
    cells, because ACTB and KRT8 do not vary with maturity -- a control that
    cannot move with the thing in question cannot certify it.

    So this pairs genes WITHIN each patient. `MS4A12` is a mature-colonocyte
    marker with no silencing story attached to it: if GUCA2A were specifically
    silenced it should fall further than MS4A12; if the mature LABEL is
    admitting less-mature cells in the diseased arm, both fall together and the
    difference is zero. That is the week-0 falsification rule in a second place.

    TWO THINGS CHANGED AFTER THE FIRST ADENOMA RUN.

    *Every pair, not just the target's.* `CDX2 - KRT8` is what decides whether
    intestinal identity is retained, and it was never computed; the claim rested
    on CDX2's delta looking small beside GUCA2A's. On the committed detection
    statistic that contrast excludes zero -- it says the opposite -- and on both
    abundance-free statistics it contains zero. A conclusion that flips with the
    statistic needs the statistic named, which is what `SPECIFICITY_STATISTICS`
    does.

    *Three statistics per pair.* See `SPECIFICITY_STATISTICS`. `cloglog` is what
    the reading rests on; the others are reported beside it so a disagreement is
    visible in the table rather than available to whoever quotes it.
    """
    rng = np.random.default_rng(seed)
    rows = []
    for rung, block in deltas.groupby("granularity_rung", observed=True):
        present = [g for g in SPECIFICITY_PANEL if g in set(block["gene"])]
        for name, spec in SPECIFICITY_STATISTICS.items():
            column = spec["column"]
            if column not in block.columns:
                continue
            wide = block.pivot_table(
                index="patient_id", columns="gene", values=column
            )
            for i, gene in enumerate(present):
                for other in present[:i] + present[i + 1:]:
                    if gene not in wide.columns or other not in wide.columns:
                        continue
                    paired = (wide[gene] - wide[other]).dropna().to_numpy(dtype=float)
                    if paired.size < 3:
                        continue
                    draws = rng.choice(
                        paired, size=(N_BOOTSTRAP, paired.size), replace=True
                    ).mean(axis=1)
                    lo, hi = (float(x) for x in np.percentile(draws, [2.5, 97.5]))
                    rows.append({
                        "granularity_rung": rung,
                        "statistic": name,
                        "standing": spec["standing"],
                        "contrast": f"{gene} - {other}",
                        "gene": gene,
                        "other": other,
                        "role_of_gene": GENE_ROLES.get(gene, "?"),
                        "role_of_other": GENE_ROLES.get(other, "?"),
                        "n_patients": int(paired.size),
                        "mean_difference": float(paired.mean()),
                        "ci_low": lo, "ci_high": hi,
                        "excludes_zero": bool(lo > 0 or hi < 0),
                        "load_bearing": name == LOAD_BEARING_STATISTIC,
                    })
    return pd.DataFrame(rows)


def contrast_matrix(
    contrasts: pd.DataFrame, rung: str, statistic: str = LOAD_BEARING_STATISTIC
) -> pd.DataFrame:
    """The pairwise table as a matrix: rows minus columns, ``*`` excludes zero.

    A block of genes that are mutually indistinguishable and jointly separated
    from another block is a TIER, and that shape is legible in a matrix and
    close to invisible in a list of thirty rows. The adenoma reading's whole
    conclusion is which genes fall into which block.
    """
    block = contrasts[
        (contrasts["granularity_rung"] == rung)
        & (contrasts["statistic"] == statistic)
    ]
    genes = [g for g in SPECIFICITY_PANEL if g in set(block["gene"])]
    out = pd.DataFrame(".", index=genes, columns=genes, dtype=object)
    for row in block.itertuples():
        if row.gene in out.index and row.other in out.columns:
            mark = "*" if row.excludes_zero else " "
            out.loc[row.gene, row.other] = f"{row.mean_difference:+.3f}{mark}"
    return out


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
    parser.add_argument("--arms", choices=tuple(ARM_MAPS), default="carcinoma",
                        help="`adenoma` scores polyp against the patient's own "
                             "normal. Default is the committed carcinoma contrast.")
    parser.add_argument(
        "--collect-fractions", action="store_true",
        help="also emit the mature FRACTION per arm, which the decomposition "
             "needs and the coexpression reading does not. Off by default so "
             "the committed coexpression tables keep their exact shape.",
    )
    parser.add_argument("--rungs", nargs="+", default=[RUNG],
                        help="granularity rungs to score. `--rungs lineage best4` "
                             "for the adenoma reading.")
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

    if args.arms == "adenoma":
        log.info("\n%s\nADENOMA READ BAR, recorded before the numbers exist\n%s",
                 "=" * 72, "=" * 72)
        for line in ADENOMA_READ_BAR["how_this_will_be_read"]:
            log.info("  %s", line)
        log.info("  %s", ADENOMA_READ_BAR["not_a_meta_analysis"])

    frames, reports = [], []
    for study in studies:
        for rung in args.rungs:
            label = f"{study} / {rung}" if len(args.rungs) > 1 else study
            log.info("\n=== %s ===", label)
            try:
                deltas, report = study_deltas(
                    args.atlas, obs, study, seed=args.seed,
                    max_patients=args.max_patients, rung=rung, reading=args.arms,
                collect_fractions=args.collect_fractions,
                )
            except SliceError as exc:
                log.error("  REFUSED: %s", exc)
                reports.append({"study_id": study, "granularity_rung": rung,
                                "error": str(exc)})
                continue
            if not deltas.empty:
                deltas = deltas.assign(granularity_rung=rung, reading=args.arms)
                holds, premise = premise_holds(deltas, seed=args.seed)
                report["premise_holds"] = bool(holds)
                report["premise_reading"] = premise
                log.info("  premise: %s", premise)
                frames.append(deltas)
            reports.append(report)

    if not frames:
        log.error("\nno study produced any scored patient.")
        return 3

    deltas = pd.concat(frames, ignore_index=True)
    # The log fold-change scale, alongside the detection rate rather than
    # instead of it. Detection deltas are not comparable BETWEEN genes -- see
    # src/reference/detection_scale.py -- and every cross-gene contrast this
    # job reports is computed on this column as well.
    deltas["delta_cloglog"] = delta_cloglog(deltas)

    # SUMMARISE PER RUNG. `summarise` groups by (study, gene, statistic) and
    # knows nothing about rungs, so concatenating two rungs and calling it once
    # averages two different populations into one row -- the first adenoma run
    # reported n_patients = 64, which is 44 at lineage plus 20 at best4, and
    # those are not the same 64 patients or the same cells.
    summary = pd.concat(
        [summarise(block, seed=args.seed).assign(granularity_rung=rung)
         for rung, block in deltas.groupby("granularity_rung")],
        ignore_index=True,
    )

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

    stem = "icbi_coexpression" if args.arms == "carcinoma" else f"icbi_{args.arms}"
    contrasts = specificity(deltas, seed=args.seed)
    if not contrasts.empty:
        log.info("\n%s\nIS THE TARGET'S FALL SPECIFIC TO IT?\n%s", "=" * 72, "=" * 72)
        log.info(
            "Reported on %d statistics; %r is load-bearing and the others are "
            "beside it\nso a disagreement is visible in the table. Matrices "
            "below are the load-bearing one.",
            len(SPECIFICITY_STATISTICS), LOAD_BEARING_STATISTIC,
        )
        for rung in sorted(contrasts["granularity_rung"].unique()):
            log.info("\n-- %s, %s --", rung, LOAD_BEARING_STATISTIC)
            log.info("%s", contrast_matrix(contrasts, rung).to_string())
        log.info(
            "\nRows minus columns; a trailing * excludes zero. A gene against a "
            "CONTROL\nexcluding zero says it falls more than housekeeping; "
            "against an IDENTITY marker\ncontaining zero says it is "
            "indistinguishable from the mature programme, and no\n"
            "gene-specific claim follows. Read the ROLE columns, not the gene "
            "names."
        )

    # THE NULL THE GRADIENT HAS TO BEAT. Fitted per rung, on the detection
    # deltas, because it is the detection deltas whose ordering gets quoted.
    null_frames, null_verdicts = [], {}
    for rung, block in deltas.groupby("granularity_rung", observed=True):
        table, verdict = uniform_thinning_null(block)
        null_verdicts[str(rung)] = verdict
        if not table.empty:
            null_frames.append(table.assign(granularity_rung=rung))
        log.info("\n%s\nCOULD ONE UNIFORM THINNING PRODUCE THIS GRADIENT? "
                 "[%s]\n%s", "=" * 72, rung, "=" * 72)
        log.info("  %s: %s", verdict["verdict"], verdict["detail"])
        if not table.empty:
            log.info("%s", table.to_string(index=False))
    thinning = (pd.concat(null_frames, ignore_index=True)
                if null_frames else pd.DataFrame())

    for frame, name in ((deltas, stem), (summary, f"{stem}_summary"),
                        (contrasts, f"{stem}_specificity"),
                        (thinning, f"{stem}_thinning_null")):
        if frame.empty:
            continue
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
                "reading": args.arms,
                "rungs": list(args.rungs),
                "arm_map": ARM_MAPS[args.arms],
                **({"adenoma_read_bar": ADENOMA_READ_BAR} if args.arms == "adenoma" else {}),
                "premise_reading": {
                    r["study_id"]: r.get("premise_reading", "not scored")
                    for r in reports if "study_id" in r
                },
                "specificity_statistics": SPECIFICITY_STATISTICS,
                "load_bearing_statistic": LOAD_BEARING_STATISTIC,
                "specificity_panel": list(SPECIFICITY_PANEL),
                "cross_gene_rule": (
                    "Detection deltas are NOT comparable between genes: the "
                    "sensitivity of a proportion depends on its baseline, and "
                    "this panel spans 0.36 to 0.98. Every cross-gene statement "
                    "is made on the load-bearing log fold-change scale, with "
                    "log2 CP10K corroborating. delta_detect is reported for "
                    "comparability with earlier tables and is not read across "
                    "genes."
                ),
                "cloglog_boundary_rule": (
                    "a rate is turned back into cell counts and given the "
                    "Jeffreys correction (k + 1/2)/(n + 1) before the "
                    "transform, so the correction scales with how many cells "
                    "the estimate rests on rather than being a fixed epsilon"
                ),
                "uniform_thinning_null": null_verdicts,
                "uniform_thinning_r2_ceiling": UNIFORM_THINNING_R2_CEILING,
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
