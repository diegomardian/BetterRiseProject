"""Can the adenoma decomposition be replicated on Becker? The gate, before the work.

    python -m src.reference.jobs.becker_feasibility --object <path> --inspect
    python -m src.reference.jobs.becker_feasibility --object <path>

Pre-registered in ``docs/prereg_becker_replication.md`` §3. **Read it before
reading any number this produces** — the per-gene gate and its four
consequences are fixed there, including the one where the replication cannot be
run at all.

WHY A GATE AND NOT JUST THE ANALYSIS. B1's risk is not statistical, it is
physical: GUCA2A and MS4A12 are cytoplasmic transcripts and snRNA-seq samples
nuclei. Their baseline detection in Chen_2021's normal-arm mature cells is 0.437
and **0.363**, the two lowest on the panel. A protocol that halves cytoplasmic
detection puts MS4A12 near 0.18 and quartering puts it near 0.09. So the whole
avenue can fail for a reason that has nothing to do with the biology, and it can
fail cheaply — this job is hours, the replication is not.

**Run this before deleting the ICBI atlas to make room.** If the gate fails,
32 GB was freed and a 25-minute re-fetch incurred for nothing.

TWO THINGS THIS JOB REFUSES TO ASSUME, both of which have cost this repository
real time before.

*The arm vocabulary.* ``ADENOMA_TISSUE_MAP`` was written against the ICBI
atlas's words, and reading a label rather than the patient grouping once put
Chen_2021's usable pairs at **zero** when the true number was 44 — its reference
samples say ``healthy normal`` and that was read as "a different donor". Becker
is a different deposit with its own words. **This job reports the observed
``sample_type`` vocabulary and does not map it**; the mapping is a decision for
a human who has seen the list.

*The identifier space.* S matrices are Ensembl, Lee's GEO matrices are symbols,
and the ICBI atlas is Ensembl in ``/var/_index`` with symbols in a separate
column. Four times the symptom was an **empty intersection reported as a
finding**. So this job reports how many panel genes match under each naming and
**refuses to report a detection of zero without saying which space it looked
in**.

``--inspect`` does both of those and nothing else. **Run it first.**

THE FORMAT QUESTION IS SETTLED, and not the way §6 feared. Verified 2026-09-06
against the GEO listing: ``GSE201348_RAW.tar``, 1.2 GB, **72 standard 10x
triplets** — not Seurat objects. What that costs instead is that **the tar
carries no metadata at all**. A filename gives a GSM, a donor and a sample and
says nothing about whether it is a polyp or unaffected mucosa, so
``GSE201348_series_matrix.txt.gz`` is a required second input and ``--tar``
refuses to run without it. The mapping it reads
(``Polyp``/``Unaffected``/``CRC``) is fixed in Amendment 1 of the prereg.
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
from src.common.provenance import DEFAULT_SEED
from src.reference.jobs.coexpression_silencing import DETECTION_MIN_UMI, GENE_ROLES

log = logging.getLogger(__name__)

#: Pre-registered gate, `docs/prereg_becker_replication.md` §3. A gene enters
#: the replication only if BOTH hold in the reference arm.
MIN_DETECTION = 0.10
MIN_PATIENT_SHARE_NONZERO = 0.75

#: What Chen_2021 measured, for the side-by-side. Not a threshold — the gate is
#: absolute — but a reader needs to see how far the nuclear protocol moved each
#: gene, and a gene that drops by 4x while clearing the floor is still news.
CHEN_BASELINE: dict[str, float] = {
    "ACTB": 0.984, "KRT8": 0.958, "EPCAM": 0.900,
    "CDX2": 0.822, "GUCA2A": 0.437, "MS4A12": 0.363,
}

#: Failing this gene ends the replication: the primary claim is GUCA2A's four
#: cross-block contrasts (prereg §1).
CRITICAL_GENE = "GUCA2A"

#: Mature-colonocyte markers for the label §3's gate is actually about.
#:
#: NONE OF THESE IS ON THE PANEL, and that is required rather than tidy.
#: Invariant 2 forbids a target gene appearing in a label; labelling mature
#: cells by GUCA2A and then measuring GUCA2A in them would return the threshold
#: it was given. The whole panel is excluded, not just the target, so the
#: controls stay usable as controls inside the label.
#:
#: Canonical differentiated / absorptive colonocyte markers. Committed here
#: BEFORE the labelled gate is run, so the set cannot be tuned to the answer.
MATURE_MARKERS: tuple[str, ...] = (
    "CA1", "CA2", "AQP8", "SLC26A3", "KRT20", "CEACAM7",
)

#: A nucleus is called mature when it detects at least this many of them. Two,
#: not one: at these detection rates a single marker is mostly noise, and
#: requiring all six would select depth rather than identity.
MIN_MARKERS_FOR_MATURE = 2

#: Below this many labelled nuclei the gate is not taken. A detection rate over
#: a handful of cells is not a detection rate, and the failure mode is reporting
#: one anyway.
MIN_MATURE_CELLS = 200

#: THE DEPTH TRAP, named because it is the way this analysis fails.
#: "Detects >= 2 markers" is correlated with sequencing depth, and a deeper
#: nucleus detects EVERYTHING more -- including GUCA2A. So a mature-cell
#: restriction can lift the target over the floor while measuring nothing but
#: library size. The controls are scored in the SAME cells for exactly this
#: reason: if GUCA2A's enrichment sits inside the band the controls describe,
#: the enrichment is depth and the pass is hollow. `enrichment_audit` computes
#: it and the job prints the verdict either way.
DEPTH_AUDIT_CONTROLS: tuple[str, ...] = ("ACTB", "KRT8", "EPCAM")


class FeasibilityError(ValueError):
    """The object cannot be read as a cohort this gate applies to."""


def inspect_deposit(tar: Path, series_matrix: Path) -> dict:
    """Report what is in the GSE201348 deposit. Assume nothing, map nothing.

    THE FIRST THING TO RUN. It needs BOTH files because the tar carries no
    metadata at all: a filename gives a GSM, a donor and a sample, and says
    nothing about whether the sample is a polyp or unaffected mucosa. That lives
    only in the series matrix.

    Reports the arm counts the mapping produces rather than applying it
    silently, so a human sees the cohort before any analysis does.
    """
    from src.reference.becker_io import (
        DISEASE_STAGE_MAP,
        PAIRED_ARMS,
        gene_symbols,
        paired_donors,
        read_series_matrix,
        read_triplet,
        sample_files,
    )

    metadata = read_series_matrix(series_matrix)
    files = sample_files(tar)
    merged = files.merge(metadata, on="gsm", how="outer", suffixes=("", "_meta"))

    counts, barcodes, features = read_triplet(tar, files.iloc[0])
    symbols = gene_symbols(features)
    panel = set(GENE_ROLES)

    # PAIRED, not merely "has an arm". `healthy_donor` is an arm and is not a
    # reference for a paired design -- those donors have no polyps at all.
    # Counting arm.nunique()==2 would call a tumour+healthy_donor pair "paired",
    # which is the cross-donor comparison Becker Amendment 2 refuses.
    scored = metadata[metadata["arm"].isin(PAIRED_ARMS)]
    paired = paired_donors(metadata)
    report = {
        "tar": str(tar), "series_matrix": str(series_matrix),
        "n_samples_in_tar": int(len(files)),
        "n_samples_in_metadata": int(len(metadata)),
        "incomplete_triplets": files.loc[~files["complete"], "gsm"].tolist(),
        "in_tar_not_metadata": merged.loc[merged["title"].isna(), "gsm"].dropna().tolist(),
        "in_metadata_not_tar": merged.loc[merged["matrix"].isna(), "gsm"].dropna().tolist(),
        "disease_stage_counts": metadata["disease_stage"].value_counts().to_dict(),
        "arm_counts": metadata["arm"].value_counts(dropna=False).to_dict(),
        "arm_map_used": {k: v for k, v in DISEASE_STAGE_MAP.items()},
        "n_donors": int(metadata["donor"].nunique()),
        "n_donors_PAIRED": int(len(paired)),
        "paired_donors": sorted(paired),
        "paired_arms": sorted(PAIRED_ARMS),
        "donors_with_polyps_but_no_same_donor_reference": sorted(
            set(metadata.loc[metadata["arm"] == "tumour", "donor"]) - set(paired)),
        "healthy_donor_only": sorted(
            set(metadata.loc[metadata["arm"] == "healthy_donor", "donor"])
            - set(paired)),
        "samples_per_donor": scored.groupby(["donor", "arm"]).size()
                                   .unstack(fill_value=0).to_dict("index"),
        "replicate_samples": metadata.loc[metadata["replicate"].notna(),
                                          "sample_id"].unique().tolist(),
        "first_sample_shape_cells_by_genes": list(counts.shape),
        "n_genes": int(len(symbols)),
        "panel_genes_found": sorted(panel & set(symbols)),
        "panel_genes_missing": sorted(panel - set(symbols)),
        "features_column_used": "column 1 (symbol); column 0 is Ensembl",
    }
    if "fap" in metadata.columns:
        report["fap_donors"] = (metadata.drop_duplicates("donor")["fap"]
                                .value_counts().to_dict())
    return report


def inspect(path: Path) -> dict:
    """Report what is in an h5ad. Kept for a deposit that arrives as one.

    The first thing to run, and on a deposit whose format is unverified it may
    be the only thing that runs. Reports the obs vocabulary, the gene naming
    space, and whether the matrix looks like raw counts — the three things that
    have to be true before any gate is meaningful.
    """
    import anndata

    adata = anndata.read_h5ad(str(path), backed="r")
    obs = adata.obs
    var_names = [str(v) for v in adata.var_names]
    panel = set(GENE_ROLES)

    by_index = panel & set(var_names)
    by_column: dict[str, int] = {}
    for column in adata.var.columns:
        values = {str(v) for v in adata.var[column].astype(str)}
        hit = len(panel & values)
        if hit:
            by_column[column] = hit

    report = {
        "path": str(path),
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "obs_columns": list(obs.columns),
        "layers": list(adata.layers.keys()),
        "var_columns": list(adata.var.columns),
        "panel_genes_in_var_names": sorted(by_index),
        "panel_genes_by_var_column": by_column,
        "var_names_look_like": (
            "ensembl" if sum(v.startswith("ENSG") for v in var_names[:200]) > 100
            else "symbols_or_other"
        ),
    }
    for candidate in ("sample_type", "tissue", "condition", "disease",
                      "sample", "Sample", "polyp_type", "patient", "patient_id",
                      "donor", "Donor"):
        if candidate in obs.columns:
            counts = obs[candidate].astype(str).value_counts()
            report[f"vocabulary::{candidate}"] = counts.head(25).to_dict()
    return report


def gate(detection: pd.DataFrame) -> pd.DataFrame:
    """Apply §3's per-gene gate. One row per gene, verdict on each.

    ``detection`` needs ``gene``, ``detection``, ``share_patients_nonzero``.
    Per gene and not globally, because the claims are not equally exposed:
    MS4A12 failing costs the secondary claim and GUCA2A failing ends the
    replication.
    """
    needed = {"gene", "detection", "share_patients_nonzero"}
    missing = sorted(needed - set(detection.columns))
    if missing:
        raise FeasibilityError(f"the detection table is missing {missing}")

    out = detection.copy()
    out["chen_baseline"] = out["gene"].map(CHEN_BASELINE)
    # THE SCALE. `detection / chen_baseline` is a ratio of PROBABILITIES and it
    # is not comparable across genes, because detection is bounded at 1 and
    # ACTB/KRT8/EPCAM enter saturated (0.984, 0.958, 0.900). Dividing a small
    # number by a saturated one compresses the loss and makes a housekeeping
    # gene look better retained than a target. It reversed the two here: on the
    # ratio scale ACTB reads as retained 2.6x better than GUCA2A, and on the
    # detection scale GUCA2A is retained 1.1x BETTER than ACTB.
    #
    # This repository already found and fixed this once -- `docs/HANDOFF.md`
    # §6d, the cross-gene scale correction -- and the fix lived in another
    # module, so this job reintroduced it. cloglog(p) = log(mu), so a difference
    # on that scale is a log fold change in expected UMIs and IS comparable.
    out["mu_becker"] = -np.log1p(-out["detection"].clip(upper=1 - 1e-12))
    out["mu_chen"] = -np.log1p(-out["chen_baseline"].clip(upper=1 - 1e-12))
    out["log_fc_vs_chen"] = np.log(out["mu_becker"] / out["mu_chen"])
    out["fold_mu_vs_chen"] = out["mu_becker"] / out["mu_chen"]
    # Kept under a name that says what it is, so the old number stays checkable.
    out["naive_ratio_of_p"] = out["detection"] / out["chen_baseline"]
    out["clears_detection"] = out["detection"] >= MIN_DETECTION
    out["clears_patient_share"] = (
        out["share_patients_nonzero"] >= MIN_PATIENT_SHARE_NONZERO)
    out["passes"] = out["clears_detection"] & out["clears_patient_share"]
    out["role"] = out["gene"].map(GENE_ROLES)
    return out.sort_values("detection", ascending=False, ignore_index=True)


def verdict(gated: pd.DataFrame, *, mature_labelled: bool = False,
            audit: "pd.DataFrame | None" = None) -> dict:
    """§3's outcome table, taken by code rather than by a reader.

    ``mature_labelled`` says whether the detection handed in is the one §3
    actually names -- mature cells of the reference arm. When it is not, a FAIL
    is a failure of a LOWER BOUND and does not license §3's "cannot run": the
    quantity that failed is smaller than the quantity the gate is about, by an
    unmeasured amount. Emitting a pre-registered verdict from an analysis that
    is not the pre-registered one is the category error this argument exists to
    prevent, and the first version of this job made it.

    A PASS on a lower bound is safe in the other direction and is reported as a
    pass, which is the asymmetry §3 already relies on.
    """
    passing = set(gated.loc[gated["passes"], "gene"])
    failing = sorted(set(gated["gene"]) - passing)

    # A PASS BOUGHT BY DEPTH IS NOT A PASS. The mature label is built from
    # marker DETECTION, which rises with library size, so the restriction lifts
    # every gene. `enrichment_audit` asks whether the critical gene rose FURTHER
    # than the controls did in the same cells. If it did not, the floor was
    # cleared by sequencing depth and §3's design is not licensed.
    #
    # This used to log a warning beside a FULL DESIGN verdict. A check that
    # warns and does not bind is a check that reports success, which is the
    # defect this repository is named after. It binds now.
    if (mature_labelled and audit is not None and not audit.empty
            and CRITICAL_GENE in passing):
        row = audit.loc[audit["gene"] == CRITICAL_GENE]
        if len(row) and not bool(row["beyond_control_band"].iloc[0]):
            r = row.iloc[0]
            # The audit now carries every arm, so restrict the contrast to the
            # arm the gate was taken on. Without this it prints one name per
            # arm and reads as three genes.
            same_arm = (audit["arm"] == r["arm"] if "arm" in audit.columns
                        else pd.Series(True, index=audit.index))
            others = sorted(set(audit.loc[
                audit["beyond_control_band"] & same_arm
                & (audit["gene"] != CRITICAL_GENE), "gene"]))
            return {
                "verdict": "CLEARED BY DEPTH — NOT LICENSED",
                "detail": (
                    f"{CRITICAL_GENE} clears the floor "
                    f"({r['detection_whole_arm']:.3f} -> "
                    f"{r['detection_mature']:.3f}) but its enrichment "
                    f"({r['log_enrichment']:+.3f}) sits INSIDE the band the "
                    f"controls describe ({r['control_band_low']:+.3f} to "
                    f"{r['control_band_high']:+.3f}). The label selected deeper "
                    f"nuclei, not mature ones, so the gate was cleared by "
                    f"library size. The replication is NOT licensed on this "
                    f"reading."
                    + (f" Genuinely enriched, for contrast: {', '.join(others)}."
                       if others else "")
                ),
            }

    if failing and not mature_labelled:
        rows = gated.loc[gated["gene"].isin(failing)]
        shortfall = ", ".join(
            f"{r.gene} {r.detection:.3f} (needs {MIN_DETECTION / r.detection:.2f}x)"
            for r in rows.itertuples())
        return {
            "verdict": "GATE NOT RUN — LOWER BOUND FAILS",
            "detail": (
                f"{shortfall}. This is NOT §3's gate: §3 names the mature cells "
                f"of the reference arm and this is every cell in it. Mature "
                f"cells enrich for these markers, so the pre-registered "
                f"quantity is LARGER than what was measured, by an amount "
                f"nobody has measured. B1 is UNDETERMINED. To close it, label "
                f"the cells and re-run; to abandon it, say why the enrichment "
                f"cannot cover the shortfall above."
            ),
        }

    if CRITICAL_GENE not in set(gated["gene"]):
        return {
            "verdict": "CANNOT RUN",
            "detail": (
                f"{CRITICAL_GENE} is not in the object at all. Before "
                f"concluding it is not measured, check the identifier space — "
                f"this repository has reported an empty intersection as a "
                f"finding four times."
            ),
        }
    if CRITICAL_GENE not in passing:
        row = gated.loc[gated["gene"] == CRITICAL_GENE].iloc[0]
        return {
            "verdict": "CANNOT RUN",
            "detail": (
                f"{CRITICAL_GENE} fails the gate: detection "
                f"{row['detection']:.3f} against a floor of {MIN_DETECTION}, "
                f"non-zero in {row['share_patients_nonzero']:.0%} of patients "
                f"against {MIN_PATIENT_SHARE_NONZERO:.0%}. The primary claim is "
                f"its four cross-block contrasts, so the replication cannot be "
                f"run. **This is a measurement about snRNA-seq, not a negative "
                f"result about the biology** — report the detection table and "
                f"stop."
            ),
        }
    if "MS4A12" in failing:
        return {
            "verdict": "PRIMARY ONLY",
            "detail": (
                f"{CRITICAL_GENE} passes and MS4A12 does not, so the primary "
                f"claim is testable and the secondary is not. "
                f"GUCA2A − MS4A12 cannot be formed, so the 'not gene-specific' "
                f"null is NOT replicated and stays single-cohort. Genes failing: "
                f"{failing}."
            ),
        }
    if failing:
        return {
            "verdict": "REDUCED PANEL",
            "detail": (
                f"{failing} fail the gate. The comparator set is smaller and "
                f"cross-block counts are out of fewer than 8 — say so wherever "
                f"they are quoted."
            ),
        }
    return {
        "verdict": "FULL DESIGN",
        "detail": "every panel gene clears the gate; the replication runs as "
                  "pre-registered.",
    }


def label_mature(counts, marker_index: dict[str, int], *,
                 min_umi: int = DETECTION_MIN_UMI,
                 min_markers: int = MIN_MARKERS_FOR_MATURE):
    """Boolean mask: nuclei calling as mature colonocytes.

    ``marker_index`` maps :data:`MATURE_MARKERS` to columns and is built by the
    caller, so a label can never be produced without the caller having said
    which identifier space it looked in — the same rule ``detection_table``
    follows for the panel.

    Raises rather than returning an all-False mask when no marker was located.
    An empty label would flow downstream as "no mature cells", which reads as a
    biological statement and is a lookup failure.
    """
    if not marker_index:
        raise FeasibilityError(
            f"none of {list(MATURE_MARKERS)} was located, so no mature label "
            f"can be built. Check the identifier space before concluding the "
            f"markers are absent — that error has been made four times here. "
            f"An all-False mask would read as 'this tissue has no mature "
            f"colonocytes', which is a claim, not a missing lookup."
        )
    hits = None
    for column in marker_index.values():
        values = counts[:, column]
        values = (np.asarray(values.todense()).ravel()
                  if hasattr(values, "todense") else np.asarray(values))
        detected = (values.astype(float) >= min_umi).astype(int)
        hits = detected if hits is None else hits + detected
    return hits >= min_markers


def enrichment_audit(whole: pd.DataFrame, mature: pd.DataFrame) -> pd.DataFrame:
    """Did the mature label enrich the target, or just select deeper nuclei?

    Per gene, the log fold change of detection from the whole arm to the mature
    label, on the detection scale. The controls define a band; a target whose
    enrichment sits INSIDE that band has not been enriched, it has been
    resampled at greater depth, and a gate it clears on that basis is hollow.

    This is the same shape as ``coexpression_silencing``'s rule that every
    control is scored in the same cells as the target, and it exists because
    the label is built from marker DETECTION, which is a depth-correlated
    quantity by construction.
    """
    mu = lambda p: -np.log1p(-np.clip(np.asarray(p, dtype=float), 0, 1 - 1e-12))
    left = whole.set_index("gene")["detection"]
    right = mature.set_index("gene")["detection"]
    genes = [g for g in left.index if g in right.index]
    out = pd.DataFrame({
        "gene": genes,
        "detection_whole_arm": left.loc[genes].to_numpy(),
        "detection_mature": right.loc[genes].to_numpy(),
    })
    out["log_enrichment"] = np.log(mu(out["detection_mature"])
                                   / mu(out["detection_whole_arm"]))
    out["role"] = out["gene"].map(GENE_ROLES)
    band = out.loc[out["gene"].isin(DEPTH_AUDIT_CONTROLS), "log_enrichment"]
    out["control_band_low"] = float(band.min()) if len(band) else float("nan")
    out["control_band_high"] = float(band.max()) if len(band) else float("nan")
    out["beyond_control_band"] = out["log_enrichment"] > out["control_band_high"]
    return out.sort_values("log_enrichment", ascending=False, ignore_index=True)


def detection_table(
    counts, gene_index: dict[str, int], patient_id, *, min_umi: int = DETECTION_MIN_UMI
) -> pd.DataFrame:
    """Detection and per-patient non-zero share, per panel gene.

    ``gene_index`` maps a panel symbol to its column. It is built by the caller
    from whichever identifier space matched, and passed in rather than resolved
    here, so a zero can never be reported without the caller having said which
    space it looked in.
    """
    patient_id = np.asarray([str(p) for p in patient_id])
    patients = np.unique(patient_id)
    rows = []
    for gene, column in gene_index.items():
        values = counts[:, column]
        values = (np.asarray(values.todense()).ravel()
                  if hasattr(values, "todense") else np.asarray(values)).astype(float)
        detected = values >= min_umi
        nonzero_patients = sum(
            bool(detected[patient_id == p].any()) for p in patients)
        rows.append({
            "gene": gene,
            "n_cells": int(values.size),
            "n_patients": int(patients.size),
            "detection": float(detected.mean()),
            "mean_cp10k_proxy": float(values.mean()),
            "share_patients_nonzero": float(nonzero_patients / max(patients.size, 1)),
        })
    return pd.DataFrame(rows)


def _read_deposit(tar: Path, series_matrix: Path, *, pool_by: str):
    """Stack every scored sample's cells into one matrix, with its pooling key.

    ``CRC`` samples are dropped here because Amendment 1 excludes them, and the
    count of what was dropped is logged rather than left implicit. The gene
    index is taken from the FIRST sample's features and asserted identical on
    every other — CellRanger writes the same reference for a series, and a
    sample with a different one cannot share a column index. Concatenating
    across a changed reference would misalign every gene silently.
    """
    from scipy.sparse import vstack

    from src.reference.becker_io import (
        gene_symbols,
        pooling_key,
        read_series_matrix,
        read_triplet,
        sample_files,
    )

    metadata = read_series_matrix(series_matrix)
    metadata["pool_key"] = pooling_key(metadata, pool_by=pool_by)
    files = sample_files(tar).merge(metadata, on="gsm", how="inner")

    scored = files[files["arm"].notna()]
    log.info("  %d of %d samples carry an arm; %d dropped (CRC, Amendment 1)",
             len(scored), len(files), len(files) - len(scored))
    if scored.empty:
        raise FeasibilityError("no sample carries an arm after the CRC exclusion")

    blocks, keys, arms, reference = [], [], [], None
    for _, row in scored.iterrows():
        counts, _, features = read_triplet(tar, row)
        symbols = gene_symbols(features)
        if reference is None:
            reference = symbols
        elif not np.array_equal(symbols, reference):
            raise FeasibilityError(
                f"{row['gsm']} has a different gene index from the first "
                f"sample. Concatenating across a changed reference misaligns "
                f"every gene, and nothing about that raises on its own."
            )
        blocks.append(counts)
        keys.extend([row["pool_key"]] * counts.shape[0])
        arms.extend([row["arm"]] * counts.shape[0])

    gene_index = {}
    for gene in tuple(GENE_ROLES) + MATURE_MARKERS:
        hit = np.flatnonzero(reference == gene)
        if hit.size:
            gene_index[gene] = int(hit[0])
    absent = sorted(set(GENE_ROLES) - set(gene_index))
    missing_markers = sorted(set(MATURE_MARKERS) - set(gene_index))
    if missing_markers:
        log.warning("  mature markers not located: %s (of %d)",
                    missing_markers, len(MATURE_MARKERS))
    log.info("  stacked %d samples -> %d cells, pooled by %s (%d units)",
             len(blocks), sum(b.shape[0] for b in blocks), pool_by,
             len(set(keys)))
    return (vstack(blocks).tocsr(), gene_index, np.asarray(keys),
            np.asarray(arms), absent)


def lesion_inventory(series_matrix: Path) -> pd.DataFrame:
    """The durable, metadata-only Becker lesion inventory.

    This is deliberately separate from ``--inspect``: the question needs only
    GEO's small series matrix, not the 1.2 GB count tar, and the answer is an
    input inventory rather than a biological result. One lesion is one unique
    ``sample_id``; technical GSM replicates remain visible in ``n_tumour_rows``.
    """
    from src.reference.becker_io import read_series_matrix, tumour_lesion_counts

    return tumour_lesion_counts(read_series_matrix(series_matrix))


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--object", type=Path, default=None,
                        help="an h5ad, if the deposit ever arrives as one")
    parser.add_argument("--tar", type=Path, default=None,
                        help="GSE201348_RAW.tar")
    parser.add_argument("--series-matrix", type=Path, default=None,
                        help="GSE201348_series_matrix.txt.gz — REQUIRED with "
                             "--tar. The tar carries no arm labels at all.")
    parser.add_argument("--inspect", action="store_true",
                        help="report the file's structure and vocabulary, and "
                             "do nothing else. RUN THIS FIRST.")
    parser.add_argument(
        "--lesion-inventory", action="store_true",
        help="write the per-donor polyp inventory from --series-matrix alone; "
             "unique sample_id is a lesion and technical replicate rows remain visible",
    )
    parser.add_argument("--gene-column", default=None,
                        help="var column holding gene symbols, from --inspect")
    parser.add_argument("--patient-column", default=None,
                        help="obs column holding the patient id, from --inspect")
    parser.add_argument("--layer", default=None, help="counts layer, if not X")
    parser.add_argument(
        "--pool-by", choices=("donor", "lesion"), default="donor",
        help="the unit. Amendment 1: 'donor' is primary and confirmatory "
             "because it reproduces Chen_2021's shape; 'lesion' is secondary "
             "and exploratory. They are different estimands.",
    )
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    if args.inspect and args.lesion_inventory:
        raise SystemExit("choose --inspect or --lesion-inventory, not both")
    if args.tar and not args.series_matrix:
        raise SystemExit(
            "--tar needs --series-matrix. The tar carries NO metadata: a "
            "filename gives a GSM, a donor and a sample, and says nothing "
            "about whether it is a polyp or unaffected mucosa. Without the "
            "series matrix there are no arms."
        )
    if not args.tar and not args.object and not args.lesion_inventory:
        raise SystemExit("pass --tar (with --series-matrix), --object, or --lesion-inventory")
    if args.lesion_inventory and args.series_matrix is None:
        raise SystemExit("--lesion-inventory requires --series-matrix")
    for candidate in (args.tar, args.series_matrix, args.object):
        if candidate is not None and not candidate.exists():
            raise SystemExit(f"{candidate} not found")

    if args.lesion_inventory:
        inventory = lesion_inventory(args.series_matrix)
        path = write_versioned_table(
            inventory,
            "becker_lesion_inventory",
            seed=args.seed,
            results_dir=args.results_dir,
            allow_dirty=args.allow_dirty,
            notes=(
                "GSE201348 metadata inventory. One lesion is one unique "
                "sample_id; technical replicate GSM rows are retained separately."
            ),
            extra_meta={
                "source": "GSE201348_series_matrix.txt.gz",
                "source_kind": "GEO series metadata",
                "unit": "unique polyp sample_id within donor",
                "inference_unit": "donor; inventory does not promote lesions to n",
            },
        )
        log.info("wrote %s", path)
        return 0

    if args.inspect:
        report = (inspect_deposit(args.tar, args.series_matrix) if args.tar
                  else inspect(args.object))
        log.info("%s\nINSPECTION — assume nothing, map nothing\n%s",
                 "=" * 72, "=" * 72)
        for key, value in report.items():
            log.info("  %-32s %s", key, value)
        log.info(
            "\n  NEXT: the sample_type vocabulary above is what the arm mapping "
            "must be\n  built from, BY A HUMAN WHO HAS SEEN IT. Reading a label "
            "rather than the\n  patient grouping once put Chen_2021's usable "
            "pairs at zero when the true\n  number was 44."
        )
        found = (report.get("panel_genes_found")
                 or report.get("panel_genes_in_var_names")
                 or report.get("panel_genes_by_var_column"))
        if not found:
            log.error(
                "\n  NO PANEL GENE MATCHED IN ANY IDENTIFIER SPACE. Suspect the "
                "identifier\n  space before concluding the genes are absent — "
                "that error has been made\n  four times in this repository."
            )
            return 4
        return 0

    if args.object and (not args.gene_column or not args.patient_column):
        raise SystemExit(
            "--gene-column and --patient-column are required and have no "
            "defaults. Run --inspect first and read them off its output. "
            "Guessing them is how an empty intersection gets reported as a "
            "finding."
        )

    if args.tar:
        counts, gene_index, patients, cell_arms, absent = _read_deposit(
            args.tar, args.series_matrix, pool_by=args.pool_by)
        log.info("panel genes located in the features symbol column: %d of %d%s",
                 len(set(gene_index) & set(GENE_ROLES)), len(GENE_ROLES),
                 f" (absent: {absent})" if absent else "")
    else:
        import anndata

        adata = anndata.read_h5ad(str(args.object))
        symbols = adata.var[args.gene_column].astype(str).to_numpy()
        gene_index = {}
        for gene in tuple(GENE_ROLES) + MATURE_MARKERS:
            hit = np.flatnonzero(symbols == gene)
            if hit.size:
                gene_index[gene] = int(hit[0])
        absent = sorted(set(GENE_ROLES) - set(gene_index))
        log.info("panel genes located in var['%s']: %d of %d%s",
                 args.gene_column, len(set(gene_index) & set(GENE_ROLES)),
                 len(GENE_ROLES),
                 f" (absent: {absent})" if absent else "")
        counts = adata.layers[args.layer] if args.layer else adata.X
        patients = adata.obs[args.patient_column].to_numpy()
        # Invariant 1: no arm vector is not an arm vector of one value.
        cell_arms = np.full(counts.shape[0], None, dtype=object)

    panel_index = {g: i for g, i in gene_index.items() if g in GENE_ROLES}
    marker_index = {g: i for g, i in gene_index.items() if g in MATURE_MARKERS}
    if not panel_index:
        raise FeasibilityError(
            "no panel gene matched. Check the identifier space before "
            "concluding they are not measured — that error has been made four "
            "times in this repository."
        )

    # THE ARM. CHEN_BASELINE is Chen_2021's NORMAL-arm mature cells, and the
    # prereg's gate is on the normal arm. The first version of this pooled every
    # scored sample -- 43 polyps against 16 normals -- and compared that to a
    # normal-arm baseline. The polyps are the arm where a maturity marker is
    # EXPECTED to be low, so the hypothesis's own predicted effect was inside
    # the feasibility baseline, biasing the gate toward CANNOT RUN. That is not
    # a conservative choice; it is the wrong comparison.
    per_arm = []
    for arm in sorted({a for a in cell_arms.tolist() if a is not None}):
        mask = cell_arms == arm
        frame = detection_table(counts[mask], panel_index, patients[mask])
        frame.insert(0, "arm", arm)
        per_arm.append(frame)
    whole = detection_table(counts, panel_index, patients)
    whole.insert(0, "arm", "ALL_ARMS_POOLED")
    per_arm.append(whole)
    by_arm = gate(pd.concat(per_arm, ignore_index=True))
    # Arm-to-arm, on the detection scale. The reference arm's own health is not
    # a detail here: `normal` is a FAP donor's uninvolved mucosa and
    # `healthy_donor` is a different person's colon (becker_io's fifth arm), so
    # the difference between them is a property of the REFERENCE, and the
    # paired design leans on that reference entirely.
    wide = by_arm.pivot(index="gene", columns="arm", values="mu_becker")
    for a, b in (("healthy_donor", "normal"), ("normal", "tumour")):
        if a in wide.columns and b in wide.columns:
            by_arm[f"logfc_{a}_vs_{b}"] = by_arm["gene"].map(
                np.log(wide[a] / wide[b]))

    gate_arm = "normal" if (cell_arms == "normal").any() else "ALL_ARMS_POOLED"
    if gate_arm != "normal":
        log.warning("no normal arm in this object; gating on the pooled object, "
                    "which is NOT the pre-registered gate.")
    arm_mask = ((cell_arms == "normal") if gate_arm == "normal"
                else np.ones(counts.shape[0], dtype=bool))
    table = detection_table(counts[arm_mask], panel_index, patients[arm_mask])
    gated = gate(table)

    # ---------------------------------------------------------------- the gate
    # §3 names the MATURE CELLS of the reference arm. Everything above is every
    # cell in that arm, which is a lower bound. This is the pre-registered
    # quantity, and it is only taken when the markers were actually located.
    audit = pd.DataFrame()
    mature_labelled = False
    if marker_index:
        mature_mask = label_mature(counts, marker_index)
        gate_mask = arm_mask & mature_mask
        n_mature = int(gate_mask.sum())
        if n_mature < MIN_MATURE_CELLS:
            log.warning(
                "only %d nuclei carry >=%d of %s in the %s arm, below the %d "
                "needed; the labelled gate is NOT taken.",
                n_mature, MIN_MARKERS_FOR_MATURE, list(marker_index), gate_arm,
                MIN_MATURE_CELLS)
        else:
            mature_table = detection_table(counts[gate_mask], panel_index,
                                           patients[gate_mask])
            audit = enrichment_audit(table, mature_table)
            depth = np.asarray(counts.sum(axis=1)).ravel()
            log.info(
                "\n%s\nTHE MATURE LABEL — §3's actual gate\n%s",
                "=" * 72, "=" * 72)
            log.info("  %d of %d nuclei in the %s arm call as mature "
                     "(>=%d of %s)", n_mature, int(arm_mask.sum()), gate_arm,
                     MIN_MARKERS_FOR_MATURE, list(marker_index))
            log.info("  median UMIs/nucleus: %.0f in the label, %.0f in the arm "
                     "-> %.2fx depth",
                     float(np.median(depth[gate_mask])),
                     float(np.median(depth[arm_mask])),
                     float(np.median(depth[gate_mask])
                           / max(np.median(depth[arm_mask]), 1e-9)))
            log.info("\n  ENRICHMENT AUDIT — is this identity or is it depth?")
            log.info("%s", audit[["gene", "role", "detection_whole_arm",
                                  "detection_mature", "log_enrichment",
                                  "beyond_control_band"]].to_string(index=False))
            critical = audit.loc[audit["gene"] == CRITICAL_GENE]
            if len(critical) and not bool(critical["beyond_control_band"].iloc[0]):
                log.warning(
                    "\n  %s's enrichment sits INSIDE the band the controls "
                    "describe (%.3f to %.3f).\n  The label selected deeper "
                    "nuclei, not mature ones, and a floor cleared on that\n  "
                    "basis is cleared by library size. Read the gate below "
                    "with that in front of it.",
                    CRITICAL_GENE, float(critical["control_band_low"].iloc[0]),
                    float(critical["control_band_high"].iloc[0]))
            table, gated, mature_labelled = mature_table, gate(mature_table), True

            # THE SAME RULE ON EVERY ARM. Descriptive, and labelled so: only
            # the gate arm is pre-registered. It is here because it is the
            # internal control this deposit can supply for free -- the same
            # markers, the same threshold, the same protocol, different
            # tissue. A target that enriches beyond its controls in one arm and
            # not another is saying something about the arm; one that never
            # outruns its controls anywhere is saying the label is depth.
            audit["arm"] = gate_arm
            audit["exploratory"] = False
            extra = []
            for other in sorted({a for a in cell_arms.tolist() if a is not None}):
                if other == gate_arm:
                    continue
                other_mask = (cell_arms == other) & mature_mask
                if int(other_mask.sum()) < MIN_MATURE_CELLS:
                    log.info("  %s arm: %d mature nuclei, below %d — not audited",
                             other, int(other_mask.sum()), MIN_MATURE_CELLS)
                    continue
                frame = enrichment_audit(
                    detection_table(counts[cell_arms == other], panel_index,
                                    patients[cell_arms == other]),
                    detection_table(counts[other_mask], panel_index,
                                    patients[other_mask]))
                frame["arm"] = other
                frame["exploratory"] = True
                extra.append(frame)
            if extra:
                audit = pd.concat([audit] + extra, ignore_index=True)
                log.info("\n  THE SAME AUDIT ON EVERY ARM — exploratory except "
                         "'%s'", gate_arm)
                log.info("%s", audit.pivot(index="gene", columns="arm",
                                           values="log_enrichment")
                         .to_string(float_format=lambda v: f"{v:+.3f}"))
                beyond = audit.loc[audit["beyond_control_band"]]
                log.info("  beyond its own arm's control band: %s",
                         ", ".join(f"{r.gene}@{r.arm}" for r in beyond.itertuples())
                         or "nothing, in any arm")

    outcome = verdict(gated, mature_labelled=mature_labelled,
                      audit=audit if not audit.empty else None)

    log.info("\n%s\nDETECTION BY ARM — the gate reads the '%s' row\n%s",
             "=" * 72, gate_arm, "=" * 72)
    log.info("%s", by_arm.pivot(index=["gene", "role"], columns="arm",
                                values="detection").to_string())

    log.info("\n%s\nTHE PRE-REGISTERED GATE — prereg §3, %s\n%s", "=" * 72,
             (f"mature cells of the {gate_arm} arm" if mature_labelled
              else f"the {gate_arm} arm, ALL cells (a LOWER BOUND)"), "=" * 72)
    log.info("  detection >= %.2f AND non-zero in >= %.0f%% of patients",
             MIN_DETECTION, 100 * MIN_PATIENT_SHARE_NONZERO)
    log.info("%s", gated[["gene", "role", "detection", "chen_baseline",
                          "fold_mu_vs_chen", "log_fc_vs_chen",
                          "share_patients_nonzero", "passes"]]
             .to_string(index=False))
    # A criterion that cannot fail on the data in hand has not been applied to
    # it. At thousands of nuclei per donor, P(all zero) is e^-{n*p}: for the
    # LOWEST detection in this table that is already astronomically small, so
    # the patient-share arm of the gate is inert and the verdict rests entirely
    # on the detection floor. Say so rather than let two criteria be read where
    # only one discriminated.
    cells_per_patient = gated["n_cells"].max() / max(gated["n_patients"].max(), 1)
    weakest = gated["detection"].min()
    log_p_all_zero = cells_per_patient * np.log1p(-weakest)
    if (gated["share_patients_nonzero"] >= MIN_PATIENT_SHARE_NONZERO).all():
        log.info(
            "\n  NOTE: the patient-share criterion passed 6/6 at 100%% and was "
            "INERT.\n  With ~%.0f nuclei per donor, the weakest gene here "
            "(detection %.4f) has\n  P(all-zero in a donor) = e^%.0f. It could "
            "not have fired. The verdict\n  rests on the detection floor "
            "alone.", cells_per_patient, weakest, log_p_all_zero)

    log.info("\n%s\nVERDICT\n%s", "=" * 72, "=" * 72)
    log.info("  %s", outcome["verdict"])
    log.info("  %s", outcome["detail"])
    if mature_labelled:
        log.info(
            "\n  NOTE: this IS §3's quantity — mature cells of the reference "
            "arm. Read the\n  enrichment audit above before reading the "
            "verdict: the label is built from\n  marker detection, which rises "
            "with depth, so a floor cleared without the\n  target outrunning "
            "the controls was cleared by library size."
        )
    else:
        log.info(
            "\n  NOTE: this is the reference arm but NOT its mature cells, "
            "which is what the\n  prereg names. Mature cells ENRICH for these "
            "markers, so this remains a LOWER\n  bound: a gene passing here "
            "passes the real gate. A gene failing here needs\n  the labelled "
            "reading before it is called dead — and B1 is UNDETERMINED until\n"
            "  that reading exists, not refuted."
        )

    meta = {
        "prereg": "docs/prereg_becker_replication.md",
        "gate": {"min_detection": MIN_DETECTION,
                 "min_patient_share_nonzero": MIN_PATIENT_SHARE_NONZERO},
        "critical_gene": CRITICAL_GENE,
        "verdict": outcome,
        "chen_baseline": CHEN_BASELINE,
        "mature_label": {
            "markers": list(MATURE_MARKERS),
            "min_markers": MIN_MARKERS_FOR_MATURE,
            "min_cells": MIN_MATURE_CELLS,
            "none_is_on_the_panel": True,
            "why": (
                "invariant 2 — a target gene in the label returns the threshold "
                "it was given. The whole panel is excluded so the controls stay "
                "usable inside the label."
            ),
            "taken": bool(mature_labelled),
            "depth_audit_controls": list(DEPTH_AUDIT_CONTROLS),
        },
        "gene_column": args.gene_column,
        "patient_column": args.patient_column,
        "layer": args.layer or "X",
        "genes_absent_from_object": absent,
        "scope": (
            f"the {gate_arm} arm, but not its mature cells, which is what the "
            f"prereg names. Mature cells enrich for these markers, so this is "
            f"still a lower bound on the pre-registered gate. Detection for "
            f"every arm is in the by_arm table."
        ),
        "exploratory": False,
        "pre_registered": True,
    }
    tables = [(gated, "becker_feasibility"),
              (by_arm, "becker_feasibility_by_arm")]
    if not audit.empty:
        tables.append((audit, "becker_feasibility_enrichment_audit"))
    for frame, name in tables:
        log.info("wrote %s", write_versioned_table(
            frame, name, seed=args.seed, results_dir=args.results_dir,
            allow_dirty=args.allow_dirty, extra_meta=meta,
        ))
    return 0 if outcome["verdict"] != "CANNOT RUN" else 5


if __name__ == "__main__":
    sys.exit(main())
