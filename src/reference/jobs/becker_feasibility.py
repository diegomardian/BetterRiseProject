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
    out["fold_vs_chen"] = out["detection"] / out["chen_baseline"]
    out["clears_detection"] = out["detection"] >= MIN_DETECTION
    out["clears_patient_share"] = (
        out["share_patients_nonzero"] >= MIN_PATIENT_SHARE_NONZERO)
    out["passes"] = out["clears_detection"] & out["clears_patient_share"]
    out["role"] = out["gene"].map(GENE_ROLES)
    return out.sort_values("detection", ascending=False, ignore_index=True)


def verdict(gated: pd.DataFrame) -> dict:
    """§3's outcome table, taken by code rather than by a reader."""
    passing = set(gated.loc[gated["passes"], "gene"])
    failing = sorted(set(gated["gene"]) - passing)

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

    blocks, keys, reference = [], [], None
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

    gene_index = {}
    for gene in GENE_ROLES:
        hit = np.flatnonzero(reference == gene)
        if hit.size:
            gene_index[gene] = int(hit[0])
    absent = sorted(set(GENE_ROLES) - set(gene_index))
    log.info("  stacked %d samples -> %d cells, pooled by %s (%d units)",
             len(blocks), sum(b.shape[0] for b in blocks), pool_by,
             len(set(keys)))
    return vstack(blocks).tocsr(), gene_index, np.asarray(keys), absent


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

    if args.tar and not args.series_matrix:
        raise SystemExit(
            "--tar needs --series-matrix. The tar carries NO metadata: a "
            "filename gives a GSM, a donor and a sample, and says nothing "
            "about whether it is a polyp or unaffected mucosa. Without the "
            "series matrix there are no arms."
        )
    if not args.tar and not args.object:
        raise SystemExit("pass --tar (with --series-matrix) or --object")
    for candidate in (args.tar, args.series_matrix, args.object):
        if candidate is not None and not candidate.exists():
            raise SystemExit(f"{candidate} not found")

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
        counts, gene_index, patients, absent = _read_deposit(
            args.tar, args.series_matrix, pool_by=args.pool_by)
        log.info("panel genes located in the features symbol column: %d of %d%s",
                 len(gene_index), len(GENE_ROLES),
                 f" (absent: {absent})" if absent else "")
    else:
        import anndata

        adata = anndata.read_h5ad(str(args.object))
        symbols = adata.var[args.gene_column].astype(str).to_numpy()
        gene_index = {}
        for gene in GENE_ROLES:
            hit = np.flatnonzero(symbols == gene)
            if hit.size:
                gene_index[gene] = int(hit[0])
        absent = sorted(set(GENE_ROLES) - set(gene_index))
        log.info("panel genes located in var['%s']: %d of %d%s",
                 args.gene_column, len(gene_index), len(GENE_ROLES),
                 f" (absent: {absent})" if absent else "")
        counts = adata.layers[args.layer] if args.layer else adata.X
        patients = adata.obs[args.patient_column].to_numpy()

    if not gene_index:
        raise FeasibilityError(
            "no panel gene matched. Check the identifier space before "
            "concluding they are not measured — that error has been made four "
            "times in this repository."
        )

    table = detection_table(counts, gene_index, patients)
    gated = gate(table)
    outcome = verdict(gated)

    log.info("\n%s\nTHE PRE-REGISTERED GATE — prereg §3\n%s", "=" * 72, "=" * 72)
    log.info("  detection >= %.2f AND non-zero in >= %.0f%% of patients",
             MIN_DETECTION, 100 * MIN_PATIENT_SHARE_NONZERO)
    log.info("%s", gated[["gene", "role", "detection", "chen_baseline",
                          "fold_vs_chen", "share_patients_nonzero",
                          "passes"]].to_string(index=False))
    log.info("\n%s\nVERDICT\n%s", "=" * 72, "=" * 72)
    log.info("  %s", outcome["verdict"])
    log.info("  %s", outcome["detail"])
    log.info(
        "\n  NOTE: this measures the WHOLE object, not the mature cells of the "
        "normal arm\n  the prereg names. Mature cells ENRICH for these markers, "
        "so this is a LOWER\n  bound: a gene passing here passes the real gate. "
        "A gene failing here needs\n  the labelled reading before it is called "
        "dead."
    )

    meta = {
        "prereg": "docs/prereg_becker_replication.md",
        "gate": {"min_detection": MIN_DETECTION,
                 "min_patient_share_nonzero": MIN_PATIENT_SHARE_NONZERO},
        "critical_gene": CRITICAL_GENE,
        "verdict": outcome,
        "chen_baseline": CHEN_BASELINE,
        "gene_column": args.gene_column,
        "patient_column": args.patient_column,
        "layer": args.layer or "X",
        "genes_absent_from_object": absent,
        "scope": (
            "whole object, not the mature cells of the normal arm. Mature cells "
            "enrich for these markers, so this is a lower bound on the "
            "pre-registered gate."
        ),
        "exploratory": False,
        "pre_registered": True,
    }
    log.info("\nwrote %s", write_versioned_table(
        gated, "becker_feasibility", seed=args.seed,
        results_dir=args.results_dir, allow_dirty=args.allow_dirty,
        extra_meta=meta,
    ))
    return 0 if outcome["verdict"] != "CANNOT RUN" else 5


if __name__ == "__main__":
    sys.exit(main())
