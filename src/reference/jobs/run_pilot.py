#!/usr/bin/env python
"""Run the five-patient pilot through the stages that exist. W1, week 2.

    python src/reference/jobs/run_pilot.py

The pilot is the critical path: W2's harness is blocked on it, and §8.3 moves the
week-5 gate to week 7 if it slips past week 3. It is a thin vertical slice —
every stage runs, none is tuned.

Stages implemented here:

  load      the five patients, reading only their columns off disk
  qc        per-cell metrics, per-batch thresholds with rationale, applied
  ambient   contamination from impossible genes, per sample

Stages still to be written (weeks 2-4), so the slice is not yet complete:

  malignancy   inferCNV / CopyKAT
  labels       axes 1 and 2 across four granularity rungs
  signature    _select_markers, then the pilot S matrix for W2

Artifacts land in results/{date}_{sha}/ with a provenance sidecar.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from src.common.io import write_versioned_table
from src.common.panel import granularity_rungs, tier_genes
from src.common.paths import s_matrix_path
from src.common.provenance import DEFAULT_SEED, set_global_seeds
from src.reference.ambient import contamination_by_sample, soup_profile_from_cells
from src.reference.gene_index import read_gene_index
from src.reference.ingest import (
    assign_compartments,
    read_gse178341,
    read_gse178341_clusters,
    read_gse178341_metadata,
)
from src.reference.labels import (
    NON_EPITHELIAL,
    UNRESOLVED,
    assign_labels,
    axis_tie_fraction,
    describe_labels,
    label_column,
    label_depth_confounding,
    mature_cell_counts,
)
from src.reference.qc import (
    apply_qc,
    cell_qc_metrics,
    differential_retention,
    qc_summary,
    qc_thresholds,
)
from src.reference.signature import build_signature_sparse

PILOT = ["C122", "C165", "C107", "C138", "C162"]

#: Compartment label from the authors' cluster file. Safe to use for the
#: contamination mask and for the S matrix's non-epithelial columns —
#: distinguishing epithelium from immune from stroma does not depend on the
#: differentiation markers under test. NOT safe for deciding which epithelial
#: cells are mature; W1 builds those labels from the frozen axes in weeks 3-4.
#:
#: The first run used an EPCAM>=1 gate instead and returned contamination of
#: exactly 1.0 on eight samples: ambient EPCAM is everywhere, so the "epithelial"
#: mask was effectively all cells and the estimator degenerated.
EPITHELIAL_COMPARTMENT = "Epi"

#: Below this many epithelial cells a per-sample contamination estimate is noise.
MIN_CELLS_FOR_CONTAMINATION = 20

#: Shared gene index version. W3 emits bulk on the same one.
GENE_INDEX_VERSION = "1.0.0"

#: S matrix version. Bump rather than overwrite — results cite it by version.
S_MATRIX_VERSION = "0.1.0-pilot"

#: Within [500, 2000] per §2.1 error 4. Modest for a five-patient pilot.
SIGNATURE_GENES = 800


def main() -> int:
    set_global_seeds(DEFAULT_SEED)
    data = Path(os.environ.get("BRP_DATA_DIR", "data")) / "raw" / "GSE178341"
    h5 = data / "GSE178341_crc10x_full_c295v4_submit.h5"
    meta_csv = data / "GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz"
    cluster_csv = data / "GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz"
    if not h5.exists():
        raise SystemExit(f"missing {h5}")

    print("=" * 70)
    print(f"PILOT: {', '.join(PILOT)}")
    print("=" * 70)

    print("\n--- load ---")
    adata = read_gse178341(h5, patients=PILOT)
    print(f"{adata.n_obs:,} cells x {adata.n_vars:,} genes")
    metadata = read_gse178341_metadata(meta_csv)
    shared = adata.obs.index.intersection(metadata.index)
    for column in ("PROCESSING_TYPE", "MMRStatus", "MLH1Status"):
        if column in metadata.columns:
            adata.obs[column] = metadata.loc[shared, column].reindex(adata.obs.index)
    print(adata.obs.groupby(["patient_id", "tissue"], observed=True).size().to_string())

    print("\n--- qc ---")
    metrics = cell_qc_metrics(
        adata.X,
        adata.var["gene_symbol"],
        batch=adata.obs["sample_id"],
        patient_id=adata.obs["patient_id"],
        tissue=adata.obs["tissue"],
    )
    thresholds = qc_thresholds(metrics)
    passes = apply_qc(metrics, thresholds)
    summary = qc_summary(metrics, passes)
    print(f"retained {int(passes.sum()):,} of {len(passes):,} cells "
          f"({passes.mean():.1%})")
    print(summary.to_string(index=False))
    print("\nthresholds (per batch, with rationale — the week-1 deliverable):")
    print(thresholds[["batch", "metric", "lower", "upper", "n_cells", "n_failed"]]
          .to_string(index=False))

    print("\n--- differential retention: does QC cut the two arms unequally? ---")
    retention = differential_retention(metrics, passes)
    print(retention.to_string(index=False))
    if retention["flagged"].any():
        print(
            "\n!! QC removes cells at materially different rates between tumour and\n"
            "   normal in the flagged patients. The compositional term is the\n"
            "   difference between those two arms, so this biases it directly.\n"
            "   See docs/open_decisions.md #12 before trusting any composition."
        )

    print("\n--- ambient contamination (impossible genes) ---")
    symbols = list(adata.var["gene_symbol"])
    contamination = pd.DataFrame()
    if not cluster_csv.exists():
        print(f"!! missing {cluster_csv}; cannot build the epithelial mask")
    else:
        clusters = read_gse178341_clusters(cluster_csv)
        adata.obs["clTopLevel"] = clusters["clTopLevel"].reindex(adata.obs.index)
        print("compartments present:")
        print(adata.obs["clTopLevel"].value_counts().to_string())
        epithelial = (
            adata.obs["clTopLevel"].astype(str) == EPITHELIAL_COMPARTMENT
        ).to_numpy() & np.asarray(passes, dtype=bool)
        print(f"\nepithelial cells passing QC: {int(epithelial.sum()):,}")
        contamination = contamination_by_sample(
            adata.X, symbols,
            sample_id=adata.obs["sample_id"], cell_mask=epithelial,
        )
        # Sorting matters here: in a CD45-enriched sample the few epithelial
        # cells sit in immune-dominated soup, so a high estimate is expected
        # rather than a failure. And an estimate from <20 cells is noise.
        # .astype(str) before the map: PROCESSING_TYPE is categorical, and
        # filling a category column with "" raises rather than adding a level.
        sorting = (
            adata.obs.assign(_pt=adata.obs["PROCESSING_TYPE"].astype(str))
            .groupby("sample_id", observed=True)["_pt"]
            .agg(lambda s: s.mode().iat[0] if len(s) else "")
            .astype(str)
        )
        contamination["processing_type"] = (
            contamination["sample_id"].map(sorting).astype(object).fillna("")
        )
        contamination["reliable"] = contamination["n_cells"] >= MIN_CELLS_FOR_CONTAMINATION
        print(contamination.to_string(index=False))

        usable = contamination[
            contamination["reliable"] & (contamination["processing_type"] == "unsorted")
        ]
        if len(usable):
            print(
                f"\nunsorted samples with >={MIN_CELLS_FOR_CONTAMINATION} epithelial "
                f"cells (n={len(usable)}): median contamination "
                f"{usable['contamination'].median():.1%}, "
                f"max {usable['contamination'].max():.1%}"
            )

        soup = soup_profile_from_cells(adata.X, symbols).sort_values(ascending=False)
        print("\ntop 10 soup genes (pooled across the pilot):")
        print(soup.head(10).to_string())

    print("\n--- labels: 2 axes x 4 rungs ---")
    labels = pd.DataFrame()
    counts = pd.DataFrame()
    if not cluster_csv.exists():
        print("!! no cluster file; cannot assign compartments, so no labels")
    else:
        compartment = assign_compartments(clusters).reindex(adata.obs.index)
        adata.obs["compartment"] = compartment
        print(compartment.value_counts().to_string())

        # Target set for this pilot run: tier A, the compositional control. Tier A
        # does not collide with either axis, so both are usable. A run testing
        # MUC2 or TFF3 would have to pass axes=["stem_pole"] (open decision #1).
        targets = tier_genes("A")
        print(f"\ntarget set for this run: {targets}")

        # Decision #11: sorted samples have had their cell-type composition
        # manipulated, so they cannot carry a compositional estimate. This filter
        # was decided but never applied — earlier runs pooled C107's unsorted,
        # CD45pMACS and LiveMACS tumour samples into one "tumour" arm.
        unsorted = (adata.obs["PROCESSING_TYPE"].astype(str) == "unsorted").to_numpy()
        keep = np.asarray(passes, dtype=bool) & unsorted
        print(f"\nrestricted to UNSORTED samples for composition: "
              f"{int(keep.sum()):,} of {int(np.asarray(passes).sum()):,} QC-passing "
              f"cells (decision #11)")
        labels = assign_labels(
            adata.X[keep],
            adata.var["gene_symbol"],
            compartment=compartment.to_numpy()[keep],
            sample_id=adata.obs["sample_id"].to_numpy()[keep],
            target_genes=targets,
            # Cut points come from each patient's own normal, applied to their
            # tumour. Without this the mature fraction is pinned to the quantile
            # in every arm and the compositional term cannot move.
            tissue=adata.obs["tissue"].to_numpy()[keep],
            patient_id=adata.obs["patient_id"].to_numpy()[keep],
            # Depth matching: without it the maturity call partly measures
            # sequencing depth (decision #14).
            depth_quantile=0.10,
            seed=DEFAULT_SEED,
            index=adata.obs.index[keep],
        )
        # The depth target trades ties against unresolved cells: a lower target
        # keeps more cells but detects fewer marker counts, so the tied block
        # grows. Report the trade rather than hiding it behind one default.
        print("\ndepth target trade-off (tie fraction vs cells lost):")
        epi_mask = compartment.to_numpy()[keep] == "epithelial"
        epi_totals = np.asarray(adata.X[keep].sum(axis=1)).ravel()[epi_mask]
        for q in (0.05, 0.10, 0.25, 0.50):
            target = float(np.quantile(epi_totals, q))
            stats = axis_tie_fraction(
                adata.X[keep], adata.var["gene_symbol"], "stem_pole",
                target_genes=targets, epithelial=epi_mask,
                depth_target=target,
            )
            print(f"  q={q:<5} target {target:>8,.0f}  "
                  f"stem_pole tied {stats['tied_fraction']:.1%}  "
                  f"cells lost {q:.0%}")

        print("\naxis resolution — how much of each score is one tied block:")
        for axis in ("stem_pole", "opposite_lineage"):
            stats = axis_tie_fraction(
                adata.X[keep], adata.var["gene_symbol"], axis,
                target_genes=targets,
                epithelial=(compartment.to_numpy()[keep] == "epithelial"),
            )
            print(f"  {axis:<18} tied {stats['tied_fraction']:.1%} "
                  f"(largest block {stats['largest_tied_block']:,} of "
                  f"{stats['n_cells']:,}, {stats['n_distinct_scores']:,} distinct)")

        print(f"\nlabelled {len(labels):,} QC-passing cells, "
              f"{len(labels.columns)} columns")
        print(describe_labels(labels).to_string(index=False))

        depth_report = label_depth_confounding(labels, metrics[keep].reset_index(drop=True))
        print("\nis the maturity call tracking DEPTH rather than biology?")
        print("(counts_ratio < 1 means mature cells are shallower — an artifact)")
        print(depth_report.to_string(index=False))
        if depth_report["flagged"].any():
            print("\n!! flagged rows: the mature fraction is partly a depth measurement,\n"
                  "   so Delta(mature fraction) between arms of differing depth is partly\n"
                  "   an artifact. Do not quote a compositional number from those rungs.")

        counts = mature_cell_counts(
            labels,
            patient_id=adata.obs["patient_id"].to_numpy()[keep],
            tissue=adata.obs["tissue"].to_numpy()[keep],
        )
        print("\nmature-cell counts. Cut points come from each patient's NORMAL "
              "arm,\nso normal sits near the rung's target fraction and TUMOUR IS "
              "FREE TO DIFFER —\nthat difference is the compositional term.")
        print(counts.sort_values(["granularity_rung", "patient_id", "tissue"])
              .to_string(index=False))

    print("\n--- pilot S matrix (the W2 handoff) ---")
    if len(labels):
        try:
            index = read_gene_index(GENE_INDEX_VERSION)
        except Exception as exc:
            print(f"!! no gene index: {exc}")
            print("   run: python src/reference/jobs/emit_gene_index.py --version "
                  f"{GENE_INDEX_VERSION}")
        else:
            # Sparse throughout: densifying 25,959 x 43,113 in float64 is
            # 8.3 GB and will not allocate. build_signature_sparse aggregates
            # per cell type without materialising that array.
            ensembl = list(adata.var["ensembl_id"])
            count_matrix = adata.X[keep]
            for rung in granularity_rungs():
                column = labels[label_column("stem_pole", rung)].astype(str).to_numpy()
                # Compartment for non-epithelial cells; the rung's own bins
                # within epithelium. §2.1 error 3 needs stromal/immune/endothelial.
                cell_type = np.where(
                    np.isin(column, [NON_EPITHELIAL, UNRESOLVED]),
                    compartment.to_numpy()[keep],
                    column,
                )
                try:
                    s_matrix = build_signature_sparse(
                        count_matrix, ensembl, cell_type,
                        target_genes=targets, gene_index=index,
                        n_genes=SIGNATURE_GENES,
                    )
                except Exception as exc:
                    print(f"  {rung:<16} FAILED: {exc}")
                    continue
                path = s_matrix_path(rung, S_MATRIX_VERSION)
                path.parent.mkdir(parents=True, exist_ok=True)
                s_matrix.to_parquet(path)
                print(f"  {rung:<16} {s_matrix.shape[0]:,} genes x "
                      f"{s_matrix.shape[1]} columns -> {path.name}")

            print("\n  Hand W2 these plus the labelled object. Have them load it")
            print("  against the frozen schema before week 2 closes — if they can")
            print("  read it without asking a question, the interface holds.")

    print("\n" + "=" * 70)
    print("WHAT STILL BLOCKS A QUOTABLE COMPOSITIONAL NUMBER")
    print("=" * 70)
    print("""
  APPLIED in this run
    depth matching          maturity scored at a common depth (decision #14)
    unsorted samples only   sorted samples cannot carry composition (#11)
    mito cap 50%            was cutting normal 22.7 points harder (#12)
    reference thresholds    cuts from each patient's own normal arm

  NOT YET, and each can move the answer
    malignancy calls        the "tumour" arm still contains NON-MALIGNANT
                            epithelium from tumour samples. inferCNV is weeks
                            2-3 and unwritten, so today's contrast is
                            sample-of-origin, not malignant-vs-normal. This is
                            the largest remaining gap.
    ambient correction      contamination is measured (median 1.6%) but never
                            subtracted; scores run on uncorrected counts
    doublet removal         flag_doublets() is a stub. Doublets co-express
                            markers from two compartments and inflate both the
                            compartment call and the maturity score
    tiny samples            C138_T_0_0_0_c2_v2 contributes 29 cells and is
                            pooled with samples 100x larger
    1,108 cells             370,115 here against a published 371,223, still
                            unreconciled
""")

    print("\n--- pct_mito distribution, for open decision #12 ---")
    print("Pick the cap from this, not from convention. 20% is a lymphocyte number;")
    print("colonic epithelium runs higher, and much of it here is ambient.")
    if "clTopLevel" in adata.obs.columns:
        mito = metrics.copy()
        mito["compartment"] = adata.obs["clTopLevel"].astype(str).to_numpy()
        table = mito.groupby(["compartment", "tissue"], observed=True)["pct_mito"].describe(
            percentiles=[0.5, 0.75, 0.9, 0.95]
        )
        print(table[["count", "50%", "75%", "90%", "95%", "max"]].to_string())
        epi = mito[mito["compartment"] == EPITHELIAL_COMPARTMENT]
        for cap in (20, 30, 40, 50):
            kept = (epi["pct_mito"] <= cap).mean()
            by_tissue = epi.groupby("tissue", observed=True)["pct_mito"].apply(
                lambda s, c=cap: (s <= c).mean()
            )
            gap = abs(by_tissue.get("tumour", np.nan) - by_tissue.get("normal", np.nan))
            print(f"  cap {cap:>3}%: epithelium kept {kept:.1%}, tumour/normal gap {gap:.1%}")

    print("\n--- compartments available for the S matrix ---")
    print("§2.1 error 3 requires stromal, immune AND endothelial columns.")
    if "clMidwayPr" in locals().get("clusters", pd.DataFrame()).columns:
        midway = clusters["clMidwayPr"].reindex(adata.obs.index)
        print(midway.value_counts().to_string())

    print("\n--- what is in the cluster file (for the real labels) ---")
    if cluster_csv.exists():
        head = pd.read_csv(cluster_csv, nrows=3)
        print("columns:", list(head.columns))
        print(head.to_string())
    else:
        print(f"missing {cluster_csv}")

    print("\n--- writing artifacts ---")
    for frame, name in (
        (thresholds, "pilot_qc_thresholds"),
        (summary, "pilot_qc_summary"),
        (retention, "pilot_differential_retention"),
        (contamination, "pilot_contamination"),
        (describe_labels(labels) if len(labels) else labels, "pilot_label_summary"),
        (counts, "pilot_mature_cell_counts"),
    ):
        if not isinstance(frame, pd.DataFrame):
            print(f"  !! {name}: expected a DataFrame, got "
                  f"{type(frame).__name__} — skipping")
            continue
        if len(frame):
            path = write_versioned_table(
                frame, name, seed=DEFAULT_SEED, allow_dirty=True,
                notes=f"five-patient pilot: {', '.join(PILOT)}",
            )
            print(f"  {path}")

    print("\nNEXT: malignancy calls, then labels, then the pilot S matrix for W2.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
