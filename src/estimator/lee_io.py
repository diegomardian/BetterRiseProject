"""Real-file adapter for the Lee cohorts (GSE132465 SMC, GSE144735 KUL3). W4.

Kept separate from ``ingest.py`` on purpose: ``ingest.py``'s ``qc_flags()`` is
generic MAD-outlier math over an already-built metrics frame, reusable for any
cohort; this module is the GEO-file-format-specific adapter that builds that
frame (and the expression/label frames) from these two specific files.
``ingest.py``, ``labels.py`` and ``kitagawa.py`` are called here, never edited.

Real-file facts this module encodes, verified by direct inspection before
writing any of it (see the W4 plan, 2026-08-15):

- Both matrices are tab-separated, genes-as-rows, integer UMI counts, on the
  SAME 33,694-gene set but in DIFFERENT row order — every lookup here is by
  gene symbol, never by row number.
- MUC2 is absent from the gene index in both cohorts. ``opposite_lineage``
  therefore runs on 3 of its 4 markers for Lee specifically — surfaced via
  ``LeeCohort.axis_gene_coverage``, not silently absorbed.
- Only SMC01-SMC10 (10 of 24 SMC patients) have matched Normal tissue; the
  other 14 are Tumor-only. Kitagawa's ``decompose()`` is a two-arm
  comparison, so those patients cannot produce a row — excluded, loudly,
  via ``LeeCohort.excluded_patients`` and a warning at load time. KUL3's 6
  patients all have both arms (plus some have a third class, Border).
- KUL3 has a third ``Class`` value, ``Border``, that SMC does not have and
  that Kitagawa's two-arm math has no slot for. Border cells are labelled
  and carried through like any other cell, but excluded from the
  decompose_cohort-bound summary.
"""

from __future__ import annotations

import gzip
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Final, Literal

import numpy as np
import pandas as pd

from src.estimator.ingest import qc_flags
from src.estimator.labels import label_cohort

logger = logging.getLogger(__name__)

#: The 13 mitochondrial-encoded protein-coding genes, verified present in
#: both Lee matrices by direct inspection. Hardcoded rather than a runtime
#: "^MT-" regex because it is now a confirmed fact, not a runtime unknown —
#: see the module docstring. If this is ever reused on a reference where
#: none of these are found, fail loudly (below) rather than silently
#: emitting pct_mito=NaN, which would make qc_flags' mito check silently
#: pass every cell.
MITO_GENES: Final[tuple[str, ...]] = (
    "MT-ATP6", "MT-ATP8", "MT-CO1", "MT-CO2", "MT-CO3", "MT-CYB",
    "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6",
)  # fmt: skip

#: Literal accession strings — matches data/manifest.csv and CLAUDE.md's own
#: wording ("GSE132465 (SMC) and GSE144735 (KUL3) ... two separate study_id
#: values").
STUDY_IDS: Final[dict[str, str]] = {"smc": "GSE132465", "kul3": "GSE144735"}

_CLASS_TO_TISSUE: Final[dict[str, str]] = {
    "Tumor": "tumour",
    "Normal": "normal",
    "Border": "border",
}

_ANNOTATION_COLUMNS = {"Index", "Patient", "Class", "Sample", "Cell_type", "Cell_subtype"}


def load_annotation(path: Path, *, study_id: str) -> pd.DataFrame:
    """Per-cell metadata: patient_id, study_id, tissue, author labels.

    ``Class`` is mapped to ``tissue`` via ``_CLASS_TO_TISSUE`` — raises on any
    value outside {Tumor, Normal, Border} so a GEO vocabulary change fails
    loudly here rather than silently producing an unmapped tissue downstream.
    """
    annotation = pd.read_csv(path, sep="\t")
    missing = _ANNOTATION_COLUMNS - set(annotation.columns)
    if missing:
        raise ValueError(f"annotation file is missing column(s): {sorted(missing)}")

    bad_classes = set(annotation["Class"].unique()) - set(_CLASS_TO_TISSUE)
    if bad_classes:
        raise ValueError(
            f"unrecognized Class value(s) {sorted(bad_classes)} in {path.name} — "
            f"expected only {sorted(_CLASS_TO_TISSUE)}"
        )

    return pd.DataFrame(
        {
            "patient_id": annotation["Patient"].to_numpy(),
            "study_id": study_id,
            "tissue": annotation["Class"].map(_CLASS_TO_TISSUE).to_numpy(),
            "author_cell_type": annotation["Cell_type"].to_numpy(),
            "author_cell_subtype": annotation["Cell_subtype"].to_numpy(),
        },
        index=annotation["Index"],
    )


def stream_matrix_stats(
    path: Path,
    *,
    genes_of_interest: list[str],
    mito_genes: tuple[str, ...] = MITO_GENES,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, list[str]]]:
    """Single streaming pass over a genes-as-rows raw UMI matrix.

    Parses one gene row at a time with ``numpy.fromstring`` rather than
    ``pandas.read_csv(chunksize=...)``. This isn't a style preference:
    pandas' block-manager has real fixed overhead per chunk that scales with
    *column* count, and these matrices are 27k-64k columns wide — a
    chunked-pandas version of this function was benchmarked at ~74ms/row
    (~42 minutes projected for SMC's 33,694 rows); this version is
    ~4.5ms/row (~2.5 minutes), a 16x difference, because it never builds a
    DataFrame during the hot loop at all. ``gzip.open(path, "rt")`` streams
    decompression, so the file is never fully decompressed to disk or memory,
    and only a few small accumulator arrays (length n_cells) plus the
    ``genes_of_interest`` rows are held at once — nowhere near the ~8.6GB a
    dense load of SMC would need.

    ``np.fromstring(..., sep="\\t")`` is deprecated for text parsing in
    favour of slower alternatives (``np.array(line.split(), dtype=...)`` was
    ~3x slower in the same benchmark); the warning is suppressed deliberately
    here, not accidentally silenced project-wide.

    Returns ``(metrics, expression, coverage)``:
      - ``metrics``: index = cell barcode, columns n_counts/n_genes/pct_mito.
      - ``expression``: cells (rows) x genes_of_interest (columns) — the
        orientation ``labels.py`` expects — raw UMI counts, NOT yet
        normalised (see ``load_lee_cohort`` for the CP10K step).
      - ``coverage``: ``{"requested": [...], "found": [...], "missing":
        [...]}`` — which of ``genes_of_interest`` actually turned up as a
        row in this file. Surfaced as data, not just a log line, so a gap
        like MUC2's survives into whatever write-up reads this later.
    """
    import warnings

    wanted_genes = set(genes_of_interest)
    wanted_mito = set(mito_genes)
    found_genes: set[str] = set()
    found_mito: set[str] = set()
    collected_genes: list[str] = []
    collected_rows: list[np.ndarray] = []

    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        header_line = fh.readline().rstrip("\n")
        cell_ids = header_line.split("\t")[1:]
        n_cells = len(cell_ids)
        counts_sum = np.zeros(n_cells, dtype=np.int64)
        genes_detected = np.zeros(n_cells, dtype=np.int64)
        mito_sum = np.zeros(n_cells, dtype=np.int64)
        n_rows = 0

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=DeprecationWarning)
            for line in fh:
                tab_idx = line.index("\t")
                gene = line[:tab_idx]
                values = np.fromstring(line[tab_idx + 1 :], dtype=np.int32, sep="\t")
                if len(values) != n_cells:
                    raise ValueError(
                        f"{path.name}: row {n_rows + 1} ({gene!r}) has {len(values)} "
                        f"values, expected {n_cells} to match the header"
                    )
                n_rows += 1

                counts_sum += values
                genes_detected += values > 0

                if gene in wanted_mito:
                    found_mito.add(gene)
                    mito_sum += values
                if gene in wanted_genes:
                    found_genes.add(gene)
                    collected_genes.append(gene)
                    collected_rows.append(values)

    if n_rows == 0:
        raise ValueError(f"{path} contains no data rows")
    if not found_mito:
        raise ValueError(
            f"none of the {len(wanted_mito)} expected mitochondrial genes were found "
            f"in {path.name} — pct_mito would silently become NaN and disable that "
            f"whole QC check. This reference's gene symbols differ from what "
            f"MITO_GENES was verified against; update it deliberately, don't ignore this."
        )

    with np.errstate(invalid="ignore", divide="ignore"):
        pct_mito = np.where(counts_sum > 0, 100.0 * mito_sum / counts_sum, 0.0)

    metrics = pd.DataFrame(
        {"n_counts": counts_sum, "n_genes": genes_detected, "pct_mito": pct_mito},
        index=cell_ids,
    )

    expression = (
        pd.DataFrame(collected_rows, index=collected_genes, columns=cell_ids).T
        if collected_rows
        else pd.DataFrame(index=cell_ids)
    )

    coverage = {
        "requested": sorted(wanted_genes),
        "found": sorted(found_genes),
        "missing": sorted(wanted_genes - found_genes),
    }
    return metrics, expression, coverage


@dataclass(frozen=True)
class LeeCohort:
    """One cohort (SMC or KUL3), QC-filtered, labelled, ready to summarise."""

    study_id: str
    cells: pd.DataFrame  # patient_id, tissue, author_cell_type/subtype
    expression: pd.DataFrame  # cells x genes-of-interest, CP10K, QC-passed
    labels: pd.DataFrame  # mature__{axis}__{rung} columns, QC-passed
    axis_gene_coverage: dict[str, list[str]]
    excluded_patients: list[str] = field(default_factory=list)
    n_border_cells: int = 0


_FILES: Final[dict[str, dict[str, str]]] = {
    "smc": {
        "annotation": "GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz",
        "matrix": "GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz",
    },
    "kul3": {
        "annotation": "GSE144735_processed_KUL3_CRC_10X_annotation.txt.gz",
        "matrix": "GSE144735_processed_KUL3_CRC_10X_raw_UMI_count_matrix.txt.gz",
    },
}


def load_lee_cohort(
    which: Literal["smc", "kul3"],
    *,
    target_genes: list[str],
    axes: tuple[str, ...] = ("stem_pole", "opposite_lineage"),
    rungs: tuple[str, ...] | None = None,
    raw_dir: Path | None = None,
) -> LeeCohort:
    """Load, QC, normalise and label one real Lee cohort end to end.

    Pipeline: ``load_annotation`` -> ``stream_matrix_stats`` ->
    ``ingest.qc_flags`` (unmodified) to drop failing cells -> CP10K-normalise
    the retained expression (linear scale — ``kitagawa.decompose()``'s
    ``total = compositional + intrinsic + interaction`` identity is additive
    and needs a linear, depth-normalised scale; raw counts confound "how
    much a cell makes" with sequencing depth, and log values would break the
    identity outright) -> exclude tumor-only patients, logged -> label with
    ``labels.label_cohort`` (unmodified, leakage-guarded against
    ``target_genes``) -> assemble.

    ``target_genes`` are the panel genes under test — passed straight through
    to ``label_cohort``'s leakage guard, so a run testing MUC2/TFF3 correctly
    loses the ``opposite_lineage`` axis for that run (docs/open_decisions.md
    #1); it is not this function's job to work around that.

    ``raw_dir`` defaults to ``data/raw/lee/`` relative to the repo root.
    """
    from src.common.paths import REPO_ROOT

    raw_dir = raw_dir or (REPO_ROOT / "data" / "raw" / "lee")
    study_id = STUDY_IDS[which]
    files = _FILES[which]

    annotation = load_annotation(raw_dir / files["annotation"], study_id=study_id)

    from src.common.panel import axis_genes

    axis_marker_genes = {g for axis in axes for g in axis_genes(axis)}
    genes_of_interest = sorted(set(target_genes) | axis_marker_genes)

    metrics, raw_expression, coverage = stream_matrix_stats(
        raw_dir / files["matrix"],
        genes_of_interest=genes_of_interest,
    )
    for gene in coverage["missing"]:
        logger.warning(
            "%s: %r requested but not found in the gene index — any axis using "
            "it as a marker runs on the remaining markers only",
            study_id,
            gene,
        )

    metrics = metrics.assign(study_id=study_id)
    fail = qc_flags(metrics)
    passed_cells = metrics.index[~fail]

    cells = annotation.loc[annotation.index.intersection(passed_cells)]
    expression_counts = raw_expression.loc[raw_expression.index.intersection(passed_cells)]

    n_border_cells = int((cells["tissue"] == "border").sum())

    has_both_arms = (
        cells.loc[cells["tissue"].isin(["normal", "tumour"])]
        .groupby("patient_id")["tissue"]
        .nunique()
    )
    complete_patients = set(has_both_arms.index[has_both_arms == 2])
    all_patients = set(cells["patient_id"].unique())
    excluded_patients = sorted(all_patients - complete_patients)
    if excluded_patients:
        logger.warning(
            "excluded %d/%d %s patients lacking matched normal+tumour tissue: %s",
            len(excluded_patients),
            len(all_patients),
            study_id,
            ", ".join(excluded_patients),
        )

    keep = cells["patient_id"].isin(complete_patients)
    cells = cells.loc[keep]
    expression_counts = expression_counts.loc[expression_counts.index.intersection(cells.index)]

    # CP10K: linear, depth-normalised — see the docstring for why this scale.
    library_size = metrics.loc[expression_counts.index, "n_counts"].replace(0, pd.NA)
    expression = expression_counts.div(library_size, axis=0).mul(1e4).fillna(0.0)

    labels = label_cohort(
        expression,
        patient_id=cells.loc[expression.index, "patient_id"].tolist(),
        tissue=cells.loc[expression.index, "tissue"].tolist(),
        target_genes=target_genes,
        axes=axes,
        rungs=rungs,
    )
    labels.index = expression.index

    return LeeCohort(
        study_id=study_id,
        cells=cells,
        expression=expression,
        labels=labels,
        axis_gene_coverage=coverage,
        excluded_patients=excluded_patients,
        n_border_cells=n_border_cells,
    )


def build_gene_rung_axis_summary(
    cohort: LeeCohort,
    *,
    genes: list[str],
    rungs: tuple[str, ...] | None = None,
    axes: tuple[str, ...] = ("stem_pole", "opposite_lineage"),
) -> pd.DataFrame:
    """One row per (patient_id, study_id, gene, granularity_rung,
    labeling_axis): frac_mature_{normal,tumour}, mean_{normal,tumour} (CP10K,
    among mature cells only), n_cells_mature (tumour side) — exactly the
    input shape ``kitagawa.decompose_cohort()`` expects, unmodified.

    Border-tissue cells are excluded here (not earlier) — they are labelled
    and available in ``cohort.labels``/``cohort.cells`` for anyone who wants
    them, but Kitagawa's two-arm ``decompose()`` has no slot for a third
    tissue class.
    """
    from src.common.panel import granularity_rungs as default_rungs

    rungs = rungs or tuple(default_rungs())
    two_arm = cohort.cells["tissue"].isin(["normal", "tumour"])
    cell_meta = cohort.cells.loc[two_arm, ["patient_id", "tissue"]]

    rows: list[dict] = []
    for gene in genes:
        if gene not in cohort.expression.columns:
            continue
        gene_expr = cohort.expression.loc[cell_meta.index, gene]
        for axis in axes:
            for rung in rungs:
                col = f"mature__{axis}__{rung}"
                if col not in cohort.labels.columns:
                    continue
                mature = cohort.labels.loc[cell_meta.index, col]
                frame = pd.DataFrame(
                    {
                        "patient_id": cell_meta["patient_id"],
                        "tissue": cell_meta["tissue"],
                        "mature": mature,
                        "expression": gene_expr,
                    }
                )
                for patient_id, group in frame.groupby("patient_id"):
                    pivot = {}
                    ok = True
                    for tissue in ("normal", "tumour"):
                        sub = group.loc[group["tissue"] == tissue]
                        if len(sub) == 0:
                            ok = False
                            break
                        pivot[f"frac_mature_{tissue}"] = float(sub["mature"].mean())
                        mature_sub = sub.loc[sub["mature"]]
                        pivot[f"mean_{tissue}"] = (
                            float(mature_sub["expression"].mean()) if len(mature_sub) else 0.0
                        )
                        if tissue == "tumour":
                            pivot["n_cells_mature"] = int(mature_sub.shape[0])
                    if not ok:
                        continue
                    rows.append(
                        {
                            "patient_id": patient_id,
                            "study_id": cohort.study_id,
                            "gene": gene,
                            "granularity_rung": rung,
                            "labeling_axis": axis,
                            **pivot,
                        }
                    )
    return pd.DataFrame(rows)


def label_agreement(
    labels: pd.DataFrame,
    author_subtype: pd.Series,
    *,
    mature_subtype_allowlist: list[str],
) -> pd.DataFrame:
    """Per axis/rung agreement rate between ``classify_maturity``'s call and
    whether the author's ``Cell_subtype`` is in ``mature_subtype_allowlist``.

    Diagnostic only — never an input to ``decompose_cohort``. The Lee
    authors' subtype calls are carried through as inert metadata everywhere
    else in this module; using them to define maturity directly would be
    exactly the shortcut README design decision 2 (independent,
    structurally-different labelling axes) exists to avoid.

    ``mature_subtype_allowlist`` has no default on purpose. Which of the
    author's ~33 Cell_subtype values (e.g. "Mature Enterocytes type 1/2" vs.
    "CMS2"/"CMS3", which are transcriptional subtypes, not differentiation
    states) count as "mature" is a CRC-biology judgment call, not an
    engineering one — this function refuses to guess at it.
    """
    if not mature_subtype_allowlist:
        raise ValueError(
            "mature_subtype_allowlist must be a deliberately chosen, non-empty list — "
            "see this function's docstring for why there is no default"
        )
    author_mature = author_subtype.isin(mature_subtype_allowlist)
    mature_cols = [c for c in labels.columns if c.startswith("mature__")]
    rows = []
    for col in mature_cols:
        _, axis, rung = col.split("__")
        agree = (labels[col] == author_mature.reindex(labels.index)).mean()
        rows.append({"labeling_axis": axis, "granularity_rung": rung, "agreement": float(agree)})
    return pd.DataFrame(rows)
