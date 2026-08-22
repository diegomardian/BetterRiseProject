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
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Final, Literal

import numpy as np
import pandas as pd

from src.estimator.ingest import (
    COMPARTMENT_COLUMN,
    RETENTION_GAP_WARN_PTS,
    differential_retention,
    qc_flags,
)
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

#: The compartment the compositional arm is built from. Matches the value the
#: Lee deposits use in ``Cell_type``.
EPITHELIAL_COMPARTMENT: Final[str] = "Epithelial cells"

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
    #: Compartment the maturity labels were computed within, or None if every
    #: cell was labelled. Cells outside it carry pd.NA, not False — see
    #: ``load_lee_cohort``'s ``label_compartment``.
    label_compartment: str | None = EPITHELIAL_COMPARTMENT
    #: The same cells and genes as ``expression``, but **raw integer counts** —
    #: what came off the matrix before the CP10K step.
    #:
    #: ``expression`` is the right scale for ``decompose_cohort``: the Kitagawa
    #: identity is additive and needs a linear, depth-normalised value. It is
    #: the wrong input for W2's pseudobulk generator, which realises a
    #: multiplicative shift by binomial thinning and Poisson augmentation —
    #: both defined on counts. Handing it CP10K used to truncate every value
    #: below 1.0 to zero silently, destroying exactly the low-expressing cells
    #: the near-zero mature-cell edge cases are made of; the generator now
    #: refuses a non-integer matrix outright, which is safe but leaves the
    #: harness with nothing to run on.
    #:
    #: Empty by default. Pass ``keep_raw_counts=True`` to populate it.
    #:
    #: NOTE for callers building a reference matrix: this frame contains the
    #: target genes, because the generator needs them to apply the shift.
    #: Exclude them before deconvolution (CLAUDE.md invariant 2) — W2's
    #: ``bulk_recovery.reference_profiles(..., exclude_genes=[target])`` does.
    raw_counts: pd.DataFrame = field(default_factory=pd.DataFrame)


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
    extra_genes: Sequence[str] = (),
    keep_raw_counts: bool = False,
    label_compartment: str | None = EPITHELIAL_COMPARTMENT,
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

    ``extra_genes`` widens the set collected off the matrix beyond
    ``target_genes | axis markers``. Deconvolution needs 500-2000 genes for
    nu-SVR to be robust (execution_plan.md §2.1 error #4), and the default set
    is roughly twenty — so W2's harness passes a marker panel here rather than
    re-parsing the matrix itself. Streaming cost is unchanged; only the number
    of rows retained grows.

    ``keep_raw_counts`` additionally populates ``LeeCohort.raw_counts`` with
    the pre-CP10K integer counts. See the field's docstring for why the
    generator cannot use the normalised frame. Off by default so existing
    callers pay nothing.

    ``label_compartment`` restricts *labelling* to one compartment, epithelium
    by default, and is not a cosmetic filter. ``classify_maturity`` thresholds
    at a quantile of the axis score over whatever it is handed, and on Lee the
    non-epithelial cells are the majority (56k of 63k on SMC). They carry no
    LGR5/OLFM4 at all, so on the inverted ``stem_pole`` axis they score as
    maximally mature and drag the quantile into the immune mass — leaving
    essentially no epithelial cell above it. The mature *fraction* is the
    compositional term, so that is not a labelling detail, it is the estimate.
    Cells outside the compartment keep their row with ``pd.NA`` in every
    maturity column, so they stay available and cannot be silently counted as
    immature. Pass None to label everything and take responsibility for the
    threshold. docs/open_decisions.md #13, "Non-epithelial cells: caller must
    pre-filter" — this is W4 doing that rather than leaving it to the caller.

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
    genes_of_interest = sorted(set(target_genes) | axis_marker_genes | set(extra_genes))

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

    # QC needs the compartment, the patient and the arm, not just the three
    # metrics: MAD bounds are computed within (study, compartment), and the
    # retention check below is per patient per arm. docs/open_decisions.md #12.
    annotated_metrics = metrics.join(
        annotation.loc[:, ["patient_id", "tissue", "author_cell_type"]], how="inner"
    ).rename(columns={"author_cell_type": COMPARTMENT_COLUMN})
    n_unannotated = len(metrics) - len(annotated_metrics)
    if n_unannotated:
        logger.warning(
            "%s: %d/%d cells in the matrix have no annotation row and are dropped "
            "before QC — they have no compartment to be an outlier within, and no "
            "arm to belong to",
            study_id,
            n_unannotated,
            len(metrics),
        )
    metrics = annotated_metrics.assign(study_id=study_id)
    fail = qc_flags(metrics)
    passed_cells = metrics.index[~fail]

    # The check decision #12 requires before any compositional number is
    # believed. Reported, not enforced: a flagged patient is a fact about the
    # deposit that the gate memo needs to carry, not a reason to fail the load.
    retention = differential_retention(metrics, ~fail)
    flagged = retention.loc[retention["flagged"]]
    if len(flagged):
        logger.warning(
            "%s: %d/%d patients have a tumour-vs-normal epithelial QC retention "
            "gap over %.0f pts (median %.1f pts). The compositional term is a "
            "within-patient difference, so this moves it directly — carry it as a "
            "limitation. docs/open_decisions.md #12.\n%s",
            study_id,
            len(flagged),
            len(retention),
            RETENTION_GAP_WARN_PTS,
            retention["gap_pts"].median(),
            flagged.to_string(index=False),
        )

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

    if label_compartment is None:
        labelled_index = expression.index
    else:
        compartment = cells.loc[expression.index, "author_cell_type"]
        labelled_index = expression.index[compartment == label_compartment]
        if len(labelled_index) == 0:
            raise ValueError(
                f"{study_id}: no cells in compartment {label_compartment!r} survived "
                f"QC — the quantile that defines maturity would have nothing to be "
                f"computed over. Observed compartments: "
                f"{sorted(cells['author_cell_type'].unique())}"
            )
        logger.info(
            "%s: labelling within %r (%d/%d QC-passed cells); the rest keep a row "
            "with pd.NA maturity, never False",
            study_id,
            label_compartment,
            len(labelled_index),
            len(expression),
        )

    labels = label_cohort(
        expression.loc[labelled_index],
        patient_id=cells.loc[labelled_index, "patient_id"].tolist(),
        tissue=cells.loc[labelled_index, "tissue"].tolist(),
        target_genes=target_genes,
        axes=axes,
        rungs=rungs,
    )
    labels.index = labelled_index

    if len(labelled_index) != len(expression):
        # Widen back to every QC-passed cell. Maturity columns become nullable
        # boolean so an unlabelled cell reads as pd.NA — "not asked" — and not
        # as False, which would be counted as an immature epithelial cell and
        # would move the compositional term. Invariant 1's shape, one level down.
        maturity_columns = [c for c in labels.columns if c.startswith("mature__")]
        widened = labels.reindex(expression.index)
        widened["patient_id"] = cells.loc[expression.index, "patient_id"].to_numpy()
        widened["tissue"] = cells.loc[expression.index, "tissue"].to_numpy()
        for col in maturity_columns:
            widened[col] = widened[col].astype("boolean")
        labels = widened

    return LeeCohort(
        study_id=study_id,
        cells=cells,
        expression=expression,
        labels=labels,
        axis_gene_coverage=coverage,
        excluded_patients=excluded_patients,
        n_border_cells=n_border_cells,
        label_compartment=label_compartment,
        # Same rows and columns as `expression`, one step earlier in the
        # pipeline. Cast to a nullable integer dtype so a caller can tell at a
        # glance that these are counts and not a normalised scale.
        raw_counts=(
            expression_counts.astype("Int64")
            if keep_raw_counts
            else pd.DataFrame(index=expression.index)
        ),
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

    Cells with ``pd.NA`` maturity — everything outside
    ``cohort.label_compartment`` — are dropped per (axis, rung) rather than
    read as immature. Counting them as immature would put the tumour's immune
    infiltrate into the denominator of the mature fraction, which is the
    compositional term; Δ(mature fraction) would then be measuring immune
    content. A patient with no labelled cell in one arm is skipped, the same
    way a patient missing an arm entirely is.
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
                ).loc[mature.notna().to_numpy()]
                frame["mature"] = frame["mature"].astype(bool)
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
