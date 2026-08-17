"""Per-cell QC and doublet flagging for GSE178341. W1, week 1.

execution_plan.md §4 asks for two things from this step, and the second is the
one usually skipped: *"Cell counts by patient and tissue tabulated; QC
thresholds documented with rationale."* So the deliverable here is not a
filtered matrix — it is :func:`qc_thresholds`, a table saying what was cut,
where, by what rule, and why. :func:`apply_qc` then just applies it.

**Thresholds are per batch, never global.** GSE178341 spans patients and
sequencing runs with genuinely different depth distributions; a single global
cutoff over- or under-filters somebody. The unit here is the sample/batch, which
for this dataset is effectively (patient, tissue).

Conventions deliberately match W4's Lee-cohort QC in
``src/estimator/ingest.py`` — 5 MADs on the modified z-score, a 20% hard
mitochondrial cap. execution_plan.md §4 tells W4 to mirror W1's pipeline shape
and coordinate rather than share code prematurely; the reason to keep the
numbers identical is that GSE178341 and the Lee cohorts have to be comparable at
the week-5 gate, and silently diverging thresholds is a difference that would
show up as biology.

``_mad_outlier`` below is a second implementation of the same five-line textbook
statistic W4 has. That duplication is deliberate and follows the precedent
``src/common/io.py`` set for its directory-naming logic: reaching into another
workstream's private helper is what CONTRIBUTING §2 forbids, and promoting it to
``src/common/`` is only worth a ``shared/`` PR once a third caller appears.

Metric computation works on a raw count matrix rather than AnnData so it is
testable without the single-cell stack, matching ``src/reference/ingest.py``.
Matrix orientation is **cells x genes**.
"""

from __future__ import annotations

from typing import Any, Final

import numpy as np
import pandas as pd

#: Modified z-score cutoff. 5 MADs is the sc-best-practices convention and is
#: what W4 uses on Lee — keep them equal.
DEFAULT_N_MADS: Final[float] = 5.0

#: Hard cap, not per-batch. A dying cell's mitochondrial fraction is not
#: expected to be protocol-relative the way library size is.
DEFAULT_MAX_PCT_MITO: Final[float] = 20.0

#: Absolute floors, applied on top of the per-batch MAD rule. A batch that is
#: uniformly poor has a low median, and the MAD rule alone would happily keep
#: 40-gene barcodes because they are typical *for that batch*.
ABSOLUTE_MIN_GENES: Final[int] = 200
ABSOLUTE_MIN_COUNTS: Final[int] = 500

#: Human mitochondrial gene symbols. Ribosomal is reported but never filtered on
#: — high ribosomal content is a real biological state in proliferating cells,
#: which is exactly the stem-pole population the labelling axes care about.
MITO_PREFIXES: Final[tuple[str, ...]] = ("MT-", "mt-")
RIBO_PREFIXES: Final[tuple[str, ...]] = ("RPS", "RPL")

METRIC_COLUMNS: Final[tuple[str, ...]] = ("n_counts", "n_genes", "pct_mito", "pct_ribo")


class QCError(ValueError):
    """Malformed QC input — a missing column, or a metrics/label length mismatch."""


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------


def _row_sums(matrix: Any) -> np.ndarray:
    return np.asarray(matrix.sum(axis=1)).ravel()


def _genes_detected(matrix: Any) -> np.ndarray:
    """Non-zero genes per cell, for dense or sparse input."""
    if hasattr(matrix, "getnnz"):
        return matrix.getnnz(axis=1).astype(np.int64)
    return (np.asarray(matrix) > 0).sum(axis=1).astype(np.int64)


def _prefix_fraction(matrix: Any, gene_names: Any, prefixes: tuple[str, ...]) -> np.ndarray:
    """Percentage of each cell's counts falling in genes matching `prefixes`."""
    names = np.asarray(list(gene_names), dtype=object)
    mask = np.zeros(names.shape, dtype=bool)
    for prefix in prefixes:
        mask |= np.array([str(n).startswith(prefix) for n in names], dtype=bool)

    totals = _row_sums(matrix)
    if not mask.any():
        return np.zeros_like(totals, dtype=float)

    subset = matrix[:, mask]
    selected = _row_sums(subset)
    with np.errstate(divide="ignore", invalid="ignore"):
        pct = np.where(totals > 0, selected / totals * 100.0, 0.0)
    return np.asarray(pct, dtype=float)


def cell_qc_metrics(
    matrix: Any,
    gene_names: Any,
    *,
    batch: Any,
    patient_id: Any | None = None,
    tissue: Any | None = None,
) -> pd.DataFrame:
    """Per-cell QC metrics. One row per cell, in input order.

    `batch` is the unit thresholds are computed within — for GSE178341 that is
    the sample, i.e. effectively (patient, tissue). Pass `patient_id` and
    `tissue` too when you have them so the summary tables come out directly.
    """
    n_cells = matrix.shape[0]
    if len(list(gene_names)) != matrix.shape[1]:
        raise QCError(
            f"gene_names has {len(list(gene_names))} entries for "
            f"{matrix.shape[1]} columns. Matrix must be cells x genes."
        )

    frame = pd.DataFrame(
        {
            "batch": pd.Series(list(batch), dtype="object"),
            "n_counts": _row_sums(matrix).astype(np.int64),
            "n_genes": _genes_detected(matrix),
            "pct_mito": _prefix_fraction(matrix, gene_names, MITO_PREFIXES),
            "pct_ribo": _prefix_fraction(matrix, gene_names, RIBO_PREFIXES),
        }
    )
    if len(frame) != n_cells:
        raise QCError(f"batch has {len(frame)} entries for {n_cells} cells")
    if patient_id is not None:
        frame.insert(0, "patient_id", list(patient_id))
    if tissue is not None:
        frame.insert(1, "tissue", list(tissue))
    return frame


# ---------------------------------------------------------------------------
# Thresholds — the actual week-1 deliverable
# ---------------------------------------------------------------------------


def _mad_outlier_bounds(values: pd.Series, n_mads: float) -> tuple[float, float]:
    """Modified z-score bounds (Iglewicz & Hoaglin 1993).

    Robust to the heavy right tail that library-size and gene-count
    distributions have; a mean/SD cutoff gets dragged around by exactly the
    cells the filter exists to catch. A zero MAD (a degenerate batch) yields
    infinite bounds — no cell is cut on a rule that cannot discriminate.
    """
    median = float(values.median())
    mad = float((values - median).abs().median())
    if mad == 0:
        return (-np.inf, np.inf)
    spread = n_mads * mad / 0.6745
    return (median - spread, median + spread)


def qc_thresholds(
    metrics: pd.DataFrame,
    *,
    n_mads: float = DEFAULT_N_MADS,
    max_pct_mito: float = DEFAULT_MAX_PCT_MITO,
    min_genes: int = ABSOLUTE_MIN_GENES,
    min_counts: int = ABSOLUTE_MIN_COUNTS,
) -> pd.DataFrame:
    """Per-batch thresholds with a written rationale. **This is the deliverable.**

    Returns one row per (batch, metric) with the bounds actually applied, the
    rule that produced them, why that rule, and how many cells it costs. Commit
    it next to the run — "QC thresholds documented with rationale" means this
    table exists, not that someone remembers what they typed.
    """
    missing = {"batch", *METRIC_COLUMNS} - set(metrics.columns)
    if missing:
        raise QCError(f"metrics is missing column(s): {sorted(missing)}")

    rows: list[dict[str, Any]] = []
    for batch, group in metrics.groupby("batch", observed=True, sort=True):
        for metric, floor in (("n_counts", min_counts), ("n_genes", min_genes)):
            mad_low, mad_high = _mad_outlier_bounds(group[metric], n_mads)
            lower = max(mad_low, float(floor))
            failed = (group[metric] < lower) | (group[metric] > mad_high)
            rows.append(
                {
                    "batch": batch,
                    "metric": metric,
                    "lower": lower,
                    "upper": mad_high,
                    "method": (
                        f"max({n_mads:g}-MAD modified z within batch, "
                        f"absolute floor {floor})"
                    ),
                    "rationale": (
                        "Per-batch because depth distributions differ by sample; a global "
                        "cutoff over-filters shallow batches. MAD rather than SD because the "
                        "tail is heavy. Absolute floor because a uniformly poor batch has a "
                        "low median and the MAD rule alone would keep near-empty barcodes."
                    ),
                    "n_cells": int(len(group)),
                    "n_failed": int(failed.sum()),
                }
            )

        mito_failed = group["pct_mito"] > max_pct_mito
        rows.append(
            {
                "batch": batch,
                "metric": "pct_mito",
                "lower": -np.inf,
                "upper": float(max_pct_mito),
                "method": f"hard cap at {max_pct_mito:g}%",
                "rationale": (
                    "Hard, not per-batch: a dying cell's mitochondrial fraction is not "
                    "protocol-relative the way library size is. Matches W4's Lee-cohort cap "
                    "so the two cohorts stay comparable at the gate."
                ),
                "n_cells": int(len(group)),
                "n_failed": int(mito_failed.sum()),
            }
        )
    return pd.DataFrame(rows)


def apply_qc(metrics: pd.DataFrame, thresholds: pd.DataFrame) -> pd.Series:
    """Boolean pass/fail per cell, from a thresholds table. Index matches `metrics`.

    Deliberately dumb: every judgement call lives in :func:`qc_thresholds`, so
    what actually got applied is auditable from the committed table rather than
    reconstructed from code.
    """
    missing = {"batch", *METRIC_COLUMNS} - set(metrics.columns)
    if missing:
        raise QCError(f"metrics is missing column(s): {sorted(missing)}")

    passes = pd.Series(True, index=metrics.index)
    for _, rule in thresholds.iterrows():
        in_batch = metrics["batch"] == rule["batch"]
        if not in_batch.any():
            continue
        values = metrics.loc[in_batch, rule["metric"]]
        ok = (values >= rule["lower"]) & (values <= rule["upper"])
        passes.loc[in_batch] &= ok
    return passes


def qc_summary(metrics: pd.DataFrame, passes: pd.Series) -> pd.DataFrame:
    """Cells before and after QC, per patient x tissue. The week-1 tabulation.

    Read the `retained` column before moving on. A patient who loses most of
    their normal cells here has quietly lost their compositional reference, and
    that is a week-1 finding, not a week-4 one.
    """
    if len(passes) != len(metrics):
        raise QCError(f"passes has {len(passes)} entries for {len(metrics)} cells")

    keys = [k for k in ("patient_id", "tissue") if k in metrics.columns] or ["batch"]
    frame = metrics.loc[:, keys].copy()
    frame["passed"] = np.asarray(passes, dtype=bool)

    summary = frame.groupby(keys, observed=True).agg(
        n_cells=("passed", "size"), n_passed=("passed", "sum")
    )
    summary["n_failed"] = summary["n_cells"] - summary["n_passed"]
    summary["retained"] = summary["n_passed"] / summary["n_cells"]
    return summary.reset_index()


# ---------------------------------------------------------------------------
# Doublets
# ---------------------------------------------------------------------------


def differential_retention(
    metrics: pd.DataFrame, passes: pd.Series, *, warn_at: float = 0.10
) -> pd.DataFrame:
    """Per-patient QC retention in tumour vs normal. **Read this before filtering.**

    QC is not neutral in this project. The compositional term is
    Delta(mature fraction) between a patient's tumour and their own normal, so if
    QC removes cells at different rates in the two arms, it moves that difference
    directly — and mature colonocytes are exactly the fragile, high-mitochondrial
    population a mitochondrial cap removes first.

    A patient whose normal loses 20 points more than their tumour has had their
    normal mature fraction understated, which inflates the apparent compositional
    loss. That is a bias in the direction of the prior hypothesis, which is the
    worst kind (README, "the bias points the wrong way").

    Returns one row per patient with both retentions, their difference, and a
    ``flagged`` column for gaps beyond `warn_at`.
    """
    if len(passes) != len(metrics):
        raise QCError(f"passes has {len(passes)} entries for {len(metrics)} cells")
    if "tissue" not in metrics.columns or "patient_id" not in metrics.columns:
        raise QCError("metrics needs patient_id and tissue columns")

    frame = metrics.loc[:, ["patient_id", "tissue"]].copy()
    frame["passed"] = np.asarray(passes, dtype=bool)
    wide = (
        frame.groupby(["patient_id", "tissue"], observed=True)["passed"]
        .mean()
        .unstack("tissue")
    )
    for column in ("tumour", "normal"):
        if column not in wide.columns:
            wide[column] = np.nan

    out = wide.loc[:, ["tumour", "normal"]].copy()
    out.columns = ["retained_tumour", "retained_normal"]
    out["difference"] = out["retained_tumour"] - out["retained_normal"]
    out["flagged"] = out["difference"].abs() > warn_at
    return out.reset_index()


def flag_doublets(*args: Any, **kwargs: Any) -> Any:
    """scDblFinder against real GSE178341 count matrices. W1, week 1.

    Run **per sample, never across samples** — a doublet is two cells in one
    droplet, which is a within-run event, and pooling makes the simulated
    doublet distribution meaningless.

    Unimplemented on purpose, the same pattern as ``_select_markers`` in
    ``signature.py`` and ``flag_doublets`` in ``src/estimator/ingest.py``:
    calling scDblFinder needs the real matrices in hand, and the cutoff is a
    judgement call over real score distributions rather than a formula.
    ``bioconductor-scdblfinder`` is pinned in ``env/w1_reference.yml``.

    When you implement it: add the per-sample doublet rate to the
    :func:`qc_thresholds` table so it is documented alongside everything else.
    """
    raise NotImplementedError(
        "W1 — doublet detection needs the real GSE178341 count matrices. "
        "Run scDblFinder per sample; see src/reference/README.md, week 1."
    )
