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
#:
#: **50, not the conventional 20.** Measured on the pilot (results/, 2026-08-17):
#: colonic epithelium runs at a median of 29.8% in normal tissue and 21.1% in
#: tumour, against 4-11% for every immune and stromal compartment. A 20% cap
#: therefore sits *below the median for normal epithelium*: it removes 59% of
#: all epithelial cells and opens a 22.7-point tumour/normal retention gap, which
#: lands directly on Delta(mature fraction) and inflates apparent compositional
#: loss — a bias pointing at the prior hypothesis.
#:
#: Two further facts settle the value. GSE178341 is **already filtered at 50%**
#: upstream (observed max 49.976 normal / 49.988 tumour), so a 50 cap is a no-op
#: rather than an opinion, and double-filtering would only compound the first
#: cut. And the mitochondrial content is genuine, not ambient: contamination is
#: ~2.7% of counts and MT genes are ~18% of the soup, so ambient contributes
#: about 0.5% of a cell's counts — nowhere near a 30% observed fraction.
#:
#: W4 uses 20 on the Lee cohorts. **They must match or justify diverging**, or
#: the two cohorts are not comparable at the gate — see docs/open_decisions.md #12.
DEFAULT_MAX_PCT_MITO: Final[float] = 50.0

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


#: Scrublet's expected doublet rate. 10x loads ~0.8% per 1,000 cells recovered,
#: and these samples run 1k-16k, so a single global value is already an
#: approximation. It is a PRIOR, not a target: Scrublet uses it to place the
#: simulated-doublet distribution, and the called rate can land either side.
DEFAULT_EXPECTED_DOUBLET_RATE: Final[float] = 0.06

#: Below this many cells, the simulated-doublet distribution is too sparse to
#: threshold and Scrublet's automatic cutoff becomes arbitrary. Such samples are
#: reported with a null rate rather than a made-up one.
MIN_CELLS_FOR_DOUBLET_CALL: Final[int] = 100


def flag_doublets(
    counts: Any,
    sample_id: Any,
    *,
    expected_rate: float = DEFAULT_EXPECTED_DOUBLET_RATE,
    min_cells: int = MIN_CELLS_FOR_DOUBLET_CALL,
    threshold: float | None = None,
    seed: int = 20260101,
) -> pd.DataFrame:
    """Scrublet per sample. W1, week 1.

    **Per sample, never across samples.** A doublet is two cells captured in one
    droplet, which is a within-run event. Scrublet works by simulating doublets
    from observed transcriptomes and asking which real cells look like them, so
    pooling samples simulates chimeras that no droplet could have produced and
    the resulting distribution means nothing.

    Returns one row per cell — `sample_id`, `doublet_score`, `predicted_doublet`,
    `threshold`, `callable` — in input order, so it can be assigned straight onto
    an obs frame. Samples below `min_cells` get a score of NaN and
    ``callable=False`` rather than a fabricated call: too few cells makes the
    simulated distribution too sparse to threshold, and a made-up boundary there
    would silently delete real cells.

    On `expected_rate`: this is Scrublet's *prior*, not a quota. It positions the
    simulated distribution; the called rate comes out of the data and may differ
    substantially. Do not tune it until the reported rate looks like what you
    wanted — that is fitting the QC to the hypothesis.

    **Expect a low rate on this deposit.** Pelka et al. filtered doublets before
    deposition, and GEO ran dropletUtils upstream. A per-sample rate of 1-2% is
    the likely and perfectly good answer: the week-1 deliverable is the rate
    being *documented*, and "already handled upstream, here is the evidence" is a
    result. A high rate would be the surprising outcome and would want explaining
    before anything is deleted.

    ``threshold`` overrides Scrublet's automatic call for every sample. Use it
    only after looking at the score histograms — the automatic threshold is a
    minimum between two modes, and it degrades quietly when the distribution is
    unimodal, which is exactly what pre-filtered data produces.
    """
    samples = np.asarray([str(s) for s in sample_id], dtype=object)
    if samples.shape[0] != counts.shape[0]:
        raise QCError(
            f"sample_id has {samples.shape[0]} entries for {counts.shape[0]} cells"
        )

    score = np.full(samples.shape[0], np.nan, dtype=float)
    called = np.zeros(samples.shape[0], dtype=bool)
    cutoff = np.full(samples.shape[0], np.nan, dtype=float)
    usable = np.zeros(samples.shape[0], dtype=bool)

    for sample in pd.unique(samples):
        rows = np.flatnonzero(samples == sample)
        if rows.size < min_cells:
            print(
                f"note: {sample} has {rows.size} cells, below {min_cells} — "
                f"no doublet call. Reported as not callable, not as zero doublets."
            )
            continue
        # Imported here, not at the top: every guard above runs without the
        # dependency, so shape errors surface in any environment and the
        # pure-pandas reporting below stays testable where scrublet is absent.
        # scrublet is heavy (it compiles annoy) and single-workstream, so it is
        # pinned in env/w1_reference.yml rather than the CI dev extra — same
        # treatment as scanpy.
        try:
            import scrublet
        except ImportError as exc:
            raise QCError(
                "scrublet is not installed. It is pinned in "
                "env/w1_reference.yml (scrublet=0.2.3). If it will not import "
                "against this numpy, use scDblFinder instead — "
                "bioconductor-scdblfinder is pinned in the same file and the "
                "plan sanctions either."
            ) from exc

        subset = counts[rows]
        detector = scrublet.Scrublet(
            subset, expected_doublet_rate=expected_rate, random_state=seed
        )
        try:
            scores, predicted = detector.scrub_doublets(verbose=False)
        except Exception as exc:  # noqa: BLE001 - Scrublet raises broadly
            print(f"note: {sample} — Scrublet failed ({exc}); not callable.")
            continue
        if scores is None:
            print(f"note: {sample} — Scrublet returned no scores; not callable.")
            continue

        score[rows] = np.asarray(scores, dtype=float)
        usable[rows] = True
        if threshold is not None:
            cutoff[rows] = threshold
            called[rows] = score[rows] >= threshold
        else:
            automatic = getattr(detector, "threshold_", np.nan)
            cutoff[rows] = automatic
            # Scrublet returns None for `predicted` when it cannot find a
            # minimum between two modes — common on pre-filtered data, where the
            # distribution is unimodal. Keep the scores, refuse the call.
            if predicted is None or not np.isfinite(automatic):
                usable[rows] = False
                print(
                    f"note: {sample} — no automatic threshold (the score "
                    f"distribution is probably unimodal, which is what "
                    f"pre-filtered data looks like). Scores kept, call withheld."
                )
            else:
                called[rows] = np.asarray(predicted, dtype=bool)

    return pd.DataFrame(
        {
            "sample_id": samples,
            "doublet_score": score,
            "predicted_doublet": called,
            "threshold": cutoff,
            "callable": usable,
        }
    )


def doublet_rate_by_sample(calls: pd.DataFrame) -> pd.DataFrame:
    """Per-sample doublet rate, for the QC deliverable.

    The week-1 requirement is "QC thresholds documented with rationale", and a
    doublet threshold is one of those. Commit this next to
    :func:`qc_thresholds`.
    """
    required = {"sample_id", "predicted_doublet", "callable", "threshold"}
    missing = required - set(calls.columns)
    if missing:
        raise QCError(f"calls is missing column(s): {sorted(missing)}")

    out = (
        calls.groupby("sample_id", observed=True)
        .agg(
            n_cells=("predicted_doublet", "size"),
            n_doublets=("predicted_doublet", "sum"),
            threshold=("threshold", "first"),
            callable=("callable", "any"),
        )
        .reset_index()
    )
    out["doublet_rate"] = np.where(
        out["callable"], out["n_doublets"] / out["n_cells"], np.nan
    )
    return out


def doublet_compartment_enrichment(
    calls: pd.DataFrame, compartment: Any
) -> pd.DataFrame:
    """Are called doublets concentrated in one compartment? **The check that
    decides whether doublets can distort the maturity call.**

    A doublet's transcriptome is the sum of two cells, so an epithelial-immune
    doublet carries both compartments' markers. Two consequences, and only the
    second is about counts:

    1. It can be assigned to the wrong compartment, moving the denominator of
       the mature fraction.
    2. **It scores as more mature than either cell was.** The maturity call on
       axis 1 is a detection gate — "no stem marker detected" — and a doublet has
       roughly twice the counts, so it clears detection more easily and lands on
       the immature side. Doublets therefore push the mature fraction *down*, and
       they are not distributed evenly between tumour and normal.

    A rate that is uniform across compartments is reassuring. A rate concentrated
    in epithelium is not, whatever its absolute size.
    """
    if "predicted_doublet" not in calls.columns:
        raise QCError("calls needs a predicted_doublet column")
    values = np.asarray([str(c) for c in compartment], dtype=object)
    if values.shape[0] != len(calls):
        raise QCError(
            f"compartment has {values.shape[0]} entries for {len(calls)} cells"
        )
    frame = calls.loc[:, ["predicted_doublet", "callable"]].copy()
    frame["compartment"] = values
    frame = frame[frame["callable"]]
    if frame.empty:
        raise QCError("no callable cells; every sample was below min_cells")

    out = (
        frame.groupby("compartment", observed=True)
        .agg(n_cells=("predicted_doublet", "size"),
             n_doublets=("predicted_doublet", "sum"))
        .reset_index()
    )
    out["doublet_rate"] = out["n_doublets"] / out["n_cells"]
    overall = float(frame["predicted_doublet"].mean())
    out["overall_rate"] = overall
    # Ratio rather than difference: a 2% rate against a 1% baseline matters as
    # much as 20% against 10%, and the absolute gap hides the first.
    out["enrichment"] = out["doublet_rate"] / overall if overall else np.nan
    out["flagged"] = out["enrichment"] > 1.5
    return out.sort_values("enrichment", ascending=False).reset_index(drop=True)
