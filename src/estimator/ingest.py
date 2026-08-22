"""Lee cohort ingest: per-study QC, doublet flagging, ambient correction. W4.

Per-study thresholds, not global (execution_plan.md W1 task table — W4 mirrors
the same pipeline SHAPE for Lee, coordinating with W1 rather than sharing code
prematurely, src/estimator/README.md). MAD-based outlier flags are computed
independently within each ``study_id``, since GSE132465 (SMC) and GSE144735
(KUL3) come from different protocols — a single global cutoff would over- or
under-filter one of them.

...AND WITHIN EACH COMPARTMENT. Measured 2026-08-22, docs/open_decisions.md #12
------------------------------------------------------------------------------
Per study was not enough, and the shortfall was not cosmetic. A MAD outlier rule
assumes one unimodal population; a whole-cohort matrix is six of them, and on
SMC the epithelial compartment runs **3.9x deeper** than the median immune one
(18,724 vs 4,861 median UMIs). Pooled across compartments, the upper MAD bound
therefore fires on epithelial cells for being epithelial:

===============================  ==================  ==================
epithelium retained (SMC)        MAD pooled          MAD per compartment
===============================  ==================  ==================
normal arm                       88.5%               100%
tumour arm                       **62.0%**           100%
median per-patient gap           **-29.6 pts**       0.0 pts
patients with a >10 pt gap       **9 of 10**         0
===============================  ==================  ==================

The gap points the wrong way. Mature colonocytes are the deepest epithelial
cells there are, so cutting the tumour arm's deepest cells understates the
tumour mature fraction, which **inflates the apparent compositional loss** —
a bias toward the prior hypothesis, produced entirely by a QC parameter.
That is the failure mode W1 raised for the mitochondrial cap on GSE178341,
reaching Lee through a different filter.

Grouping is by (study, compartment) and deliberately NOT by tissue. Grouping by
tissue would equalise retention by construction and hide the artifact rather
than remove it — the same trap as the per-sample quantile bug in
docs/open_decisions.md #13. Compartment labels come from the deposit's own
annotation, which is mildly circular (they were clustered post-QC); it is the
standard per-cell-type QC recommendation and it is visible here rather than
assumed.

THE MITOCHONDRIAL CAP IS A NO-OP ON LEE
---------------------------------------
Decision #12 asked W4 to set ``max_pct_mito`` from Lee's own per-compartment
distribution rather than inherit W1's 50.0 or keep 20.0 by convention. Measured:
both deposits are **already filtered at 20%** upstream — observed maxima 19.995
(SMC) and 19.994 (KUL3) — so the cap fails **zero** cells at any value >= 20 and
GEO does not ship the droplets that would let W4 revisit it. Lee's epithelial
medians (8.6% SMC, 9.7% KUL3) are not comparable to GSE178341's 29.8% for the
same reason: they are conditioned on the authors' cut. The default stays 20.0
because that is what the data already had done to it; the number is inherited,
not chosen, and no value W4 picks changes a single cell. See
docs/open_decisions.md #12 and #11 — the same shape of finding as GSE178341
shipping no unfiltered droplets.

Doublet scoring (scDblFinder or Scrublet) and ambient correction (SoupX) are
external-tool calls against real count matrices — not meaningfully fakeable
without the actual Lee files in hand, and getting the thresholds right is a
judgment call, not a formula. Their entry points are stubbed here, the same
pattern as ``src/reference/signature.py``'s ``_select_markers``: the scaffold
owns the guard rails (per-study QC, documenting the rationale), the judgment
call over real data is yours once GSE132465/GSE144735 are downloaded.
"""

from __future__ import annotations

from typing import Final

import numpy as np
import pandas as pd

DEFAULT_N_MADS: Final[float] = 5.0

#: Inherited from the Lee deposits, which are already cut at 20% upstream — not
#: chosen, and a no-op at any value >= 20. See the module docstring and
#: docs/open_decisions.md #12. W1 uses 50.0 on GSE178341 for a measured reason;
#: the two cohorts differ because the deposits differ, not because one number
#: went unchecked.
DEFAULT_MAX_PCT_MITO: Final[float] = 20.0

#: Column carrying the coarse cell compartment MAD bounds are computed within.
#: On Lee this is the deposit's ``Cell_type`` (Epithelial cells, T cells, ...).
COMPARTMENT_COLUMN: Final[str] = "compartment"

_REQUIRED_METRIC_COLUMNS = {
    "study_id",
    "n_genes",
    "n_counts",
    "pct_mito",
    COMPARTMENT_COLUMN,
}


def _mad_outlier(values: pd.Series, n_mads: float) -> pd.Series:
    """Modified z-score outlier flag (Iglewicz & Hoaglin 1993), robust to the
    heavy right tail library-size and gene-count distributions usually have —
    a plain mean/SD cutoff would be dragged around by exactly the cells the
    filter exists to catch."""
    median = values.median()
    mad = (values - median).abs().median()
    if mad == 0:
        return pd.Series(False, index=values.index)
    modified_z = 0.6745 * (values - median) / mad
    return modified_z.abs() > n_mads


def qc_flags(
    metrics: pd.DataFrame,
    *,
    n_mads: float = DEFAULT_N_MADS,
    max_pct_mito: float = DEFAULT_MAX_PCT_MITO,
) -> pd.Series:
    """Per-cell QC-fail flag. ``metrics`` needs columns ``study_id``,
    ``n_genes``, ``n_counts``, ``pct_mito`` and ``compartment`` — one row per
    cell.

    A cell fails if ``n_genes`` or ``n_counts`` is a MAD-outlier *within its own
    study and compartment* (``n_mads``, default 5 — the sc-best-practices
    convention), or ``pct_mito`` exceeds ``max_pct_mito``. The mitochondrial cap
    is a single hard threshold, not per-study — a dying cell's mitochondrial
    fraction isn't expected to be protocol-relative the way library size is.

    ``compartment`` is required rather than optional. Defaulting it away would
    restore the pooled grouping that cut SMC's tumour epithelium 29.6 points
    harder than its normal arm, in the direction of the prior hypothesis, with
    nothing in the output to show for it — the module docstring has the
    measurement. A caller without compartment labels genuinely cannot run this
    check safely and should say so, not pass a constant.

    execution_plan.md wants "QC thresholds documented with rationale," not just
    applied — log the per-study fail counts this produces, and run
    :func:`differential_retention` before believing any compositional number.
    """
    missing = _REQUIRED_METRIC_COLUMNS - set(metrics.columns)
    if missing:
        raise ValueError(f"metrics is missing column(s): {sorted(missing)}")

    fail = pd.Series(False, index=metrics.index)
    for _, group in metrics.groupby(["study_id", COMPARTMENT_COLUMN], observed=True):
        idx = group.index
        fail.loc[idx] |= _mad_outlier(group["n_genes"], n_mads)
        fail.loc[idx] |= _mad_outlier(group["n_counts"], n_mads)
    fail |= metrics["pct_mito"] > max_pct_mito
    return fail


#: A tumour/normal retention gap beyond this many percentage points is flagged.
#: W1 uses the same 10-point threshold in ``src/reference/qc.py``; the two
#: cohorts are checked against one number so the gate can compare them.
RETENTION_GAP_WARN_PTS: Final[float] = 10.0


def differential_retention(
    metrics: pd.DataFrame,
    passes: pd.Series,
    *,
    compartment: str | None = "Epithelial cells",
    warn_at: float = RETENTION_GAP_WARN_PTS,
) -> pd.DataFrame:
    """Per-patient QC retention, tumour vs normal. **Read this before filtering.**

    QC is not neutral in this project. The compositional term is Δ(mature
    fraction) between a patient's tumour and their own normal, so a filter that
    removes cells at different rates in the two arms moves that difference
    directly, and it does so without appearing anywhere in the result.

    Both directions are bad and they are bad asymmetrically. A patient whose
    *normal* loses more has their normal mature fraction understated, inflating
    the apparent compositional loss — W1's mitochondrial-cap finding on
    GSE178341. A patient whose *tumour* loses more, which is what pooled MAD did
    on Lee, understates the tumour mature fraction and inflates it too, because
    the cells a depth filter takes first are the deep mature colonocytes on
    either side. Neither is a wash.

    ``compartment`` restricts the check to the compartment the compositional arm
    is actually built from; pass None to check every cell. Returns one row per
    patient with both retentions, their difference in percentage points, and a
    ``flagged`` column — the counterpart of W1's ``qc.differential_retention``,
    written here rather than imported because CONTRIBUTING §2 puts
    ``src/reference/`` off-limits to W4 and it is not exported.
    """
    if len(passes) != len(metrics):
        raise ValueError(f"passes has {len(passes)} entries for {len(metrics)} cells")
    for col in ("patient_id", "tissue"):
        if col not in metrics.columns:
            raise ValueError(f"metrics needs a {col!r} column for the retention check")

    frame = metrics.copy()
    frame["passed"] = np.asarray(passes, dtype=bool)
    if compartment is not None:
        if COMPARTMENT_COLUMN not in frame.columns:
            raise ValueError(
                f"metrics needs a {COMPARTMENT_COLUMN!r} column to restrict the "
                f"check to {compartment!r}; pass compartment=None to check all cells"
            )
        frame = frame.loc[frame[COMPARTMENT_COLUMN] == compartment]

    wide = (
        frame.loc[frame["tissue"].isin(["normal", "tumour"])]
        .groupby(["patient_id", "tissue"], observed=True)["passed"]
        .mean()
        .unstack("tissue")
        .mul(100.0)
    )
    for arm in ("normal", "tumour"):
        if arm not in wide.columns:
            wide[arm] = np.nan
    wide = wide.loc[:, ["normal", "tumour"]].dropna()
    wide["gap_pts"] = wide["tumour"] - wide["normal"]
    wide["flagged"] = wide["gap_pts"].abs() > warn_at
    return wide.reset_index()


def flag_doublets(*args, **kwargs):
    """scDblFinder or Scrublet against real Lee count matrices. W4, weeks 1-2.

    Needs the actual GSE132465/GSE144735 matrices to run and to tune — not
    scaffolded further than this signature. See src/estimator/README.md.
    """
    raise NotImplementedError("W4 — doublet detection needs real Lee count matrices.")


def correct_ambient(*args, **kwargs):
    """SoupX against a real raw+filtered matrix pair. W4, weeks 1-2.

    Same pipeline shape as W1 — coordinate on it, do not share code
    prematurely (src/estimator/README.md). Needs real per-sample soup
    profiles to fit; ``r-soupx`` is pinned in env/w4_estimator.yml for this.
    """
    raise NotImplementedError("W4 — ambient correction needs real Lee count matrices.")
