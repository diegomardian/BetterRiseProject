"""Lee cohort ingest: per-study QC, doublet flagging, ambient correction. W4.

Per-study thresholds, not global (execution_plan.md W1 task table — W4 mirrors
the same pipeline SHAPE for Lee, coordinating with W1 rather than sharing code
prematurely, src/estimator/README.md). MAD-based outlier flags are computed
independently within each ``study_id``, since GSE132465 (SMC) and GSE144735
(KUL3) come from different protocols — a single global cutoff would over- or
under-filter one of them.

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

import pandas as pd

DEFAULT_N_MADS: Final[float] = 5.0
DEFAULT_MAX_PCT_MITO: Final[float] = 20.0

_REQUIRED_METRIC_COLUMNS = {"study_id", "n_genes", "n_counts", "pct_mito"}


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
    ``n_genes``, ``n_counts``, ``pct_mito`` — one row per cell.

    A cell fails if ``n_genes`` or ``n_counts`` is a MAD-outlier *within its
    own study* (``n_mads``, default 5 — the sc-best-practices convention), or
    ``pct_mito`` exceeds ``max_pct_mito``. The mitochondrial cap is a single
    hard threshold, not per-study — a dying cell's mitochondrial fraction
    isn't expected to be protocol-relative the way library size is.

    execution_plan.md wants "QC thresholds documented with rationale," not
    just applied — log the per-study fail counts this produces before moving
    on, the same way W1's task table asks for a tabulation, not just a filter.
    """
    missing = _REQUIRED_METRIC_COLUMNS - set(metrics.columns)
    if missing:
        raise ValueError(f"metrics is missing column(s): {sorted(missing)}")

    fail = pd.Series(False, index=metrics.index)
    for _, group in metrics.groupby("study_id"):
        idx = group.index
        fail.loc[idx] |= _mad_outlier(group["n_genes"], n_mads)
        fail.loc[idx] |= _mad_outlier(group["n_counts"], n_mads)
    fail |= metrics["pct_mito"] > max_pct_mito
    return fail


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
