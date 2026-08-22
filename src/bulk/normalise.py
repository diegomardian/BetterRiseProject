"""Both scales, kept apart on purpose. W3.

Deconvolution needs linear non-log space; survival modelling wants log. The
brief's constraint is that neither ever silently substitutes for the other, and
that there is **an assertion at every entry point that checks the scale of the
input**. That is what this module is: the conversions are four lines each, the
guards are the reason the file exists.

W2 hit the same class of bug from the other direction — ``lee_io`` handed CP10K
floats to a generator expecting raw counts, and ``astype(int64)`` silently
turned every value below 1.0 into a zero (docs/open_decisions.md #8). The fix
there was to refuse the wrong scale loudly. Same principle here.

Orientation is **samples x genes** throughout, matching the AnnData convention
W1's ingest already uses, so nobody has to remember which arm transposes.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

#: log2(CPM + 1) of a full-depth bulk library tops out around 20 (a gene at
#: ~1e6 CPM). Anything far above that is linear data mislabelled as log.
MAX_PLAUSIBLE_LOG2_CPM = 32.0

#: A linear bulk matrix has a long right tail — highly expressed genes reach
#: thousands of TPM. If nothing exceeds this, it is almost certainly logged.
MIN_PLAUSIBLE_LINEAR_MAX = 100.0

#: TPM sums to this per sample by construction. Allowed drift before we call it
#: a different quantity: 0.5%, which covers float accumulation over ~60k genes
#: but not a matrix that was subset without renormalising.
TPM_TOTAL = 1e6
TPM_TOTAL_RTOL = 5e-3


class ScaleError(AssertionError):
    """A matrix was passed on the wrong scale. Do not coerce; find the caller."""


def _finite_values(matrix: pd.DataFrame) -> np.ndarray:
    values = matrix.to_numpy(dtype=float).ravel()
    return values[np.isfinite(values)]


def _describe(matrix: pd.DataFrame) -> str:
    values = _finite_values(matrix)
    if values.size == 0:
        return "empty"
    return (
        f"min {values.min():.4g}, median {np.median(values):.4g}, "
        f"max {values.max():.4g}, shape {matrix.shape}"
    )


# ---------------------------------------------------------------------------
# Guards — call these at every entry point
# ---------------------------------------------------------------------------


def assert_counts(matrix: pd.DataFrame, *, context: str = "the matrix") -> None:
    """Raise unless `matrix` holds raw non-negative integer counts."""
    values = _finite_values(matrix)
    if values.size == 0:
        raise ScaleError(f"{context} is empty or all-NaN.")
    if values.min() < 0:
        raise ScaleError(
            f"{context} has negative values (min {values.min():.4g}). "
            f"That is scaled or centred data, not counts."
        )
    if not np.allclose(values, np.rint(values)):
        bad = values[~np.isclose(values, np.rint(values))][:3]
        raise ScaleError(
            f"{context} has non-integer values (e.g. {bad}). Counts are integers; "
            f"this is already normalised. TPM and log2-CPM both derive from raw "
            f"counts, so normalising twice is a silent error."
        )


def assert_linear_scale(matrix: pd.DataFrame, *, context: str = "the matrix") -> None:
    """Raise if `matrix` looks log-transformed. Deconvolution needs linear space.

    Deconvolution solves a mixture equation that is linear in cell fractions.
    Handing it log values does not fail — it returns fractions that are wrong in
    a plausible-looking way, which is worse.
    """
    values = _finite_values(matrix)
    if values.size == 0:
        raise ScaleError(f"{context} is empty or all-NaN.")
    if values.min() < 0:
        raise ScaleError(
            f"{context} has negative values (min {values.min():.4g}) and cannot be "
            f"a linear expression scale."
        )
    if values.max() < MIN_PLAUSIBLE_LINEAR_MAX:
        raise ScaleError(
            f"{context} looks log-transformed, not linear: {_describe(matrix)}. "
            f"A linear bulk matrix has highly expressed genes in the thousands. "
            f"Deconvolution is linear in cell fractions — log input returns "
            f"wrong fractions that look plausible. Pass the TPM matrix."
        )


def assert_log_scale(matrix: pd.DataFrame, *, context: str = "the matrix") -> None:
    """Raise if `matrix` looks linear. Survival models want log2-CPM."""
    values = _finite_values(matrix)
    if values.size == 0:
        raise ScaleError(f"{context} is empty or all-NaN.")
    if values.max() > MAX_PLAUSIBLE_LOG2_CPM:
        raise ScaleError(
            f"{context} looks linear, not log2: {_describe(matrix)}. log2(CPM+1) "
            f"cannot exceed ~20 for a real library. Pass the log2-CPM matrix."
        )


def assert_tpm(matrix: pd.DataFrame, *, context: str = "the matrix") -> None:
    """Raise unless every row sums to 1e6. Catches a subset that was not renormalised.

    This is the specific failure the brief's "document the exact filtering rule"
    is guarding against: GDC computes TPM over its full gene set, so restricting
    the matrix to the shared index leaves rows that no longer sum to 1e6 and are
    therefore no longer comparable across samples.
    """
    assert_linear_scale(matrix, context=context)
    totals = matrix.sum(axis=1).to_numpy(dtype=float)
    off = ~np.isclose(totals, TPM_TOTAL, rtol=TPM_TOTAL_RTOL)
    if off.any():
        worst = totals[off]
        raise ScaleError(
            f"{context}: {int(off.sum())} of {len(totals)} sample(s) do not sum to "
            f"1e6 (e.g. {worst[:3].round(1)}). A TPM matrix subset to a gene index "
            f"must be renormalised — call renormalise_tpm()."
        )


# ---------------------------------------------------------------------------
# Conversions
# ---------------------------------------------------------------------------


def counts_to_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    """Counts -> counts per million. Linear."""
    assert_counts(counts, context="counts_to_cpm input")
    totals = counts.sum(axis=1).replace(0, np.nan)
    return counts.div(totals, axis=0) * 1e6


def counts_to_log2_cpm(counts: pd.DataFrame, *, prior_count: float = 1.0) -> pd.DataFrame:
    """Counts -> log2(CPM + prior_count).

    The prior is 1.0 rather than a small epsilon deliberately: it keeps zeros at
    exactly 0.0 and compresses the variance of low-count genes, which is the
    behaviour a Cox model wants. Stated here because "log-CPM" names a family,
    not a transform, and the offset changes the numbers.
    """
    if prior_count <= 0:
        raise ValueError(f"prior_count must be positive, got {prior_count}")
    out = np.log2(counts_to_cpm(counts) + prior_count)
    assert_log_scale(out, context="counts_to_log2_cpm output")
    return out


def renormalise_tpm(tpm: pd.DataFrame) -> pd.DataFrame:
    """Rescale a gene-subset TPM matrix so each sample sums to 1e6 again.

    Correct without gene lengths: TPM is already length-normalised, so a subset
    is renormalised by rescaling, not recomputed. This is why the pipeline keeps
    GDC's ``tpm_unstranded`` rather than reconstructing TPM from counts — the
    effective lengths GDC used are not published with the counts, and guessing
    them would introduce a discrepancy nobody could later diagnose.
    """
    assert_linear_scale(tpm, context="renormalise_tpm input")
    totals = tpm.sum(axis=1).replace(0, np.nan)
    out = tpm.div(totals, axis=0) * TPM_TOTAL
    assert_tpm(out, context="renormalise_tpm output")
    return out
