"""Bulk -> cell-type fractions, for Stage 4's variance question.

FRACTIONS ONLY. CLAUDE.md invariant 6: bulk recovers fractions at r ~ 0.92, and
cell-type-specific expression recovered from bulk comes back attenuated x0.6-0.8
in the direction of this project's prior hypothesis. Nothing here returns an
expression estimate, and nothing here should be extended to.

WHAT THIS MODULE EXISTS TO PREVENT. The committed reference matrices are built
as the mean of ``log1p(CP10K)`` -- ``build_signature``'s "mean of log, not log
of mean", which is the right choice for marker SELECTION and the wrong scale for
a linear mixture. Bulk TPM/CPM is a linear mixture of linear per-cell-type
profiles. Deconvolving it against a log-scale reference is a misspecification,
and ON THE COMMITTED MATRICES its failure mode is not noise:

    linear bulk, committed log lineage reference, NNLS
        -> the `differentiated` fraction comes back EXACTLY 0.0 on 200 of 200
           synthetic samples whose true mean is 0.29; linearised, the same
           matrix returns a mean of 0.306 with no exact zeros

That column is Stage 4's entire predictor. A constant predictor gives every gene
an R-squared of zero, which reads as "fraction explains nothing" -- a result, in
the pre-registered direction, produced by the instrument being broken.

BE PRECISE ABOUT WHY, because the general claim is false. A log reference does
not zero the predictor on just any matrix: on a synthetic reference where each
cell type owns a private block of marker genes, the same mismatch still recovers
the mature fraction fine, and `tests/test_bulk_deconvolution.py` asserts that so
the narrower claim cannot quietly widen. It takes the mismatch AND a reference
whose columns sit close together, and the committed one does -- log compression
puts `differentiated` at cosine 0.982 to `epithelial_unscored`, against 0.929
linearised. Dropping that column from the log matrix reduces the collapse from
200/200 to 121/200 without fixing it, so the scale is the primary cause and the
near-collinearity compounds it.

The pre-committed instrument gate does not catch it. The gate reads the
non-epithelial aggregate against ABSOLUTE purity, and that aggregate is recovered
at r = 0.881 in the same run where the mature fraction is a constant. Epithelial
versus everything-else is a coarse, high-contrast split that survives a scale
error the epithelial-internal split does not. **The gate passes for a quantity
the analysis does not use.**

So this module refuses the mismatch up front rather than deconvolving and hoping
the gate notices, and it refuses a constant fraction column afterwards under
invariant 1 -- a quantity that is not estimable is `None`, never `0.0`, and a
column of exact zeros for every patient is that coercion at cohort scale.

    from src.bulk.deconvolution import load_reference, deconvolve_cohort
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from src.harness.deconvolve.base import Deconvolver, available_methods
from src.harness.deconvolve.nnls import NNLSDeconvolver
from src.harness.deconvolve.nusvr import NuSVRDeconvolver

log = logging.getLogger(__name__)


class DeconvolutionError(RuntimeError):
    """The inputs cannot support a fraction estimate, or the output is not one."""


#: How a reference matrix's values are scaled. Not cosmetic: it decides whether
#: the linear mixture model the deconvolvers assume is the model that generated
#: the bulk.
LOG1P_CP10K = "log1p_cp10k"
LINEAR_CP10K = "linear_cp10k"
SCALES = (LOG1P_CP10K, LINEAR_CP10K)

#: Columns that are not epithelium, at every rung. The instrument gate's
#: aggregate. Named here rather than inferred, so a new rung with a new
#: non-epithelial column fails loudly instead of being silently counted as
#: epithelial.
NON_EPITHELIAL = ("immune", "stromal", "endothelial")

#: Epithelial columns carrying no maturity call. Part of neither the mature
#: fraction nor the non-epithelial aggregate -- they are the reason those two do
#: not sum to one, and folding them into either would be a choice nobody made.
UNSCORED = ("epithelial_unscored",)

#: Below this, a recovered fraction column is not varying enough to regress on.
#: A column that is CONSTANT is refused outright (invariant 1); this is the
#: softer neighbouring case, reported rather than raised.
MIN_FRACTION_SD = 0.01


def default_methods() -> list[Deconvolver]:
    """The two adapters that run in the pinned env. Never averaged.

    Returned as a list so the driver reports both side by side. A mean of two
    deconvolvers is a number no method produced, and the failure this module
    exists to catch is method-specific -- NNLS collapses under a scale mismatch
    where nu-SVR degrades gracefully, so averaging them would half-hide it.
    """
    return [NNLSDeconvolver(), NuSVRDeconvolver()]


def mature_column(rung: str) -> str | None:
    """The column this rung calls mature, or None where the rung has no split.

    Read off ``RUNG_SPECS`` rather than hard-coded: the bins are ordered least-
    to most-mature, so the mature one is the last. The ``epithelial`` rung has a
    single bin and therefore no maturity call at all -- it returns None, and a
    caller that treats None as a column name gets a TypeError rather than a
    silently wrong fraction.
    """
    from src.reference.labels import RUNG_SPECS

    bins = RUNG_SPECS[rung].bins
    return bins[-1] if len(bins) > 1 else None


# ---------------------------------------------------------------------------
# The reference matrix


@dataclass(frozen=True)
class Reference:
    """A signature matrix, its scale, and where it came from."""

    matrix: pd.DataFrame          # genes x cell types, gene ids on the index
    rung: str
    scale: str
    source: Path
    #: Set when `matrix` was derived rather than read -- see `linearise`.
    derived_from_scale: str | None = None

    @property
    def cell_types(self) -> list[str]:
        return list(self.matrix.columns)

    def describe(self) -> str:
        via = f", linearised from {self.derived_from_scale}" if self.derived_from_scale else ""
        return (
            f"{self.rung} rung, {self.matrix.shape[0]} genes x "
            f"{self.matrix.shape[1]} types, scale={self.scale}{via}, "
            f"from {self.source.name}"
        )


def load_reference(
    path: str | Path,
    *,
    rung: str,
    scale: str = LOG1P_CP10K,
    targets: list[str] | None = None,
) -> Reference:
    """Read an ``S_matrix_{rung}_{version}.parquet`` and check what it is.

    The committed matrices carry the gene id in a ``gene`` COLUMN and a
    positional 0..n-1 index. Promoting it is not tidying: a positional index
    silently aligns to whatever the bulk's positional index happens to be, which
    is the leakage class fixed in 23b1d83, and every join downstream depends on
    the promotion having happened here.
    """
    source = Path(path)
    if scale not in SCALES:
        raise DeconvolutionError(f"scale must be one of {SCALES}, got {scale!r}")
    matrix = pd.read_parquet(source)

    if "gene" not in matrix.columns:
        raise DeconvolutionError(
            f"{source.name} has no `gene` column, so its rows cannot be aligned "
            f"to bulk by identity. A positional index would align by ORDER, "
            f"which is the defect fixed in 23b1d83."
        )
    matrix = matrix.set_index("gene")
    if not matrix.index.is_unique:
        dupes = matrix.index[matrix.index.duplicated()].unique().tolist()
        raise DeconvolutionError(f"{source.name} has duplicate genes: {dupes[:5]}")
    if matrix.isna().to_numpy().any():
        raise DeconvolutionError(f"{source.name} contains NaN")

    numeric = matrix.select_dtypes("number")
    if numeric.shape[1] != matrix.shape[1]:
        non_numeric = sorted(set(matrix.columns) - set(numeric.columns))
        raise DeconvolutionError(f"{source.name} has non-numeric columns {non_numeric}")
    if (numeric.to_numpy() < 0).any():
        raise DeconvolutionError(f"{source.name} has negative expression")

    if targets:
        leaked = sorted(set(targets) & set(matrix.index))
        if leaked:
            raise DeconvolutionError(
                f"{source.name} contains target genes {leaked}. CLAUDE.md "
                f"invariant 2: a silenced mature cell must not be readable as "
                f"an absent mature cell."
            )

    mature = mature_column(rung)
    if mature is not None and mature not in matrix.columns:
        raise DeconvolutionError(
            f"{source.name} has no {mature!r} column, which is the {rung} rung's "
            f"mature bin. Columns: {sorted(matrix.columns)}"
        )
    return Reference(matrix=matrix, rung=rung, scale=scale, source=source)


def linearise(reference: Reference) -> Reference:
    """``expm1`` a log1p reference, so it can meet a linear bulk.

    APPROXIMATE, AND THE DIRECTION OF THE ERROR IS KNOWN. The committed matrix
    holds the mean of ``log1p(x)``, so ``expm1`` of it is a geometric rather
    than an arithmetic mean of ``CP10K + 1``. By Jensen it is biased LOW, and
    the bias grows with a gene's dispersion across cells of its type -- so the
    most variable genes are understated most.

    This is a repair for a reference that already exists, not the right way to
    build one. A reference built linearly in W1 needs the cells, which live on
    the cluster; until that exists, the sidecar of anything using this must say
    it was derived here. Refusing to run at all would leave the scale defect
    undetected in the committed matrix, which is worse.
    """
    if reference.scale == LINEAR_CP10K:
        return reference
    if reference.scale != LOG1P_CP10K:
        raise DeconvolutionError(f"cannot linearise a {reference.scale!r} reference")
    return Reference(
        matrix=np.expm1(reference.matrix),
        rung=reference.rung,
        scale=LINEAR_CP10K,
        source=reference.source,
        derived_from_scale=LOG1P_CP10K,
    )


def assert_scales_agree(reference: Reference, bulk_scale: str) -> None:
    """Require BOTH sides linear. Matching is not the criterion -- linearity is.

    Every deconvolver here solves ``bulk ~ S @ f``, and that identity holds only
    in linear space. It is tempting to read the requirement as "put the two on
    the same scale", which would make a log reference against log bulk
    acceptable. IT IS NOT: ``log(S @ f)`` is not ``log(S) @ f``, so a log/log
    pair is scale-matched and still not a mixture model. Two misspecifications
    that agree with each other are not a specification.

    So this refuses anything that is not linear on both sides, and the error
    says which side is wrong.

    THIS IS ALSO THE CHECK THE INSTRUMENT GATE CANNOT DO. The pre-committed gate
    correlates the non-epithelial aggregate against ABSOLUTE purity, and that
    aggregate survives the mismatch: in the run where NNLS returns a
    `differentiated` fraction of exactly 0.0 for every sample, the gate reads
    r = 0.881 and passes. Epithelium versus not-epithelium is a coarse split;
    the epithelial-INTERNAL split, which is the only part Stage 4 uses, is not.

    A gate that passes for a quantity the analysis does not use is not a gate
    for the analysis, so this is refused here instead -- before any number is
    produced, rather than after one is certified.
    """
    if bulk_scale not in SCALES:
        raise DeconvolutionError(f"bulk_scale must be one of {SCALES}, got {bulk_scale!r}")
    wrong = [
        name for name, scale in (("reference", reference.scale), ("bulk", bulk_scale))
        if scale != LINEAR_CP10K
    ]
    if wrong:
        raise DeconvolutionError(
            f"the {' and '.join(wrong)} {'are' if len(wrong) > 1 else 'is'} not "
            f"on the linear scale (reference={reference.scale!r}, "
            f"bulk={bulk_scale!r}). A linear mixture model needs BOTH linear, "
            f"and matching them on the log scale does not help: log(S @ f) is "
            f"not log(S) @ f, so a log/log pair is still not a mixture model.\n"
            f"For the reference, call `linearise()` -- and read its docstring, "
            f"because it is an approximation with a known direction -- or "
            f"rebuild it linearly in W1. For the bulk, use TPM or CPM rather "
            f"than log2-CPM; the log matrix is the Stage 4 OUTCOME, not its "
            f"deconvolution input.\n"
            f"This does not merely add noise. On the committed lineage matrix "
            f"NNLS returns the mature fraction as exactly 0.0 on every sample "
            f"-- and the pre-committed instrument gate passes anyway, because "
            f"the aggregate it reads survives."
        )


# ---------------------------------------------------------------------------
# Running it


def align(bulk: pd.DataFrame, reference: Reference) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Restrict both sides to the genes they share, aligned BY IDENTITY.

    ``bulk`` is samples x genes. Returns ``(bulk_aligned, signature_aligned)``
    on the same gene order, which is the join CLAUDE.md's repo layout promises
    ("integration is a join, not a negotiation").
    """
    shared = reference.matrix.index.intersection(bulk.columns)
    if len(shared) == 0:
        raise DeconvolutionError(
            "bulk and reference share no genes. Both must be on the 1.0.0 gene "
            "index -- check whether one carries symbols and the other Ensembl ids."
        )
    if len(shared) < 100:
        raise DeconvolutionError(
            f"only {len(shared)} shared genes. nu-SVR's robustness is a "
            f"high-dimensionality property (execution_plan.md §2.1 error #4); "
            f"below a few hundred genes the fractions are not trustworthy."
        )
    ordered = reference.matrix.index[reference.matrix.index.isin(shared)]
    return bulk.loc[:, ordered], reference.matrix.loc[ordered]


def deconvolve_cohort(
    bulk: pd.DataFrame,
    reference: Reference,
    *,
    bulk_scale: str = LINEAR_CP10K,
    methods: list[Deconvolver] | None = None,
) -> tuple[pd.DataFrame, dict[str, str]]:
    """Fractions for every sample under every available method.

    Returns ``(long_frame, skipped)`` where the frame is one row per
    (sample, method, cell_type) and ``skipped`` maps a method name to why it did
    not run -- named rather than absent, so a missing method is visible in the
    results instead of silently changing what was compared.
    """
    assert_scales_agree(reference, bulk_scale)
    bulk_aligned, signature = align(bulk, reference)

    usable, skipped = available_methods(methods or default_methods())
    if not usable:
        raise DeconvolutionError(f"no deconvolution method is available: {skipped}")

    rows: list[dict] = []
    values = bulk_aligned.to_numpy(dtype=float)
    for method in usable:
        for sample_id, vector in zip(bulk_aligned.index, values, strict=True):
            fractions = method.fit_predict(vector, signature)
            rows.extend(
                {
                    "sample_id": str(sample_id),
                    "method": method.name,
                    "cell_type": str(cell_type),
                    "fraction": float(value),
                }
                for cell_type, value in fractions.items()
            )
    return pd.DataFrame(rows), skipped


def summarise_fractions(long: pd.DataFrame, rung: str) -> pd.DataFrame:
    """One row per (sample, method): the two aggregates Stage 4 reads.

    ``mature_colonocyte_fraction`` is the predictor; ``non_epithelial_fraction``
    is what the instrument gate checks against ABSOLUTE purity. They are emitted
    side by side precisely so a reader can see the case this module exists for,
    where the second is fine and the first is not.
    """
    mature = mature_column(rung)
    wide = long.pivot_table(
        index=["sample_id", "method"], columns="cell_type", values="fraction"
    )
    present_non_epi = [c for c in NON_EPITHELIAL if c in wide.columns]
    if not present_non_epi:
        raise DeconvolutionError(
            f"none of {NON_EPITHELIAL} is a column of the {rung} reference, so "
            f"the instrument gate's aggregate cannot be formed. "
            f"Columns: {sorted(wide.columns)}"
        )

    out = pd.DataFrame(index=wide.index)
    out["non_epithelial_fraction"] = wide[present_non_epi].sum(axis=1)
    # Invariant 1 in its positive form: a rung with no maturity call does not
    # get a mature fraction of 0.0, it gets None, and the estimability column
    # says which of those it is.
    if mature is None:
        out["mature_colonocyte_fraction"] = np.nan
        out["mature_cell_type"] = None
        out["estimability"] = "not_estimable"
        out["estimability_reason"] = f"the {rung} rung has a single bin and makes no maturity call"
    else:
        out["mature_colonocyte_fraction"] = wide[mature]
        out["mature_cell_type"] = mature
        out["estimability"] = "estimated"
        out["estimability_reason"] = ""
    unscored = [c for c in UNSCORED if c in wide.columns]
    out["epithelial_unscored_fraction"] = wide[unscored].sum(axis=1) if unscored else 0.0
    out["granularity_rung"] = rung
    return out.reset_index()


# ---------------------------------------------------------------------------
# Invariant 1, applied to the predictor itself


@dataclass(frozen=True)
class PredictorCheck:
    """Whether a recovered fraction is an estimate at all."""

    method: str
    rung: str
    n_samples: int
    n_exact_zero: int
    sd: float
    is_constant: bool
    verdict: str
    detail: str

    @property
    def usable(self) -> bool:
        return self.verdict == "usable"


def check_predictor(summary: pd.DataFrame, *, rung: str, method: str) -> PredictorCheck:
    """Is ``mature_colonocyte_fraction`` a variable, or a coerced zero?

    CLAUDE.md invariant 1 says an unestimable quantity is ``None`` and never
    ``0.0``. A column that is 0.0 for every patient in the cohort is that
    coercion at cohort scale, and it arrives looking like data rather than like
    a missing value -- so nothing downstream can tell it apart from a real
    finding of "no mature cells anywhere".

    The distinction this makes, which the R-squared cannot: a CONSTANT predictor
    gives every gene R-squared 0, and 0 is exactly what the pre-registered
    primary arm expects for GUCA2A. Reported without this check, total
    instrument failure and the predicted result are the same number.
    """
    rows = summary[summary["method"] == method]
    if rows.empty:
        raise DeconvolutionError(f"no rows for method {method!r}")
    values = rows["mature_colonocyte_fraction"].to_numpy(dtype=float)
    finite = values[np.isfinite(values)]

    n_zero = int((finite == 0.0).sum())
    sd = float(finite.std()) if finite.size else float("nan")
    constant = bool(finite.size and np.ptp(finite) == 0.0)

    if not finite.size:
        verdict = "not_estimable"
        detail = f"the {rung} rung makes no maturity call, so there is no predictor"
    elif constant and finite[0] == 0.0:
        verdict = "refused"
        detail = (
            f"{method} returned mature_colonocyte_fraction == 0.0 on all "
            f"{finite.size} samples. That is invariant 1's coercion at cohort "
            f"scale: not-estimable written as zero. It is the documented "
            f"signature of a log-scale reference against linear bulk, and the "
            f"instrument gate does NOT catch it -- the non-epithelial aggregate "
            f"is recovered fine in the same run."
        )
    elif constant:
        verdict = "refused"
        detail = (
            f"{method} returned a constant mature_colonocyte_fraction of "
            f"{finite[0]:.4f}. A constant predictor gives every gene an "
            f"R-squared of 0, which is indistinguishable from the pre-registered "
            f"expectation for GUCA2A."
        )
    elif sd < MIN_FRACTION_SD:
        verdict = "degenerate"
        detail = (
            f"{method}'s mature_colonocyte_fraction varies by sd={sd:.4f}, below "
            f"the {MIN_FRACTION_SD} floor. Not refused, because it is not "
            f"constant, but every R-squared from it is near-noise and must be "
            f"reported with this number beside it."
        )
    else:
        verdict = "usable"
        detail = f"sd={sd:.4f} over {finite.size} samples, {n_zero} exact zeros"

    return PredictorCheck(
        method=method, rung=rung, n_samples=int(finite.size), n_exact_zero=n_zero,
        sd=sd, is_constant=constant, verdict=verdict, detail=detail,
    )


def require_usable_predictor(checks: list[PredictorCheck]) -> None:
    """Raise unless at least one method produced a predictor worth regressing.

    Not "all methods" -- the failure is method-specific by design, and NNLS
    collapsing while nu-SVR survives is a reportable fact rather than a reason
    to stop. But if EVERY method returns a constant, there is nothing to
    regress and the Stage 4 result is the instrument, not the biology.
    """
    if any(c.usable for c in checks):
        return
    raise DeconvolutionError(
        "no method produced a usable mature-colonocyte fraction:\n"
        + "\n".join(f"  {c.method}/{c.rung}: {c.verdict} -- {c.detail}" for c in checks)
    )
