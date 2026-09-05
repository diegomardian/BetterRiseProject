"""Stage 4's pre-committed instrument check, and what it does and does not cover.

The locked pre-specification says:

    The deconvolved non-epithelial (stromal + immune + endothelial) fraction
    correlates with (1 - ABSOLUTE purity) at Pearson r >= 0.5 on the samples
    with an ABSOLUTE call.
    on_failure: STOP. Report the instrument failure as the Stage 4 result.

That is implemented here exactly as written, because it is locked and this
module does not get to reinterpret it. Two things about it are worth stating
plainly rather than discovering later.

**ABSOLUTE, and only ABSOLUTE.** The committed purity table carries three
methods: `absolute` (556 barcodes, called from copy number), `aran_cpe` (623,
consensus) and `estimate_affy_extrapolated` (675, derived from EXPRESSION). The
largest is the one that must not be used. The gate's entire logical force is
that purity is measured independently of the expression the fractions came from;
correlate an expression-derived fraction against an expression-derived purity
and it passes because both are reading the same signal. `expression_derived` is
a column in that table and this module filters on it as well as on the method
name, so a future method that is expression-derived is excluded by default
rather than by having been thought of.

**It passes for a quantity Stage 4 does not use.** The gate reads the
non-epithelial aggregate -- epithelium against everything else, a coarse and
high-contrast split. Stage 4's predictor is the epithelial-INTERNAL split, and
that one is far more fragile: in a run where the committed log-scale reference
drives `mature_colonocyte_fraction` to exactly 0.0 on every sample, this gate
still reads r = 0.881 and passes. See `src/bulk/deconvolution.py`. The gate is
necessary and it is not sufficient, and `check_predictor` is the other half.

Saying so here is not an amendment to the locked spec. The spec's own
`estimability` section already requires that an unestimable fraction be `None`
rather than `0.0` (invariant 1); refusing a constant-zero predictor enforces
that clause rather than adding one.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)


class InstrumentError(RuntimeError):
    """The instrument check cannot be evaluated on these inputs."""


#: The one purity call independent of expression. See the module docstring.
INDEPENDENT_PURITY_METHOD = "absolute"


@dataclass(frozen=True)
class InstrumentResult:
    """The gate's verdict, with everything needed to audit it."""

    method: str
    granularity_rung: str
    n_samples: int
    r: float
    threshold: float
    purity_method: str
    passed: bool
    detail: str

    def as_row(self) -> dict:
        return {
            "method": self.method,
            "granularity_rung": self.granularity_rung,
            "purity_method": self.purity_method,
            "n_samples": self.n_samples,
            "pearson_r": self.r,
            "threshold": self.threshold,
            "passed": self.passed,
            "detail": self.detail,
        }


def independent_purity(purity: pd.DataFrame) -> pd.DataFrame:
    """ABSOLUTE calls only, and only those not derived from expression.

    Two filters where one would do, deliberately. The method name is the spec's
    wording; `expression_derived` is the property the spec's *rationale* depends
    on, and a table gaining a new independent-sounding method should not
    silently join the gate.
    """
    for column in ("method", "purity", "expression_derived", "barcode"):
        if column not in purity.columns:
            raise InstrumentError(
                f"the purity table has no {column!r} column, so the gate cannot "
                f"establish that its purity call is independent of expression. "
                f"Columns: {sorted(purity.columns)}"
            )
    rows = purity[
        (purity["method"] == INDEPENDENT_PURITY_METHOD)
        & (~purity["expression_derived"].astype(bool))
    ]
    if rows.empty:
        available = sorted(purity["method"].unique())
        raise InstrumentError(
            f"no {INDEPENDENT_PURITY_METHOD!r} purity calls that are not "
            f"expression-derived. Available methods: {available}. The gate is "
            f"only a check because purity is measured independently of the "
            f"expression the fractions came from -- it cannot fall back to an "
            f"expression-derived call."
        )
    return rows.dropna(subset=["purity"])


def sample_key(barcode: str) -> str:
    """TCGA barcode -> patient-level key, so DNA and RNA aliquots can meet.

    ABSOLUTE is called from DNA and the expression is RNA, so the two never
    share an aliquot barcode. Delegates to the purity module's own key so the
    two cannot drift.
    """
    from src.bulk.purity import sample_key as _key

    return _key(barcode)


def run_instrument_check(
    fractions: pd.DataFrame,
    purity: pd.DataFrame,
    *,
    method: str,
    rung: str,
    threshold: float = 0.5,
) -> InstrumentResult:
    """Pearson r between non-epithelial fraction and (1 - ABSOLUTE purity).

    ``fractions`` is ``summarise_fractions``' output, one row per
    (sample_id, method). Samples are matched on the patient-level key rather
    than the aliquot barcode.
    """
    absolute = independent_purity(purity)
    rows = fractions[fractions["method"] == method]
    if rows.empty:
        raise InstrumentError(f"no fractions for method {method!r}")

    left = rows.assign(_key=rows["sample_id"].map(sample_key))
    right = absolute.assign(_key=absolute["barcode"].map(sample_key))
    # One purity call per patient. A patient with two aliquots would otherwise
    # weight twice, and the duplicate is an artefact of aliquoting rather than
    # of anything the gate is measuring.
    right = right.groupby("_key", as_index=False)["purity"].median()

    merged = left.merge(right, on="_key", how="inner").dropna(
        subset=["non_epithelial_fraction", "purity"]
    )
    n = len(merged)
    if n < 3:
        raise InstrumentError(
            f"only {n} sample(s) carry both a deconvolved fraction and an "
            f"ABSOLUTE call. The gate cannot be evaluated, which is not the "
            f"same as it failing -- report it as unevaluable."
        )

    x = merged["non_epithelial_fraction"].to_numpy(dtype=float)
    y = 1.0 - merged["purity"].to_numpy(dtype=float)
    if np.ptp(x) == 0.0 or np.ptp(y) == 0.0:
        constant = "the non-epithelial fraction" if np.ptp(x) == 0 else "1 - purity"
        return InstrumentResult(
            method=method, granularity_rung=rung, n_samples=n, r=float("nan"),
            threshold=threshold, purity_method=INDEPENDENT_PURITY_METHOD,
            passed=False,
            detail=(
                f"{constant} is constant over all {n} samples, so the "
                f"correlation is undefined. Undefined is not a pass: a "
                f"comparison against NaN is False in numpy either way, and "
                f"reporting it as a failure to correlate would be reporting a "
                f"number that does not exist."
            ),
        )

    r = float(np.corrcoef(x, y)[0, 1])
    passed = bool(r >= threshold)
    return InstrumentResult(
        method=method, granularity_rung=rung, n_samples=n, r=r, threshold=threshold,
        purity_method=INDEPENDENT_PURITY_METHOD, passed=passed,
        detail=(
            f"r = {r:.3f} against a threshold of {threshold} on {n} samples with "
            f"an ABSOLUTE call. "
            + (
                "The instrument recovers a fraction we already know from an "
                "orthogonal assay. This says nothing about the epithelial-"
                "internal split Stage 4 regresses on -- see check_predictor."
                if passed else
                "STOP per the locked pre-specification: report the instrument "
                "failure as the Stage 4 result and report no R-squared."
            )
        ),
    )


def gate_verdict(results: list[InstrumentResult]) -> tuple[bool, str]:
    """Does Stage 4 proceed? Per-method, because the failure is method-specific.

    Returns ``(proceed, message)``. Proceeding needs at least one method whose
    gate passed; that method's fractions are the only ones any R-squared may be
    computed from, and a method that failed the gate does not get quietly
    averaged in with one that passed.
    """
    if not results:
        return False, "the instrument check was never run"
    passing = [r for r in results if r.passed]
    lines = [f"  {r.method}/{r.granularity_rung}: {r.detail}" for r in results]
    if not passing:
        return False, (
            "STOP. No method passed the pre-committed instrument check, so the "
            "Stage 4 result is the instrument failure and no R-squared is "
            "reported (locked prespec, instrument_checks.positive_control."
            "on_failure):\n" + "\n".join(lines)
        )
    return True, (
        f"{len(passing)} of {len(results)} (method, rung) combinations passed:\n"
        + "\n".join(lines)
    )
