"""Forcing inputs for the Crowell reproduction comparison.

The failure this guards is a reproduction that reports success because it never
compared anything: a mismatch in a gene's mean, its block count, or its interval
verdict must flip `reproduces` to False, and a gene missing from either side
must raise rather than drop out.
"""

from __future__ import annotations

import pandas as pd
import pytest

from src.reference.crowell_reproduction import (
    ReproductionError,
    compare_summaries,
)


def _summary(*, guca=(-1.307, False), n_blocks=7) -> pd.DataFrame:
    return pd.DataFrame([
        {"gene": "GUCA2A", "role": "target", "n_blocks": n_blocks,
         "mean_did": guca[0], "ci_low": guca[0] - 0.5, "ci_high": guca[0] + 0.5,
         "interval": "student_t", "excludes_zero": guca[1]},
        {"gene": "CDX2", "role": "identity", "n_blocks": n_blocks,
         "mean_did": 0.327, "ci_low": -0.150, "ci_high": 0.805,
         "interval": "student_t", "excludes_zero": False},
    ])


def test_an_exact_reproduction_is_flagged_true():
    out = compare_summaries(_summary(), _summary())
    assert out["reproduces"].all()
    assert out["mean_did_diff"].max() == 0.0


def test_a_changed_mean_fails_the_reproduction():
    reproduced = _summary(guca=(-1.2, True))
    committed = _summary(guca=(-1.307, True))
    out = compare_summaries(reproduced, committed).set_index("gene")
    assert not out.loc["GUCA2A", "reproduces"]
    assert out.loc["CDX2", "reproduces"]
    assert out.loc["GUCA2A", "mean_did_diff"] == pytest.approx(0.107, abs=1e-6)


def test_a_changed_block_count_fails_even_with_the_same_mean():
    out = compare_summaries(_summary(n_blocks=6), _summary(n_blocks=7))
    assert not out["reproduces"].any()


def test_a_changed_interval_verdict_fails():
    reproduced = _summary(guca=(-1.307, False))
    committed = _summary(guca=(-1.307, True))  # excludes zero changed
    out = compare_summaries(reproduced, committed).set_index("gene")
    assert not out.loc["GUCA2A", "reproduces"]


def test_a_gene_on_only_one_side_raises():
    """A reproduction must cover the same genes. Dropping one is how a missing
    result reads as a passing check."""
    reproduced = _summary().iloc[:1]     # only GUCA2A
    with pytest.raises(ReproductionError, match="only one summary"):
        compare_summaries(reproduced, _summary())


def test_missing_columns_are_refused():
    with pytest.raises(ReproductionError, match="missing"):
        compare_summaries(pd.DataFrame({"gene": ["GUCA2A"]}), _summary())
