"""Forcing inputs for the adenoma claim-sensitivity consolidation.

The claims under test are the ones a robustness table can get quietly wrong:
that a `None` intrinsic term is dropped rather than read as zero, that both
denominators are actually carried (and can differ), that a log-ratio undefined at
a zero compositional term is counted rather than floored, and that the cross-block
count halves ordered pairs the way the prose does.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.reference.claim_sensitivity import (
    SensitivityError,
    contrasts_by_denominator,
    cross_block_survival,
    estimability_attrition,
)


def _split(*, intrinsic_zero_last=None, weighting="normal", rung="lineage",
           denominator_shift=0.0) -> pd.DataFrame:
    """Four patients x two genes (a target and a control), one weighting."""
    rows = []
    base = {"GUCA2A": (4.0, 1.0), "KRT8": (2.0, 1.0)}
    for pid in range(4):
        for gene, (intr, comp) in base.items():
            intrinsic = intr + 0.1 * pid + denominator_shift
            if intrinsic_zero_last is not None and gene == "GUCA2A" and pid == 3:
                intrinsic = intrinsic_zero_last
            rows.append({
                "patient_id": f"P{pid}", "gene": gene,
                "granularity_rung": rung, "weighting": weighting,
                "intrinsic": intrinsic, "compositional": comp,
                "estimability": "ok",
            })
    return pd.DataFrame(rows)


def _compositional(split: pd.DataFrame, rung="lineage") -> pd.DataFrame:
    frame = split[["patient_id", "gene", "granularity_rung"]].copy()
    frame["granularity_rung"] = rung
    frame["compositional_estimability"] = "ok"
    return frame


def test_both_denominators_and_all_statistics_are_carried():
    split = _split()
    out = contrasts_by_denominator(split, split, seed=1)
    assert set(out["denominator"]) == {"resolved", "all_epithelial"}
    assert set(out["statistic"]) == {"log_ratio", "share_abs", "share_signed",
                                     "ratio"}
    assert out["load_bearing"].sum() == (out["statistic"] == "log_ratio").sum()
    # The target-control pair is the only pair, and it appears in both orders.
    assert set(out["cross_block"]) == {True}


def test_denominators_can_differ_and_are_not_averaged():
    primary = _split(denominator_shift=0.0)
    shifted = _split(denominator_shift=5.0)
    out = contrasts_by_denominator(primary, shifted, seed=1)
    lineage = out[(out.granularity_rung == "lineage")
                  & (out.statistic == "share_abs")]
    resolved = lineage[lineage.denominator == "resolved"]["centre"].to_numpy()
    epithelial = lineage[lineage.denominator == "all_epithelial"]["centre"].to_numpy()
    assert not np.allclose(resolved, epithelial)


def test_none_intrinsic_is_dropped_and_counted_not_zeroed():
    """Invariant 1 in the sensitivity table.

    One patient's intrinsic term is `None`. They must leave the contrast (n drops
    from 4 to 3) and be counted in the attrition table, never enter a numerator
    as a zero that would pull the centre down.
    """
    split = _split(intrinsic_zero_last=None)
    split.loc[(split.gene == "GUCA2A") & (split.patient_id == "P3"), "intrinsic"] = None
    out = contrasts_by_denominator(split, split, seed=1)
    guca = out[(out.contrast == "GUCA2A - KRT8")
               & (out.denominator == "resolved")
               & (out.statistic == "share_abs")].iloc[0]
    assert guca["n_patients"] == 3
    assert guca["n_patients_cohort"] == 4
    assert guca["n_missing"] == 1

    attrition = estimability_attrition(split, split, _compositional(split))
    guca_att = attrition[(attrition.gene == "GUCA2A")
                         & (attrition.denominator == "resolved")].iloc[0]
    assert guca_att["n_intrinsic_missing"] == 1
    assert guca_att["n_intrinsic_present"] == 3


def test_log_ratio_is_undefined_at_a_zero_compositional_term():
    split = _split()
    split.loc[split.granularity_rung == "lineage", "compositional"] = 0.0
    attrition = estimability_attrition(
        split, split, _compositional(split)
    )
    resolved = attrition[attrition.denominator == "resolved"]
    assert resolved["n_log_ratio_undefined"].sum() == len(split)
    assert resolved["n_log_ratio_defined"].sum() == 0


def test_cross_block_survival_halves_ordered_pairs():
    split = _split()
    sensitivity = contrasts_by_denominator(split, split, seed=1)
    summary = cross_block_survival(sensitivity)
    cell = summary[(summary.granularity_rung == "lineage")
                   & (summary.weighting == "normal")
                   & (summary.statistic == "share_abs")
                   & (summary.denominator == "resolved")].iloc[0]
    # Two ordered pairs, so one unordered contrast; GUCA2A has a much larger
    # share than KRT8 on every patient, so it excludes zero.
    assert cell["n_contrasts"] == 1
    assert cell["n_excluding_zero"] == 1
    assert cell["share_excluding_zero"] == 1.0


def test_missing_columns_are_refused():
    bad = _split().drop(columns=["compositional"])
    with pytest.raises(SensitivityError, match="missing"):
        contrasts_by_denominator(bad, bad, seed=1)
