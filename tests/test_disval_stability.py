"""The DIS/VAL split, against the branches §5 fixed before the numbers.

The verdict function is where the pre-registration stops being prose. Its three
outcomes are not symmetric and the asymmetry is the whole design: a half-miss is
weak evidence because the interval is 1.80x wider at n=15, and a sign reversal
is strong because width does not flip signs. A verdict function that treated
them alike would let a power artefact read as a refutation.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.reference.jobs.adenoma_decomposition_scales import LOAD_BEARING
from src.reference.jobs.disval_stability import (
    PRIMARY_HALVES,
    StabilityError,
    assign_halves,
    verdict,
)

CONTRASTS = ["GUCA2A - ACTB", "GUCA2A - CDX2", "GUCA2A - EPCAM", "GUCA2A - KRT8"]


def _table(discovery, validation, *, means=None) -> pd.DataFrame:
    """`discovery`/`validation` are how many of the 4 contrasts exclude zero."""
    means = means or {}
    rows = []
    for half, n_holding in zip(PRIMARY_HALVES, (discovery, validation), strict=True):
        for i, contrast in enumerate(CONTRASTS):
            rows.append({
                "half": half, "statistic": LOAD_BEARING, "contrast": contrast,
                "mean": means.get((half, contrast), 1.0),
                "excludes_zero": i < n_holding,
            })
    return pd.DataFrame(rows)


def test_all_four_in_both_halves_excludes_batch_drivenness():
    out = verdict(_table(4, 4))
    assert out["verdict"] == "BATCH-DRIVENNESS EXCLUDED"
    assert "does not license removing that qualifier" in out["detail"], (
        "the single-cohort caveat survives even the best outcome"
    )


def test_a_half_miss_is_ambiguous_by_design_and_not_a_failure():
    """THE ASYMMETRY. At n=15 a real effect can miss on width alone."""
    out = verdict(_table(1, 4))
    assert out["verdict"] == "AMBIGUOUS AT THIS N"
    assert "must not be read as a failure" in out["detail"]
    assert "1.80x" in out["detail"]


def test_a_sign_reversal_is_strong_evidence_and_is_reported_as_such():
    """Width does not flip signs, so this is not a power artefact."""
    flipped = _table(4, 4, means={(PRIMARY_HALVES[1], "GUCA2A - CDX2"): -1.0})
    out = verdict(flipped)
    assert out["verdict"] == "SIGN REVERSAL"
    assert "withdraws to the" in out["detail"]


def test_a_sign_flip_on_an_UNRESOLVED_contrast_is_not_a_reversal():
    """A contrast whose interval contains zero has no established sign.

    Otherwise noise around zero in a small half would trip the strongest
    branch, which is the opposite of what the asymmetry is for.
    """
    noisy = _table(4, 3, means={(PRIMARY_HALVES[1], "GUCA2A - KRT8"): -0.02})
    assert verdict(noisy)["verdict"] != "SIGN REVERSAL"


def test_an_empty_load_bearing_table_says_so_rather_than_passing():
    assert verdict(pd.DataFrame(columns=["half", "statistic", "contrast",
                                         "mean", "excludes_zero"]))["verdict"] \
        == "NOT COMPUTED"


# ---------------------------------------------------------------------------
# The split itself
# ---------------------------------------------------------------------------


def test_a_patient_in_two_collections_is_excluded_from_both():
    """Assigning it to either puts the same patient on both sides."""
    split = pd.DataFrame({
        "patient_id": ["P1", "P2", "P3"],
        "collection": ["VUMC_HTAN_discovery", "VUMC_HTAN_validation",
                       "VUMC_HTAN_discovery|VUMC_HTAN_validation"],
    })
    out = assign_halves(split, ["P1", "P2", "P3"]).set_index("patient_id")
    assert not out.loc["P1", "shared"]
    assert out.loc["P3", "shared"]


def test_a_scored_patient_with_no_collection_stops_the_run():
    """An identifier-space mismatch, not a patient without a collection."""
    split = pd.DataFrame({"patient_id": ["P1"], "collection": ["VUMC_HTAN_discovery"]})
    with pytest.raises(StabilityError, match="identifier-space"):
        assign_halves(split, ["P1", "P_MISSING"])


def test_the_prefixed_identifier_form_also_resolves():
    """The atlas writes 'Chen_2021_Cell.HTA11_866', the decomposition 'HTA11_866'."""
    split = pd.DataFrame({"patient_id": ["Chen_2021_Cell.HTA11_866"],
                          "collection": ["VUMC_HTAN_discovery"]})
    out = assign_halves(split, ["HTA11_866"])
    assert out.loc[0, "collection"] == "VUMC_HTAN_discovery"


def test_the_width_penalty_quoted_in_the_prereg_is_real():
    """§6's 1.80x, re-derived rather than trusted."""
    from scipy import stats

    at = lambda n: stats.t.ppf(0.975, n - 1) / np.sqrt(n)  # noqa: E731
    assert at(15) / at(43) == pytest.approx(1.80, abs=0.01)
    assert at(13) / at(43) == pytest.approx(1.96, abs=0.01)
