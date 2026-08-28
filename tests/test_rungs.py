"""Can the estimator tell two granularity rungs apart?

W1's pilot found lineage and crypt_position are the same partition on axis 1.
Two explanations: the partitions really are identical, or the estimator cannot
resolve a difference that is there. These tests exclude the second, so the
finding can be reported as a statement about the labelling.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.harness.pseudobulk import generate_pseudobulk
from src.harness.rungs import (
    estimate_under_rung,
    estimator_can_separate,
    rung_mature_mask,
    rung_separation,
    separation_summary,
)

TARGET = "GUCA2A"
#: A nested structure: the coarse rung pools crypt_top and crypt_mid into
#: "mature"; the fine rung calls only crypt_top mature. The two subpopulations
#: express the target very differently, so the rungs genuinely disagree.
CRYPT_TOP, CRYPT_MID = "crypt_top", "crypt_mid"
TYPES = [CRYPT_TOP, CRYPT_MID, "stem", "stromal", "immune"]
GENES = [TARGET, "ACTB", "F1", "F2", "F3", "F4"]

PROFILE = {
    CRYPT_TOP: [80.0, 50.0, 5.0, 5.0, 5.0, 5.0],   # high target
    CRYPT_MID: [8.0, 50.0, 5.0, 5.0, 5.0, 5.0],    # low target — pooling dilutes
    "stem": [1.0, 50.0, 5.0, 5.0, 5.0, 5.0],
    "stromal": [1.0, 55.0, 6.0, 4.0, 5.0, 5.0],
    "immune": [1.0, 45.0, 4.0, 6.0, 5.0, 5.0],
}

COARSE = [CRYPT_TOP, CRYPT_MID]   # lineage-like: pools the two
FINE = [CRYPT_TOP]                # crypt_position-like: the top only


@pytest.fixture(scope="module")
def cohort():
    rng = np.random.default_rng(4)
    rows, ctypes, patients = [], [], []
    for p in range(10):
        for t in TYPES:
            for _ in range(60):
                rows.append(rng.poisson(PROFILE[t]))
                ctypes.append(t)
                patients.append(f"P{p:02d}")
    return np.array(rows), np.array(ctypes), np.array(patients)


def _comp(top, mid):
    rest = (1.0 - top - mid) / 3
    return {CRYPT_TOP: top, CRYPT_MID: mid, "stem": rest, "stromal": rest, "immune": rest}


@pytest.fixture(scope="module")
def sample(cohort):
    counts, ctypes, patients = cohort
    return generate_pseudobulk(
        counts, ctypes, patients, GENES,
        composition_normal=_comp(0.20, 0.20),
        composition_tumour=_comp(0.05, 0.20),   # the top compartment is depleted
        shift={TARGET: 0.5},
        held_out_patients=["P08", "P09"],
        n_cells=2000, seed=1, mature_label=CRYPT_TOP,
    )


# ---------------------------------------------------------------------------
# masks and single-rung estimates
# ---------------------------------------------------------------------------


def test_mask_selects_exactly_the_named_types(sample):
    coarse = rung_mature_mask(sample, "tumour", COARSE)
    fine = rung_mature_mask(sample, "tumour", FINE)
    assert coarse.sum() > fine.sum()
    assert np.all(fine <= coarse)  # fine is nested inside coarse


def test_mask_rejects_an_arm_that_does_not_exist(sample):
    with pytest.raises(KeyError, match="drawn_cell_type"):
        rung_mature_mask(sample, "border", COARSE)


def test_estimate_under_rung_reports_the_compartment_it_used(sample):
    fine = estimate_under_rung(sample, TARGET, FINE)
    coarse = estimate_under_rung(sample, TARGET, COARSE)
    assert coarse["n_cells_mature"] > fine["n_cells_mature"]
    # Pooling a low-expressing subpopulation must drag the per-cell mean down.
    assert coarse["mean_normal"] < fine["mean_normal"]


def test_estimate_under_rung_rejects_an_unknown_gene(sample):
    with pytest.raises(KeyError, match="drawn expression"):
        estimate_under_rung(sample, "NOTAGENE", FINE)


# ---------------------------------------------------------------------------
# THE question: does the estimator separate rungs that genuinely differ?
# ---------------------------------------------------------------------------


def test_estimator_separates_rungs_that_genuinely_differ(sample):
    """Excludes "the estimator cannot resolve it" as an explanation for W1's
    observed degeneracy."""
    sep = rung_separation(sample, TARGET, {"lineage": COARSE, "crypt_position": FINE})
    assert len(sep) == 2
    assert estimator_can_separate(sep)

    summary = separation_summary(sep)
    assert summary["relative"].iloc[0] > 0.10
    assert not summary["identical_partition"].iloc[0]


def test_identical_partitions_give_identical_estimates(sample):
    """The degeneracy case. Two rungs naming the same cell types must agree
    exactly — so when W1 sees agreement, the labels explain it fully."""
    sep = rung_separation(sample, TARGET, {"lineage": COARSE, "crypt_position": list(COARSE)})
    summary = separation_summary(sep)
    assert summary["absolute"].iloc[0] == pytest.approx(0.0)
    assert summary["identical_partition"].iloc[0]
    assert not estimator_can_separate(sep)


def test_separation_needs_two_rungs(sample):
    with pytest.raises(ValueError, match="at least two rungs"):
        rung_separation(sample, TARGET, {"lineage": COARSE})


def test_separation_is_none_when_a_rung_is_not_estimable(cohort):
    """Not comparable, rather than compared against zero. Invariant 1 reaches
    the rung comparison too."""
    counts, ctypes, patients = cohort
    starved = generate_pseudobulk(
        counts, ctypes, patients, GENES,
        composition_normal=_comp(0.20, 0.20),
        composition_tumour=_comp(0.0, 0.20),   # no crypt_top left at all
        shift={TARGET: 0.5},
        held_out_patients=["P08", "P09"],
        n_cells=2000, seed=2, mature_label=CRYPT_TOP,
    )
    sep = rung_separation(starved, TARGET, {"lineage": COARSE, "crypt_position": FINE})
    fine = sep[sep["rung"] == "crypt_position"].iloc[0]
    assert fine["estimability"] == "not_estimable"
    # <NA> in the frame, None from the dict API. What matters is that neither
    # is 0.0 — an unestimable intrinsic term reported as zero is the single
    # most likely route to a wrong conclusion in this project.
    assert fine["intrinsic"] is pd.NA
    assert estimate_under_rung(starved, TARGET, FINE)["intrinsic"] is None
    # And the compositional term survives — the row is not dropped.
    assert not pd.isna(fine["compositional"])

    summary = separation_summary(sep)
    assert summary["absolute"].iloc[0] is None
    assert not estimator_can_separate(sep)


def test_compositional_term_also_moves_with_the_rung(sample):
    """The granularity choice is not only an intrinsic-term story — §6.2 asks
    whether the whole split swings."""
    sep = rung_separation(sample, TARGET, {"lineage": COARSE, "crypt_position": FINE})
    summary = separation_summary(sep, term="compositional")
    assert summary["relative"].iloc[0] > 0.10
