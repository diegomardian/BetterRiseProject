"""The committed Lee driver — depth matching semantics. W4.

The Lee tables the project quotes were produced by a driver that was never
committed. This file pins the three choices that driver has to make the same
way W1's `build_decomposition_summary` makes them, because a meta-analysis
across the two cohorts is only meaningful if both arms were built the same way:

  1. matched WITHIN patient,
  2. AFTER labelling,
  3. on the SCORED epithelium, with an unmatchable patient contributing nothing.

If a test here goes red, the two cohorts stopped being comparable.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.common.panel import granularity_rungs
from src.estimator.lee_io import LeeCohort
from src.estimator.run_lee_decomposition import (
    _scorable,
    _subset,
    depth_matched_index,
)
from src.reference.labels import (
    NON_EPITHELIAL,
    RUNG_SPECS,
    TRANSCRIPT_AXES,
    UNRESOLVED,
    label_column,
)

SEED = 20260829


def _cohort(
    *,
    n_per_arm=200,
    n_patients=2,
    normal_depth=3_000.0,
    tumour_depth=7_200.0,
    sigma=0.8,
    unresolved_frac=0.0,
    one_armed: tuple[str, ...] = (),
    seed=0,
):
    """A synthetic Lee cohort whose arms differ in depth by construction.

    Depth is LOGNORMAL, which is what sequencing depth actually looks like and
    — more to the point — what makes the arms OVERLAP. `match_arm_depth` bins
    on pooled quantiles and keeps min(n_normal, n_tumour) per bin, so two arms
    with disjoint support have nothing to match and it correctly retains
    nothing. A Gaussian fixture at this separation is disjoint in practice and
    tests the refusal path instead of the matching path; `sigma` is what keeps
    the supports overlapping. Real Lee arms sit ~2.4x apart in median depth
    with heavy overlap, which is what these defaults reproduce.

    Cells marked unresolved carry UNRESOLVED in every label column and must
    never be matched, scored, or counted as immature.
    """
    rng = np.random.default_rng(seed)
    rows, depths, labels_first = [], [], []
    for p in range(n_patients):
        patient = f"P{p}"
        arms = ("tumour",) if patient in one_armed else ("normal", "tumour")
        for arm in arms:
            median = normal_depth if arm == "normal" else tumour_depth
            for _ in range(n_per_arm):
                rows.append({"patient_id": patient, "tissue": arm})
                depths.append(float(rng.lognormal(np.log(median), sigma)))
                labels_first.append(
                    UNRESOLVED if rng.random() < unresolved_frac else "mature_colonocyte"
                )

    index = pd.Index([f"c{i}" for i in range(len(rows))], name="cell")
    cells = pd.DataFrame(rows, index=index)
    cells["study_id"] = "GSE132465"
    cells["n_counts"] = np.clip(depths, 1.0, None)

    labels = pd.DataFrame(index=index)
    for axis in TRANSCRIPT_AXES:
        for rung in granularity_rungs():
            labels[label_column(axis, rung)] = labels_first
            labels[f"mature__{axis}__{rung}"] = pd.Series(
                [v != UNRESOLVED for v in labels_first], index=index, dtype="boolean"
            )
    labels["patient_id"] = cells["patient_id"].to_numpy()
    labels["tissue"] = cells["tissue"].to_numpy()

    return LeeCohort(
        study_id="GSE132465",
        cells=cells,
        expression=pd.DataFrame(index=index),
        labels=labels,
        axis_gene_coverage={"found": [], "missing": [], "requested": []},
    )


def _ratio(cohort, index=None):
    cells = cohort.cells if index is None else cohort.cells.loc[index]
    med = cells.groupby("tissue")["n_counts"].median()
    return float(med["tumour"] / med["normal"])


# ---------------------------------------------------------------------------
# The positive control, and the negative control that makes it mean something
# ---------------------------------------------------------------------------


def test_depth_matching_equalises_the_arms_rather_than_merely_thinning_them():
    """~2.4x apart before, ~1x after. Thinning both arms keeps the ratio; only
    matching removes it. This is the whole point of decision #24.1."""
    cohort = _cohort(normal_depth=3_000.0, tumour_depth=7_200.0)
    assert _ratio(cohort) == pytest.approx(2.4, rel=0.3), "fixture no longer confounded"

    keep, counts = depth_matched_index(cohort, seed=SEED)
    assert _ratio(cohort, keep) == pytest.approx(1.0, abs=0.25)
    assert counts["n_scorable_after_matching"] < counts["n_scorable_before_matching"]


def test_matching_leaves_an_already_balanced_cohort_close_to_intact():
    """The negative control. If matching also gutted a cohort that needed
    nothing, the cell loss above would not be evidence of anything."""
    cohort = _cohort(normal_depth=5_000.0, tumour_depth=5_000.0)
    _keep, counts = depth_matched_index(cohort, seed=SEED)
    retained = counts["n_scorable_after_matching"] / counts["n_scorable_before_matching"]
    assert retained > 0.8


def test_arms_with_disjoint_depth_cannot_be_matched_and_say_so():
    """Not a cohort-size failure -- a support failure, and the one case where
    'retained nothing' is the correct answer rather than a bug. Two arms that
    never overlap in depth have no common distribution to match onto, so any
    comparison between them is confounded and cannot be repaired by
    subsampling."""
    cohort = _cohort(normal_depth=2_000.0, tumour_depth=80_000.0, sigma=0.05)
    with pytest.raises(SystemExit, match="retained 0 of"):
        depth_matched_index(cohort, seed=SEED)


# ---------------------------------------------------------------------------
# The three semantics that have to match W1
# ---------------------------------------------------------------------------


def test_a_patient_without_two_scorable_arms_contributes_nothing():
    """An unmatchable patient is not a matched patient. Letting it through is
    how a depth confound survives a matched read."""
    cohort = _cohort(n_patients=3, one_armed=("P2",))
    keep, _counts = depth_matched_index(cohort, seed=SEED)
    assert "P2" not in set(cohort.cells.loc[keep, "patient_id"])
    assert {"P0", "P1"} == set(cohort.cells.loc[keep, "patient_id"])


def test_matching_is_within_patient_so_each_patient_keeps_both_arms():
    cohort = _cohort(n_patients=3)
    keep, _ = depth_matched_index(cohort, seed=SEED)
    kept = cohort.cells.loc[keep]
    for patient, block in kept.groupby("patient_id"):
        assert set(block["tissue"]) == {"normal", "tumour"}, patient


def test_unresolved_cells_are_never_matched_because_they_are_not_scored():
    cohort = _cohort(unresolved_frac=0.5)
    keep, counts = depth_matched_index(cohort, seed=SEED)
    first = cohort.labels.loc[keep, label_column(TRANSCRIPT_AXES[0], granularity_rungs()[0])]
    assert (first.astype(str) != UNRESOLVED).all()
    assert counts["n_scorable_before_matching"] < len(cohort.cells)


def test_the_scorable_mask_does_not_depend_on_axis_or_rung():
    """Matching on one axis's mask would give each rung a different subsample
    and stop the rungs being comparable to each other."""
    cohort = _cohort(unresolved_frac=0.3)
    base = _scorable(cohort.labels)
    for axis in TRANSCRIPT_AXES:
        for rung in granularity_rungs():
            values = cohort.labels[label_column(axis, rung)].astype(str).to_numpy()
            mask = (values != UNRESOLVED) & (values != NON_EPITHELIAL)
            assert np.array_equal(base, mask), f"{axis}/{rung} disagrees"


# ---------------------------------------------------------------------------
# Refusals, reproducibility, and the subset helper
# ---------------------------------------------------------------------------


def test_matching_refuses_loudly_when_every_bin_empties_an_arm():
    """A cohort smaller than the bin count cannot be matched. It must say so,
    not return an empty set that surfaces later as 'no summary rows'."""
    cohort = _cohort(n_per_arm=3, n_patients=1)
    with pytest.raises(SystemExit, match="retained 0 of"):
        depth_matched_index(cohort, seed=SEED)


def test_matching_refuses_when_no_patient_has_two_arms():
    cohort = _cohort(n_patients=1, one_armed=("P0",))
    with pytest.raises(SystemExit, match="two scorable arms"):
        depth_matched_index(cohort, seed=SEED)


def test_matching_is_reproducible_given_the_same_seed():
    cohort = _cohort()
    a, _ = depth_matched_index(cohort, seed=SEED)
    b, _ = depth_matched_index(cohort, seed=SEED)
    assert a.equals(b)


def test_subset_keeps_every_frame_in_step():
    cohort = _cohort()
    keep, _ = depth_matched_index(cohort, seed=SEED)
    subset = _subset(cohort, keep)
    assert subset.cells.index.equals(keep)
    assert subset.labels.index.equals(keep)
    assert subset.study_id == cohort.study_id


def test_the_mature_call_is_not_recomputed_by_subsetting():
    """Matching happens AFTER labelling. Subsetting must carry the existing
    calls, never re-derive them -- re-deriving would move each patient's cut
    points and turn the matched read into a different measurement."""
    cohort = _cohort()
    keep, _ = depth_matched_index(cohort, seed=SEED)
    subset = _subset(cohort, keep)
    col = f"mature__{TRANSCRIPT_AXES[0]}__{granularity_rungs()[0]}"
    pd.testing.assert_series_equal(subset.labels[col], cohort.labels.loc[keep, col])


def test_rung_specs_still_names_a_mature_value_for_every_rung():
    """The fixture writes 'mature_colonocyte' as a scored call. If RUNG_SPECS
    stopped agreeing, these tests would pass while measuring nothing."""
    for rung in granularity_rungs():
        assert RUNG_SPECS[rung].mature
