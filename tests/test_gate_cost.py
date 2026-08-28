"""Re-costing the gate at the cohort sizes that exist. W2, handoff §5 task 3.

Promised to W4 twice and never done. The numbers here are the answer, pinned so
that "we re-costed it" stays true rather than becoming folklore.
"""

from __future__ import annotations

import pandas as pd
import pytest

from src.harness.gate_cost import (
    COHORT_SIZES,
    cohort_ci_width_by_n,
    effective_decision_line,
    g4_operating_characteristic,
    g4_over_rungs,
    g4_pass_probability,
    largest_clean_pass,
    matched_and_unmatched,
    widening_vs_plan,
)
from src.harness.positivity import (
    NON_IDENTIFIABILITY_HEADLINE_FRACTION,
    gate_g4_verdict,
    wilson_interval,
)

# ---------------------------------------------------------------------------
# G4's verdict, at the sizes that exist
# ---------------------------------------------------------------------------


def test_the_pre_committed_line_has_not_been_moved():
    """Re-costing a criterion is allowed. Moving it because n is small is not."""
    assert NON_IDENTIFIABILITY_HEADLINE_FRACTION == 0.50


def test_smc_at_ten_patients_cannot_resolve_g4_at_all():
    """The headline of the re-costing.

    A PASS the cohort can defend needs 1 or fewer of 10 patients below the
    threshold. At 2 the interval already straddles the 50% line, so G4 on SMC
    alone is not a question the cohort can answer.
    """
    assert largest_clean_pass(10) == 1
    assert effective_decision_line(10) == pytest.approx(0.10)


def test_gse178341_matched_still_resolves_decision_19s_verdict():
    """The good news, and the reason the re-costing had to be done rather than feared.

    Decision #19 recorded 6 of 36 matched patients below threshold — 16.7%. A
    defensible PASS at n=36 allows up to 12 (33.3%), so that verdict survives
    the shortfall from 60 with room to spare.
    """
    assert largest_clean_pass(36) == 12
    assert 6 <= largest_clean_pass(36)
    verdict = gate_g4_verdict([200] * 30 + [1] * 6, n_unmatched_patients=26)
    assert verdict["passes"] and verdict["resolvable"]


def test_the_effective_line_is_below_the_rule_at_every_real_cohort_size():
    """G4 says 50%; no cohort in this project can enforce 50%."""
    for n in COHORT_SIZES.values():
        assert effective_decision_line(n) < NON_IDENTIFIABILITY_HEADLINE_FRACTION


def test_g4_verdict_says_when_the_cohort_cannot_resolve_it():
    """A PASS and a PASS-it-could-not-have-contradicted must not read the same."""
    smc = gate_g4_verdict([200] * 8 + [1] * 2, n_unmatched_patients=0)
    assert smc["passes"] is True
    assert smc["resolvable"] is False
    assert "indeterminate" in smc["precision"]


def test_g4_verdict_keeps_reporting_the_pre_committed_consequence():
    """The added precision fields must not displace what the gate acts on."""
    failing = gate_g4_verdict([2, 3, 4, 200], n_unmatched_patients=0)
    assert not failing["passes"]
    assert "headline" in failing["consequence"]


@pytest.mark.parametrize("n", [10, 36, 60])
def test_pass_probability_is_certain_at_the_extremes(n):
    assert g4_pass_probability(n, 0.0) == pytest.approx(1.0)
    assert g4_pass_probability(n, 1.0) == pytest.approx(0.0)


def test_pass_probability_is_near_a_coin_flip_on_the_line():
    """At a true 50% the gate is a coin flip whatever n is — n does not fix that.

    More patients sharpen the verdict away from the line, not on it. Worth
    stating because "get more patients" is the reflex, and near the line it buys
    almost nothing.
    """
    for n in (10, 36, 60):
        assert 0.3 < g4_pass_probability(n, 0.50) < 0.5


def test_more_patients_sharpen_the_verdict_away_from_the_line():
    below = [g4_pass_probability(n, 0.35) for n in (10, 36, 60)]
    above = [g4_pass_probability(n, 0.60) for n in (10, 36, 60)]
    assert below == sorted(below)  # more n -> more likely to pass, correctly
    assert above == sorted(above, reverse=True)  # more n -> more likely to fail


def test_operating_characteristic_covers_every_cohort_and_fraction():
    table = g4_operating_characteristic()
    assert set(table["cohort"]) == set(COHORT_SIZES)
    assert table["p_gate_says_pass"].between(0.0, 1.0).all()


# ---------------------------------------------------------------------------
# Wilson, since the whole re-costing rests on it
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(("k", "n"), [(0, 10), (10, 10), (0, 1), (6, 36), (32, 62)])
def test_wilson_stays_inside_the_unit_interval(k, n):
    low, high = wilson_interval(k, n)
    assert 0.0 <= low <= high <= 1.0


def test_wilson_narrows_as_n_grows():
    widths = [
        wilson_interval(round(0.2 * n), n)[1] - wilson_interval(round(0.2 * n), n)[0]
        for n in (10, 36, 60, 200)
    ]
    assert widths == sorted(widths, reverse=True)


@pytest.mark.parametrize(("k", "n"), [(-1, 10), (11, 10), (5, 0)])
def test_wilson_rejects_impossible_counts(k, n):
    with pytest.raises(ValueError):
        wilson_interval(k, n)


# ---------------------------------------------------------------------------
# The cohort band on the terms
# ---------------------------------------------------------------------------


def test_the_cohort_band_widens_as_patients_are_lost():
    """Measured through W4's own ``bootstrap_over_patients``, not a formula.

    Cheap settings: the claim is the ordering, and the magnitudes live in the
    gate memo where they are quoted from a versioned table.
    """
    widths = cohort_ci_width_by_n(seed=99, n_boot=60, n_replicates=2)
    median = widening_vs_plan(widths)
    for term in median["term"].unique():
        by_n = median[median["term"] == term].sort_values("n_patients")
        assert list(by_n["ci_width"]) == sorted(by_n["ci_width"], reverse=True), term
        assert by_n.iloc[-1]["ratio_vs_plan"] == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# G4 against a real counts frame — the criterion's actual answer
# ---------------------------------------------------------------------------


def _counts_frame():
    """Two axes x two rungs x 3 patients, one of them single-armed."""
    rows = []
    for axis in ("stem_pole", "opposite_lineage"):
        for rung, mature in (("lineage", [200, 150, 4]), ("best4", [5, 6, 2])):
            for patient, n in zip(["P1", "P2", "P3"], mature, strict=True):
                for tissue in ("normal", "tumour"):
                    if patient == "P3" and tissue == "normal":
                        continue  # single-armed: must be excluded, not zeroed
                    rows.append(
                        {
                            "patient_id": patient,
                            "tissue": tissue,
                            "labeling_axis": axis,
                            "granularity_rung": rung,
                            "n_cells_mature": n if tissue == "tumour" else 500,
                            "n_cells_resolved": 900,
                        }
                    )
    return pd.DataFrame(rows)


def test_matched_is_defined_by_observed_arms_not_by_a_column():
    matched, unmatched = matched_and_unmatched(_counts_frame())
    assert matched == {"P1", "P2"}
    assert unmatched == {"P3"}


def test_g4_reports_one_verdict_per_rung_because_it_has_no_single_answer():
    """The mature population is defined by the rung, so the fraction below the
    cutpoint is too. One number would be picking a rung and calling it the
    cohort."""
    table = g4_over_rungs(_counts_frame())
    assert len(table) == 4
    assert set(table["granularity_rung"]) == {"lineage", "best4"}
    # The finest rung has the fewest mature cells, so it fails where lineage does not.
    by_rung = table.groupby("granularity_rung", observed=True)["fraction_below"].max()
    assert by_rung["best4"] > by_rung["lineage"]


def test_g4_excludes_the_single_armed_patient_from_the_fraction():
    table = g4_over_rungs(_counts_frame())
    assert (table["n_patients"] == 2).all()
    assert (table["n_unmatched_excluded"] == 1).all()


def test_g4_counts_each_patient_once():
    """A patient counted twice is invariant 5's mistake one level up — the
    fourth instance of this repo's bug family (issue #36)."""
    doubled = pd.concat([_counts_frame(), _counts_frame()], ignore_index=True)
    assert (g4_over_rungs(doubled)["n_patients"] == 2).all()


def test_g4_over_rungs_names_a_missing_column():
    with pytest.raises(ValueError, match="granularity_rung"):
        g4_over_rungs(_counts_frame().drop(columns=["granularity_rung"]))


def test_g4_refuses_a_cohort_with_no_matched_patient():
    frame = _counts_frame()
    frame = frame[frame["tissue"] == "tumour"]
    with pytest.raises(ValueError, match="both arms"):
        g4_over_rungs(frame)
