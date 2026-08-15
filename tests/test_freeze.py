"""The freeze, asserted.

CLAUDE.md invariant 3: the panel and the labelling axes are frozen. Changing
either needs a PR, two approvals and a written reason.

The gene lists below are a deliberate second copy of config/*.yaml. Editing the
config alone turns CI red — that is the friction working. Change both, in one
PR, with the reason in the PR body.
"""

from __future__ import annotations

import pytest

from src.common.panel import (
    axis_genes,
    granularity_rungs,
    load_axes,
    load_panel,
    panel_genes,
    tier_expectation,
    tier_genes,
    tier_of,
    wnt_targets,
)

# ---------------------------------------------------------------------------
# Panel — execution_plan.md §3.1
# ---------------------------------------------------------------------------

FROZEN_PANEL = {
    "A": ["GUCA2A", "GUCA2B", "OTOP2", "CA7"],
    "B": ["MLH1", "SFRP1", "SFRP2"],
    "C": ["CDX2"],
    "D": ["MS4A12"],
    "E": [
        "AQP8", "CA1", "CA2", "CA4", "CEACAM7", "SLC26A3", "FABP1",
        "PIGR", "KRT20", "LGALS4", "VIL1", "SATB2", "MUC2", "TFF3",
    ],
}

FROZEN_EXPECTATIONS = {
    "A": "compositional",
    "B": "intrinsic",
    "C": "mixed",
    "D": "neither",
    "E": "unspecified",
}


@pytest.mark.parametrize("tier", sorted(FROZEN_PANEL))
def test_panel_tier_is_frozen(tier):
    assert tier_genes(tier) == FROZEN_PANEL[tier]


@pytest.mark.parametrize("tier", sorted(FROZEN_EXPECTATIONS))
def test_tier_expectations_are_frozen(tier):
    """The falsification rule compares tiers A, B and D against these.

    They are read from config, not retyped after seeing results.
    """
    assert tier_expectation(tier) == FROZEN_EXPECTATIONS[tier]


def test_panel_has_no_extra_tiers():
    assert sorted(load_panel()["tiers"]) == ["A", "B", "C", "D", "E"]


def test_falsification_rule_has_three_distinct_expectations():
    """A/B/D must expect different answers or the rule cannot fire."""
    assert len({tier_expectation(t) for t in ("A", "B", "D")}) == 3


def test_tier_b_carries_mlh1():
    """Tier B is load-bearing. A compositional control alone validates only
    the half nobody doubted."""
    assert "MLH1" in tier_genes("B")
    assert tier_of("MLH1") == "B"


def test_panel_genes_are_deduplicated():
    all_genes = panel_genes()
    assert len(all_genes) == len(set(all_genes))


# ---------------------------------------------------------------------------
# Wnt — CLAUDE.md invariant 8
# ---------------------------------------------------------------------------


def test_wnt_score_is_targets_not_ctnnb1():
    """CTNNB1 / TCF7L2 transcript level is not Wnt activity — activation is
    post-translational."""
    core = wnt_targets()
    assert core == ["AXIN2", "NKD1", "RNF43", "NOTUM", "TCF7"]
    assert "CTNNB1" not in core
    assert "TCF7L2" not in core


def test_stem_overlapping_wnt_genes_are_opt_in():
    """ASCL2 and LGR5 are the stem axis. Including them by default makes
    'high mature fraction + high Wnt' near-contradictory (§2.3)."""
    assert "ASCL2" not in wnt_targets()
    assert "LGR5" not in wnt_targets()
    assert set(wnt_targets(include_stem_overlapping=True)) >= {"ASCL2", "LGR5"}


# ---------------------------------------------------------------------------
# Labelling axes — execution_plan.md §3.2
# ---------------------------------------------------------------------------

FROZEN_AXES = {
    "stem_pole": ["LGR5", "ASCL2", "MKI67", "OLFM4", "SMOC2"],
    "opposite_lineage": ["MUC2", "TFF3", "SPDEF", "ITLN1"],
    "chromatin": [],
    "spatial": [],
}


@pytest.mark.parametrize("axis", sorted(FROZEN_AXES))
def test_axis_is_frozen(axis):
    assert axis_genes(axis) == FROZEN_AXES[axis]


def test_a_non_transcript_axis_exists():
    """Axis 3 is the strongest defence against label leakage precisely because
    it is not made of transcripts."""
    axes = load_axes()["axes"]
    assert any(not spec["transcript_based"] for spec in axes.values())


def test_granularity_rungs_are_frozen():
    """Four rungs. The split is reported as a curve across them — a single
    point estimate would present a modelling choice as a measurement."""
    assert granularity_rungs() == ["epithelial", "lineage", "crypt_position", "best4"]


def test_axis_and_schema_vocabularies_agree():
    from src.schema import GRANULARITY_RUNGS, LABELING_AXES

    assert tuple(granularity_rungs()) == GRANULARITY_RUNGS
    assert tuple(load_axes()["axes"]) == LABELING_AXES


# ---------------------------------------------------------------------------
# The known collision — docs/open_decisions.md #1
# ---------------------------------------------------------------------------


def test_known_target_label_collision_is_still_the_documented_one():
    """MUC2 and TFF3 are both tier-E targets and axis-2 labels.

    This test does not forbid the collision — it pins it, so that if a future
    panel or axis edit introduces a NEW one, someone has to look at this test
    and think about invariant 2 rather than discovering it in week 9.
    """
    collisions = set(panel_genes()) & set(axis_genes("opposite_lineage"))
    assert collisions == {"MUC2", "TFF3"}, (
        f"target/label collision set changed to {sorted(collisions)}. "
        f"Invariant 2: target genes never appear in labels. Resolve it in "
        f"docs/open_decisions.md before updating this test."
    )


def test_stem_pole_axis_does_not_collide_with_the_panel():
    assert not set(panel_genes()) & set(axis_genes("stem_pole"))
