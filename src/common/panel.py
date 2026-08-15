"""Loader for the frozen panel and labelling axes.

Read these from here. Do not retype gene lists into your own module — a second
copy is a second thing to forget to update, and the freeze exists precisely so
that nothing gets quietly updated.
"""

from __future__ import annotations

from functools import lru_cache
from typing import Any

import yaml

from src.common.paths import CONFIG_DIR

PANEL_PATH = CONFIG_DIR / "panel.yaml"
AXES_PATH = CONFIG_DIR / "labeling_axes.yaml"


@lru_cache(maxsize=1)
def load_panel() -> dict[str, Any]:
    return yaml.safe_load(PANEL_PATH.read_text(encoding="utf-8"))


@lru_cache(maxsize=1)
def load_axes() -> dict[str, Any]:
    return yaml.safe_load(AXES_PATH.read_text(encoding="utf-8"))


def tier_genes(tier: str) -> list[str]:
    """Genes in one tier, e.g. ``tier_genes("B") -> ["MLH1", "SFRP1", "SFRP2"]``."""
    tiers = load_panel()["tiers"]
    key = tier.upper()
    if key not in tiers:
        raise KeyError(f"unknown tier {tier!r}; frozen tiers are {sorted(tiers)}")
    return list(tiers[key]["genes"])


def panel_genes(tiers: tuple[str, ...] | None = None) -> list[str]:
    """All panel genes, or just the named tiers, deduplicated, order preserved."""
    wanted = tuple(t.upper() for t in tiers) if tiers else tuple(load_panel()["tiers"])
    seen: dict[str, None] = {}
    for t in wanted:
        for g in tier_genes(t):
            seen.setdefault(g, None)
    return list(seen)


def tier_of(gene: str) -> str | None:
    """Which tier a gene belongs to, or None if it is not on the panel."""
    for tier, spec in load_panel()["tiers"].items():
        if gene in spec["genes"]:
            return tier
    return None


def tier_expectation(tier: str) -> str:
    """The pre-registered expectation: compositional / intrinsic / mixed / neither.

    The falsification rule (execution_plan.md §3.1) compares tiers A, B and D
    against these. Read it from here so the expectation cannot drift after
    seeing results.
    """
    return load_panel()["tiers"][tier.upper()]["expectation"]


def wnt_targets(include_stem_overlapping: bool = False) -> list[str]:
    """Wnt-target signature. CLAUDE.md invariant 8.

    Defaults to the core five (AXIN2, NKD1, RNF43, NOTUM, TCF7). ASCL2 and LGR5
    are excluded unless explicitly requested, because they are the stem axis and
    including them makes "high mature fraction + high Wnt" near-contradictory.
    """
    spec = load_panel()["wnt_targets"]
    genes = list(spec["core"])
    if include_stem_overlapping:
        genes += list(spec["stem_overlapping"])
    return genes


def axis_genes(axis: str) -> list[str]:
    """Marker genes for a labelling axis. Empty for the non-transcript axes."""
    axes = load_axes()["axes"]
    if axis not in axes:
        raise KeyError(f"unknown axis {axis!r}; frozen axes are {sorted(axes)}")
    return list(axes[axis]["genes"])


def granularity_rungs() -> list[str]:
    return [r["name"] for r in load_axes()["granularity_rungs"]]
