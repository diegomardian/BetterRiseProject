"""Is the maturity call measuring maturity, or sequencing depth? W2.

Written after the first real-data run, where it turned out to be depth
(issue #44). On Lee/SMC, normal epithelium is called 71.0% mature and tumour
25.1% — a 46-point apparent compositional loss in the hypothesised direction —
because normal epithelium was sequenced 4.3x shallower, dropout put 32% of
epithelial cells at zero for every stem marker, and the maturity axis is
*inverted*, so absence of the marker is the top of the scale.

WHY THIS LIVES IN THE HARNESS
-----------------------------
The harness is the only place ground truth exists, and this is a question about
whether a measurement measures what it claims. It is not W2's job to fix another
workstream's labeller; it *is* W2's job to make the failure detectable by
something other than one person happening to look.

Nothing here is a gate criterion. G1-G4 are pre-registered and this is not among
them, so these are diagnostics with a stated reading, to be quoted in the gate
memo beside the decomposition rather than used to pass or fail anything. Adding a
gate after seeing a result is the move this project refuses.

WHAT A CLEAN LABELLER LOOKS LIKE HERE
-------------------------------------
`src/reference/labels.py` (W1's) thins marker counts to a common depth and marks
cells below the target ``unresolved`` — *"not scored, not counted as immature"* —
then depth-matches the population that sets the cut points too. Run this against
labels built that way and ``maturity_tracks_depth`` should be false. That is the
comparison this module exists to make cheap.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Final

import numpy as np
import pandas as pd

#: Arms this far apart in median depth cannot be compared on a dropout-sensitive
#: score without saying so. 1.5x is deliberately lax — Lee/SMC is 4.3x, and a bar
#: that only catches the egregious case is the one worth having when the
#: alternative is no bar at all.
DEPTH_RATIO_TOLERANCE: Final = 1.5

#: |Spearman| between depth and the maturity call, *within* an arm. Within, not
#: pooled: pooling lets a genuine between-arm difference in maturity masquerade
#: as a depth effect and vice versa. This is the same "compute the statistic
#: inside the group" rule that decision #19 turned on.
MATURITY_DEPTH_RHO_TOLERANCE: Final = 0.20

#: Deciles for the dose-response table. Dropout is monotone in depth, so the
#: shape across bins is more convincing than any single correlation.
N_DEPTH_BINS: Final = 10


def mature_share_by_depth(
    depth: np.ndarray, is_mature: np.ndarray, *, n_bins: int = N_DEPTH_BINS
) -> pd.DataFrame:
    """Share called mature in each depth bin. The dose-response table.

    A labeller that measures biology should show no trend here. One that is
    reading dropout shows a monotone one — on Lee/SMC the share ran from 80.8%
    in the shallowest decile to 8.5% in the deepest.
    """
    depth = np.asarray(depth, dtype=float)
    is_mature = np.asarray(is_mature, dtype=bool)
    if depth.shape != is_mature.shape:
        raise ValueError(f"depth has {depth.shape}, is_mature has {is_mature.shape}")
    if depth.size == 0:
        raise ValueError("no cells")

    bins = pd.qcut(pd.Series(depth), n_bins, labels=False, duplicates="drop")
    frame = pd.DataFrame({"bin": bins, "depth": depth, "is_mature": is_mature})
    return (
        frame.groupby("bin", observed=True)
        .agg(
            n_cells=("is_mature", "size"),
            median_depth=("depth", "median"),
            mature_share=("is_mature", "mean"),
        )
        .reset_index()
    )


def max_attainable_rho(prevalence: float) -> float:
    """The largest |rho| a BINARY label can reach, at this prevalence.

    ``sqrt(3p(1-p))``. Spearman between a continuous variable and a binary one is
    bounded by the binary's own variance, so a rare label cannot correlate
    strongly with anything however completely it is determined:

    ======  ==============
    p       max |rho|
    ======  ==============
    0.50    0.866
    0.05    0.377
    0.0137  0.200
    0.0014  0.066
    ======  ==============

    **Below p = 1.37% the 0.20 tolerance is mathematically unreachable**, so
    ``maturity_tracks_depth`` cannot fire there no matter how completely depth
    determines the label. On Lee/SMC the ``best4`` rung's tumour arm has 8 mature
    cells in 5,564, capping |rho| at 0.066 — its "clean" reading was not evidence
    of cleanliness, it was the ceiling.

    It also makes raw |rho| incomparable ACROSS rungs whose prevalence differs by
    orders of magnitude. Compare ``rho_vs_ceiling`` instead, which is what
    :func:`depth_confound_report` reports alongside.
    """
    if not 0.0 <= prevalence <= 1.0:
        raise ValueError(f"prevalence={prevalence} outside [0, 1]")
    return float(np.sqrt(3.0 * prevalence * (1.0 - prevalence)))


def _spearman(x: np.ndarray, y: np.ndarray) -> float:
    """Rank correlation, without pulling scipy into a diagnostic."""
    if len(x) < 3 or len(np.unique(y)) < 2 or len(np.unique(x)) < 2:
        return float("nan")
    rx = pd.Series(x).rank().to_numpy()
    ry = pd.Series(y).rank().to_numpy()
    return float(np.corrcoef(rx, ry)[0, 1])


def depth_confound_report(
    depth: np.ndarray,
    is_mature: np.ndarray,
    arm: np.ndarray,
    *,
    depth_ratio_tolerance: float = DEPTH_RATIO_TOLERANCE,
    rho_tolerance: float = MATURITY_DEPTH_RHO_TOLERANCE,
) -> dict[str, object]:
    """Whether a maturity call can be told apart from sequencing depth.

    ``arm`` is the tumour/normal label per cell. Two things have to go wrong
    together for the compositional term to be manufactured, and the report keeps
    them separate because the fixes differ:

    - **``maturity_tracks_depth``** — within an arm, deeper cells get a
      different call. That is the labeller reading dropout.
    - **``arms_are_depth_matched``** — the two arms were sequenced comparably.
      If they were not, any depth sensitivity is converted into a *difference*
      between arms, which is exactly the compositional term.

    Either alone is a caveat. **Both together means the compositional term and
    the depth imbalance are not separable in this data**, and a decomposition
    from it should not be quoted without saying so.

    Returns the numbers, not a verdict on the science. ``confounded`` is a
    conjunction of the two stated conditions and nothing more.
    """
    depth = np.asarray(depth, dtype=float)
    is_mature = np.asarray(is_mature, dtype=bool)
    arm = np.asarray(arm)
    if not (depth.shape == is_mature.shape == arm.shape):
        raise ValueError(
            f"depth {depth.shape}, is_mature {is_mature.shape}, arm {arm.shape} "
            f"must be the same length — one entry per cell"
        )
    if depth.size == 0:
        raise ValueError("no cells")
    if np.any(depth < 0):
        raise ValueError("depth has negative entries")

    arms = sorted(set(arm.tolist()))
    per_arm: dict[str, dict[str, float]] = {}
    for a in arms:
        m = arm == a
        per_arm[str(a)] = {
            "n_cells": int(m.sum()),
            "median_depth": float(np.median(depth[m])) if m.any() else float("nan"),
            "mature_share": float(is_mature[m].mean()) if m.any() else float("nan"),
            "rho_depth_vs_mature": _spearman(depth[m], is_mature[m].astype(float)),
            "max_attainable_rho": max_attainable_rho(float(is_mature[m].mean()))
            if m.any()
            else float("nan"),
        }

    medians = [v["median_depth"] for v in per_arm.values() if np.isfinite(v["median_depth"])]
    if len(medians) >= 2 and min(medians) > 0:
        depth_ratio = max(medians) / min(medians)
    else:
        depth_ratio = float("nan")

    # Each arm's correlation belongs with THAT ARM's bound. Taking the two
    # maxima independently — max(rhos) over here, max(ceilings) over there —
    # lets the statistic come from one arm and its ceiling from the other, so a
    # rare arm's unreachable reading gets normalised by a common arm's headroom
    # and disappears. That mispairing is the one this diagnostic exists to
    # catch, and it was live in this function until it got caught.
    for v in per_arm.values():
        rho_a, ceil_a = abs(v["rho_depth_vs_mature"]), v["max_attainable_rho"]
        v["rho_vs_ceiling"] = (
            float(rho_a / ceil_a)
            if np.isfinite(rho_a) and np.isfinite(ceil_a) and ceil_a > 0
            else float("nan")
        )
        v["tolerance_is_reachable"] = bool(
            np.isfinite(ceil_a) and ceil_a >= rho_tolerance
        )

    # An arm at prevalence exactly 0 or 1 has no correlation to report. That is
    # UNDEFINED, not clean, and it does not get to vote by being dropped: it is
    # counted and named so the caller can abstain rather than read the survivors
    # as the whole picture.
    defined = {
        a: v for a, v in per_arm.items() if np.isfinite(v["rho_depth_vs_mature"])
    }
    undefined_arms = sorted(set(per_arm) - set(defined))

    if defined:
        worst_arm = max(defined, key=lambda a: abs(defined[a]["rho_depth_vs_mature"]))
        worst = defined[worst_arm]
        worst_rho = abs(worst["rho_depth_vs_mature"])
        ceiling = worst["max_attainable_rho"]
        rho_vs_ceiling = worst["rho_vs_ceiling"]
    else:
        worst_arm = None
        worst_rho = ceiling = rho_vs_ceiling = float("nan")

    # A rare label cannot reach the tolerance however completely depth drives
    # it, so report whether the test was capable of firing — for the arm the
    # verdict actually rests on, and separately for every arm.
    reachable = bool(np.isfinite(ceiling) and ceiling >= rho_tolerance)
    reachable_arms = sum(v["tolerance_is_reachable"] for v in per_arm.values())

    tracks = bool(np.isfinite(worst_rho) and worst_rho >= rho_tolerance)
    matched = bool(np.isfinite(depth_ratio) and depth_ratio <= depth_ratio_tolerance)
    return {
        "per_arm": per_arm,
        "depth_ratio_between_arms": depth_ratio,
        "worst_within_arm_rho": worst_rho,
        "worst_arm": worst_arm,
        "max_attainable_rho": ceiling,
        "rho_vs_ceiling": rho_vs_ceiling,
        "tolerance_is_reachable": reachable,
        "n_arms_tolerance_reachable": int(reachable_arms),
        "n_arms": len(per_arm),
        "arms_with_undefined_rho": undefined_arms,
        "maturity_tracks_depth": tracks,
        "arms_are_depth_matched": matched,
        "confounded": bool(tracks and not matched),
        "reading": _reading(
            tracks, matched, depth_ratio, worst_rho, reachable, undefined_arms
        ),
    }


def _reading(
    tracks: bool,
    matched: bool,
    ratio: float,
    rho: float,
    reachable: bool = True,
    undefined_arms: Sequence[str] | None = None,
) -> str:
    if undefined_arms is None:
        undefined_arms = []
    if not np.isfinite(rho):
        return (
            f"UNDEFINED: no arm carries a correlation — every arm sits at "
            f"prevalence exactly 0 or 1 ({', '.join(map(str, undefined_arms))}). "
            f"A label with no variance has no correlation to bound. This is an "
            f"abstention, not a clean reading."
        )
    caveat = (
        f" {len(undefined_arms)} arm(s) had no correlation to report "
        f"({', '.join(map(str, undefined_arms))}) and are excluded rather than "
        f"scored clean."
        if undefined_arms
        else ""
    )
    if not reachable and not tracks:
        return (
            f"NOT TESTABLE: the mature label is too rare here for |rho| to reach "
            f"the tolerance at all (ceiling sqrt(3p(1-p))). |rho|={rho:.2f} is "
            f"not evidence of cleanliness — the test could only ever have "
            f"returned clean. Read rho_vs_ceiling instead." + caveat
        )
    if tracks and not matched:
        return (
            f"CONFOUNDED: the maturity call tracks depth within an arm "
            f"(|rho|={rho:.2f}) and the arms differ {ratio:.1f}x in median depth. "
            f"The compositional term and the depth imbalance are not separable "
            f"in this data." + caveat
        )
    if tracks:
        return (
            f"The maturity call tracks depth (|rho|={rho:.2f}), but the arms are "
            f"depth-matched ({ratio:.1f}x), so it does not convert into a "
            f"between-arm difference. Report it; it is not driving the "
            f"compositional term." + caveat
        )
    if not matched:
        return (
            f"The arms differ {ratio:.1f}x in median depth, but the maturity call "
            f"does not track depth within an arm (|rho|={rho:.2f}). The imbalance "
            f"is real and this particular route from it to the compositional term "
            f"is closed." + caveat
        )
    return (
        f"Clean on both counts: arms within {ratio:.1f}x on median depth and the "
        f"maturity call does not track depth (|rho|={rho:.2f})." + caveat
    )


# ---------------------------------------------------------------------------
# Matching the arms rather than only flooring them
# ---------------------------------------------------------------------------


def match_arm_depth(
    depth: np.ndarray,
    arm: np.ndarray,
    *,
    seed: int,
    n_bins: int = 20,
) -> np.ndarray:
    """Subsample cells so both arms share a depth distribution. Returns a mask.

    WHY A FLOOR IS NOT ENOUGH
    -------------------------
    A depth floor drops the shallowest cells and thins the rest to a common
    target, which removes the *mechanism* by which depth reaches the maturity
    call. It does not make the arms comparable: on Lee/SMC, after W1's floor, the
    arms still differ 2.36x in median depth. Any residual depth sensitivity is
    then still converted into a between-arm difference — which is the
    compositional term.

    This equalises the distributions by construction. Depth is binned on pooled
    quantiles and each bin keeps ``min(n_normal, n_tumour)`` cells from each arm,
    sampled without replacement. Afterwards the arms have the same depth
    histogram, so a surviving difference in the maturity call cannot be a depth
    difference.

    IT IS EXPENSIVE AND THAT IS THE POINT
    -------------------------------------
    Matching discards most of the larger arm. If the compositional gap survives
    on a fraction of the data, that is a stronger statement than the same gap on
    all of it; if it disappears, the floor was hiding a confound rather than
    removing one. Either way the answer is not an artefact of choosing a
    threshold, because there is no threshold here to choose.

    The caller must report ``mask.sum()`` alongside anything computed on it. A
    matched subsample is a smaller cohort and its intervals are wider.
    """
    depth = np.asarray(depth, dtype=float)
    arm = np.asarray(arm)
    if depth.shape != arm.shape:
        raise ValueError(f"depth has {depth.shape}, arm has {arm.shape}")
    if depth.size == 0:
        raise ValueError("no cells")
    arms = sorted(set(arm.tolist()))
    if len(arms) != 2:
        raise ValueError(f"need exactly two arms to match, got {arms}")

    rng = np.random.default_rng(seed)
    bins = pd.qcut(pd.Series(depth), n_bins, labels=False, duplicates="drop").to_numpy()
    keep = np.zeros(len(depth), dtype=bool)
    for b in np.unique(bins[~pd.isna(bins)]):
        in_bin = bins == b
        per_arm = {a: np.flatnonzero(in_bin & (arm == a)) for a in arms}
        take = min(len(idx) for idx in per_arm.values())
        if take == 0:
            continue
        for idx in per_arm.values():
            keep[rng.choice(idx, size=take, replace=False)] = True
    return keep
