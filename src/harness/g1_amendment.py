"""Can the proposed G1 fail? W2's adjudication of prereg Amendment 2. Issue #37.

W1 cannot ratify its own pre-registration — that is the one thing an amendment to
a pre-registered statistic cannot survive — so W1 asked the team, and the part of
the team with ground truth is the harness. This module is W2's answer.

WHAT IS BEING ADJUDICATED
-------------------------
Decision #17 specified G1 as Spearman between gene abundance and apparent loss,
where loss is a raw Δ of per-cell means. Amendment 2 says that statistic is
arithmetically broken two ways — an absolute difference scales with what it
differences, so the correlation measures the units; and tier D, which carries
G1's falsification logic, holds one gene, over which a correlation is undefined.
It proposes the MA construction (M = log₂ tumour/normal, A = the mean of the
logs) genome-wide, with panel genes located by percentile of M within their own
abundance bin, and three thresholds committed before any real M or A exists.

WHAT THIS MODULE ADDS THAT THE AMENDMENT DOES NOT
-------------------------------------------------
The amendment demonstrates the old statistic fails on a null. That establishes
the old one is broken. It does **not** establish the new one works, and those are
different claims — this repository has now shipped four checks that could not
fail, and replacing one with another would be the fifth.

So the question here is the opposite one: **can the proposed G1 both pass and
fail, for the right reasons?** Four worlds with known truth, three of which must
fail and one of which must pass, plus the number nobody had computed — how often
the gate passes by chance when nothing is going on.

Simulation only. Nothing here is a G1 result, and running it does not clear G1;
`src/reference/checks.py` returns ``not_estimable`` until the team ratifies, and
that is the correct state.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from typing import Final

import numpy as np
import pandas as pd

#: Amendment 2 §4, committed 2026-08-25 before any M or A existed on real data.
#: These are W1's numbers. W2 does not get to adjust them — testing whether a
#: pre-committed threshold discriminates is legitimate; tuning it here is not.
TIER_A_MAX_MEDIAN: Final = 0.20
TIER_D_MIN_PERCENTILE: Final = 0.50
MIN_SEPARATION: Final = 0.30

#: Amendment 2 §4: 20 equal-size bins by A, ≈1,960 genes each on the 39,236-gene
#: shared index. Binned rather than residuals-from-a-fit because Spearman
#: produces no fitted values, so "residual from the trend" has no definition
#: until someone picks a functional form — a free parameter nobody pre-registered.
N_ABUNDANCE_BINS: Final = 20

#: Where the panel genes sit on the abundance range, as quantiles. Fixed and
#: spread out so that no world gets a lucky draw and the four worlds are
#: compared on the same footing. Real panel genes have real abundances; this is
#: a simulation of the *rule*, not of the panel.
TIER_A_QUANTILES: Final = (0.10, 0.35, 0.60, 0.85)
TIER_D_QUANTILE: Final = 0.475

FoldChangeModel = Callable[[np.ndarray, np.ndarray, int, np.random.Generator], np.ndarray]


@dataclass(frozen=True)
class G1World:
    """A simulated world with a known answer, and the pass rate G1 owes it.

    A **band**, not a single verdict, because two of these worlds are genuinely
    stochastic in the gate's own terms and asserting "it passed on my seed" is
    the n=1 conclusion this project keeps retracting. The band is what the
    ratification test checks.
    """

    name: str
    fold_change: FoldChangeModel
    expected_pass_rate: tuple[float, float]
    why: str


# ---------------------------------------------------------------------------
# The construction
# ---------------------------------------------------------------------------


def simulate_arms(
    abundance: np.ndarray,
    fold: np.ndarray,
    *,
    n_cells: int,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    """Per-gene mean counts in each arm, Poisson around ``abundance``.

    Counting noise is the point: the MA transform takes logs, and a construction
    without sampling noise would make every ratio exact and hide whether the
    percentile rule survives the noise a real count matrix has.

    The mean of ``n_cells`` iid Poisson(λ) is Poisson(n_cells·λ)/n_cells, so the
    per-cell matrix is never materialised. That is an identity, not an
    approximation, and it is the difference between 60 million variates per arm
    and 20 thousand — which is what makes replicating the worlds affordable
    rather than a single lucky draw.
    """
    if n_cells < 1:
        raise ValueError(f"n_cells={n_cells} must be at least 1")
    normal = rng.poisson(abundance * n_cells) / n_cells
    tumour = rng.poisson(abundance * fold * n_cells) / n_cells
    return normal, tumour


def ma_transform(
    normal: np.ndarray, tumour: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(A, M, keep)`` — Amendment 2 §1's construction.

    ``keep`` is the genes with a non-zero mean in **both** arms. The amendment is
    explicit that a log ratio against zero is not a number and that a pseudocount
    there invents the very quantity being measured, so those genes are dropped
    rather than floored. The mask is returned so the caller can map percentiles
    back onto the full gene index instead of silently re-indexing.
    """
    keep = (normal > 0) & (tumour > 0)
    m = np.log2(tumour[keep] / normal[keep])
    a = 0.5 * (np.log2(tumour[keep]) + np.log2(normal[keep]))
    return a, m, keep


def within_bin_percentile(
    a: np.ndarray, m: np.ndarray, *, n_bins: int = N_ABUNDANCE_BINS
) -> np.ndarray:
    """Percentile of each gene's M among genes of comparable abundance.

    **Low means more lost.** Equal-size bins by A, then the average-rank
    percentile of M inside each bin, so the result is uniform on [0, 1] within
    every bin whatever the shape of the abundance–loss trend. That is what makes
    the rule assumption-free: it does not need the trend to be linear, monotone
    or anything else, only that "genes of comparable abundance" is a meaningful
    comparison set.
    """
    if a.shape != m.shape:
        raise ValueError(f"A has {a.shape} and M has {m.shape}")
    if len(a) < n_bins:
        raise ValueError(
            f"{len(a)} genes cannot be split into {n_bins} bins. The amendment "
            f"assumes a genome-wide index (≈39,236 genes); this rule is not "
            f"meaningful on a panel-sized set."
        )
    order = np.argsort(a, kind="stable")
    bin_of = np.empty(len(a), dtype=int)
    bin_of[order] = np.arange(len(a)) * n_bins // len(a)

    percentile = np.empty(len(a), dtype=float)
    for b in range(n_bins):
        selected = bin_of == b
        values = m[selected]
        ranks = values.argsort().argsort().astype(float)
        percentile[selected] = (ranks + 0.5) / len(values)
    return percentile


# ---------------------------------------------------------------------------
# The verdict
# ---------------------------------------------------------------------------


def amendment2_verdict(
    tier_a_percentiles: np.ndarray, tier_d_percentile: float
) -> dict[str, object]:
    """Amendment 2 §4's three pre-committed thresholds. PASS only if all clear.

    G1 **fails** if tier A's median percentile exceeds 0.20 (the compositional
    targets are not lost more than four in five genes of comparable abundance),
    or MS4A12 sits below 0.50 (the gene chosen to be retained is more lost than a
    typical gene of its abundance), or the two are separated by less than 0.30
    (`config/panel.yaml`'s falsification rule in its own terms).

    NaN IS REFUSED, NOT COMPARED
    ----------------------------
    ``scipy.stats.spearmanr`` on one observation returns ``nan``, and it does not
    raise. Under decision #17 as written, tier D's correlation is that ``nan``,
    and ``abs(nan) > 0.5`` is ``False`` — so the rule "fail if |ρ| > 0.5" would
    have **silently passed** the one tier that carries the falsification logic.
    That is this repository's recurring defect, in the gate itself: not a check
    that reports a wrong number, a check that reports success because it could
    not be computed. This function raises instead.
    """
    tier_a = np.asarray(tier_a_percentiles, dtype=float)
    if tier_a.size == 0:
        raise ValueError("tier A has no percentiles; the verdict is not computable")
    if not np.all(np.isfinite(tier_a)) or not np.isfinite(tier_d_percentile):
        raise ValueError(
            "a percentile is nan or infinite, so the threshold comparison would "
            "be False and G1 would read as PASS without having been evaluated. "
            "scipy returns nan for a one-observation correlation rather than "
            "raising; that is how decision #17's tier D would have passed "
            "silently. Fix the input — do not compare it."
        )

    median = float(np.median(tier_a))
    separation = float(tier_d_percentile) - median
    failures = []
    if median > TIER_A_MAX_MEDIAN:
        failures.append(f"tier A median {median:.3f} > {TIER_A_MAX_MEDIAN}")
    if tier_d_percentile < TIER_D_MIN_PERCENTILE:
        failures.append(f"MS4A12 {tier_d_percentile:.3f} < {TIER_D_MIN_PERCENTILE}")
    if separation < MIN_SEPARATION:
        failures.append(f"separation {separation:.3f} < {MIN_SEPARATION}")
    return {
        "verdict": "FAIL" if failures else "PASS",
        "tier_a_median_pct": median,
        "tier_d_pct": float(tier_d_percentile),
        "separation": separation,
        "failures": failures,
    }


# ---------------------------------------------------------------------------
# The four worlds
# ---------------------------------------------------------------------------


def _uniform_loss(
    abundance: np.ndarray, panel: np.ndarray, n_genes: int, rng: np.random.Generator
) -> np.ndarray:
    return np.full(n_genes, 0.7)


def _isolated_tier_a_loss(
    abundance: np.ndarray, panel: np.ndarray, n_genes: int, rng: np.random.Generator
) -> np.ndarray:
    """Tier A gone, everything else — MS4A12 included — untouched."""
    fold = np.ones(n_genes)
    fold[panel[:-1]] = 0.15
    fold[panel[-1]] = 1.0
    return fold


def _broad_loss_tier_d_retained(
    abundance: np.ndarray, panel: np.ndarray, n_genes: int, rng: np.random.Generator
) -> np.ndarray:
    """Broad differentiation loss, tier A hardest, MS4A12 genuinely maintained.

    This is the project's own premise — marker loss is broad — and it is the
    world in which MS4A12 has something to be retained *against*.
    """
    fold = np.exp(rng.normal(-0.55, 0.45, n_genes))
    fold[panel[:-1]] = 0.15
    fold[panel[-1]] = 1.0
    return fold


def _pure_soup(
    abundance: np.ndarray, panel: np.ndarray, n_genes: int, rng: np.random.Generator
) -> np.ndarray:
    # Loss is a function of abundance and nothing else — no gene-specific
    # biology anywhere. This is what an ambient artefact looks like, and it is
    # the case G1 exists to catch.
    logs = np.log10(abundance)
    return 0.5 + 0.45 * (logs - logs.min()) / np.ptp(logs)


def _tiers_drift_together(
    abundance: np.ndarray, panel: np.ndarray, n_genes: int, rng: np.random.Generator
) -> np.ndarray:
    fold = np.ones(n_genes)
    fold[panel] = 0.15  # MS4A12 lost as hard as tier A
    return fold


G1_WORLDS: Final = (
    G1World(
        name="uniform_loss",
        fold_change=_uniform_loss,
        expected_pass_rate=(0.00, 0.15),
        why="every gene loses 30%; tier A is not lost more than its abundance "
        "peers, so the project's premise is unsupported and PASS would be wrong",
    ),
    G1World(
        name="broad_loss_tier_d_retained",
        fold_change=_broad_loss_tier_d_retained,
        expected_pass_rate=(0.95, 1.00),
        why="the world the project claims — broad loss, tier A hardest, MS4A12 "
        "maintained against it. A gate that cannot pass here is as useless as "
        "one that cannot fail",
    ),
    G1World(
        name="isolated_tier_a_loss",
        fold_change=_isolated_tier_a_loss,
        expected_pass_rate=(0.35, 0.65),
        why="tier A gone and NOTHING else moves. Also a true signal, and the "
        "gate is a coin flip on it — see AMENDMENT_2_POWER_CAVEAT",
    ),
    G1World(
        name="pure_soup",
        fold_change=_pure_soup,
        expected_pass_rate=(0.00, 0.10),
        why="loss is a function of abundance and nothing else — the ambient "
        "artefact G1 exists to catch",
    ),
    G1World(
        name="tiers_drift_together",
        fold_change=_tiers_drift_together,
        expected_pass_rate=(0.00, 0.02),
        why="config/panel.yaml's falsification rule: if the tiers return the "
        "same answer, no biological claim may be made",
    ),
)

#: The one finding W2 would not have got from reading the amendment.
AMENDMENT_2_POWER_CAVEAT: Final = """
Threshold 2 asks MS4A12's within-bin percentile to be >= 0.50. A gene that is
unchanged AGAINST AN UNCHANGED BACKGROUND sits at 0.50 by definition, so the
test is a fair coin: measured P(PASS) = 0.497 over 300 replicates of
`isolated_tier_a_loss`, a world in which the project's claim is TRUE. Against a
broadly-lost background MS4A12 sits at 0.890 +/- 0.015 and P(PASS) = 1.000.

So G1's power is not a property of the gate alone — it depends on there being a
broad loss background for the retained control to stand out against. That
premise is the project's own, and it is probably satisfied, but it was unstated,
and if it fails the gate rejects a true signal half the time for reasons that
have nothing to do with ambient RNA.

Amendment 2 made tier D computable. It did not make n=1 powerful, and those are
different achievements.

Not a reason to reject the amendment: the old statistic returned rho ~ -1 on a
null and could not be computed on tier D at all. This is strictly better. It is
a reason to state the premise and to read a G1 FAIL with the caveat attached.
"""


# ---------------------------------------------------------------------------
# Running them
# ---------------------------------------------------------------------------


def run_world(
    world: G1World,
    *,
    seed: int,
    replicate: int = 0,
    n_genes: int = 20_000,
    n_cells: int = 3_000,
    n_bins: int = N_ABUNDANCE_BINS,
) -> dict[str, object]:
    """One draw of one world. Returns a row, including both statistics.

    ``rho_decision_17`` is the superseded statistic — cohort abundance against a
    raw Δ — carried alongside so the two can be read against the same truth
    rather than against each other's write-ups.
    """
    from scipy.stats import spearmanr

    rng = np.random.default_rng(seed + replicate)
    abundance = 10 ** rng.uniform(-2, 2, n_genes)

    # Panel genes at fixed abundance quantiles: tier A first, tier D last.
    ranked = np.argsort(abundance, kind="stable")
    panel = np.array([ranked[int(q * (n_genes - 1))] for q in (*TIER_A_QUANTILES, TIER_D_QUANTILE)])

    fold = world.fold_change(abundance, panel, n_genes, rng)
    normal, tumour = simulate_arms(abundance, fold, n_cells=n_cells, rng=rng)
    a, m, keep = ma_transform(normal, tumour)

    percentile = np.full(n_genes, np.nan)
    percentile[np.flatnonzero(keep)] = within_bin_percentile(a, m, n_bins=n_bins)

    result = amendment2_verdict(percentile[panel[:-1]], float(percentile[panel[-1]]))
    return {
        "world": world.name,
        "replicate": replicate,
        "verdict": result["verdict"],
        "expected_pass_low": world.expected_pass_rate[0],
        "expected_pass_high": world.expected_pass_rate[1],
        "rho_decision_17": float(spearmanr(0.5 * (tumour + normal), tumour - normal).statistic),
        "rho_g1a": float(spearmanr(a, m).statistic),
        "tier_a_median_pct": result["tier_a_median_pct"],
        "tier_d_pct": result["tier_d_pct"],
        "separation": result["separation"],
        "n_genes_kept": int(keep.sum()),
        "n_genes": n_genes,
        "seed": seed + replicate,
    }


def results_table(*, seed: int, n_replicates: int = 60, **kwargs: object) -> pd.DataFrame:
    """:func:`run_all_worlds` in ``G1_AMENDMENT_COLUMNS`` order, ready to write.

    The gate memo quotes these numbers, so they belong in versioned parquet with
    a sha and a seed like every other result (CLAUDE.md invariant 10) rather than
    living only in a docstring.
    """
    from src.harness.results import G1_AMENDMENT_COLUMNS

    return run_all_worlds(seed=seed, n_replicates=n_replicates, **kwargs).loc[
        :, list(G1_AMENDMENT_COLUMNS)
    ]


def run_all_worlds(*, seed: int, n_replicates: int = 20, **kwargs: object) -> pd.DataFrame:
    """Every world, every replicate. One row each.

    Replicated rather than run once because a single draw of the null world
    lands wherever the panel quantiles happen to fall, and "it failed on my
    seed" is the kind of n=1 conclusion this project keeps having to retract.
    """
    rows = [
        run_world(world, seed=seed, replicate=r, **kwargs)  # type: ignore[arg-type]
        for world in G1_WORLDS
        for r in range(n_replicates)
    ]
    return pd.DataFrame(rows)


def null_pass_rate(*, seed: int, n_trials: int = 100_000) -> dict[str, float]:
    """How often the three thresholds pass when nothing is going on.

    Under a true null the panel genes are ordinary genes, so their within-bin
    percentiles are uniform on [0, 1] and independent. This is then arithmetic
    about the thresholds rather than a claim about colorectal cancer, which is
    why it is computed on uniforms directly and not simulated through counts.

    The per-threshold rates matter more than the joint one. A threshold that
    passes half the time under the null is a coin flip, whatever it is called.
    """
    rng = np.random.default_rng(seed)
    tier_a_median = np.median(rng.uniform(size=(n_trials, len(TIER_A_QUANTILES))), axis=1)
    tier_d = rng.uniform(size=n_trials)
    one = tier_a_median <= TIER_A_MAX_MEDIAN
    two = tier_d >= TIER_D_MIN_PERCENTILE
    three = (tier_d - tier_a_median) >= MIN_SEPARATION
    return {
        "p_threshold_1_tier_a": float(one.mean()),
        "p_threshold_2_tier_d": float(two.mean()),
        "p_threshold_3_separation": float(three.mean()),
        "p_g1_passes_on_a_null": float((one & two & three).mean()),
        "n_trials": float(n_trials),
    }


# ---------------------------------------------------------------------------
# The world W2's ratification could not represent — issue #46
# ---------------------------------------------------------------------------
#
# The five worlds above vary per-gene fold change against a FLAT background.
# None of them has a mature/immature compartment, so in none of them does a
# gene's mean over all epithelium differ from its mean within mature cells.
# That is exactly the axis along which G1 turned out to be broken, which is why
# the ratification passed it: a harness that cannot represent a confound cannot
# rule it out.
#
# This is not a proposed repair. Threshold 2 is frozen by having been seen to
# fail, and W2 — having ratified the amendment — is the last party who should
# propose a replacement. It encodes the NEGATIVE result so it cannot be
# rediscovered as news.

#: Mature-restricted genes are expressed in mature cells and ~nowhere else.
#: Tier A and MS4A12 are both in this class — panel.yaml calls MS4A12
#: "colonocyte-restricted yet frequently maintained" in one line, and that
#: conjunction is the whole problem.
RESTRICTED_SHARE: Final = 0.05


def simulate_compartment_world(
    *,
    seed: int,
    n_genes: int = 20_000,
    n_cells: int = 3_000,
    frac_mature_normal: float = 0.50,
    frac_mature_tumour: float = 0.15,
    tier_a_silencing: float = 0.15,
) -> dict[str, object]:
    """A cohort with a mature compartment, so depletion and silencing differ.

    ``tier_a_silencing`` is the multiplicative shift applied to tier A **inside
    surviving mature cells**. Set it to 1.0 for a composition-only world and
    below 1.0 for the world where the project's claim is true.

    Returns the within-bin percentiles the amendment's thresholds read, plus the
    verdict, so a caller can ask what G1 says about a world with a known answer.
    """
    rng = np.random.default_rng(seed)
    abundance = 10 ** rng.uniform(-2, 2, n_genes)

    restricted = rng.random(n_genes) < RESTRICTED_SHARE
    ranked = np.argsort(abundance, kind="stable")
    panel = np.array(
        [ranked[int(q * (n_genes - 1))] for q in (*TIER_A_QUANTILES, TIER_D_QUANTILE)]
    )
    restricted[panel] = True  # tier A and MS4A12 are all mature-restricted

    # Per-cell mean in each compartment. Restricted genes are absent from
    # immature cells; unrestricted genes are the same in both.
    mature_normal = abundance.copy()
    immature = np.where(restricted, 0.0, abundance)

    mature_tumour = mature_normal.copy()
    mature_tumour[panel[:-1]] *= tier_a_silencing  # tier A silenced, MS4A12 not

    normal_mean = frac_mature_normal * mature_normal + (1 - frac_mature_normal) * immature
    tumour_mean = frac_mature_tumour * mature_tumour + (1 - frac_mature_tumour) * immature

    normal, tumour = (
        rng.poisson(normal_mean * n_cells) / n_cells,
        rng.poisson(tumour_mean * n_cells) / n_cells,
    )
    a, m, keep = ma_transform(normal, tumour)
    percentile = np.full(n_genes, np.nan)
    percentile[np.flatnonzero(keep)] = within_bin_percentile(a, m)

    tier_a = percentile[panel[:-1]]
    tier_d = float(percentile[panel[-1]])
    finite = np.isfinite(tier_a) & True
    result = (
        amendment2_verdict(tier_a[finite], tier_d)
        if finite.any() and np.isfinite(tier_d)
        else {"verdict": "not_estimable", "failures": ["a panel gene dropped out"]}
    )
    return {
        "hypothesis_is_true": tier_a_silencing < 1.0,
        "tier_a_median_pct": float(np.nanmedian(tier_a)),
        "tier_d_pct": tier_d,
        "verdict": result["verdict"],
        "failures": result.get("failures", []),
    }
