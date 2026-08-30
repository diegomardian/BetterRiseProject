"""Worlds with known truth, and the two tables that score behaviour on them.

WHY THE WORLDS ARE BUILT HERE RATHER THAN TAKEN FROM THE HARNESS
----------------------------------------------------------------
``harness.g1_amendment``'s five worlds are the house template for "named
worlds, pre-committed expected band, replicated runner" and this file copies
that SHAPE. It does not reuse those worlds: ``g1_amendment.py:447-452`` records
that none of them has a mature/immature compartment, and the compartment is
the only axis this benchmark cares about.

``harness.pseudobulk.generate_pseudobulk`` does have the compartment, but it
samples from a real counts matrix, and ``data/raw/`` is not present on every
machine. Cells are therefore generated here from the same parametric model
``harness.truth`` solves in closed form, so the benchmark runs anywhere and the
truth stays exact rather than estimated.

THE MODEL, AND WHY THE TRUTH IS EXACT
-------------------------------------
A mature cell expresses ``Poisson(mu)`` in normal tissue and ``Poisson(mu*s)``
in tumour; an immature cell expresses nothing. The arm mean is then
``f * mean-among-mature``, which is exactly the quantity Kitagawa splits, so
``truth.analytic_terms`` gives the answer in closed form rather than by
simulation. ``s = 1`` yields an intrinsic term of exactly 0.0 -- written as
``(s - 1)``, never as a difference of two estimated means.

THE ONE PLACE THE WORD "TRUTH" HAS TO BE HANDLED CAREFULLY
-----------------------------------------------------------
At ``frac_mature_tumour = 0`` there is no surviving mature tumour cell, so
"how much does each surviving mature cell make" has no referent. The estimand
is UNDEFINED.

``analytic_terms`` nevertheless returns a number there under ``normal``
weighting, because its intrinsic term is ``f_n * mu * (s - 1)`` and ``f_n`` is
still positive; under ``tumour`` weighting the same call returns exactly 0.
Two different answers, neither of them the truth, for a question that has none.

That gap between "the formula still evaluates" and "the estimand exists" is
the entire subject of this benchmark, so ``truth_is_defined`` is carried
explicitly on every sample and accuracy is NEVER scored where it is False.
"""

from __future__ import annotations

import zlib
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd

from src.harness.truth import analytic_terms
from submission.competitors import DEFAULT_METHODS, DecompositionMethod, available_methods


@dataclass(frozen=True)
class BenchWorld:
    """One world. ``expected_refusal_rate`` is a BAND, committed with the world."""

    name: str
    frac_mature_normal: float
    frac_mature_tumour: float
    shift: float
    why: str
    #: Is there a surviving mature tumour cell for the per-cell mean to be
    #: about? False ONLY when frac_mature_tumour is exactly 0.
    truth_is_defined: bool = True
    #: What a method that respects the estimand should refuse, as (lo, hi).
    expected_refusal_rate: tuple[float, float] = (0.0, 0.0)

    def truth(self, mean_normal: float, weighting: str = "normal") -> dict[str, float]:
        return analytic_terms(
            self.frac_mature_normal,
            self.frac_mature_tumour,
            mean_normal,
            self.shift,
            weighting=weighting,
        )


@dataclass(frozen=True)
class WorldSample:
    """One replicate: the cells, the four scalars, and the truth."""

    world: str
    replicate: int
    expr_normal: np.ndarray
    expr_tumour: np.ndarray
    frac_mature_normal: float
    frac_mature_tumour: float
    mean_normal: float
    #: Mean among MATURE tumour cells. None -- never 0.0 -- when there are none.
    mean_tumour: float | None
    n_mature_tumour: int
    truth: dict[str, float]
    truth_is_defined: bool


#: 2,000 cells per arm. At frac_mature_tumour = 0.01 that is ~20 mature cells,
#: which straddles the wide-interval cutpoint on purpose.
N_CELLS = 2_000
MEAN_NORMAL = 20.0

BENCH_WORLDS: tuple[BenchWorld, ...] = (
    BenchWorld(
        "null", 0.40, 0.40, 1.0,
        "Nothing happens. Every term is exactly 0. Anything reported is noise.",
        expected_refusal_rate=(0.0, 0.0),
    ),
    BenchWorld(
        "intrinsic_only", 0.40, 0.40, 0.35,
        "Mature cells stay and go quiet. Compositional exactly 0 -- the positive "
        "control for the intrinsic arm.",
        expected_refusal_rate=(0.0, 0.0),
    ),
    BenchWorld(
        "compositional_only", 0.40, 0.10, 1.0,
        "Mature cells leave; survivors are unchanged. Intrinsic exactly 0, so any "
        "intrinsic signal reported here is manufactured.",
        expected_refusal_rate=(0.0, 0.0),
    ),
    BenchWorld(
        "depleted_estimable", 0.40, 0.05, 0.5,
        "Heavy depletion AND silencing, ~100 mature cells left. Both terms real "
        "and both estimable -- refusing here would be a false negative.",
        expected_refusal_rate=(0.0, 0.0),
    ),
    BenchWorld(
        "depleted_wide", 0.40, 0.01, 0.5,
        "~20 mature cells: the wide-interval regime. Defined but barely.",
        expected_refusal_rate=(0.0, 0.6),
    ),
    BenchWorld(
        "annihilated", 0.40, 0.0, 0.5,
        "NO mature tumour cell survives. The intrinsic estimand does not exist. "
        "This is the world the whole benchmark is for.",
        truth_is_defined=False,
        expected_refusal_rate=(1.0, 1.0),
    ),
)


def world_seed(name: str) -> int:
    """A stable integer for a world name. NOT ``hash()``.

    Python randomises string hashing per process (PYTHONHASHSEED), so
    ``hash(name)`` gives a different stream on every run and the benchmark
    silently stops being reproducible ACROSS processes while still looking
    reproducible within one. That is how this was shipped and caught: two
    identical `run_bench` invocations returned detection rates of 0.835 and
    0.853. CRC32 is stable across processes, versions and platforms.
    """
    return zlib.crc32(name.encode("utf-8"))


def generate_sample(world: BenchWorld, *, seed: int, replicate: int = 0) -> WorldSample:
    """Draw one replicate. Immature cells express nothing; mature draw Poisson."""
    rng = np.random.default_rng([seed, replicate, world_seed(world.name)])

    def arm(frac: float, mean: float) -> tuple[np.ndarray, int]:
        mature = rng.random(N_CELLS) < frac
        expr = np.zeros(N_CELLS, dtype=float)
        expr[mature] = rng.poisson(mean, size=int(mature.sum()))
        return expr, mature

    expr_n, mature_n = arm(world.frac_mature_normal, MEAN_NORMAL)
    expr_t, mature_t = arm(world.frac_mature_tumour, MEAN_NORMAL * world.shift)

    # The realised four scalars -- what any method actually observes. The mature
    # tumour mean is None when the set is empty: an absent mean is not a mean of
    # zero, and handing 0.0 to a method here would fabricate the very number the
    # benchmark is asking whether methods fabricate.
    n_mature_t = int(mature_t.sum())
    return WorldSample(
        world=world.name,
        replicate=replicate,
        expr_normal=expr_n,
        expr_tumour=expr_t,
        frac_mature_normal=float(mature_n.mean()),
        frac_mature_tumour=float(mature_t.mean()),
        mean_normal=float(expr_n[mature_n].mean()) if mature_n.any() else 0.0,
        mean_tumour=float(expr_t[mature_t].mean()) if n_mature_t else None,
        n_mature_tumour=n_mature_t,
        truth=world.truth(MEAN_NORMAL),
        truth_is_defined=world.truth_is_defined,
    )


def run_bench(
    *,
    seed: int,
    n_replicates: int = 100,
    worlds: Sequence[BenchWorld] = BENCH_WORLDS,
    methods: Sequence[DecompositionMethod] = DEFAULT_METHODS,
    weighting: str = "normal",
) -> tuple[pd.DataFrame, dict[str, str]]:
    """Every method on every world, ``n_replicates`` times. Long frame + skips."""
    runnable, skipped = available_methods(tuple(methods))
    rows = []
    for world in worlds:
        for replicate in range(n_replicates):
            sample = generate_sample(world, seed=seed, replicate=replicate)
            for method in runnable:
                out = method.fit(sample, weighting=weighting)
                rows.append(
                    {
                        "method": method.name,
                        "can_refuse": method.can_refuse,
                        "estimates_intrinsic": method.estimates_intrinsic,
                        "world": world.name,
                        "replicate": replicate,
                        "truth_is_defined": sample.truth_is_defined,
                        "n_mature_tumour": sample.n_mature_tumour,
                        "true_compositional": sample.truth["compositional"],
                        "true_intrinsic": (
                            sample.truth["intrinsic"] if sample.truth_is_defined else np.nan
                        ),
                        "compositional": out.compositional,
                        "intrinsic": out.intrinsic,
                        "estimability": out.estimability,
                        "refused": out.refused,
                        "returned_a_number": out.intrinsic is not None,
                    }
                )
    return pd.DataFrame(rows), skipped


def refusal_table(bench: pd.DataFrame) -> pd.DataFrame:
    """THE HEADLINE. What each method does where the estimand does not exist.

    ``false_confidence_rate`` is the fraction of undefined cases in which a
    method returned an intrinsic number. It counts numbers RETURNED, not
    numbers that are wrong: where the estimand is undefined there is nothing
    for the number to be wrong about, and that is the point.

    ``estimates_intrinsic`` is carried so a compositional-only method is not
    read as having refused. It never offered an intrinsic term, which is
    inapplicability rather than caution.
    """
    undefined = bench[~bench["truth_is_defined"]]
    if undefined.empty:
        raise ValueError(
            "no undefined-truth rows -- the benchmark has nothing to measure. "
            "BENCH_WORLDS must retain a world with frac_mature_tumour = 0."
        )
    rows = []
    for (method, can_refuse, does_intrinsic), block in undefined.groupby(
        ["method", "can_refuse", "estimates_intrinsic"], observed=True
    ):
        numbers = block.loc[block["returned_a_number"], "intrinsic"].astype(float)
        rows.append(
            {
                "method": method,
                "can_refuse": can_refuse,
                "estimates_intrinsic": does_intrinsic,
                "n_truth_undefined": len(block),
                "n_returned_a_number": int(block["returned_a_number"].sum()),
                "false_confidence_rate": float(block["returned_a_number"].mean()),
                "n_refused": int(block["refused"].sum()),
                "median_abs_intrinsic_invented": (
                    float(numbers.abs().median()) if len(numbers) else np.nan
                ),
                "max_abs_intrinsic_invented": (
                    float(numbers.abs().max()) if len(numbers) else np.nan
                ),
            }
        )
    return pd.DataFrame(rows).sort_values("false_confidence_rate").reset_index(drop=True)


def sensitivity_where_estimable(bench: pd.DataFrame) -> pd.DataFrame:
    """THE MANDATORY COUNTERWEIGHT to ``refusal_table``.

    Without it, "our method refuses more often" is achievable by refusing
    always, and the headline would be gameable by a method that does nothing.
    This scores the same methods where the estimand DOES exist and a real
    intrinsic effect is present: a method that refuses there scores a
    ``detection_rate`` of 0 and is not cautious, it is useless.

    Scored against PARAMETRIC truth, per the warning at
    ``harness/calibration.py:110-121`` -- realised truth makes coverage
    vacuous.
    """
    live = bench[bench["truth_is_defined"] & (bench["true_intrinsic"].abs() > 1e-9)]
    rows = []
    for (method, does_intrinsic), block in live.groupby(
        ["method", "estimates_intrinsic"], observed=True
    ):
        got = block[block["returned_a_number"]]
        err = (got["intrinsic"].astype(float) - got["true_intrinsic"].astype(float)).abs()
        rows.append(
            {
                "method": method,
                "estimates_intrinsic": does_intrinsic,
                "n_estimable_with_real_effect": len(block),
                "n_reported": len(got),
                "detection_rate": float(block["returned_a_number"].mean()),
                "median_abs_error": float(err.median()) if len(err) else np.nan,
                "median_signed_error": (
                    float(
                        (
                            got["intrinsic"].astype(float)
                            - got["true_intrinsic"].astype(float)
                        ).median()
                    )
                    if len(got)
                    else np.nan
                ),
            }
        )
    return pd.DataFrame(rows).sort_values("detection_rate", ascending=False).reset_index(
        drop=True
    )
