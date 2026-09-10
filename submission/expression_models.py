"""Count distributions for the mature-cell arm, and why the mean is held fixed.

WHY THIS MODULE EXISTS
----------------------
``FINDINGS.md`` limitation 1 says, of the whole benchmark:

    "Fully synthetic, and the generative model is *exactly* the model Kitagawa
    assumes. This is a favourable setting for the estimator. It tests the
    **gate**, not robustness to a misspecified model."

That is an admission, not a measurement. This module turns it into one: the
same worlds, the same committed method classes, a different count
distribution. Nothing in ``competitors.py`` is changed or subclassed for it.

THE ONE DESIGN DECISION, AND IT IS LOAD-BEARING
-----------------------------------------------
Every model here has marginal mean **exactly ``mu``** among mature cells. Only
the variance and the shape change.

This is not tidiness. ``bench.BenchWorld.truth`` calls
``harness.truth.analytic_terms(f_n, f_t, mu, s)``, which is a statement about
*means*. If a misspecified arm also moved the mean -- as a zero-inflated draw
of ``NB(mu)`` with inflation ``pi`` does, to ``(1 - pi) * mu`` -- then the
truth column would be wrong and every accuracy number computed against it
would be measuring the discrepancy in the generator rather than in the method.
The finding would read as "the estimator is biased under zero inflation" when
what actually happened is that the benchmark lied about the answer.

So the zero-inflated arm draws its non-structural-zero counts at
``mu / (1 - pi)``. The marginal mean is ``mu``; the variance is inflated; the
truth stays exact. **Misspecification is confined to the second moment and the
shape, which is the axis the limitation is about.**

THE PARAMETERISATION
--------------------
``dispersion`` is the negative binomial's size parameter ``r``, so

    Var = mu + mu^2 / r

and ``r -> inf`` recovers Poisson. Smaller ``r`` is *more* overdispersed:
``r = 10`` is mild, ``r = 0.5`` is severe (variance about 3 mu^2 at mu = 20, a
coefficient of variation above 1.4 -- heavier than anything a UMI count matrix
plausibly shows, which is the point of including it).

BIT-IDENTITY WITH THE COMMITTED BENCHMARK
------------------------------------------
``POISSON.draw`` is literally ``rng.poisson(mean, size=size)`` -- the exact
call ``bench.generate_sample`` made before this module existed, in the same
position in the stream. The committed six-world benchmark therefore returns
the same numbers to the last bit after this refactor, and
``tests/test_bench_extensions.py`` asserts that against committed values
rather than trusting it.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


class ExpressionModelError(ValueError):
    """A count model whose parameters do not describe a distribution."""


@dataclass(frozen=True)
class ExpressionModel:
    """A count distribution for a mature cell, with marginal mean ``mu``.

    ``kind`` is one of ``poisson``, ``nb``, ``zip``, ``zinb``. ``dispersion``
    is the NB size ``r`` and must be absent for the Poisson kinds;
    ``zero_inflation`` is the structural-zero probability ``pi``.
    """

    name: str
    kind: str
    dispersion: float | None = None
    zero_inflation: float = 0.0

    def __post_init__(self) -> None:
        if self.kind not in ("poisson", "nb", "zip", "zinb"):
            raise ExpressionModelError(f"unknown kind {self.kind!r}")
        needs_disp = self.kind in ("nb", "zinb")
        if needs_disp and not (self.dispersion is not None and self.dispersion > 0):
            raise ExpressionModelError(
                f"{self.name}: kind {self.kind!r} needs a positive dispersion "
                f"(the NB size r), got {self.dispersion!r}"
            )
        if not needs_disp and self.dispersion is not None:
            raise ExpressionModelError(
                f"{self.name}: kind {self.kind!r} has no dispersion parameter, but "
                f"{self.dispersion!r} was given. Silently ignoring it would let a "
                f"sweep report an overdispersion arm that never overdispersed."
            )
        if not 0.0 <= self.zero_inflation < 1.0:
            raise ExpressionModelError(
                f"{self.name}: zero_inflation must be in [0, 1), got "
                f"{self.zero_inflation!r}. At pi = 1 no cell expresses anything and "
                f"the conditional mean mu/(1-pi) is undefined."
            )
        if (self.zero_inflation > 0) != (self.kind in ("zip", "zinb")):
            raise ExpressionModelError(
                f"{self.name}: kind {self.kind!r} and zero_inflation "
                f"{self.zero_inflation!r} disagree about whether this model is "
                f"zero-inflated."
            )

    # -- the moments this model claims, in closed form -----------------------

    def theoretical_variance(self, mean: float) -> float:
        """Var of one mature cell's count at marginal mean ``mean``.

        Written out rather than measured, so that the test asserting the
        sampler matches it is a test and not a tautology. For the
        zero-inflated mixture with component mean ``m = mu / (1 - pi)`` and
        component variance ``v``:

            Var = (1 - pi) * (v + m^2) - mu^2
        """
        pi = self.zero_inflation
        if pi == 0.0:
            return mean if self.kind == "poisson" else mean + mean**2 / self.dispersion
        component_mean = mean / (1.0 - pi)
        component_var = (
            component_mean
            if self.kind == "zip"
            else component_mean + component_mean**2 / self.dispersion
        )
        return (1.0 - pi) * (component_var + component_mean**2) - mean**2

    def draw(self, rng: np.random.Generator, mean: float, size: int) -> np.ndarray:
        """``size`` counts with marginal mean ``mean``.

        ``size = 0`` returns an empty array and -- for every kind here --
        consumes nothing from ``rng``. That matters for the ``annihilated``
        world, where the mature tumour set is empty by construction.
        """
        if self.kind == "poisson":
            # Byte-for-byte the call bench.generate_sample used before this
            # module existed. Do not fold this into the branch below.
            return rng.poisson(mean, size=size).astype(float)
        if self.kind == "nb":
            return self._nb(rng, mean, size).astype(float)

        pi = self.zero_inflation
        component_mean = mean / (1.0 - pi)
        if self.kind == "zip":
            counts = rng.poisson(component_mean, size=size).astype(float)
        else:
            counts = self._nb(rng, component_mean, size).astype(float)
        structural_zero = rng.random(size) < pi
        counts[structural_zero] = 0.0
        return counts

    def _nb(self, rng: np.random.Generator, mean: float, size: int) -> np.ndarray:
        """NB with mean ``mean`` and variance ``mean + mean^2 / r``.

        numpy parameterises by (r, p) with mean ``r(1-p)/p``; solving for p at
        a fixed mean gives ``p = r / (r + mean)``.
        """
        r = float(self.dispersion)
        p = r / (r + mean)
        return rng.negative_binomial(r, p, size=size)


#: The committed generator. Everything else is measured against it.
POISSON = ExpressionModel("poisson", "poisson")

#: Overdispersion arms. r = 10 mild, r = 2 strong, r = 0.5 extreme.
NB_10 = ExpressionModel("nb_disp10", "nb", dispersion=10.0)
NB_2 = ExpressionModel("nb_disp2", "nb", dispersion=2.0)
NB_05 = ExpressionModel("nb_disp0.5", "nb", dispersion=0.5)

#: Zero inflation, the other way a real count matrix departs from Poisson.
#: ``zip_pi0.3`` isolates inflation from overdispersion; ``zinb`` has both,
#: which is what a droplet protocol with ambient dropout actually looks like.
ZIP_30 = ExpressionModel("zip_pi0.3", "zip", zero_inflation=0.3)
ZINB_30_2 = ExpressionModel(
    "zinb_pi0.3_disp2", "zinb", dispersion=2.0, zero_inflation=0.3
)

#: The sweep, committed in order. Poisson first so a reader sees the baseline
#: it is being compared against on the first row of every table.
EXPRESSION_MODELS: tuple[ExpressionModel, ...] = (
    POISSON,
    NB_10,
    NB_2,
    NB_05,
    ZIP_30,
    ZINB_30_2,
)

MODELS_BY_NAME: dict[str, ExpressionModel] = {m.name: m for m in EXPRESSION_MODELS}
