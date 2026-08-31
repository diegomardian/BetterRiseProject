"""Decomposition methods, as interchangeable adapters.

WHAT THIS BENCHMARK IS FOR
--------------------------
README.md:83-86 claims, of every existing method:

    "Every existing method returns a number. None flags that the intrinsic
    estimate is meaningless in a tumour with no mature cells left."

That is an **empirical claim about other methods** and it has never been
measured. This module makes it measurable by putting every method behind one
interface, so the same worlds can be put to all of them and the answer is a
table rather than an assertion.

THE HEADLINE IS NOT ACCURACY
----------------------------
Accuracy where the estimand exists is the ordinary question and every method
here can be scored on it. The question this benchmark exists for is what a
method does where the estimand does **not** exist -- a tumour with no mature
cells, where "how much does each surviving mature cell make" has no referent.

Three behaviours have to stay distinguishable, and collapsing any two of them
would make the headline meaningless:

``refused=True``
    The method determined the estimand was undefined and said so. Only a
    method with ``can_refuse`` can do this.
``intrinsic=None`` with ``estimates_intrinsic=False``
    The method never attempts an intrinsic term at all -- a compositional-only
    method is not silent about silencing, it is inapplicable to it. That is a
    different fact from refusing and is not a point in its favour.
``intrinsic`` is a float where the truth is undefined
    **False confidence.** The number is not wrong, which would at least be
    checkable; it is unfounded. This is what the claim above is about.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np

from src.estimator.kitagawa import decompose
from src.harness.positivity import CUTPOINTS, classify_estimability

#: Below this many mature tumour cells the intrinsic estimand is treated as
#: undefined. Not this benchmark's number -- W2's, from
#: ``harness.positivity``, so the method under test is scored against the
#: project's own committed rule rather than one chosen to flatter it.
NOT_ESTIMABLE_BELOW = CUTPOINTS.wide


@dataclass(frozen=True)
class MethodOutput:
    """One method's answer for one sample.

    ``intrinsic is None`` is not a small intrinsic term and must never be read
    as 0.0 -- CLAUDE.md invariant 1, which is the rule this whole benchmark is
    testing other methods against.
    """

    compositional: float | None
    intrinsic: float | None
    #: The method's own estimability verdict, or None if it has no vocabulary
    #: for one. None here means "could not have said", not "did not say".
    estimability: str | None = None
    #: True only when the method actively determined the estimand is undefined.
    refused: bool = False
    refusal_reason: str | None = None


class DecompositionMethod(ABC):
    """One decomposition method behind a fixed interface."""

    name: str = "abstract"
    #: Can this method express "undefined" at all? A method that structurally
    #: cannot is a different kind of object from one that can and chose not
    #: to, and only recording the difference makes the README claim testable.
    can_refuse: bool = False
    #: Does it attempt an intrinsic term? A compositional-only method does not,
    #: and must not be counted as having refused one.
    estimates_intrinsic: bool = True
    #: External runtime (R, Docker, a token). Reported, never silently skipped.
    requires_external: bool = False

    @abstractmethod
    def fit(self, sample, *, weighting: str = "normal") -> MethodOutput: ...

    def is_available(self) -> tuple[bool, str]:
        return True, "available"


def _mean_or_none(values: np.ndarray) -> float | None:
    """Mean of a possibly-empty set. None, never 0.0, when empty."""
    return float(values.mean()) if values.size else None


# ---------------------------------------------------------------------------
# The project's method, and the ablation that isolates what its gate buys
# ---------------------------------------------------------------------------


class KitagawaPositivityMethod(DecompositionMethod):
    """This project: Kitagawa standardisation behind the positivity gate."""

    name = "kitagawa+positivity"
    can_refuse = True

    def fit(self, sample, *, weighting: str = "normal") -> MethodOutput:
        verdict = classify_estimability(sample.n_mature_tumour)
        compositional = (sample.frac_mature_tumour - sample.frac_mature_normal) * (
            sample.mean_normal
        )
        if verdict == "not_estimable":
            return MethodOutput(
                compositional=compositional,
                intrinsic=None,
                estimability=verdict,
                refused=True,
                refusal_reason=(
                    f"{sample.n_mature_tumour} mature tumour cells is below the "
                    f"cutpoint of {NOT_ESTIMABLE_BELOW}; the per-cell mean has no "
                    "referent, so the intrinsic term is undefined -- not zero"
                ),
            )
        d = decompose(
            sample.frac_mature_normal,
            sample.frac_mature_tumour,
            sample.mean_normal,
            sample.mean_tumour,
            n_cells_mature=sample.n_mature_tumour,
            weighting=weighting,
        )
        return MethodOutput(d.compositional, d.intrinsic, estimability=verdict)


class KitagawaNoGateMethod(DecompositionMethod):
    """THE ABLATION. Identical arithmetic, gate removed.

    Every other comparison in this benchmark differs from ours in more than one
    way at once, so none of them can isolate what the estimability gate is
    responsible for. This one differs in exactly one way. If our method's
    advantage does not survive against this, it is not the gate doing the work.

    When no mature tumour cell exists it does what an ungated implementation
    does: takes the mean of an empty set as 0.0 and reports the result. That is
    the invariant-1 violation, written out.
    """

    name = "kitagawa-no-gate"
    can_refuse = False

    def fit(self, sample, *, weighting: str = "normal") -> MethodOutput:
        mean_tumour = sample.mean_tumour if sample.mean_tumour is not None else 0.0
        d = decompose(
            sample.frac_mature_normal,
            sample.frac_mature_tumour,
            sample.mean_normal,
            mean_tumour,
            n_cells_mature=max(sample.n_mature_tumour, 0),
            weighting=weighting,
        )
        return MethodOutput(d.compositional, d.intrinsic)


class VarianceGateKitagawaMethod(DecompositionMethod):
    """Kitagawa behind a WIDTH gate rather than a COUNT gate.

    The natural alternative primitive, and the fairest competitor in this file
    after the ablation. Instead of asking "are there enough mature tumour
    cells", it asks "is the interval on the intrinsic term narrow enough to
    resolve the effect the design says it must detect" -- and abstains when it
    is not.

    Both the effect size and the confidence level come from the SAME design
    document as the count cutpoint (a halving of per-cell output, s = 0.5), so
    neither gate is tuned against the other on this benchmark.

    THE FINITENESS GUARD IS THE WHOLE POINT. With no mature tumour cell the
    standard error is undefined, and ``nan > threshold`` is ``False`` in numpy
    and in pandas -- so a width gate written the obvious way SILENTLY DECLINES
    TO ABSTAIN in exactly the world the estimand does not exist. That is the
    same NaN coercion this project found inside its own calibration routine.
    We guard it explicitly here, which makes this competitor stronger than the
    version most people would write, and the comparison correspondingly fairer.
    """

    name = "kitagawa+variance-gate"
    can_refuse = True

    #: Detect a halving of per-cell output. The design document's own s = 0.5.
    DETECTABLE_SHIFT = 0.5
    #: Two-sided 95%, matching the interval the estimator already reports.
    Z = 1.959963984540054

    @staticmethod
    def _mature_moments(expr: np.ndarray, n_mature: int) -> tuple[float, float]:
        """Mean and variance among mature cells, without needing the mask.

        Immature cells express exactly 0.0, and a zero contributes nothing to
        either the sum or the sum of squares, so both moments over the mature
        subset are recoverable from the full arm given ``n_mature``.
        """
        total = float(expr.sum())
        total_sq = float((expr**2).sum())
        mean = total / n_mature
        var = max(total_sq / n_mature - mean**2, 0.0)
        return mean, var

    def fit(self, sample, *, weighting: str = "normal") -> MethodOutput:
        n_t = int(sample.n_mature_tumour)
        n_n = int(round(sample.frac_mature_normal * sample.expr_normal.size))

        compositional = (sample.frac_mature_tumour - sample.frac_mature_normal) * (
            sample.mean_normal
        )

        if n_t < 1 or n_n < 1:
            return MethodOutput(
                compositional=compositional,
                intrinsic=None,
                estimability="not_estimable",
                refused=True,
                refusal_reason=(
                    "no mature cell in one arm, so the standard error of the "
                    "per-cell mean is undefined -- not large. Reached only "
                    "because the finiteness of the SE is checked before it is "
                    "compared against the threshold."
                ),
            )

        mean_t, var_t = self._mature_moments(sample.expr_tumour, n_t)
        mean_n, var_n = self._mature_moments(sample.expr_normal, n_n)

        # SE of the intrinsic term f_N (m_T - m_N), treating f_N as fixed.
        se = sample.frac_mature_normal * float(np.sqrt(var_t / n_t + var_n / n_n))
        half_width = self.Z * se
        detectable = abs(
            sample.frac_mature_normal * (self.DETECTABLE_SHIFT - 1.0) * sample.mean_normal
        )

        if not np.isfinite(half_width) or half_width > detectable:
            return MethodOutput(
                compositional=compositional,
                intrinsic=None,
                estimability="not_estimable",
                refused=True,
                refusal_reason=(
                    f"half-width {half_width:.4g} exceeds the detectable effect "
                    f"{detectable:.4g}; the interval cannot resolve a halving"
                ),
            )

        d = decompose(
            sample.frac_mature_normal,
            sample.frac_mature_tumour,
            sample.mean_normal,
            sample.mean_tumour,
            n_cells_mature=n_t,
            weighting=weighting,
        )
        return MethodOutput(d.compositional, d.intrinsic, estimability="ok")


# ---------------------------------------------------------------------------
# What the literature actually does
# ---------------------------------------------------------------------------


class NaiveDeltaMeanMethod(DecompositionMethod):
    """"Difference the mature-cell means." The most common informal approach.

    No standardisation, so it never separates the two mechanisms: the whole
    change is reported as expression. Its compositional term is absent by
    construction rather than estimated as zero, so it is reported as None.
    """

    name = "naive-delta-mean"
    can_refuse = False

    def fit(self, sample, *, weighting: str = "normal") -> MethodOutput:
        mean_tumour = sample.mean_tumour if sample.mean_tumour is not None else 0.0
        return MethodOutput(compositional=None, intrinsic=mean_tumour - sample.mean_normal)


class PseudobulkDEMethod(DecompositionMethod):
    """Pseudobulk differential expression -- the "cell-intrinsic DEG" framing.

    Sums over ALL cells, mature and not, so a purely compositional change reads
    as an expression change. It is never undefined, because summing an empty
    mature set is still a well-formed sum over the arm. That is precisely why
    it cannot notice that the intrinsic question stopped having an answer.
    """

    name = "pseudobulk-de"
    can_refuse = False

    def fit(self, sample, *, weighting: str = "normal") -> MethodOutput:
        return MethodOutput(
            compositional=None,
            intrinsic=float(sample.expr_tumour.mean() - sample.expr_normal.mean()),
        )


class CompositionOnlyMethod(DecompositionMethod):
    """Milo / scCODA-shaped: abundance testing, no expression arm.

    A THIRD failure mode, and the one most easily mistaken for a virtue. It
    never invents an intrinsic number -- but not because it detected that the
    estimand was undefined. It has no intrinsic arm at all. Scoring it as
    "refused" would credit inapplicability as caution.
    """

    name = "composition-only"
    can_refuse = False
    estimates_intrinsic = False

    def fit(self, sample, *, weighting: str = "normal") -> MethodOutput:
        delta = sample.frac_mature_tumour - sample.frac_mature_normal
        return MethodOutput(compositional=delta * sample.mean_normal, intrinsic=None)


class CacoaMethod(DecompositionMethod):
    """cacoa (kharchenkolab), via an Rscript subprocess. NOT INSTALLED.

    Shipped unrunnable on purpose. ``env/w4_estimator.yml`` pins ``r-devtools``
    with the comment "cacoa installs from GitHub" and nothing has installed it,
    so the honest report is a named skip with a reason rather than silence.
    "cacoa was not run because it is not installed" and "cacoa was not run" are
    materially different statements and only one of them is reportable
    (``harness/deconvolve/base.py``'s ethos, applied here).
    """

    name = "cacoa"
    can_refuse = False
    requires_external = True

    def is_available(self) -> tuple[bool, str]:
        return False, (
            "cacoa is not installed; env/w4_estimator.yml pins r-devtools and cacoa "
            "installs from GitHub (kharchenkolab/cacoa). Not run, not skipped silently."
        )

    def fit(self, sample, *, weighting: str = "normal") -> MethodOutput:
        raise RuntimeError("cacoa is not available; check is_available() before fit()")


#: Everything runnable without R, Docker or a token. CI-safe by construction:
#: numpy and the project's own modules only.
DEFAULT_METHODS: tuple[DecompositionMethod, ...] = (
    KitagawaPositivityMethod(),
    KitagawaNoGateMethod(),
    VarianceGateKitagawaMethod(),
    NaiveDeltaMeanMethod(),
    PseudobulkDEMethod(),
    CompositionOnlyMethod(),
    CacoaMethod(),
)


def available_methods(
    methods: tuple[DecompositionMethod, ...] = DEFAULT_METHODS,
) -> tuple[list[DecompositionMethod], dict[str, str]]:
    """Split into runnable methods and {name: why not}. Never drops silently."""
    runnable, skipped = [], {}
    for method in methods:
        ok, why = method.is_available()
        (runnable.append(method) if ok else skipped.update({method.name: why}))
    return runnable, skipped
