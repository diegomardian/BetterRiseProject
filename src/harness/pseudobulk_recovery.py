"""Does the reference RECIPE cost the instrument its margin? A1.

WHAT THIS IS AND IS NOT. This is a decision experiment, not a result. It cannot
say anything about colorectal cancer. Its entire output is one number -- the
recovery gap between a reference built as ``mean(CP10K)`` and one built as
``mean(log1p(CP10K))`` -- and that number decides one thing:

    gap >= ~0.02   the Stage 4 gate failure (0.462, 0.479 against 0.5) is an
                   artifact of the reference recipe, and the W1 linear rebuild
                   is justified
    gap negligible the recipe is not the problem, the gate failure is a
                   bulk-deconvolution limit, and the variance arm is
                   unreachable by this route

Both branches are instrument-level. Neither is a biological finding, and a
confirmed variance arm downstream would still be consistency with the
compositional account rather than mechanism (invariant 6, and the locked
prespec says so in its own words).

WHY PSEUDOBULK GROUND TRUTH BEATS THE GATE'S OWN. The Stage 4 gate correlates
the deconvolved non-epithelial fraction against (1 - ABSOLUTE purity), and
those are not the same quantity -- purity is the malignant-cell share, not the
epithelial share, so the gate carries a definitional ceiling below 1.0 that
nobody has measured. Mixing cells into fractions WE set removes that confound
entirely: the truth is the truth, and any shortfall is the instrument.

TWO LEGS, because one of them is a confound check.

    SMC -> KUL3   cross-cohort, like TCGA. This is the leg that matters.
    SMC -> SMC    same cohort, held-out patients. If this recovers well and the
                  cross-cohort leg collapses, the loss is BATCH rather than
                  recipe. Batch hits both recipes roughly equally so the gap
                  should survive it, but reading a gap without knowing which
                  regime you are in is how a batch effect gets published as a
                  recipe effect.

WHAT WILL NOT TRANSFER. Absolute correlations here will beat TCGA's, because
pseudobulk has no ambient contamination, no dropout beyond the cells' own, and
no library preparation between the cells and the mixture. **The gap transfers;
the absolutes do not.** An r of 0.9 here is not a prediction that the gate
passes at 0.5 on real bulk.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)


class RecoveryError(RuntimeError):
    """The recovery experiment cannot be run on these inputs."""


#: Lee's `Cell_type` vocabulary -> the compartment axis the gate reads. Lee has
#: no endothelial label; its `Stromal cells` covers that compartment, which is
#: why this collapses to epithelial / non-epithelial rather than pretending to
#: the S matrix's six types.
EPITHELIAL_LABELS: tuple[str, ...] = ("Epithelial cells",)

#: Reference recipes. `log1p` is `build_signature`'s committed choice -- right
#: for marker SELECTION, and the thing under test for the mixture model.
RECIPES: tuple[str, ...] = ("linear", "log1p")


def build_reference(
    expression: pd.DataFrame,
    cell_type: pd.Series,
    *,
    recipe: str,
) -> pd.DataFrame:
    """Per-compartment profile on one scale. Genes x compartments.

    ``expression`` is CP10K, linear, cells x genes -- what ``load_lee_cohort``
    emits. ``linear`` takes the arithmetic mean, which is the quantity a linear
    mixture actually sums. ``log1p`` reproduces ``build_signature``: the mean of
    ``log1p``, then ``expm1`` back, which is a geometric mean of ``CP10K + 1``
    and biased low by Jensen, worst for the most dispersed genes.

    Both recipes see the SAME cells and the SAME genes. That is the whole
    design: any difference downstream is the recipe and nothing else.
    """
    if recipe not in RECIPES:
        raise RecoveryError(f"recipe must be one of {RECIPES}, got {recipe!r}")
    labels = cell_type.reindex(expression.index)
    if labels.isna().any():
        raise RecoveryError("some cells carry no compartment label")

    if recipe == "linear":
        profile = expression.groupby(labels, observed=True).mean()
    else:
        profile = np.expm1(np.log1p(expression).groupby(labels, observed=True).mean())
    return profile.T.astype(float)


def make_pseudobulk(
    expression: pd.DataFrame,
    cell_type: pd.Series,
    *,
    n_samples: int,
    n_cells: int,
    seed: int,
    rng_alpha: float = 1.0,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Mix cells into known compartment fractions. Returns (bulk, truth).

    The mixture is a mean over sampled cells on the LINEAR CP10K scale, which
    is what a bulk library is: every cell contributes its transcripts and the
    sequencer counts the pool. Fractions are drawn per sample so the truth
    varies across the cohort -- a constant truth would make every correlation
    undefined, which is the degenerate case this whole project keeps finding.
    """
    labels = cell_type.reindex(expression.index)
    compartments = sorted(labels.dropna().unique())
    if len(compartments) < 2:
        raise RecoveryError(f"need at least two compartments, got {compartments}")

    rng = np.random.default_rng(seed)
    by_compartment = {c: np.flatnonzero((labels == c).to_numpy()) for c in compartments}
    thin = [c for c, idx in by_compartment.items() if len(idx) < 20]
    if thin:
        raise RecoveryError(f"compartments with under 20 cells: {thin}")

    values = expression.to_numpy(dtype=float)
    truth_rows, bulk_rows, names = [], [], []
    for i in range(n_samples):
        fractions = rng.dirichlet([rng_alpha] * len(compartments))
        counts = rng.multinomial(n_cells, fractions)
        picked = np.concatenate([
            rng.choice(by_compartment[c], size=k, replace=True)
            for c, k in zip(compartments, counts, strict=True) if k > 0
        ])
        bulk_rows.append(values[picked].mean(axis=0))
        # The REALISED fraction, not the requested one. The multinomial draw
        # does not land on the Dirichlet exactly, and scoring against the
        # request rather than the draw is the defect this project documents at
        # length (the recovery curve, WMHS section 2).
        truth_rows.append(counts / counts.sum())
        names.append(f"PB{i:04d}")

    bulk = pd.DataFrame(bulk_rows, index=names, columns=expression.columns)
    truth = pd.DataFrame(truth_rows, index=names, columns=compartments)
    return bulk, truth


def non_epithelial(frame: pd.DataFrame) -> pd.Series:
    """The gate's quantity: everything that is not epithelium.

    Refuses a non-numeric column rather than summing it. A pivoted result still
    carries its `method` label, and "everything that is not epithelium" happily
    includes a string -- which fails loudly here but would sum silently if the
    label happened to be numeric.
    """
    columns = [c for c in frame.columns if c not in EPITHELIAL_LABELS]
    if not columns:
        raise RecoveryError(f"no non-epithelial columns in {list(frame.columns)}")
    subset = frame[columns]
    bad = [c for c in columns if not pd.api.types.is_numeric_dtype(subset[c])]
    if bad:
        raise RecoveryError(
            f"non-numeric column(s) {bad} would be summed into the "
            f"non-epithelial fraction. Restrict to compartment columns first."
        )
    return subset.sum(axis=1)


@dataclass(frozen=True)
class RecoveryResult:
    """One (leg, recipe, method) cell of the experiment."""

    leg: str
    recipe: str
    method: str
    n_samples: int
    r_non_epithelial: float
    r_epithelial: float
    rmse_epithelial: float
    epithelial_exact_zero: int

    def as_row(self) -> dict:
        return {
            "leg": self.leg, "recipe": self.recipe, "method": self.method,
            "n_samples": self.n_samples,
            "r_non_epithelial": self.r_non_epithelial,
            "r_epithelial": self.r_epithelial,
            "rmse_epithelial": self.rmse_epithelial,
            "epithelial_exact_zero": self.epithelial_exact_zero,
        }


def _correlate(a: np.ndarray, b: np.ndarray) -> float:
    if np.ptp(a) == 0 or np.ptp(b) == 0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def recover(
    bulk: pd.DataFrame,
    truth: pd.DataFrame,
    reference: pd.DataFrame,
    *,
    leg: str,
    recipe: str,
    methods=None,
) -> list[RecoveryResult]:
    """Deconvolve and score against the realised truth."""
    from src.bulk.deconvolution import LINEAR_CP10K, Reference, deconvolve_cohort
    from src.harness.deconvolve.nnls import NNLSDeconvolver
    from src.harness.deconvolve.nusvr import NuSVRDeconvolver

    # Both recipes are declared LINEAR: the log1p one has already been expm1'd
    # back, which is exactly the repair `--linearise-reference` applies to the
    # committed matrices. That is what is under test, not a raw log reference.
    wrapped = Reference(
        matrix=reference, rung="compartment", scale=LINEAR_CP10K,
        source=__import__("pathlib").Path(f"{leg}:{recipe}"),
    )
    long, _ = deconvolve_cohort(
        bulk, wrapped, bulk_scale=LINEAR_CP10K,
        methods=list(methods or [NNLSDeconvolver(), NuSVRDeconvolver()]),
    )
    wide = long.pivot_table(
        index=["sample_id", "method"], columns="cell_type", values="fraction"
    ).reset_index()

    missing = [c for c in truth.columns if c not in wide.columns]
    if missing:
        raise RecoveryError(
            f"the deconvolution returned no column for {missing}, so those "
            f"compartments cannot be scored against their own truth."
        )
    out = []
    for method, group in wide.groupby("method"):
        # Restrict to the truth's compartments: the pivot still carries the
        # `method` label, and it is not a cell type.
        got = group.set_index("sample_id").reindex(truth.index)[list(truth.columns)]
        epi = EPITHELIAL_LABELS[0]
        got_epi = got[epi].to_numpy(dtype=float) if epi in got else np.zeros(len(truth))
        true_epi = truth[epi].to_numpy(dtype=float)
        out.append(RecoveryResult(
            leg=leg, recipe=recipe, method=str(method), n_samples=len(truth),
            r_non_epithelial=_correlate(
                non_epithelial(got).to_numpy(dtype=float),
                non_epithelial(truth).to_numpy(dtype=float),
            ),
            r_epithelial=_correlate(got_epi, true_epi),
            rmse_epithelial=float(np.sqrt(np.mean((got_epi - true_epi) ** 2))),
            epithelial_exact_zero=int((got_epi == 0.0).sum()),
        ))
    return out


def recipe_gap(results: list[RecoveryResult]) -> pd.DataFrame:
    """linear minus log1p, per (leg, method). The number the experiment is for."""
    frame = pd.DataFrame([r.as_row() for r in results])
    wide = frame.pivot_table(
        index=["leg", "method"], columns="recipe",
        values=["r_non_epithelial", "r_epithelial", "epithelial_exact_zero"],
    )
    out = pd.DataFrame(index=wide.index)
    for metric in ("r_non_epithelial", "r_epithelial", "epithelial_exact_zero"):
        if (metric, "linear") in wide and (metric, "log1p") in wide:
            out[f"{metric}_linear"] = wide[(metric, "linear")]
            out[f"{metric}_log1p"] = wide[(metric, "log1p")]
            out[f"{metric}_gap"] = wide[(metric, "linear")] - wide[(metric, "log1p")]
    return out.reset_index()
