"""How much of the signal survives ambient contamination. W2, handoff §5 task 7.

G1 asks whether ambient correction eliminates the intrinsic signal. Since
[open_decisions #11] the correction the criterion assumed cannot be run — GSE178341
ships no unfiltered droplets, so CellBender has nothing to learn from and SoupX
runs degraded. The gate is therefore choosing between *degraded correction* and
*no correction*, and the memo's §0 has been saying that "degraded" is a word and
not a number.

This turns it into a number. It does not tell the gate whether this cohort's
correction worked — the harness starts from corrected cells and cannot know that.
It tells the gate **what a given level of residual contamination does to the
decomposition**, which is the question that makes the choice decidable.

THE MECHANISM, WHICH IS WHY THE ANSWER IS NOT SYMMETRIC
------------------------------------------------------
Ambient RNA is a shared background of lysed-cell transcripts added to every
barcode. Its composition is the *sample's own* average expression — so the soup
in a normal sample is rich in mature-colonocyte transcripts, because that sample
is full of mature colonocytes, and the soup in a depleted tumour is not.

Now take a tumour whose loss is **purely compositional**: mature cells are gone,
the survivors are untouched, and the true intrinsic term is exactly zero. Under
contamination each arm's mature cells are pulled toward *their own arm's* soup.
The normal arm is pulled up and the tumour arm is pulled up less, because there
is less of the gene in the tumour's soup. The mature-cell means therefore
separate **where the truth says they should not**, and a compositional-only world
acquires an apparent intrinsic term.

That is the exact failure mode invariant 6 and gate G1 exist to guard, and it
runs in the direction of the project's prior hypothesis' opposite — it
manufactures *intrinsic* signal out of *compositional* truth. Measuring which way
it leans, and how far, is the point.

WHAT THIS DOES NOT MODEL, STATED SO IT IS NOT ASSUMED
-----------------------------------------------------
Contamination here perturbs **expression**, not **labels**. In a real pipeline
enough soup also moves cells across a maturity threshold, which would change the
fractions as well as the means. That is a labelling question, it belongs to W1's
axes, and pretending to simulate it here would produce a number with W1's
uncertainty hidden inside W2's. The fractions are held at their true values, so
every effect reported below is a *lower bound* on ambient's total damage.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Final

import numpy as np
import pandas as pd

from src.harness.pseudobulk import PseudobulkSample, generate_pseudobulk

#: Contamination fractions swept. Dense below 0.15 because that is where W1's
#: measured cohort sits and where the exclusion threshold was drawn; the tail is
#: carried so the shape of the curve is visible rather than extrapolated.
DEFAULT_AMBIENT_GRID: Final = (0.0, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30)

#: Decision #16, pre-committed 2026-08-23: measure ambient, do not correct it,
#: and exclude a sample above 10%. The sweep reports what that threshold buys.
W1_EXCLUSION_THRESHOLD: Final = 0.10

TERMS: Final = ("compositional", "intrinsic", "interaction")


@dataclass(frozen=True)
class AmbientRegime:
    """A world with a known decomposition, so contamination has a truth to damage."""

    name: str
    composition_normal: Mapping[str, float]
    composition_tumour: Mapping[str, float]
    shift: Mapping[str, float]
    why: str

    def parametric_zero_terms(self, target_gene: str, mature_label: str) -> set[str]:
        """Which terms this world sets to **exactly** zero, by construction.

        Read off the design, never off the realised numbers. The realised
        intrinsic term in a compositional-only world is zero *plus sampling
        noise* — around -0.05 against a compositional term of -17 — so a
        tolerance test calls it non-zero and then forms a retention ratio with a
        near-zero denominator. That produces a confident 2.4x where the honest
        statement is "the truth here is zero and 0.79 appeared".

        This is the same mistake as scoring coverage against realised truth
        instead of parametric truth (handoff §3), one module along.
        """
        zero = set()
        if self.shift.get(target_gene, 1.0) == 1.0:
            zero.add("intrinsic")
        if self.composition_normal.get(mature_label) == self.composition_tumour.get(
            mature_label
        ):
            zero.add("compositional")
        if len(zero) == 2:  # nothing moves, so the cross term cannot either
            zero.add("interaction")
        elif "intrinsic" in zero or "compositional" in zero:
            # The interaction is a product of the two changes; if either factor
            # is exactly zero, so is the cross term.
            zero.add("interaction")
        return zero


def default_regimes(
    target_gene: str,
    *,
    mature_label: str = "mature_colonocyte",
    other_types: Sequence[str] = ("stem", "stromal", "immune"),
) -> tuple[AmbientRegime, ...]:
    """The three worlds worth separating: compositional only, intrinsic only, both.

    ``compositional_only`` is the important one. Its true intrinsic term is
    **exactly zero**, so anything the estimator reports there under contamination
    was manufactured, and its size is directly comparable against the
    compositional term it was manufactured from.
    """

    def composition(mature: float) -> dict[str, float]:
        rest = (1.0 - mature) / len(other_types)
        return {mature_label: mature, **{t: rest for t in other_types}}

    return (
        AmbientRegime(
            name="compositional_only",
            composition_normal=composition(0.40),
            composition_tumour=composition(0.10),
            shift={target_gene: 1.0},
            why="mature cells gone, survivors untouched. True intrinsic is EXACTLY "
            "zero, so whatever appears under contamination is manufactured",
        ),
        AmbientRegime(
            name="intrinsic_only",
            composition_normal=composition(0.40),
            composition_tumour=composition(0.40),
            shift={target_gene: 0.35},
            why="mature cells present and silenced. True compositional is exactly "
            "zero; this is the arm the project's contribution rests on",
        ),
        AmbientRegime(
            name="both",
            composition_normal=composition(0.40),
            composition_tumour=composition(0.15),
            shift={target_gene: 0.50},
            why="the realistic mixture, where the interaction term is also non-zero",
        ),
    )


# ---------------------------------------------------------------------------
# The contamination itself
# ---------------------------------------------------------------------------


def soup_rate(sample: PseudobulkSample, gene: str, arm: str) -> float:
    """The ambient profile's per-cell rate for one gene, in one arm.

    The soup is the sample's own average expression, because that is what lysed
    cells release: every drawn cell in the arm, not only the mature ones. This is
    the single line that makes the two arms' soups differ, and therefore the
    single line the whole finding rests on.
    """
    values = sample.drawn_expression[gene][arm]
    return float(np.mean(values)) if len(values) else 0.0


def contaminate(
    values: np.ndarray,
    *,
    fraction: float,
    rate: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """``(1-c)`` of each cell's own counts, plus ``c`` of the soup. On counts.

    Binomial thinning for the retained part and a Poisson draw for the ambient
    part, matching ``pseudobulk._apply_shift``'s idiom rather than scaling floats
    — a count matrix that stops being integers stops being able to go back
    through the generator, and the near-zero cells are exactly where rounding
    would do its damage.
    """
    if not 0.0 <= fraction <= 1.0:
        raise ValueError(f"fraction={fraction} outside [0, 1]")
    if rate < 0:
        raise ValueError(f"rate={rate} is negative")
    if fraction == 0.0:
        return np.asarray(values).copy()
    retained = rng.binomial(np.asarray(values).astype(np.int64), 1.0 - fraction)
    ambient = rng.poisson(np.full(len(retained), fraction * rate))
    return retained + ambient


def contaminated_stats(
    sample: PseudobulkSample,
    gene: str,
    *,
    fraction: float,
    rng: np.random.Generator,
) -> dict[str, float]:
    """The four scalars ``decompose()`` takes, after contamination.

    The **fractions are the true ones** — see the module docstring on what is not
    modelled. Only the mature-cell means move, which is why every number this
    produces is a lower bound on ambient's damage.
    """
    stats = dict(sample.truth.realised_stats[gene])
    for arm, key in (("normal", "mean_normal"), ("tumour", "mean_tumour")):
        mature = sample.drawn_expression[gene][arm][sample.drawn_is_mature[arm]]
        if len(mature) == 0:
            continue
        dirty = contaminate(
            mature, fraction=fraction, rate=soup_rate(sample, gene, arm), rng=rng
        )
        stats[key] = float(np.mean(dirty))
    return stats


# ---------------------------------------------------------------------------
# The sweep
# ---------------------------------------------------------------------------


def ambient_sweep(
    counts: np.ndarray,
    cell_type: Sequence[str],
    patient_id: Sequence[str],
    genes: Sequence[str],
    *,
    target_gene: str,
    seed: int,
    held_out_patients: Sequence[str],
    regimes: Sequence[AmbientRegime] | None = None,
    grid: Sequence[float] = DEFAULT_AMBIENT_GRID,
    n_cells: int = 2000,
    n_replicates: int = 20,
    weighting: str = "normal",
    mature_label: str = "mature_colonocyte",
) -> pd.DataFrame:
    """One row per (regime, ambient fraction, replicate, term).

    Each replicate generates a **clean** sample, decomposes it to get the clean
    terms, then contaminates the same drawn cells at every grid point. Sharing
    the sample across the grid is deliberate: the difference between two grid
    points is then contamination alone, not contamination plus a fresh draw.
    """
    from src.estimator.kitagawa import decompose

    regimes = tuple(regimes or default_regimes(target_gene, mature_label=mature_label))
    rows: list[dict] = []

    for regime in regimes:
        for replicate in range(n_replicates):
            sample = generate_pseudobulk(
                counts,
                cell_type,
                patient_id,
                genes,
                composition_normal=dict(regime.composition_normal),
                composition_tumour=dict(regime.composition_tumour),
                shift=dict(regime.shift),
                held_out_patients=held_out_patients,
                n_cells=n_cells,
                seed=seed + replicate,
                mature_label=mature_label,
            )
            clean = sample.truth.realised_stats[target_gene]
            clean_terms = decompose(
                clean["frac_mature_normal"],
                clean["frac_mature_tumour"],
                clean["mean_normal"],
                clean["mean_tumour"],
                n_cells_mature=sample.truth.n_cells_mature,
                weighting=weighting,
            )
            zero_terms = regime.parametric_zero_terms(target_gene, mature_label)
            for fraction in grid:
                rng = np.random.default_rng(
                    seed + replicate * 1000 + int(round(fraction * 1000))
                )
                dirty = contaminated_stats(
                    sample, target_gene, fraction=fraction, rng=rng
                )
                dirty_terms = decompose(
                    dirty["frac_mature_normal"],
                    dirty["frac_mature_tumour"],
                    dirty["mean_normal"],
                    dirty["mean_tumour"],
                    n_cells_mature=sample.truth.n_cells_mature,
                    weighting=weighting,
                )
                for term in TERMS:
                    rows.append(
                        {
                            "regime": regime.name,
                            "ambient_fraction": float(fraction),
                            "replicate": replicate,
                            "weighting": weighting,
                            "term": term,
                            "truth_is_zero": term in zero_terms,
                            "value_clean": float(getattr(clean_terms, term)),
                            "value_contaminated": float(getattr(dirty_terms, term)),
                            "n_cells_mature": int(sample.truth.n_cells_mature),
                            "seed": seed + replicate,
                        }
                    )
    return pd.DataFrame(rows)


#: Below this, a term's clean value is treated as "the truth says zero here" and
#: a retention ratio is not formed. Dividing by a near-zero denominator is how a
#: sweep produces a spectacular number that means nothing.
ZERO_TERM_TOLERANCE: Final = 1e-9


def summarise_ambient(sweep: pd.DataFrame) -> pd.DataFrame:
    """Retention where the truth is non-zero, manufacture where it is zero.

    Two different questions, and a single "recovery" column would conflate them:

    - **retention** — ``contaminated / clean`` for a term the truth says exists.
      1.0 is untouched; below 1.0 the term is being attenuated.
    - **manufactured** — the term's absolute contaminated value where the truth
      says **zero**, expressed against the largest true term in the same regime
      so it is readable as "this much signal, out of nothing".

    Both are medians over replicates. The mean is the wrong summary here: a
    single replicate whose clean term lands near zero produces an enormous ratio
    and would set the mean on its own.
    """
    if "truth_is_zero" not in sweep.columns:
        raise ValueError(
            "sweep has no truth_is_zero column. It comes from the regime's "
            "design, and inferring it from the realised values instead is the "
            "bug this column exists to prevent — see "
            "AmbientRegime.parametric_zero_terms."
        )
    # The scale a manufactured term is read against is the largest term the
    # regime's design says is REAL — not the largest observed, which would let a
    # manufactured term normalise itself.
    real = sweep[~sweep["truth_is_zero"].astype(bool)]
    scale = (
        real.assign(magnitude=real["value_clean"].abs())
        .groupby(["regime", "replicate"], observed=True)["magnitude"]
        .max()
        .rename("regime_scale")
    )
    frame = sweep.merge(scale, on=["regime", "replicate"], how="left")

    is_zero = frame["truth_is_zero"].astype(bool)
    frame["retention"] = np.where(
        is_zero, np.nan, frame["value_contaminated"] / frame["value_clean"]
    )
    frame["manufactured"] = np.where(
        is_zero & (frame["regime_scale"] > ZERO_TERM_TOLERANCE),
        frame["value_contaminated"].abs() / frame["regime_scale"],
        np.nan,
    )
    out = (
        frame.groupby(["regime", "term", "ambient_fraction"], observed=True)
        .agg(
            n_replicates=("replicate", "nunique"),
            truth_is_zero=("truth_is_zero", "first"),
            median_clean=("value_clean", "median"),
            median_contaminated=("value_contaminated", "median"),
            median_retention=("retention", "median"),
            median_manufactured=("manufactured", "median"),
        )
        .reset_index()
    )
    return out.sort_values(["regime", "term", "ambient_fraction"]).reset_index(drop=True)


def cost_at_threshold(
    summary: pd.DataFrame, *, fraction: float = W1_EXCLUSION_THRESHOLD
) -> pd.DataFrame:
    """The single slice the gate needs: what decision #16's 10% cap still admits.

    A sample at exactly the exclusion threshold is *kept*, so this is the
    worst case the cohort tolerates by design rather than a hypothetical.
    """
    return summary[np.isclose(summary["ambient_fraction"], fraction)].reset_index(drop=True)
