"""Does the estimator distinguish granularity rungs? W2.

WHY THIS EXISTS
---------------
README design decision 3 and execution_plan.md §6.2 promise the split reported
"as a curve across four resolutions", with divergence between rungs as a
contribution. W1's pilot then found that on axis 1, ``lineage`` and
``crypt_position`` are **the same partition** in all ten arms, and ``best4`` is
unusable at sensitivity 0.04. A four-point curve with two identical points and
one missing one is not a curve.

Before anyone concludes the granularity knob does nothing *in the biology*, the
estimator has to be cleared. There are two very different explanations for two
rungs agreeing:

1. **The partitions are identical**, so any estimator must agree. Then the
   finding is about the labels, and it is real.
2. **The estimator cannot resolve a difference that is there.** Then the finding
   is about us, and reporting it as biology would be wrong.

This module separates those. It re-estimates *the same generated sample* under
different partitions of the same cells — no regeneration, so nothing but the
partition changes — and reports whether the recovered terms move when the truth
says they should.

WHAT A RESULT HERE MEANS
------------------------
``rung_separation`` on a cohort where two rungs genuinely differ should recover
a difference close to the true one. If it does, explanation 2 is excluded and
W1's degeneracy stands as a statement about the labelling. If it does not, the
curve was never measurable with this estimator and that is a methods finding
that belongs in the paper rather than in a footnote.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd

from src.estimator.kitagawa import decompose
from src.harness.positivity import classify_estimability
from src.harness.pseudobulk import PseudobulkSample


def rung_mature_mask(
    sample: PseudobulkSample, arm: str, mature_types: Sequence[str]
) -> np.ndarray:
    """Which drawn cells count as mature under a given rung's definition.

    A rung is just a set of cell-type labels that the rung calls mature. The
    coarse rung names more types than the fine one; identical sets mean
    identical partitions, which is exactly the degeneracy under investigation.
    """
    labels = sample.drawn_cell_type.get(arm)
    if labels is None:
        raise KeyError(
            f"sample carries no drawn_cell_type for arm {arm!r}; it was generated "
            f"before that field existed"
        )
    return np.isin(labels, list(mature_types))


def estimate_under_rung(
    sample: PseudobulkSample,
    gene: str,
    mature_types: Sequence[str],
    *,
    weighting: str = "normal",
) -> dict[str, float | str | None]:
    """Decompose ``sample`` treating ``mature_types`` as the mature compartment.

    Goes through W4's ``decompose()`` rather than re-deriving the arithmetic, so
    a change to the estimator shows up here.
    """
    values = sample.drawn_expression.get(gene)
    if values is None:
        raise KeyError(f"sample carries no drawn expression for {gene!r}")

    masks = {arm: rung_mature_mask(sample, arm, mature_types) for arm in ("normal", "tumour")}
    fracs = {arm: float(mask.mean()) for arm, mask in masks.items()}
    means = {
        arm: float(values[arm][masks[arm]].mean()) if masks[arm].any() else 0.0
        for arm in ("normal", "tumour")
    }
    n_mature = int(masks["tumour"].sum())
    estimability = classify_estimability(n_mature)

    d = decompose(
        fracs["normal"], fracs["tumour"], means["normal"], means["tumour"],
        n_cells_mature=n_mature, weighting=weighting,
    )
    return {
        "frac_mature_normal": fracs["normal"],
        "frac_mature_tumour": fracs["tumour"],
        "mean_normal": means["normal"],
        "mean_tumour": means["tumour"],
        "n_cells_mature": n_mature,
        "estimability": estimability,
        "compositional": d.compositional,
        # Undefined, not zero, when the compartment is too thin — invariant 1.
        "intrinsic": None if estimability == "not_estimable" else d.intrinsic,
        "interaction": d.interaction,
    }


def rung_separation(
    sample: PseudobulkSample,
    gene: str,
    rungs: Mapping[str, Sequence[str]],
    *,
    weighting: str = "normal",
) -> pd.DataFrame:
    """One row per rung: the decomposition under that rung's mature definition.

    ``rungs`` maps a rung name to the cell types it calls mature, e.g.::

        {"lineage": ["mature_colonocyte", "crypt_top"],   # coarse: pools both
         "crypt_position": ["crypt_top"]}                  # fine: crypt top only

    Identical value sets produce identical rows, by construction — that is the
    degeneracy case and the test suite pins it.
    """
    if len(rungs) < 2:
        raise ValueError("rung_separation needs at least two rungs to compare")
    rows = []
    for name, mature_types in rungs.items():
        row = {"rung": name, "gene": gene, "weighting": weighting}
        row.update(estimate_under_rung(sample, gene, mature_types, weighting=weighting))
        rows.append(row)
    out = pd.DataFrame(rows)
    # Nullable Float64, so an unestimable intrinsic term stays <NA> rather than
    # being coerced to a bare NaN alongside real zeros. Same convention the
    # frozen schema uses, and for the same reason: None is not 0.0.
    for col in ("compositional", "intrinsic", "interaction"):
        out[col] = out[col].astype("Float64")
    return out


def separation_summary(separation: pd.DataFrame, term: str = "intrinsic") -> pd.DataFrame:
    """Pairwise: how far apart are the rungs' estimates of ``term``?

    ``relative`` scales the absolute gap by the larger of the two magnitudes, so
    it reads as "what fraction of the effect does the granularity choice move".
    Near 0 means the rungs agree; near or above 1 means the choice of rung is
    doing as much work as the biology, which is the divergence §6.2 calls the
    contribution.

    A pair where either side is ``not_estimable`` yields ``None`` rather than a
    number — the terms are not comparable when one of them does not exist.
    """
    rows = []
    records = separation.to_dict("records")
    for i, a in enumerate(records):
        for b in records[i + 1 :]:
            va, vb = a[term], b[term]
            if va is None or vb is None or pd.isna(va) or pd.isna(vb):
                absolute = relative = None
            else:
                absolute = abs(va - vb)
                scale = max(abs(va), abs(vb))
                relative = absolute / scale if scale else 0.0
            rows.append(
                {
                    "term": term,
                    "rung_a": a["rung"],
                    "rung_b": b["rung"],
                    "value_a": va,
                    "value_b": vb,
                    "absolute": absolute,
                    "relative": relative,
                    "identical_partition": (
                        a["n_cells_mature"] == b["n_cells_mature"]
                        and a["frac_mature_tumour"] == b["frac_mature_tumour"]
                    ),
                }
            )
    return pd.DataFrame(rows)


def estimator_can_separate(
    separation: pd.DataFrame,
    *,
    term: str = "intrinsic",
    threshold: float = 0.10,
) -> bool:
    """True if any pair of rungs differs by more than ``threshold``, relatively.

    Use it on a cohort constructed so the rungs *do* differ. A False there means
    the estimator cannot resolve a difference known to be present, and W1's
    observed degeneracy cannot be attributed to the labels until that is fixed.
    """
    summary = separation_summary(separation, term=term)
    usable = summary["relative"].dropna()
    return bool(len(usable) and usable.max() > threshold)
