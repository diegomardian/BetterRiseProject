"""W3.4 — batch and technical structure. Measure it; do not correct it.

CLAUDE.md invariant 4 and the brief both say the same thing: **document batch
structure, test its confounding with stage and MMR, but do not correct the
expression matrix.** Between-dataset variation is where the compositional signal
lives. Nothing in this module writes an expression matrix, and that is on
purpose.

WHAT COUNTS AS A BATCH IN TCGA
------------------------------
There is no explicit batch column. What exists is encoded in the aliquot
barcode: tissue source site (which hospital), plate (which 96-well plate the
aliquot was processed on), sequencing centre, analyte and vial. **Plate is the
usual batch proxy** — it is the finest processing unit that groups samples
handled together — and TSS is the usual site proxy. Both come from
``src/bulk/gdc.py``'s barcode parser, so if that is wrong this whole analysis is
wrong; it has its own tests for that reason.

THE TRAP THIS MODULE IS BUILT AROUND
------------------------------------
A categorical factor with many levels explains variance *by construction*. TSS
has ~50 levels and plate ~100 across 624 tumours; fit either as a one-way ANOVA
and you get an impressive R-squared that means nothing. Three guards:

1. **Adjusted R-squared**, which penalises by the number of levels.
2. **A permutation null.** Shuffle the factor labels, keeping the level sizes,
   and recompute. What a factor explains *above its own null* is the only
   number worth reading, and it is reported as ``excess_over_null``.
3. **Level counts are always printed next to the statistic**, so a reader can
   see what they are being asked to believe.

The same problem afflicts the confounding tests: a chi-squared on a 100x4 table
with most cells empty has no valid asymptotic null. Those use permutation
p-values throughout, never the chi-squared approximation.

WHAT THE OUTPUT IS NOT
----------------------
The variance figures do not sum to 100% and are not a partition. The factors are
correlated with each other — plates sit inside sites, sites correlate with
stage — so each figure is a marginal association, and adding them up would
double-count. A real variance partition would need a mixed model with all
factors jointly, which is more machinery than a "document it" task warrants.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

#: Permutation replicates. 999 gives a smallest attainable p of 0.001 and runs
#: in seconds at this cohort size.
N_PERMUTATIONS = 999

#: PCs retained for the variance analysis: enough to cover this much of the
#: total variance. The tail is mostly noise and including it dilutes every
#: factor equally, which flatters nothing but wastes time.
PC_VARIANCE_TARGET = 0.80

#: Below this many samples in a level, the level is pooled into "other". A plate
#: with two samples contributes a perfectly-fit group and inflates R-squared.
MIN_LEVEL_SIZE = 5


class BatchError(RuntimeError):
    """A batch-structure analysis could not be run."""


# ---------------------------------------------------------------------------
# Categorical association
# ---------------------------------------------------------------------------


def cramers_v(a: pd.Series, b: pd.Series, *, bias_correct: bool = True) -> float:
    """Cramer's V between two categorical series, on complete pairs only.

    Bias-corrected by default (Bergsma 2013). Uncorrected V is badly upward-biased
    on sparse tables, which is exactly what a plate-by-stage table is: the
    correction is not optional here, it is the difference between "confounded"
    and "not".
    """
    pair = pd.concat([a.rename("a"), b.rename("b")], axis=1).dropna()
    if len(pair) < 2 or pair["a"].nunique() < 2 or pair["b"].nunique() < 2:
        return float("nan")

    table = pd.crosstab(pair["a"], pair["b"]).to_numpy(dtype=float)
    n = table.sum()
    row = table.sum(axis=1, keepdims=True)
    col = table.sum(axis=0, keepdims=True)
    expected = row @ col / n
    chi2 = float(((table - expected) ** 2 / expected).sum())
    r, k = table.shape

    if not bias_correct:
        return float(np.sqrt(chi2 / (n * (min(r, k) - 1))))

    phi2 = max(0.0, chi2 / n - (r - 1) * (k - 1) / (n - 1))
    r_hat = r - (r - 1) ** 2 / (n - 1)
    k_hat = k - (k - 1) ** 2 / (n - 1)
    denominator = min(r_hat - 1, k_hat - 1)
    if denominator <= 0:
        return float("nan")
    return float(np.sqrt(phi2 / denominator))


def permutation_p(
    a: pd.Series,
    b: pd.Series,
    *,
    n_permutations: int = N_PERMUTATIONS,
    seed: int,
) -> tuple[float, float, float]:
    """Permutation test for association. Returns (observed V, null mean V, p).

    Permutation rather than a chi-squared p-value because these tables are
    sparse enough that the asymptotic null does not hold. The p-value counts
    permutations at least as extreme, with the +1 correction so it can never be
    exactly zero.
    """
    pair = pd.concat([a.rename("a"), b.rename("b")], axis=1).dropna()
    observed = cramers_v(pair["a"], pair["b"])
    if not np.isfinite(observed):
        return observed, float("nan"), float("nan")

    rng = np.random.default_rng(seed)
    values = pair["b"].to_numpy()
    null = np.empty(n_permutations)
    for i in range(n_permutations):
        null[i] = cramers_v(pair["a"], pd.Series(rng.permutation(values), index=pair.index))
    null = null[np.isfinite(null)]
    if null.size == 0:
        return observed, float("nan"), float("nan")
    p = (1.0 + float((null >= observed).sum())) / (1.0 + null.size)
    return observed, float(null.mean()), p


def confounding_table(
    annotations: pd.DataFrame,
    technical: tuple[str, ...],
    clinical: tuple[str, ...],
    *,
    seed: int,
    n_permutations: int = N_PERMUTATIONS,
) -> pd.DataFrame:
    """Every technical x clinical pair, with a permutation p-value.

    ``excess_over_null`` is observed V minus the mean V under permutation. It is
    the honest effect size: a 100-level plate factor has a large null V against
    anything, and the excess is what is left after subtracting that.
    """
    rows = []
    for tech in technical:
        for clin in clinical:
            if tech not in annotations or clin not in annotations:
                raise BatchError(f"missing column: {tech!r} or {clin!r}")
            observed, null_mean, p = permutation_p(
                annotations[tech],
                annotations[clin],
                n_permutations=n_permutations,
                seed=seed,
            )
            pair = annotations[[tech, clin]].dropna()
            rows.append(
                {
                    "technical_factor": tech,
                    "clinical_factor": clin,
                    "n": int(len(pair)),
                    "n_levels_technical": int(pair[tech].nunique()),
                    "n_levels_clinical": int(pair[clin].nunique()),
                    "cramers_v": round(observed, 4) if np.isfinite(observed) else None,
                    "null_mean_v": round(null_mean, 4) if np.isfinite(null_mean) else None,
                    "excess_over_null": round(observed - null_mean, 4)
                    if np.isfinite(observed) and np.isfinite(null_mean)
                    else None,
                    "permutation_p": round(p, 4) if np.isfinite(p) else None,
                }
            )
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Variance explained
# ---------------------------------------------------------------------------


def pool_rare_levels(values: pd.Series, *, min_size: int = MIN_LEVEL_SIZE) -> pd.Series:
    """Collapse levels with fewer than ``min_size`` members into ``"other"``.

    A level containing one sample is fitted perfectly and adds a full degree of
    freedom of spurious explanation. Pooling them is standard and the pooled
    count is always reported alongside.
    """
    counts = values.value_counts()
    rare = set(counts[counts < min_size].index)
    return values.where(~values.isin(rare), other="other")


def _one_way_r2(y: np.ndarray, groups: np.ndarray) -> tuple[float, float, int]:
    """R-squared and adjusted R-squared of a one-way ANOVA. Returns (r2, adj, k)."""
    codes, _ = pd.factorize(groups)
    k = int(codes.max()) + 1
    n = y.size
    if k < 2 or n <= k:
        return float("nan"), float("nan"), k

    grand = y.mean()
    ss_total = float(((y - grand) ** 2).sum())
    if ss_total <= 0:
        return float("nan"), float("nan"), k

    sums = np.bincount(codes, weights=y, minlength=k)
    sizes = np.bincount(codes, minlength=k).astype(float)
    means = np.divide(sums, sizes, out=np.zeros_like(sums), where=sizes > 0)
    ss_between = float((sizes * (means - grand) ** 2).sum())

    r2 = ss_between / ss_total
    adjusted = 1.0 - (1.0 - r2) * (n - 1) / (n - k)
    return float(r2), float(adjusted), k


def variance_explained(
    expression: pd.DataFrame,
    factors: pd.DataFrame,
    *,
    seed: int,
    variance_target: float = PC_VARIANCE_TARGET,
    n_permutations: int = 100,
    min_level_size: int = MIN_LEVEL_SIZE,
) -> tuple[pd.DataFrame, dict]:
    """PVCA-style: how much expression variance each factor is associated with.

    Principal components are computed once on the sample-by-gene matrix, then
    each retained PC is regressed on each factor and the R-squared values are
    averaged weighted by the PC's share of variance. Reported per factor:

    - ``r2`` weighted raw R-squared
    - ``adjusted_r2`` penalised for the number of levels
    - ``null_r2`` the same statistic under label permutation
    - ``excess_over_null`` ``r2 - null_r2`` — **read this one**

    Not a partition. Factors are correlated and the figures do not sum to 1.
    """
    if expression.empty:
        raise BatchError("expression matrix is empty")
    shared = expression.index.intersection(factors.index)
    if len(shared) < 10:
        raise BatchError(f"only {len(shared)} samples shared between expression and factors")
    matrix = expression.loc[shared].to_numpy(dtype=float)
    factors = factors.loc[shared]

    # Centre genes; drop the constant ones, which contribute nothing and make
    # the covariance singular.
    matrix = matrix[:, matrix.std(axis=0) > 0]
    matrix = matrix - matrix.mean(axis=0, keepdims=True)

    # Economy SVD: samples (n~600) is the small dimension, so this is cheap.
    # PC scores are U*S; computing the SVD once and reusing it avoids paying
    # for a second decomposition of a 600 x 60,000 matrix.
    u, singular, _ = np.linalg.svd(matrix, full_matrices=False)
    variances = singular**2
    ratios = variances / variances.sum()
    keep = int(np.searchsorted(np.cumsum(ratios), variance_target) + 1)
    keep = max(2, min(keep, len(ratios)))

    scores = u[:, :keep] * singular[:keep]
    weights = ratios[:keep] / ratios[:keep].sum()

    meta = {
        "n_samples": int(len(shared)),
        "n_genes": int(matrix.shape[1]),
        "n_pcs_retained": keep,
        "variance_covered": round(float(np.cumsum(ratios)[keep - 1]), 4),
        "min_level_size": min_level_size,
        "n_permutations": n_permutations,
    }

    rng = np.random.default_rng(seed)
    rows = []
    for name in factors.columns:
        raw = factors[name].dropna()
        if raw.nunique() < 2:
            continue
        pooled = pool_rare_levels(raw.astype(str), min_size=min_level_size)
        mask = factors.index.isin(pooled.index)
        groups = pooled.reindex(factors.index[mask]).to_numpy()
        sub_scores = scores[mask]

        r2 = adj = 0.0
        for j in range(keep):
            this_r2, this_adj, _ = _one_way_r2(sub_scores[:, j], groups)
            if np.isfinite(this_r2):
                r2 += weights[j] * this_r2
                adj += weights[j] * this_adj

        null = np.empty(n_permutations)
        for b in range(n_permutations):
            shuffled = rng.permutation(groups)
            acc = 0.0
            for j in range(keep):
                this_r2, _, _ = _one_way_r2(sub_scores[:, j], shuffled)
                if np.isfinite(this_r2):
                    acc += weights[j] * this_r2
            null[b] = acc

        rows.append(
            {
                "factor": name,
                "n": int(mask.sum()),
                "n_levels": int(pd.Series(groups).nunique()),
                "n_levels_before_pooling": int(raw.nunique()),
                "r2": round(r2, 4),
                "adjusted_r2": round(adj, 4),
                "null_r2": round(float(null.mean()), 4),
                "excess_over_null": round(r2 - float(null.mean()), 4),
                "permutation_p": round((1.0 + float((null >= r2).sum())) / (1.0 + null.size), 4),
            }
        )

    table = pd.DataFrame(rows).sort_values("excess_over_null", ascending=False)
    return table.reset_index(drop=True), meta
