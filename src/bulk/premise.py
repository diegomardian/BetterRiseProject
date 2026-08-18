"""W3.2 — is bulk GUCA2A loss actually bimodal? The premise check.

THE POINT
---------
The project's classification framing rests on an assumption nobody has tested:
that bulk GUCA2A loss separates tumours into two groups. Executive-Brief error
#6, "assuming bulk GUCA2A is negligible" (execution_plan.md §2.1). If the
distribution is **continuous rather than bimodal, the two-type classification
dissolves into a regression** and the project changes shape — in week 2, not
week 12.

This module is deliberately unglamorous: two tests, run the same way on every
gene and every stratum, reported side by side whether or not they agree.

TWO TESTS, BOTH REPORTED
------------------------
1. **Hartigan's dip test** — non-parametric. Null is unimodality. A small
   p-value says "not unimodal" without claiming the alternative is two Gaussians.
2. **1- vs 2-component Gaussian mixture, compared by BIC.** Parametric, and it
   *will* prefer two components for a skewed unimodal distribution often enough
   to matter.

The brief says report both and do not pick whichever agrees. They answer
different questions and disagreement is informative: dip says "lumpy", BIC says
"two Gaussians fit better than one". Only the first is evidence of two groups.

THREE THINGS THAT FAKE BIMODALITY, AND WHAT IS DONE ABOUT THEM
--------------------------------------------------------------
**Zero inflation.** log2(CPM+1) maps every undetected sample to exactly 0.0.
A spike at zero plus a broad hump is bimodal to any test, and means only that
the gene is undetected in some samples. :func:`assess` reports ``zero_fraction``
and re-runs both tests on the non-zero samples as a sensitivity. If bimodality
survives only with the zeros in, it is a detection floor, not biology.

**Purity.** Bulk expression is diluted by stromal and immune content. Two
apparent groups may be two purity regimes. **Nothing here is adjusted for
purity, so every result is provisional until W3.3 lands** — this is flagged in
the output as a column, not left to the write-up.

**Sample size.** READ normal-adjacent is n=10. The dip test on n=10 is not
meaningful; :data:`MIN_N_FOR_TESTS` gates it and the result is reported as
``insufficient_n`` rather than as a number that looks like an answer.

INVARIANT 2
-----------
This module reads panel genes as **outcome variables** — their distributions are
the thing being measured. That is allowed and is not what invariant 2 forbids.
Nothing here uses a panel gene to define a group, a label, or a reference.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd

#: Below this, neither test says anything trustworthy. READ normal-adjacent is
#: n=10, which is why this gate exists rather than being assumed away.
MIN_N_FOR_TESTS = 20

#: BIC difference below which the 1- and 2-component fits are indistinguishable.
#: The conventional "positive evidence" threshold; below it, preferring two
#: components is noise-fitting.
BIC_DECISIVE = 10.0

#: Dip-test significance. Not corrected for multiplicity — this is a premise
#: check across a handful of genes, and the verdict is read off effect size and
#: the figure, not off a threshold.
DIP_ALPHA = 0.05


@dataclass
class Assessment:
    """One gene, one stratum, both tests, plus everything needed to distrust it."""

    gene: str
    stratum: str
    n: int
    zero_fraction: float
    median: float
    iqr: float
    dip_statistic: float | None
    dip_pvalue: float | None
    bic_1: float | None
    bic_2: float | None
    bic_delta: float | None  # bic_1 - bic_2; positive favours two components
    gmm_weights: str | None
    gmm_means: str | None
    verdict: str
    purity_adjusted: bool = False  # W3.3 has not landed; see module docstring

    def to_row(self) -> dict:
        return asdict(self)


def _dip(values: np.ndarray) -> tuple[float | None, float | None]:
    """Hartigan's dip test. Null is unimodality."""
    import diptest

    stat, pval = diptest.diptest(values)
    return float(stat), float(pval)


def _gmm_bic(values: np.ndarray, *, seed: int) -> tuple[float, float, list, list]:
    """BIC for 1- and 2-component Gaussian mixtures on 1-D data.

    Two non-default choices, both to keep this off a BLAS that is not reliably
    callable on every developer machine:

    ``covariance_type="diag"`` rather than ``"full"``. For univariate data the
    two are mathematically identical — there is one variance either way — but
    ``"full"`` routes through a Cholesky factorisation in LAPACK.

    ``init_params="random_from_data"`` rather than the default ``"kmeans"``.
    sklearn's KMeans calls ``threadpoolctl`` to introspect the loaded BLAS, and
    on a machine where that introspection fails the process dies with a native
    error rather than an exception. Random init costs nothing here: the data is
    one-dimensional and ``n_init=10`` covers the initialisation sensitivity that
    kmeans init exists to avoid.
    """
    from sklearn.mixture import GaussianMixture

    x = values.reshape(-1, 1)
    fits = {}
    for k in (1, 2):
        gm = GaussianMixture(
            n_components=k,
            covariance_type="diag",
            init_params="random_from_data",
            random_state=seed,
            n_init=10,
        ).fit(x)
        fits[k] = gm
    two = fits[2]
    order = np.argsort(two.means_.ravel())
    return (
        float(fits[1].bic(x)),
        float(two.bic(x)),
        [round(float(w), 4) for w in two.weights_[order]],
        [round(float(m), 4) for m in two.means_.ravel()[order]],
    )


def _verdict(
    n: int, dip_p: float | None, bic_delta: float | None, zero_fraction: float
) -> str:
    """A short machine-readable call. The prose verdict is the deliverable."""
    if n < MIN_N_FOR_TESTS:
        return "insufficient_n"
    dip_says_lumpy = dip_p is not None and dip_p < DIP_ALPHA
    bic_says_two = bic_delta is not None and bic_delta > BIC_DECISIVE
    if zero_fraction > 0.10 and (dip_says_lumpy or bic_says_two):
        return "bimodal_but_zero_inflated"
    if dip_says_lumpy and bic_says_two:
        return "bimodal"
    if bic_says_two and not dip_says_lumpy:
        return "two_gaussians_fit_better_but_dip_says_unimodal"
    if dip_says_lumpy and not bic_says_two:
        return "dip_says_lumpy_but_one_gaussian_suffices"
    return "continuous"


def assess(values: pd.Series, *, gene: str, stratum: str, seed: int) -> list[Assessment]:
    """Both tests on one gene x stratum, plus the drop-the-zeros sensitivity.

    Returns two rows: the full sample, and the non-zero subset. Comparing them is
    how a detection floor is told apart from two biological groups.
    """
    out: list[Assessment] = []
    for label, x in (
        (stratum, values.to_numpy(dtype=float)),
        (f"{stratum}|nonzero", values.loc[values > 0].to_numpy(dtype=float)),
    ):
        x = x[np.isfinite(x)]
        n = int(x.size)
        zero_fraction = float((x == 0).mean()) if n else float("nan")

        dip_stat = dip_p = bic1 = bic2 = delta = None
        weights = means = None
        if n >= MIN_N_FOR_TESTS and np.ptp(x) > 0:
            dip_stat, dip_p = _dip(x)
            bic1, bic2, w, m = _gmm_bic(x, seed=seed)
            delta = bic1 - bic2
            weights, means = str(w), str(m)

        out.append(
            Assessment(
                gene=gene,
                stratum=label,
                n=n,
                zero_fraction=zero_fraction,
                median=float(np.median(x)) if n else float("nan"),
                iqr=float(np.subtract(*np.percentile(x, [75, 25]))) if n else float("nan"),
                dip_statistic=dip_stat,
                dip_pvalue=dip_p,
                bic_1=bic1,
                bic_2=bic2,
                bic_delta=delta,
                gmm_weights=weights,
                gmm_means=means,
                verdict=_verdict(n, dip_p, delta, zero_fraction),
            )
        )
    return out


def strata(manifest: pd.DataFrame, samples: pd.Index) -> dict[str, np.ndarray]:
    """Boolean masks for the strata the brief asks for.

    COAD tumours are the primary analysis; the COAD/READ split is the stratified
    repeat. Normal-adjacent is included because the tumour distribution means
    little without knowing what the gene looks like when the tissue is intact —
    but note n=41 and n=10, so those two are description, not inference.
    """
    m = manifest.loc[samples]
    tumour = (m["sample_type"] == "01").to_numpy()
    coad = (m["project"] == "TCGA-COAD").to_numpy()
    read = (m["project"] == "TCGA-READ").to_numpy()
    return {
        "COAD_tumour": tumour & coad,
        "READ_tumour": tumour & read,
        "COAD+READ_tumour": tumour,
        "COAD_normal": ~tumour & coad,
        "READ_normal": ~tumour & read,
    }


def run_premise_check(
    expression: pd.DataFrame,
    manifest: pd.DataFrame,
    gene_ids: dict[str, str],
    *,
    seed: int,
) -> pd.DataFrame:
    """Every gene x stratum, both tests, one tidy frame."""
    masks = strata(manifest, expression.index)
    rows: list[dict] = []
    for symbol, ensembl_id in gene_ids.items():
        if ensembl_id not in expression.columns:
            raise KeyError(f"{symbol} ({ensembl_id}) is not in the expression matrix")
        series = expression[ensembl_id]
        for stratum, mask in masks.items():
            for a in assess(series.loc[mask], gene=symbol, stratum=stratum, seed=seed):
                rows.append(a.to_row())
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# The figure
# ---------------------------------------------------------------------------


def plot_distributions(
    expression: pd.DataFrame,
    manifest: pd.DataFrame,
    gene_ids: dict[str, str],
    out_path: str | Path,
) -> Path:
    """Tumour vs normal-adjacent, log scale, one row per gene.

    The figure is the deliverable the team reads, so it carries the caveats:
    n in the legend, and the zero fraction annotated, because a spike at the
    left edge is the single most misleading feature of these distributions.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    masks = strata(manifest, expression.index)
    genes = list(gene_ids)
    fig, axes = plt.subplots(len(genes), 2, figsize=(11, 3.2 * len(genes)), squeeze=False)

    for row, symbol in enumerate(genes):
        values = expression[gene_ids[symbol]]
        for col, (tum_key, norm_key, title) in enumerate(
            [
                ("COAD_tumour", "COAD_normal", "COAD"),
                ("READ_tumour", "READ_normal", "READ"),
            ]
        ):
            ax = axes[row][col]
            tum = values.loc[masks[tum_key]]
            norm = values.loc[masks[norm_key]]
            bins = np.linspace(0, max(float(values.max()), 1.0), 45)
            ax.hist(tum, bins=bins, alpha=0.65, label=f"tumour (n={len(tum)})", density=True)
            ax.hist(
                norm, bins=bins, alpha=0.55, label=f"normal-adj (n={len(norm)})", density=True
            )
            zf = float((tum == 0).mean())
            ax.set_title(f"{symbol} — {title}", fontsize=11)
            ax.set_xlabel("log2(CPM + 1)")
            ax.set_ylabel("density")
            ax.legend(fontsize=8, frameon=False)
            ax.annotate(
                f"tumour zeros: {zf:.1%}",
                xy=(0.98, 0.72),
                xycoords="axes fraction",
                ha="right",
                fontsize=8,
                color="#666666",
            )

    fig.suptitle(
        "W3.2 premise check — PROVISIONAL, not adjusted for tumour purity (W3.3)",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


# ---------------------------------------------------------------------------
# W3.2 re-run, conditioned on purity — the scheduled follow-up
# ---------------------------------------------------------------------------
#
# The first pass was marked provisional because bulk expression is diluted by
# stromal and immune content, and two apparent groups could be two purity
# regimes. The finding was "no groups", so the question inverts: could purity
# be *masking* structure rather than creating it? Adjusting for it is how you
# find out, and it is cheap now that W3.3 exists.
#
# ABSOLUTE is the primary covariate — it is called from copy number, not
# expression, so it is not circular with the outcome. ESTIMATE is the
# sensitivity analysis, and it covers the 18% of samples ABSOLUTE misses.


def residualise_on_purity(values: pd.Series, purity: pd.Series) -> pd.Series:
    """Residuals of ``values`` regressed on ``purity``, recentred on the mean.

    Samples without a purity call are dropped rather than imputed — imputing a
    covariate to preserve n is how a coverage gap becomes an invisible
    assumption.

    Recentring changes nothing statistically (the dip statistic is invariant to
    location and scale) and keeps the residuals on a readable expression-like
    axis, so the same histogram code works on both.
    """
    joined = pd.concat([values.rename("y"), purity.rename("p")], axis=1).dropna()
    if len(joined) < 3 or joined["p"].nunique() < 2:
        return pd.Series(dtype=float)
    slope, intercept = np.polyfit(joined["p"], joined["y"], 1)
    resid = joined["y"] - (intercept + slope * joined["p"])
    return resid + joined["y"].mean()


def purity_tertiles(purity: pd.Series) -> dict[str, np.ndarray]:
    """Low/mid/high purity masks, as a non-parametric alternative to regressing.

    Residualising assumes the expression/purity relationship is linear.
    Stratifying assumes nothing about its shape. If the two disagree, the
    linearity assumption is doing work and should be looked at.
    """
    clean = purity.dropna()
    if len(clean) < 30:
        return {}
    lo, hi = clean.quantile([1 / 3, 2 / 3])
    return {
        "purity_low": (purity <= lo).to_numpy(),
        "purity_mid": ((purity > lo) & (purity <= hi)).to_numpy(),
        "purity_high": (purity > hi).to_numpy(),
    }


def purity_conditioned_check(
    expression: pd.DataFrame,
    manifest: pd.DataFrame,
    gene_ids: dict[str, str],
    purity: pd.Series,
    *,
    method: str,
    seed: int,
) -> pd.DataFrame:
    """Re-run both tests on purity-residualised expression, and within tertiles.

    ``purity`` is indexed by barcode. ``method`` is recorded on every row so an
    ABSOLUTE-based result is never mistaken for an ESTIMATE-based one.
    """
    tumour = manifest.reindex(expression.index)["sample_type"] == "01"
    rows: list[dict] = []

    for symbol, ensembl_id in gene_ids.items():
        if ensembl_id not in expression.columns:
            raise KeyError(f"{symbol} ({ensembl_id}) is not in the expression matrix")
        values = expression.loc[tumour.to_numpy(), ensembl_id]
        aligned = purity.reindex(values.index)

        resid = residualise_on_purity(values, aligned)
        for a in assess(resid, gene=symbol, stratum="tumour|purity_residual", seed=seed):
            rows.append({**a.to_row(), "purity_method": method})

        sub_purity = aligned.loc[values.index]
        for name, mask in purity_tertiles(sub_purity).items():
            for a in assess(
                values.loc[mask], gene=symbol, stratum=f"tumour|{name}", seed=seed
            ):
                rows.append({**a.to_row(), "purity_method": method})

    out = pd.DataFrame(rows)
    # These rows ARE purity-adjusted; the base Assessment default says otherwise.
    out["purity_adjusted"] = True
    return out


def purity_association(
    expression: pd.DataFrame,
    manifest: pd.DataFrame,
    gene_ids: dict[str, str],
    purity: pd.Series,
    *,
    method: str,
) -> pd.DataFrame:
    """How much of each gene's variance purity explains. The number that says
    whether conditioning on it was ever going to matter."""
    tumour = manifest.reindex(expression.index)["sample_type"] == "01"
    rows = []
    for symbol, ensembl_id in gene_ids.items():
        values = expression.loc[tumour.to_numpy(), ensembl_id]
        joined = pd.concat(
            [values.rename("y"), purity.reindex(values.index).rename("p")], axis=1
        ).dropna()
        r = float(joined["y"].corr(joined["p"])) if len(joined) >= 3 else float("nan")
        rows.append(
            {
                "gene": symbol,
                "purity_method": method,
                "n": int(len(joined)),
                "pearson_r": round(r, 4),
                "r_squared": round(r * r, 4),
                "spearman_rho": round(float(joined["y"].corr(joined["p"], method="spearman")), 4)
                if len(joined) >= 3
                else float("nan"),
            }
        )
    return pd.DataFrame(rows)
