"""W3.2 — the premise check must be able to find bimodality when it is there.

A test that only shows "continuous data looks continuous" proves nothing: a
broken detector agrees. So the load-bearing tests here are the positive
controls — synthetic two-group data that both tests must flag — and the
zero-inflation control, which is the way a real distribution most easily fakes
two groups.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.bulk.premise import (
    BIC_DECISIVE,
    DIP_ALPHA,
    MIN_N_FOR_TESTS,
    assess,
    purity_association,
    purity_conditioned_check,
    purity_tertiles,
    residualise_on_purity,
    run_premise_check,
    strata,
)

RNG = np.random.default_rng(20260817)
SEED = 20260817


def _unimodal(n=400, skew=False):
    x = RNG.lognormal(1.0, 0.6, n) if skew else RNG.normal(8.0, 1.0, n)
    return pd.Series(x)


def _bimodal(n=400, gap=6.0):
    half = n // 2
    return pd.Series(np.r_[RNG.normal(1.0, 0.7, half), RNG.normal(1.0 + gap, 0.7, half)])


def _full(rows, stratum="s"):
    """The whole-sample row (assess returns [full, nonzero])."""
    return next(r for r in rows if r.stratum == stratum)


# ---------------------------------------------------------------------------
# Positive controls — the detector must fire
# ---------------------------------------------------------------------------


def test_dip_detects_genuinely_two_groups():
    r = _full(assess(_bimodal(), gene="X", stratum="s", seed=SEED))
    assert r.dip_pvalue < DIP_ALPHA
    assert r.bic_delta > BIC_DECISIVE
    assert r.verdict == "bimodal"


def test_gmm_recovers_the_two_means():
    r = _full(assess(_bimodal(gap=6.0), gene="X", stratum="s", seed=SEED))
    lo, hi = eval(r.gmm_means)  # noqa: S307 — our own round-tripped list
    assert hi - lo == pytest.approx(6.0, abs=0.6)


def test_a_clean_unimodal_sample_is_called_continuous():
    r = _full(assess(_unimodal(), gene="X", stratum="s", seed=SEED))
    assert r.dip_pvalue > DIP_ALPHA
    assert r.verdict == "continuous"


# ---------------------------------------------------------------------------
# The two ways this analysis lies
# ---------------------------------------------------------------------------


def test_a_skewed_unimodal_sample_can_fool_bic_but_not_dip():
    """Why the brief says report both and do not pick whichever agrees.

    BIC compares 'two Gaussians' against 'one Gaussian'. For a skewed unimodal
    distribution two Gaussians genuinely do fit better — and that is not
    evidence of two groups. The dip test is the one that answers the question
    actually being asked.
    """
    r = _full(assess(_unimodal(n=600, skew=True), gene="X", stratum="s", seed=SEED))
    assert r.dip_pvalue > DIP_ALPHA
    assert r.verdict in {
        "continuous",
        "two_gaussians_fit_better_but_dip_says_unimodal",
    }
    if r.bic_delta > BIC_DECISIVE:
        assert r.verdict == "two_gaussians_fit_better_but_dip_says_unimodal"


def test_zero_inflation_is_flagged_not_reported_as_biology():
    """A spike at zero plus a hump is bimodal to any test, and means only that
    the gene was undetected in some samples."""
    values = pd.Series(np.r_[np.zeros(150), RNG.normal(8.0, 0.8, 350)])
    rows = assess(values, gene="X", stratum="s", seed=SEED)
    full, nonzero = rows[0], rows[1]

    assert full.zero_fraction == pytest.approx(0.30, abs=0.01)
    assert full.verdict == "bimodal_but_zero_inflated"
    # Drop the zeros and the second "mode" goes with them.
    assert nonzero.zero_fraction == 0.0
    assert nonzero.verdict == "continuous"


def test_the_nonzero_sensitivity_row_is_always_emitted():
    rows = assess(_unimodal(), gene="X", stratum="s", seed=SEED)
    assert [r.stratum for r in rows] == ["s", "s|nonzero"]


# ---------------------------------------------------------------------------
# Guards
# ---------------------------------------------------------------------------


def test_small_strata_report_insufficient_n_rather_than_a_number():
    """READ normal-adjacent is n=10. A dip p-value there would look like an
    answer and would not be one."""
    r = _full(assess(pd.Series(RNG.normal(8, 1, 10)), gene="X", stratum="s", seed=SEED))
    assert r.n < MIN_N_FOR_TESTS
    assert r.verdict == "insufficient_n"
    assert r.dip_pvalue is None and r.bic_delta is None


def test_constant_values_do_not_crash_the_mixture():
    r = _full(assess(pd.Series(np.full(50, 3.0)), gene="X", stratum="s", seed=SEED))
    assert r.dip_pvalue is None
    assert r.verdict == "continuous"


def test_results_are_deterministic_under_a_fixed_seed():
    values = _bimodal()
    a = _full(assess(values, gene="X", stratum="s", seed=SEED))
    b = _full(assess(values, gene="X", stratum="s", seed=SEED))
    assert a.bic_delta == b.bic_delta
    assert a.gmm_means == b.gmm_means


def test_every_row_records_that_purity_was_not_adjusted():
    """W3.3 has not landed. The caveat is a column, not a sentence in a note
    someone may not read."""
    rows = assess(_unimodal(), gene="X", stratum="s", seed=SEED)
    assert all(r.purity_adjusted is False for r in rows)


# ---------------------------------------------------------------------------
# Strata and the top-level runner
# ---------------------------------------------------------------------------


def _manifest():
    return pd.DataFrame(
        {
            "barcode": ["b1", "b2", "b3", "b4"],
            "sample_type": ["01", "11", "01", "11"],
            "project": ["TCGA-COAD", "TCGA-COAD", "TCGA-READ", "TCGA-READ"],
        }
    ).set_index("barcode")


def test_strata_split_tumour_normal_and_project():
    idx = pd.Index(["b1", "b2", "b3", "b4"])
    masks = strata(_manifest(), idx)
    assert masks["COAD_tumour"].tolist() == [True, False, False, False]
    assert masks["READ_normal"].tolist() == [False, False, False, True]
    assert masks["COAD+READ_tumour"].sum() == 2


def test_run_premise_check_covers_every_gene_and_stratum():
    idx = pd.Index(["b1", "b2", "b3", "b4"])
    expr = pd.DataFrame({"ENSG1": [1.0, 8.0, 2.0, 9.0]}, index=idx)
    out = run_premise_check(expr, _manifest(), {"GUCA2A": "ENSG1"}, seed=SEED)
    assert set(out["gene"]) == {"GUCA2A"}
    assert len(out) == 10  # 5 strata x (full + nonzero)
    assert (out["verdict"] == "insufficient_n").all()  # every stratum is tiny


def test_run_premise_check_fails_loudly_on_a_missing_gene():
    idx = pd.Index(["b1", "b2", "b3", "b4"])
    expr = pd.DataFrame({"ENSG1": [1.0, 8.0, 2.0, 9.0]}, index=idx)
    with pytest.raises(KeyError, match="not in the expression matrix"):
        run_premise_check(expr, _manifest(), {"GUCA2A": "ENSG_absent"}, seed=SEED)


# ---------------------------------------------------------------------------
# The purity-conditioned re-run (W3.2 second pass, after W3.3)
# ---------------------------------------------------------------------------


def _purity_manifest(barcodes):
    return pd.DataFrame(
        {
            "barcode": barcodes,
            "sample_type": ["01"] * len(barcodes),
            "project": ["TCGA-COAD"] * len(barcodes),
        }
    ).set_index("barcode")


def test_residualising_removes_the_purity_relationship():
    purity = pd.Series(RNG.uniform(0.2, 0.95, 300))
    values = 3.0 + 5.0 * purity + RNG.normal(0, 0.4, 300)
    resid = residualise_on_purity(values, purity)
    assert abs(float(resid.corr(purity))) < 0.05
    # Recentred, not zero-centred, so the axis stays readable.
    assert float(resid.mean()) == pytest.approx(float(values.mean()), abs=1e-9)


def test_samples_without_a_purity_call_are_dropped_not_imputed():
    """An imputed covariate turns a coverage gap into an invisible assumption."""
    purity = pd.Series([0.5, 0.6, np.nan, 0.8, 0.9, np.nan])
    values = pd.Series([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    resid = residualise_on_purity(values, purity)
    assert len(resid) == 4
    assert list(resid.index) == [0, 1, 3, 4]


def test_purity_driven_bimodality_is_removed_by_conditioning():
    """THE positive control for this pass.

    The brief's worry was that apparent bimodality could be a purity artifact.
    Build exactly that — two purity regimes, one underlying biology — confirm
    both tests call it bimodal, then confirm residualising collapses it. If this
    fails, the re-run cannot detect the artifact it exists to rule out.
    """
    half = 200
    purity = pd.Series(np.r_[RNG.normal(0.30, 0.04, half), RNG.normal(0.85, 0.04, half)])
    values = 1.0 + 9.0 * purity + RNG.normal(0, 0.25, 2 * half)

    before = _full(assess(values, gene="X", stratum="s", seed=SEED))
    assert before.dip_pvalue < DIP_ALPHA, "fixture no longer reproduces the artifact"

    after = _full(
        assess(residualise_on_purity(values, purity), gene="X", stratum="s", seed=SEED)
    )
    assert after.dip_pvalue > DIP_ALPHA
    assert after.verdict in {"continuous", "two_gaussians_fit_better_but_dip_says_unimodal"}


def test_genuine_bimodality_survives_conditioning():
    """The converse, and just as necessary: conditioning must not erase real
    structure that has nothing to do with purity."""
    purity = pd.Series(RNG.uniform(0.3, 0.9, 400))
    values = _bimodal(400).to_numpy() + 0.5 * purity
    values = pd.Series(values)
    after = _full(
        assess(residualise_on_purity(values, purity), gene="X", stratum="s", seed=SEED)
    )
    assert after.dip_pvalue < DIP_ALPHA
    assert after.verdict == "bimodal"


def test_purity_tertiles_split_roughly_evenly():
    purity = pd.Series(RNG.uniform(0, 1, 300))
    masks = purity_tertiles(purity)
    assert set(masks) == {"purity_low", "purity_mid", "purity_high"}
    sizes = [int(m.sum()) for m in masks.values()]
    assert sum(sizes) == 300
    assert max(sizes) - min(sizes) <= 2


def test_purity_tertiles_refuse_a_tiny_sample():
    assert purity_tertiles(pd.Series(RNG.uniform(0, 1, 10))) == {}


def test_conditioned_rows_are_marked_purity_adjusted():
    """The base Assessment defaults purity_adjusted=False. These rows are the
    exception and must say so, or the two passes are indistinguishable in the
    parquet."""
    barcodes = [f"TCGA-AA-{i:04d}-01A-01R-1410-07" for i in range(120)]
    expr = pd.DataFrame({"ENSG1": RNG.normal(8, 1, 120)}, index=barcodes)
    purity = pd.Series(RNG.uniform(0.3, 0.9, 120), index=barcodes)
    out = purity_conditioned_check(
        expr, _purity_manifest(barcodes), {"GUCA2A": "ENSG1"}, purity,
        method="absolute", seed=SEED,
    )
    assert out["purity_adjusted"].all()
    assert (out["purity_method"] == "absolute").all()
    assert "tumour|purity_residual" in set(out["stratum"])


def test_purity_association_reports_variance_explained():
    barcodes = [f"TCGA-AA-{i:04d}-01A-01R-1410-07" for i in range(200)]
    purity = pd.Series(RNG.uniform(0.2, 0.95, 200), index=barcodes)
    expr = pd.DataFrame({"ENSG1": 2.0 + 4.0 * purity.to_numpy()}, index=barcodes)
    out = purity_association(
        expr, _purity_manifest(barcodes), {"GUCA2A": "ENSG1"}, purity, method="absolute"
    )
    assert out.loc[0, "r_squared"] == pytest.approx(1.0, abs=1e-6)
    assert out.loc[0, "n"] == 200
