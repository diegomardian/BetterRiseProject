"""The Wnt test, against data whose Wnt-to-target relationship we set.

The pipeline's statistical core is `wnt_score` -> `partial_spearman` ->
`summarise`. These drive it with a per-cell latent Wnt level and a known
relationship to the target gene, and check three things in the order that
matters:

1. it FINDS an injected suppression,
2. it does NOT find one that was never injected,
3. it does not find one that is really RESIDUAL MATURITY — the confound the
   partial correlation exists to remove, and the most likely false positive in
   the design.

Test 3 is the load-bearing one. A version of this analysis without the maturity
conditioner passes 1 and 2 and fails 3, which is exactly how it would have
reported the labeller's own gradient as a mechanism.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.reference.jobs.wnt_mechanism import summarise
from src.reference.wnt_score import (
    SIGNATURE,
    WntScoreError,
    assert_no_signature_leakage,
    partial_spearman,
    wnt_score,
    wnt_stem_verdict,
)

N_CELLS = 600


def _cells(*, wnt_effect: float, maturity_effect: float, seed: int = 3):
    """A synthetic mature population with a stated causal structure.

    ``wnt_effect``      how strongly the latent Wnt level suppresses the target
    ``maturity_effect`` how strongly maturity RAISES it, with maturity also
                        suppressing Wnt — the confound
    """
    rng = np.random.default_rng(seed)
    maturity = rng.normal(0, 1, N_CELLS)
    wnt_latent = -maturity_effect * maturity + rng.normal(0, 1, N_CELLS)
    depth = rng.lognormal(8.6, 0.35, N_CELLS)

    genes = list(SIGNATURE) + ["GUCA2A", "ACTB"]
    rates = {}
    for g in SIGNATURE:
        rates[g] = 1.5 * np.exp(0.8 * wnt_latent)
    rates["GUCA2A"] = 8.0 * np.exp(
        -wnt_effect * wnt_latent + maturity_effect * maturity)
    rates["ACTB"] = np.full(N_CELLS, 25.0)

    counts = np.column_stack([
        rng.poisson(rates[g] * depth / 1e4 * 1e4 / 1e4 * (depth / depth.mean()))
        for g in genes
    ]).astype(float)
    return counts, genes, depth, maturity


def _rho(counts, genes, depth, maturity, gene, *, condition: bool) -> float:
    score, _ = wnt_score(counts, genes, depth=depth)
    values = counts[:, genes.index(gene)] / depth * 1e4
    conditioners = (np.column_stack([maturity, np.log(depth)]) if condition
                    else np.zeros((len(depth), 1)))
    return partial_spearman(score, values, conditioners)


# ---------------------------------------------------------------------------
# Does it find what is there, and not what is not
# ---------------------------------------------------------------------------


def test_an_injected_wnt_suppression_is_recovered():
    counts, genes, depth, maturity = _cells(wnt_effect=0.6, maturity_effect=0.0)
    rho = _rho(counts, genes, depth, maturity, "GUCA2A", condition=True)
    assert rho < -0.15, f"injected suppression not recovered: rho={rho:+.3f}"


def test_no_association_is_found_where_none_was_injected():
    counts, genes, depth, maturity = _cells(wnt_effect=0.0, maturity_effect=0.0)
    rho = _rho(counts, genes, depth, maturity, "GUCA2A", condition=True)
    assert abs(rho) < 0.15, f"fired on nothing: rho={rho:+.3f}"


def test_the_housekeeping_control_stays_flat_under_a_real_effect():
    """ACTB has no Wnt dependence by construction, so it is the floor."""
    counts, genes, depth, maturity = _cells(wnt_effect=0.6, maturity_effect=0.0)
    assert abs(_rho(counts, genes, depth, maturity, "ACTB", condition=True)) < 0.15


# ---------------------------------------------------------------------------
# THE ONE THAT MATTERS: maturity, not Wnt
# ---------------------------------------------------------------------------


def test_residual_maturity_is_not_reported_as_wnt():
    """THE INPUT THAT FORCES THE CONDITIONER TO EARN ITS PLACE.

    Here Wnt does NOT touch the target at all. Maturity raises the target and
    suppresses Wnt, so an unconditioned correlation sees a strong negative
    association that is entirely the labeller's own gradient. Conditioning on
    maturity must remove it.

    An analysis without this conditioner passes every other test in this file
    and reports this fixture as a Wnt mechanism.
    """
    counts, genes, depth, maturity = _cells(wnt_effect=0.0, maturity_effect=1.2)

    naive = _rho(counts, genes, depth, maturity, "GUCA2A", condition=False)
    conditioned = _rho(counts, genes, depth, maturity, "GUCA2A", condition=True)

    assert naive < -0.25, (
        "the fixture must actually contain the confound, or this test is not "
        f"exercising anything: naive rho={naive:+.3f}"
    )
    assert abs(conditioned) < 0.15, (
        f"conditioning failed to remove residual maturity: {naive:+.3f} -> "
        f"{conditioned:+.3f}. Without this the labeller's own gradient is "
        f"reported as a Wnt mechanism."
    )


# ---------------------------------------------------------------------------
# The guards
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("gene", ["ASCL2", "LGR5", "CTNNB1", "TCF7L2"])
def test_invariant_8s_forbidden_genes_are_refused(gene):
    """Each is named in invariant 8 and each would break the test differently."""
    with pytest.raises(WntScoreError, match="invariant 8"):
        assert_no_signature_leakage([], [], signature=(*SIGNATURE, gene))


def test_a_signature_containing_a_label_marker_is_refused():
    """Invariant 2's failure mode with predictor and label swapped."""
    with pytest.raises(WntScoreError, match="labelling-axis markers"):
        assert_no_signature_leakage([], ["OLFM4"], signature=(*SIGNATURE, "OLFM4"))


def test_a_signature_containing_a_scored_gene_is_refused():
    with pytest.raises(WntScoreError, match="frozen panel"):
        assert_no_signature_leakage(["GUCA2A"], [], signature=(*SIGNATURE, "GUCA2A"))


def test_the_real_panel_and_axes_pass_the_guard():
    """The other half — it must not refuse everything."""
    from src.common.panel import load_axes, panel_genes

    markers = [g for a in load_axes()["axes"].values() for g in (a.get("genes") or [])]
    assert_no_signature_leakage(panel_genes(), markers)


def test_a_signature_of_one_gene_is_refused():
    rng = np.random.default_rng(0)
    counts = rng.poisson(3.0, (50, 2)).astype(float)
    with pytest.raises(WntScoreError, match="signature of one gene"):
        wnt_score(counts, ["AXIN2", "ACTB"], depth=np.full(50, 5000.0))


def test_a_wnt_score_that_is_really_maturity_withholds_the_reading():
    assert wnt_stem_verdict(0.85)["verdict"] == "WNT SCORE TRACKS MATURITY"
    assert "WITHHELD" in wnt_stem_verdict(0.85)["detail"]
    assert wnt_stem_verdict(0.20)["verdict"] == "SEPARABLE"
    assert wnt_stem_verdict(float("nan"))["verdict"] == "UNDEFINED"


def test_too_few_cells_returns_undefined_rather_than_a_number():
    rng = np.random.default_rng(0)
    x, y = rng.normal(size=10), rng.normal(size=10)
    assert np.isnan(partial_spearman(x, y, np.zeros((10, 1))))


def test_the_summary_uses_patients_as_the_unit():
    rows = pd.DataFrame({
        "arm": "tumour", "gene": "GUCA2A", "role": "target",
        "patient_id": [f"P{i}" for i in range(12)],
        "partial_rho": np.linspace(-0.4, -0.2, 12),
        "unconditioned_rho": np.linspace(-0.6, -0.4, 12),
    })
    out = summarise(rows, seed=1)
    assert set(out["statistic"]) == {"partial_rho", "unconditioned_rho"}
    assert (out["n_patients"] == 12).all()
    assert out.loc[out.statistic == "partial_rho", "excludes_zero"].all()
