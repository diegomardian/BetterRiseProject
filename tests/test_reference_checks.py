"""G1's pre-committed thresholds. W1.

The tests that matter here are the ones asserting an unestimable tier does
**not** pass the gate. Decision #17 fixed the thresholds before any G1 number
existed; the failure mode this file guards is a tier that cannot be measured
quietly clearing the gate by absence.
"""

from __future__ import annotations

import pandas as pd
import pytest

from src.common.panel import panel_genes
from src.reference.checks import (
    G1_TIERS,
    MIN_TIER_SEPARATION,
    TIER_D_MAX_RHO,
    G1Error,
    compare_statistics,
    g1_verdict,
    loss_correlation,
    retention_correlation,
    tier_correlations,
)


def _corr(rho_a: float, rho_b: float, rho_d: float) -> pd.DataFrame:
    """A correlations frame with every tier estimated."""
    return pd.DataFrame({
        "tier": list(G1_TIERS),
        "n_genes": [8, 8, 8],
        "rho": [rho_a, rho_b, rho_d],
        "p_value": [0.01, 0.01, 0.01],
        "estimability": ["estimated"] * 3,
        "reason": ["", "", ""],
    })


def _panel_frame() -> pd.DataFrame:
    genes = panel_genes(("A", "B", "D"))
    return pd.DataFrame({
        "gene": genes,
        "abundance": [float(i + 1) for i in range(len(genes))],
        "loss": [float(i) * 0.1 for i in range(len(genes))],
        "retention": [1.0 - i * 0.05 for i in range(len(genes))],
    })


class TestTierCorrelations:
    def test_always_reports_every_tier(self):
        out = tier_correlations(_panel_frame(), value_column="loss")
        assert list(out["tier"]) == list(G1_TIERS)

    def test_thin_tier_is_not_estimable_not_zero(self):
        out = tier_correlations(_panel_frame(), value_column="loss")
        thin = out[out["n_genes"] < 4]
        assert len(thin), "the frozen panel has tiers below the Spearman floor"
        assert (thin["estimability"] == "not_estimable").all()
        assert not (thin["rho"] == 0.0).any()

    def test_missing_column_raises(self):
        with pytest.raises(G1Error, match="missing column"):
            tier_correlations(_panel_frame(), value_column="nope")

    def test_empty_raises(self):
        with pytest.raises(G1Error, match="empty"):
            tier_correlations(pd.DataFrame(columns=["gene", "abundance", "loss"]),
                              value_column="loss")

    def test_records_which_statistic_it_ran(self):
        assert (loss_correlation(_panel_frame())["statistic"] == "loss").all()
        out = retention_correlation(_panel_frame())
        assert (out["statistic"] == "retention").all()


class TestVerdict:
    def test_unestimable_tier_does_not_pass(self):
        """The one that matters. Absence is not a pass."""
        corr = _corr(0.9, 0.8, 0.0)
        corr.loc[corr["tier"] == "D", "estimability"] = "not_estimable"
        corr.loc[corr["tier"] == "D", "rho"] = None
        assert g1_verdict(corr)["verdict"] == "not_estimable"

    def test_frozen_panel_cannot_currently_evaluate_g1(self):
        """Tiers B and D are below the floor, so the gate is not evaluable.

        This is a statement about the frozen panel, not about this code:
        decision #17 assumed n≈8 per tier and tier D holds one gene.
        """
        verdict = g1_verdict(loss_correlation(_panel_frame()))
        assert verdict["verdict"] == "not_estimable"

    def test_fails_when_tier_d_tracks_abundance(self):
        verdict = g1_verdict(_corr(0.1, 0.1, TIER_D_MAX_RHO + 0.2))
        assert verdict["verdict"] == "fail"
        assert any("tier D" in r for r in verdict["reasons"])

    def test_fails_when_all_tiers_alike(self):
        verdict = g1_verdict(_corr(0.30, 0.35, 0.32))
        assert verdict["verdict"] == "fail"
        assert any("within" in r for r in verdict["reasons"])

    def test_passes_when_d_flat_and_both_separate(self):
        verdict = g1_verdict(_corr(0.8, 0.75, 0.05))
        assert verdict["verdict"] == "pass"

    def test_gap_in_the_rule_is_surfaced_not_guessed(self):
        """D flat, A separated, B not — #17 covers neither outcome."""
        verdict = g1_verdict(_corr(0.8, 0.10, 0.05))
        assert verdict["verdict"] == "indeterminate"
        assert any("does not cover" in r for r in verdict["reasons"])

    def test_separation_boundary_is_strict(self):
        exact = 0.05 + MIN_TIER_SEPARATION
        assert g1_verdict(_corr(exact, exact, 0.05))["verdict"] != "pass"

    def test_missing_column_raises(self):
        with pytest.raises(G1Error, match="missing column"):
            g1_verdict(pd.DataFrame({"tier": list(G1_TIERS)}))


class TestCompareStatistics:
    def test_gate_follows_the_named_statistic(self):
        named = _corr(0.30, 0.35, 0.32)          # fail
        secondary = _corr(0.8, 0.75, 0.05)       # pass
        out = compare_statistics(named, secondary)
        assert out["gate_verdict"] == "fail"
        assert out["agree"] is False
        assert "disagree" in out["note"]

    def test_agreement_is_reported(self):
        both = _corr(0.8, 0.75, 0.05)
        out = compare_statistics(both, both.copy())
        assert out["agree"] is True
        assert out["gate_verdict"] == "pass"
