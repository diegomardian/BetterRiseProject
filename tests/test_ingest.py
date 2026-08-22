"""Per-study, per-compartment QC flagging. execution_plan.md W1 task table:
"per-study thresholds, not global" -- the same shape W4 mirrors for the Lee
cohorts. Per compartment as well, since docs/open_decisions.md #12: a MAD rule
assumes one unimodal population and a cohort matrix is six of them.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.estimator.ingest import (
    correct_ambient,
    differential_retention,
    flag_doublets,
    qc_flags,
)


def _metrics(seed=0):
    """Two studies with different library-size scales -- a global cutoff
    would either over-filter the low-depth study or under-filter the
    high-depth one. A few planted outliers in each.

    Each study carries two compartments with a realistic depth separation:
    epithelium runs ~3.9x deeper than immune on SMC (measured, module
    docstring), which is exactly what makes a pooled MAD misfire.
    """
    rng = np.random.default_rng(seed)

    def study(study_id, immune_genes, immune_counts, n_immune=800, n_epithelium=200):
        # Immune outnumbers epithelium ~4:1, as on SMC (45k non-epithelial vs
        # 18.5k epithelial). That is what puts the pooled median inside the
        # immune mode and makes the pooled MAD small enough to misfire.
        immune = pd.DataFrame(
            {
                "study_id": study_id,
                "compartment": "T cells",
                "n_genes": rng.normal(immune_genes, immune_genes / 8, n_immune),
                "n_counts": rng.normal(immune_counts, immune_counts / 8, n_immune),
                "pct_mito": rng.normal(5, 2, n_immune).clip(0),
            }
        )
        epithelium = pd.DataFrame(
            {
                "study_id": study_id,
                "compartment": "Epithelial cells",
                "n_genes": rng.normal(immune_genes * 3.9, immune_genes / 4, n_epithelium),
                "n_counts": rng.normal(immune_counts * 3.9, immune_counts / 4, n_epithelium),
                "pct_mito": rng.normal(9, 2, n_epithelium).clip(0),
            }
        )
        return pd.concat([immune, epithelium], ignore_index=True)

    metrics = pd.concat(
        [study("GSE132465", 1500, 5000), study("GSE144735", 3000, 9000)],
        ignore_index=True,
    ).reset_index(drop=True)
    # Plant unambiguous outliers, one per study, inside a compartment.
    metrics.loc[0, "n_genes"] = 5  # SMC immune: near-empty droplet
    metrics.loc[1000, "n_counts"] = 200000  # KUL3 immune: doublet-scale library
    metrics.loc[5, "pct_mito"] = 80  # dying cell, either study
    return metrics


def test_qc_flags_catches_planted_outliers_in_both_studies():
    flags = qc_flags(_metrics())
    assert flags.loc[0]  # near-empty SMC droplet
    assert flags.loc[1000]  # doublet-scale KUL3 library
    assert flags.loc[5]  # high-mito cell


def test_qc_flags_thresholds_are_computed_per_study_not_globally():
    """A cell with GSE132465-typical depth must not be flagged just because
    GSE144735 runs deeper -- that would be a global, not per-study, cutoff."""
    metrics = _metrics()
    smc_typical = metrics[(metrics["study_id"] == "GSE132465") & ~metrics.index.isin([0, 5])]
    flags = qc_flags(metrics)
    assert not flags.loc[smc_typical.index].any()


def test_qc_flags_requires_the_expected_columns():
    bad = _metrics().drop(columns=["pct_mito"])
    with pytest.raises(ValueError, match="missing column"):
        qc_flags(bad)


def test_qc_flags_mito_cap_is_a_single_hard_threshold():
    metrics = _metrics()
    lenient = qc_flags(metrics, max_pct_mito=100)
    strict = qc_flags(metrics, max_pct_mito=1)
    assert strict.sum() > lenient.sum()


def test_doublet_and_ambient_correction_are_explicit_stubs_pending_real_data():
    """These need actual GSE132465/GSE144735 matrices to fit -- not fakeable
    without them. Confirm they fail loudly, not silently."""
    with pytest.raises(NotImplementedError, match="real Lee count matrices"):
        flag_doublets()
    with pytest.raises(NotImplementedError, match="real Lee count matrices"):
        correct_ambient()


# ---------------------------------------------------------------------------
# docs/open_decisions.md #12 -- the grouping, and the check that catches it
# ---------------------------------------------------------------------------


def test_deep_epithelium_is_not_an_outlier_for_being_epithelium():
    """The measured artifact: pooled across compartments, the upper MAD bound
    fires on epithelial cells because immune cells set the median. On SMC that
    cut the tumour arm 29.6 points harder than the normal arm."""
    metrics = _metrics()
    epithelium = metrics.index[metrics["compartment"] == "Epithelial cells"]
    flags = qc_flags(metrics)
    assert not flags.loc[epithelium].any()

    pooled = pd.Series(False, index=metrics.index)
    for _, group in metrics.groupby("study_id"):
        for col in ("n_genes", "n_counts"):
            values = group[col]
            mad = (values - values.median()).abs().median()
            z = 0.6745 * (values - values.median()) / mad
            pooled.loc[group.index] |= z.abs() > 5.0
    assert pooled.loc[epithelium].any(), "fixture no longer reproduces the artifact"


def test_qc_flags_requires_the_compartment_rather_than_defaulting_it_away():
    """Not optional on purpose -- defaulting it restores the pooled grouping,
    and the bias it produces points at the hypothesis."""
    bad = _metrics().drop(columns=["compartment"])
    with pytest.raises(ValueError, match="compartment"):
        qc_flags(bad)


def _retention_frame():
    rows = []
    for patient, normal_pass, tumour_pass in [
        ("P1", [True] * 10, [True] * 10),  # even
        ("P2", [True] * 10, [True] * 5 + [False] * 5),  # tumour cut harder
    ]:
        for tissue, passes in (("normal", normal_pass), ("tumour", tumour_pass)):
            for ok in passes:
                rows.append(
                    {
                        "patient_id": patient,
                        "tissue": tissue,
                        "compartment": "Epithelial cells",
                        "passed": ok,
                    }
                )
    frame = pd.DataFrame(rows)
    return frame.drop(columns=["passed"]), frame["passed"]


def test_differential_retention_flags_an_uneven_arm():
    metrics, passes = _retention_frame()
    out = differential_retention(metrics, passes).set_index("patient_id")
    assert not out.loc["P1", "flagged"]
    assert out.loc["P2", "flagged"]
    assert out.loc["P2", "gap_pts"] == pytest.approx(-50.0)


def test_differential_retention_reports_the_sign_not_just_the_size():
    """Both directions bias the compositional term and they are different
    findings -- W1 measured normal losing more, Lee has tumour losing more."""
    metrics, passes = _retention_frame()
    out = differential_retention(metrics, passes).set_index("patient_id")
    assert out.loc["P2", "gap_pts"] < 0  # tumour retained less


def test_differential_retention_restricts_to_the_compositional_compartment():
    """The compositional arm is built from epithelium; a gap in the T cells is
    not the quantity decision #12 asks about."""
    metrics, passes = _retention_frame()
    metrics = metrics.copy()
    metrics["compartment"] = "T cells"
    out = differential_retention(metrics, passes)
    assert out.empty
    assert len(differential_retention(metrics, passes, compartment=None)) == 2


def test_differential_retention_rejects_a_mismatched_mask():
    metrics, passes = _retention_frame()
    with pytest.raises(ValueError, match="entries for"):
        differential_retention(metrics, passes.iloc[:-1])
