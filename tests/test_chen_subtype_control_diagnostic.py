"""The post-hoc diagnostic describes measurement properties and nothing else."""

from __future__ import annotations

import pandas as pd

from src.reference.jobs.chen_subtype_control_diagnostic import COVARIATES, describe


def _inputs(rows):
    return pd.DataFrame(
        [
            {"patient_id": p, "gene": g, "granularity_rung": "lineage",
             "depth_normal": d, "frac_mature_tumour": f}
            for p, g, d, f in rows
        ]
    )


def test_no_covariate_is_a_transcript_programme():
    """§8 bars the interesting readings; this must not smuggle one in."""
    for name in COVARIATES:
        assert any(
            name.startswith(prefix)
            for prefix in ("depth", "n_", "unresolved", "frac_mature")
        ), name


def test_a_patient_is_counted_once_not_once_per_gene():
    """These are patient-level properties repeated across gene rows."""
    rows = [("P1", g, 100.0, 0.5) for g in ("ACTB", "CDX2", "GUCA2A")]
    rows += [("P2", g, 200.0, 0.9) for g in ("ACTB", "CDX2", "GUCA2A")]
    assigned = pd.DataFrame({"patient_id": ["P1", "P2"], "arm": ["AD", "SER"]})
    out = describe(_inputs(rows), assigned).set_index("covariate")
    assert out.loc["depth_normal", "n_AD"] == 1
    assert out.loc["depth_normal", "n_SER"] == 1


def test_the_output_declares_itself_post_hoc_and_interval_free():
    rows = [("P1", "ACTB", 100.0, 0.5), ("P2", "ACTB", 200.0, 0.9)]
    assigned = pd.DataFrame({"patient_id": ["P1", "P2"], "arm": ["AD", "SER"]})
    out = describe(_inputs(rows), assigned)
    assert out["post_hoc"].all()
    assert not out["interval_reported"].any()
    assert "ci_low" not in out.columns
