"""Family G closes on the floor, from open labels, with no interval."""

from __future__ import annotations

import pandas as pd

from src.reference.jobs.chen_subtype_family_g import TRUNCATING, arms, floor_table


def _crosswalk(rows):
    return pd.DataFrame(
        [{"scRNA_biospecimen_id": s, "patient_id": p, "specimen_exact": e} for s, p, e in rows]
    )


def _mutations(rows):
    return pd.DataFrame(
        [{"sampleId": s, "gene": g, "mutation_type": t, "protein_change": c} for s, g, t, c in rows]
    )


def test_arms_are_positive_only_and_a_double_positive_joins_neither():
    cw = _crosswalk([("S1", "P1", True), ("S2", "P2", True)])
    mut = _mutations([
        ("S1", "APC", "Nonsense_Mutation", "p.X"),
        ("S1", "BRAF", "Missense_Mutation", "V600E"),
        ("S2", "BRAF", "Missense_Mutation", "V600E"),
    ])
    out = arms(mut, cw)
    assert out["APC_truncating"] == set()
    assert out["BRAF_V600E"] == {"P2"}


def test_a_non_exact_specimen_is_not_scoped_in():
    cw = _crosswalk([("S1", "P1", False)])
    mut = _mutations([("S1", "APC", "Nonsense_Mutation", "p.X")])
    assert arms(mut, cw)["APC_truncating"] == set()


def test_a_missense_apc_is_not_truncating():
    assert "Missense_Mutation" not in TRUNCATING
    cw = _crosswalk([("S1", "P1", True)])
    mut = _mutations([("S1", "APC", "Missense_Mutation", "p.X")])
    assert arms(mut, cw)["APC_truncating"] == set()


def test_an_arm_below_the_floor_reports_no_interval():
    dec = pd.DataFrame(
        [{"patient_id": f"P{i}", "granularity_rung": "lineage", "intrinsic": 1.0}
         for i in range(3)]
    )
    table = floor_table({"APC_truncating": {"P0", "P1"}, "BRAF_V600E": {"P2"}}, dec)
    assert not table["clears_floor"].any()
    assert table["mean_difference"].isna().all()
    assert set(table["estimability"]) == {"not_estimable"}
