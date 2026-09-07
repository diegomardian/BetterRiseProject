"""§4's DiD across blocks, and the interval it refuses to report.

Pre-registration: ``docs/prereg_crowell_multisection.md`` + Amendments 1–2.
"""

from __future__ import annotations

import json
import pathlib

import numpy as np
import pandas as pd
import pytest

from src.harness.meta import MIN_STUDIES
from src.reference.crowell_io import SECTION_TO_BLOCK, CrowellError
from src.reference.jobs.crowell_multisection import (
    CONTROL_GENE,
    DISCRIMINATOR_GENE,
    TARGET_GENE,
    aggregate,
    per_block_did,
    resolve_label,
    verdict,
)


def _write_section(root, run, section, adenoma, reference, seps,
                   detections=None, depths=None):
    """A per-section table in the shape the real job emits.

    `detection` and `median_counts_per_cell` are required columns: the
    aggregator records Amendment 2's dlog(mu) and depth ratio from them rather
    than leaving them to be recomputed by hand, which is how Becker's KRT8 got
    into a Crowell row.
    """
    d = root / run
    d.mkdir(parents=True, exist_ok=True)
    default_depth = {dom: 600.0 for per in seps.values() for dom in per}
    depths = {**default_depth, **(depths or {})}
    rows = []
    for g, per in seps.items():
        for dom, v in per.items():
            det = (detections or {}).get(g, {}).get(dom)
            if det is None:
                # a detection consistent with the separation, so dlog(mu) is
                # not degenerate in the fixture
                det = float(1 - np.exp(-np.exp(v) * 0.01))
            rows.append({
                "gene": g, "domain": dom, "log_separation": v,
                "role": {"KRT8": "control", "EPCAM": "epithelial",
                         "CDX2": "identity", "MS4A12": "identity",
                         "GUCA2A": "target"}[g],
                "detection": det,
                "median_counts_per_cell": depths[dom],
            })
    pd.DataFrame(rows).to_parquet(d / "crowell_feasibility.parquet")
    (d / "crowell_feasibility.meta.json").write_text(json.dumps({
        "section": f"{section}.h5ad", "adenoma_domain": adenoma,
        "reference_domain": reference}))
    return d


def test_the_did_cancels_the_false_positive_floor(tmp_path):
    """log_sep is log(mu_gene) - log(mu_floor), so the DiD is a difference of
    dlog(mu) and the floor drops out. Section 110's floor moves 0.0102 -> 0.0173
    between domains, which would otherwise sit inside every number."""
    _write_section(tmp_path, "r1", "231", "231_TVA", "231_REF", {
        "KRT8": {"231_REF": 3.0, "231_TVA": 3.6},
        TARGET_GENE: {"231_REF": 2.0, "231_TVA": 1.0}})
    got = per_block_did(tmp_path).set_index("gene")
    assert got.loc[CONTROL_GENE, "did"] == pytest.approx(0.0)
    assert got.loc[TARGET_GENE, "did"] == pytest.approx((1.0 - 2.0) - (3.6 - 3.0))


def test_two_blocks_get_no_interval(tmp_path):
    """Prereg §6: below MIN_STUDIES the per-block values go side by side and
    there is no interval. A refused interval is not an interval of zero width."""
    for run, sec, pre in (("r1", "231", "231"), ("r2", "110", "110")):
        _write_section(tmp_path, run, sec, f"{pre}_TVA", f"{pre}_REF", {
            "KRT8": {f"{pre}_REF": 3.0, f"{pre}_TVA": 3.5},
            TARGET_GENE: {f"{pre}_REF": 2.0, f"{pre}_TVA": 0.9}})
    summary = aggregate(per_block_did(tmp_path)).set_index("gene")
    assert summary.loc[TARGET_GENE, "n_blocks"] == 2 < MIN_STUDIES
    assert summary.loc[TARGET_GENE, "interval"] == "refused"
    assert summary.loc[TARGET_GENE, "ci_low"] is None
    out = verdict(summary.reset_index())
    assert out["verdict"].startswith("BELOW THE INTERVAL FLOOR")
    assert "not evidence of an effect size" in out["detail"]


def test_231_and_232_are_one_block_not_two(tmp_path):
    """Invariant 5, prereg §5. They are FOV ranges of one tissue block."""
    assert SECTION_TO_BLOCK["231"] == SECTION_TO_BLOCK["232"]
    for run, sec in (("r1", "231"), ("r2", "232")):
        _write_section(tmp_path, run, sec, f"{sec}_TVA", f"{sec}_REF", {
            "KRT8": {f"{sec}_REF": 3.0, f"{sec}_TVA": 3.5},
            TARGET_GENE: {f"{sec}_REF": 2.0, f"{sec}_TVA": 0.9}})
    summary = aggregate(per_block_did(tmp_path)).set_index("gene")
    assert summary.loc[TARGET_GENE, "n_blocks"] == 1, "one block, not two"


def test_the_discriminator_firing_beats_the_target_falling(tmp_path):
    """§7's decisive row. A target falling alongside CDX2 is a TIER moving, and
    that is not the claim — the same discriminator avenue A used."""
    for i, pre in enumerate(("110", "120", "210", "221")):
        _write_section(tmp_path, f"r{i}", pre, f"{pre}_TVA", f"{pre}_REF", {
            "KRT8": {f"{pre}_REF": 3.0, f"{pre}_TVA": 3.5},
            TARGET_GENE: {f"{pre}_REF": 2.0, f"{pre}_TVA": 0.4 + 0.02 * i},
            DISCRIMINATOR_GENE: {f"{pre}_REF": 2.0, f"{pre}_TVA": 0.5 + 0.02 * i}})
    out = verdict(aggregate(per_block_did(tmp_path)))
    assert out["verdict"] == "A TIER MOVED, NOT A GENE"
    assert "not gene-specific" in out["detail"]


def test_a_target_falling_while_the_discriminator_holds_is_the_pass(tmp_path):
    for i, pre in enumerate(("110", "120", "210", "221")):
        _write_section(tmp_path, f"r{i}", pre, f"{pre}_TVA", f"{pre}_REF", {
            "KRT8": {f"{pre}_REF": 3.0, f"{pre}_TVA": 3.5},
            TARGET_GENE: {f"{pre}_REF": 2.0, f"{pre}_TVA": 0.4 + 0.02 * i},
            DISCRIMINATOR_GENE: {f"{pre}_REF": 2.0, f"{pre}_TVA": 2.5 + 0.02 * i}})
    out = verdict(aggregate(per_block_did(tmp_path)))
    assert out["verdict"] == "TARGET FALLS AND THE DISCRIMINATOR DOES NOT"
    assert "not silencing" in out["detail"] and "not per-cell" in out["detail"]


def test_a_wide_interval_that_includes_zero_may_not_be_quoted_as_a_negative(tmp_path):
    rng = np.random.default_rng(3)
    for i, pre in enumerate(("110", "120", "210", "221")):
        jitter = float(rng.normal(0, 0.9))
        _write_section(tmp_path, f"r{i}", pre, f"{pre}_TVA", f"{pre}_REF", {
            "KRT8": {f"{pre}_REF": 3.0, f"{pre}_TVA": 3.5},
            TARGET_GENE: {f"{pre}_REF": 2.0, f"{pre}_TVA": 2.5 + jitter}})
    out = verdict(aggregate(per_block_did(tmp_path)))
    assert out["verdict"] == "NO CLAIM"
    assert "must not be" in out["detail"] and "negative" in out["detail"]


def test_a_sidecar_label_that_is_not_in_the_table_raises(tmp_path):
    """The sidecar records the CLI request and the table records the pooled
    label. A near-miss must not be guessed at — it would score another domain."""
    present = {"110_REF", "110_TVA1+110_TVA2+110_TVA3"}
    assert resolve_label("110_REF", present) == "110_REF"
    assert resolve_label(["110_TVA1", "110_TVA2", "110_TVA3"], present) == \
        "110_TVA1+110_TVA2+110_TVA3"
    assert resolve_label(None, present) is None
    with pytest.raises(CrowellError, match="Do not guess"):
        resolve_label(["110_TVA9"], present)


def test_an_unmapped_section_stops_the_run(tmp_path):
    """Blocks are the unit; a section with no block cannot be counted."""
    _write_section(tmp_path, "r1", "999", "999_TVA", "999_REF", {
        "KRT8": {"999_REF": 3.0, "999_TVA": 3.5},
        TARGET_GENE: {"999_REF": 2.0, "999_TVA": 0.9}})
    with pytest.raises(CrowellError, match="SECTION_TO_BLOCK"):
        per_block_did(tmp_path)


def test_amendment_2s_diagnosis_is_recorded_not_re_derived():
    """THE ERROR THIS COLUMN EXISTS TO PREVENT.

    Amendment 2's table was computed by hand from the per-section parquets and
    one row used Becker's normal-arm KRT8 (0.154186) in place of Crowell
    231_REF's (0.360144), giving +1.746 where the truth is +0.766. Every input
    was committed and re-derivable — that was never the gap. The gap was that
    checking it required recomputing it, and the recomputation is where the
    wrong number came from.
    """
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        root = pathlib.Path(tmp)
        _write_section(
            root, "r1", "231", "231_TVA", "231_REF",
            {"KRT8": {"231_REF": 3.281, "231_TVA": 3.923},
             TARGET_GENE: {"231_REF": 1.571, "231_TVA": 0.700}},
            detections={"KRT8": {"231_REF": 0.360144, "231_TVA": 0.617135},
                        TARGET_GENE: {"231_REF": 0.077585, "231_TVA": 0.037539}},
            depths={"231_REF": 575.0, "231_TVA": 821.0})
        got = per_block_did(root).set_index("gene")

    mu = lambda p: -np.log1p(-p)
    assert got.loc["KRT8", "dlog_mu"] == pytest.approx(
        float(np.log(mu(0.617135) / mu(0.360144))), abs=1e-6)
    assert got.loc["KRT8", "dlog_mu"] == pytest.approx(0.766, abs=1e-3)
    # and NOT the value the hand computation produced from Becker's number
    assert abs(got.loc["KRT8", "dlog_mu"] - 1.746) > 0.9

    assert got.loc["KRT8", "log_depth_ratio"] == pytest.approx(
        float(np.log(821 / 575)), abs=1e-6)
    # Amendment 2's split: the epithelial control outruns depth, the target does not
    assert bool(got.loc["KRT8", "rises_faster_than_depth"])
    assert not bool(got.loc[TARGET_GENE, "rises_faster_than_depth"])


def test_an_interval_whose_sign_turns_on_below_floor_blocks_is_indeterminate():
    """Amendment 4, fixed before the four-block aggregate existed.

    Block 210 has GUCA2A below the negative-probe floor in its lesion — the
    observed signal is smaller than what noise alone supplies, so the true rate
    is consistent with zero and the DiD is a BOUND. Including it flips the
    target's interval from including zero to excluding it. That is not a
    positive; it is indeterminate, and the check BINDS rather than warning.

    The three clean blocks must SPAN zero on their own or the branch is never
    reached — the first version of this fixture had them already excluding it.
    """
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        root = pathlib.Path(tmp)
        # DiD = (tva - ref) - (k_tva - k_ref); ref 2.0, control delta 0.5
        # -> DiDs of -0.9, -0.3, -0.6: mean -0.60, and at n=3 the interval spans 0
        for i, (pre, tva) in enumerate((("110", 1.6), ("120", 2.2), ("221", 1.9))):
            _write_section(root, f"r{i}", pre, f"{pre}_TVA", f"{pre}_REF",
                           {"KRT8": {f"{pre}_REF": 3.0, f"{pre}_TVA": 3.5},
                            TARGET_GENE: {f"{pre}_REF": 2.0, f"{pre}_TVA": tva}})
        # 210: target BELOW the floor in the lesion (negative separation),
        # DiD = (-0.35 - 0.6) - 0.55 = -1.50
        _write_section(root, "r9", "210", "210_TVA", "210_REF",
                       {"KRT8": {"210_REF": 3.0, "210_TVA": 3.55},
                        TARGET_GENE: {"210_REF": 0.6, "210_TVA": -0.35}})
        per_block = per_block_did(root)
        summary = aggregate(per_block)

    assert bool(per_block.set_index(["gene", "block"]).loc[
        (TARGET_GENE, SECTION_TO_BLOCK["210"]), "below_floor_adenoma"])
    row = summary.set_index("gene").loc[TARGET_GENE]
    assert row["n_blocks_below_floor"] == 1
    assert bool(row["excludes_zero"]), "the floored block flips it"
    assert not bool(row["excludes_zero_excluding_floored"]), "without it, zero is in"

    out = verdict(summary)
    assert out["verdict"].startswith("INDETERMINATE")
    assert "BOUND and not a point" in out["detail"]
    assert "Neither interval is the answer" in out["detail"]


def test_below_floor_blocks_that_do_not_change_the_sign_do_not_trigger_it():
    """The other half: the rule must not fire on every floored block, or it is
    an exclusion rule wearing a sensitivity's clothes."""
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        root = pathlib.Path(tmp)
        # three blocks that already exclude zero on their own
        for i, (pre, tva) in enumerate((("110", 1.3), ("120", 1.55), ("221", 1.05))):
            _write_section(root, f"r{i}", pre, f"{pre}_TVA", f"{pre}_REF",
                           {"KRT8": {f"{pre}_REF": 3.0, f"{pre}_TVA": 3.5},
                            TARGET_GENE: {f"{pre}_REF": 2.0, f"{pre}_TVA": tva}})
        _write_section(root, "r9", "210", "210_TVA", "210_REF",
                       {"KRT8": {"210_REF": 3.0, "210_TVA": 3.55},
                        TARGET_GENE: {"210_REF": 0.6, "210_TVA": -0.35}})
        summary = aggregate(per_block_did(root))

    row = summary.set_index("gene").loc[TARGET_GENE]
    assert row["n_blocks_below_floor"] == 1
    assert bool(row["excludes_zero"]) == bool(row["excludes_zero_excluding_floored"])
    assert not verdict(summary)["verdict"].startswith("INDETERMINATE")


def test_a_value_just_above_the_floor_is_reported_and_does_not_move_the_flag():
    """Amendment 5. Block 221's CDX2 reference separates at +0.0033 — a signal
    0.3% above noise — and CDX2 is §7's decisive row.

    The flag stays a hard `< 0`. Adding a tolerance the first time a value lands
    near the edge is the tuning this design exists to prevent. What is added is
    `min_separation`, so a reader can see it.
    """
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        root = pathlib.Path(tmp)
        _write_section(root, "r1", "221", "221_TVA", "221_REF",
                       {"KRT8": {"221_REF": 2.804, "221_TVA": 3.710},
                        DISCRIMINATOR_GENE: {"221_REF": 0.0033, "221_TVA": 1.808},
                        TARGET_GENE: {"221_REF": 0.951, "221_TVA": 0.439}})
        per_block = per_block_did(root).set_index("gene")

    # above the floor, so the hard flag does NOT fire
    assert not bool(per_block.loc[DISCRIMINATOR_GENE, "target_below_floor"])
    # but how close it sat is on the record
    assert per_block.loc[DISCRIMINATOR_GENE, "min_separation"] == pytest.approx(0.0033)
    assert per_block.loc[TARGET_GENE, "min_separation"] == pytest.approx(0.439)
    # and the discriminator's change is the largest, almost all of it the floor
    assert per_block.loc[DISCRIMINATOR_GENE, "did"] > per_block.loc[TARGET_GENE, "did"]
