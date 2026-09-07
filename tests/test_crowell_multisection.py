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
                   detections=None, depths=None, n_cells=None, floors=None):
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
    default_cells = {dom: 500 for per in seps.values() for dom in per}
    n_cells = {**default_cells, **(n_cells or {})}
    default_floors = {dom: 0.01 for per in seps.values() for dom in per}
    floors = {**default_floors, **(floors or {})}
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
                "n_cells": n_cells[dom],
                "detection": det,
                "floor_per_probe_mean": floors[dom],
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

    def mu(p):
        return -np.log1p(-p)
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


def test_the_discriminators_conclusion_turning_on_a_bound_is_also_indeterminate():
    """Amendment 6. §7 rests on TWO conclusions and Amendment 4 protected one.

    Block 222 has CDX2 below the floor in its reference, and a below-floor
    reference inflates that block's DiD upward — toward the direction §7's pass
    needs. The rule is symmetric whether or not it fires.
    """
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        root = pathlib.Path(tmp)
        # target falls cleanly in all four; discriminator only clears zero
        # BECAUSE of the block where its reference is below the floor
        for i, pre in enumerate(("110", "120", "210")):
            _write_section(root, f"r{i}", pre, f"{pre}_TVA", f"{pre}_REF",
                           {"KRT8": {f"{pre}_REF": 3.0, f"{pre}_TVA": 3.5},
                            TARGET_GENE: {f"{pre}_REF": 2.0, f"{pre}_TVA": 0.4},
                            DISCRIMINATOR_GENE: {f"{pre}_REF": 2.0,
                                                 f"{pre}_TVA": 2.55 + 0.02 * i}})
        _write_section(root, "r9", "222", "222_TVA", "222_REF",
                       {"KRT8": {"222_REF": 3.0, "222_TVA": 3.5},
                        TARGET_GENE: {"222_REF": 2.0, "222_TVA": 0.4},
                        DISCRIMINATOR_GENE: {"222_REF": -0.113, "222_TVA": 2.4}})
        summary = aggregate(per_block_did(root))

    row = summary.set_index("gene").loc[DISCRIMINATOR_GENE]
    assert row["n_blocks_below_floor"] == 1
    assert bool(row["excludes_zero"]) != bool(row["excludes_zero_excluding_floored"])
    out = verdict(summary)
    assert out["verdict"] == ("INDETERMINATE — THE DISCRIMINATOR TURNS ON "
                             "BELOW-FLOOR BLOCKS")
    assert "rests on the discriminator as much as on the target" in out["detail"]


def test_a_discriminator_clearing_zero_the_other_way_gets_its_own_verdict():
    """Amendment 7. The final branch is reached in TWO states — the
    discriminator not clearing zero, and clearing it in the opposite direction
    — and the first version asserted "does not exclude zero" in both."""
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        root = pathlib.Path(tmp)
        for i, pre in enumerate(("110", "120", "210", "221")):
            _write_section(root, f"r{i}", pre, f"{pre}_TVA", f"{pre}_REF",
                           {"KRT8": {f"{pre}_REF": 3.0, f"{pre}_TVA": 3.5},
                            TARGET_GENE: {f"{pre}_REF": 2.0, f"{pre}_TVA": 0.4 + 0.02 * i},
                            # discriminator rises well clear of zero, opposite sign
                            DISCRIMINATOR_GENE: {f"{pre}_REF": 1.0,
                                                 f"{pre}_TVA": 3.5 + 0.02 * i}})
        summary = aggregate(per_block_did(root))

    disc = summary.set_index("gene").loc[DISCRIMINATOR_GENE]
    tgt = summary.set_index("gene").loc[TARGET_GENE]
    assert bool(disc["excludes_zero"]) and disc["mean_did"] > 0
    assert bool(tgt["excludes_zero"]) and tgt["mean_did"] < 0

    out = verdict(summary)
    assert out["verdict"] == "TARGET FALLS AND THE DISCRIMINATOR MOVES THE OTHER WAY"
    assert "OPPOSITE direction" in out["detail"]
    # and it must refuse to be read as a mechanism claim
    assert "not a stronger claim about mechanism" in out["detail"]
    assert "more CDX2-positive cells" in out["detail"]


def test_the_does_not_exclude_zero_branch_prints_the_interval_it_claims():
    """A verdict that states a conclusion must show the interval behind it."""
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        root = pathlib.Path(tmp)
        rng = np.random.default_rng(5)
        for i, pre in enumerate(("110", "120", "210", "221")):
            _write_section(root, f"r{i}", pre, f"{pre}_TVA", f"{pre}_REF",
                           {"KRT8": {f"{pre}_REF": 3.0, f"{pre}_TVA": 3.5},
                            TARGET_GENE: {f"{pre}_REF": 2.0, f"{pre}_TVA": 0.4 + 0.02 * i},
                            DISCRIMINATOR_GENE: {f"{pre}_REF": 2.0,
                                                 f"{pre}_TVA": 2.5 + float(rng.normal(0, 0.9))}})
        summary = aggregate(per_block_did(root))

    out = verdict(summary)
    assert out["verdict"] == "TARGET FALLS AND THE DISCRIMINATOR DOES NOT"
    assert "does not exclude zero" in out["detail"]
    assert "[" in out["detail"].split(DISCRIMINATOR_GENE)[-1], "show the interval"


# ---------------------------------------------------------------------------
# Pre-registered rules the code did not implement (found in review, 2026-09-07)
# ---------------------------------------------------------------------------


def test_rule_3_excludes_a_block_whose_control_fails_the_bar():
    """§5 RULE 3, WHICH THE CODE DID NOT ENFORCE.

    The rule is "a KRT8 separation clearing MIN_LOG_SEPARATION in BOTH" — the
    sensitivity gate. The first version checked only that KRT8 EXISTS in both
    domains. It never bit (KRT8's worst separation over six blocks is 2.804
    against a bar of 1.099), but an inclusion rule that cannot exclude is a
    check that cannot fail.
    """
    import tempfile

    from src.reference.jobs.crowell_feasibility import MIN_LOG_SEPARATION

    with tempfile.TemporaryDirectory() as tmp:
        root = pathlib.Path(tmp)
        _write_section(root, "r1", "110", "110_TVA", "110_REF",
                       {"KRT8": {"110_REF": 3.0, "110_TVA": 3.5},
                        TARGET_GENE: {"110_REF": 2.0, "110_TVA": 0.5}})
        # the control fails the bar in the reference: this block must not enter
        _write_section(root, "r2", "120", "120_TVA", "120_REF",
                       {"KRT8": {"120_REF": MIN_LOG_SEPARATION - 0.01,
                                 "120_TVA": 3.5},
                        TARGET_GENE: {"120_REF": 2.0, "120_TVA": 0.5}})
        got = per_block_did(root)

    assert set(got["block"]) == {SECTION_TO_BLOCK["110"]}
    assert SECTION_TO_BLOCK["120"] not in set(got["block"])


def test_rule_3_gates_on_the_control_and_never_on_the_target():
    """A block whose TARGET is far below the bar still enters — that is the
    expected outcome in lesion domains and excluding it would select on the
    outcome. §5 says so in terms."""
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        root = pathlib.Path(tmp)
        _write_section(root, "r1", "210", "210_TVA", "210_REF",
                       {"KRT8": {"210_REF": 3.0, "210_TVA": 3.9},
                        TARGET_GENE: {"210_REF": -0.5, "210_TVA": -0.9}})
        got = per_block_did(root)

    assert SECTION_TO_BLOCK["210"] in set(got["block"]), "gated on the control only"


def test_a_target_that_rises_gets_section_7s_third_row_not_a_pass():
    """§7's third falsifier, which the code did not implement.

    Without it a positive target falls through to "A TIER MOVED" or "TARGET
    FALLS…", both false statements about the direction. Latent — GUCA2A is
    negative in all six blocks — and a falsifier that cannot fire is not one.
    """
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        root = pathlib.Path(tmp)
        for i, pre in enumerate(("110", "120", "210", "221")):
            _write_section(root, f"r{i}", pre, f"{pre}_TVA", f"{pre}_REF",
                           {"KRT8": {f"{pre}_REF": 3.0, f"{pre}_TVA": 3.5},
                            # target RISES against its control
                            TARGET_GENE: {f"{pre}_REF": 1.0, f"{pre}_TVA": 3.6 + 0.02 * i},
                            DISCRIMINATOR_GENE: {f"{pre}_REF": 2.0, f"{pre}_TVA": 2.4}})
        summary = aggregate(per_block_did(root))

    out = verdict(summary)
    assert out["verdict"].startswith("TARGET RISES")
    assert "opposite of §7's prediction" in out["detail"]
    assert "report it and stop" in out["detail"].lower()


def test_the_control_gene_is_pinned_and_changing_it_changes_the_estimand():
    """M2. The only test touching CONTROL_GENE asserted its DiD is zero, which
    is true of ANY control by construction — so swapping KRT8 for EPCAM would
    have changed the estimand with a green suite.

    ACTB is absent from this deposit (feasibility Amendment 1 §3), so KRT8 is
    the only `control`-role gene present, and Amendment 2 turns on it being an
    epithelial keratin. It is not interchangeable with EPCAM, which is
    `epithelial` by role and is reported, never used as a control.
    """
    from src.reference.jobs.coexpression_silencing import GENE_ROLES

    assert CONTROL_GENE == "KRT8"
    assert GENE_ROLES[CONTROL_GENE] == "control"
    assert GENE_ROLES["EPCAM"] == "epithelial", "not a control, by the frozen roles"
    assert "ACTB" not in (CONTROL_GENE, DISCRIMINATOR_GENE, TARGET_GENE)
    assert TARGET_GENE == "GUCA2A" and DISCRIMINATOR_GENE == "CDX2"


def test_carcinoma_is_secondary_and_pools_same_class_domains_at_cell_level(tmp_path):
    """§5 secondary contrast: no invented seventh CRC, no post-pool average.

    Section 120 has two CRC regions. The pooled detection is the cell-count
    weighted fraction before log separation, not an average of two log
    separations. Section 222 has no CRC label and must therefore not contribute
    a made-up carcinoma observation.
    """
    floor = 0.01

    def separation(p):
        return float(np.log((-np.log1p(-p)) / (-np.log1p(-floor))))

    detections = {
        "KRT8": {"120_REF": 0.50, "120_TVA": 0.55,
                 "120_CRC1": 0.40, "120_CRC2": 0.80},
        TARGET_GENE: {"120_REF": 0.30, "120_TVA": 0.20,
                      "120_CRC1": 0.20, "120_CRC2": 0.10},
        DISCRIMINATOR_GENE: {"120_REF": 0.30, "120_TVA": 0.32,
                             "120_CRC1": 0.20, "120_CRC2": 0.40},
    }
    seps = {
        gene: {domain: separation(value) for domain, value in values.items()}
        for gene, values in detections.items()
    }
    _write_section(
        tmp_path, "r1", "120", "120_TVA", "120_REF", seps,
        detections=detections,
        n_cells={"120_REF": 500, "120_TVA": 500, "120_CRC1": 300, "120_CRC2": 900},
    )
    _write_section(
        tmp_path, "r2", "222", "222_TVA", "222_REF",
        {"KRT8": {"222_REF": 3.0, "222_TVA": 3.5},
         TARGET_GENE: {"222_REF": 2.0, "222_TVA": 0.5}},
    )

    got = per_block_did(tmp_path, lesion_class="carcinoma")

    assert set(got["block"]) == {SECTION_TO_BLOCK["120"]}
    target = got.set_index("gene").loc[TARGET_GENE]
    pooled_target = (300 * 0.20 + 900 * 0.10) / 1200
    pooled_control = (300 * 0.40 + 900 * 0.80) / 1200
    expected = (
        separation(pooled_target) - separation(0.30)
        - (separation(pooled_control) - separation(0.50))
    )
    assert target["detection_carcinoma"] == pytest.approx(pooled_target)
    assert target["did"] == pytest.approx(expected)
    assert target["carcinoma_domain"] == "120_CRC1+120_CRC2"
    assert not bool(target["depth_summary_available"])


def test_the_carcinoma_verdict_refuses_comparison_with_the_adenoma_numbers():
    """The predictable misreading. GUCA2A is −1.307 in the adenoma contrast and
    −1.061 in the carcinoma one, and someone will read that as progression.

    They share a reference domain, so they are correlated rather than
    independent; they run over different blocks (222 has no carcinoma domain)
    and different n. The verdict has to say so, because the two tables sit next
    to each other in the same results directory.
    """
    from src.reference.jobs.crowell_multisection import secondary_verdict

    out = secondary_verdict()
    assert out["verdict"].startswith("SECONDARY")
    assert "NOT comparable to the adenoma numbers" in out["detail"]
    assert "same reference" in out["detail"].lower()
    assert "pre-registered" in out["detail"]
