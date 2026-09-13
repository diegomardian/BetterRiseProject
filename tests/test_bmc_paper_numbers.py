"""Every number in the BMC Results section, re-derived from its table.

BMC requires that no figure be transcribed. This is the BMC counterpart of
``tests/test_paper_numbers.py`` for the WMHS source: each assertion is paired
with a literal from ``paper/bmc/sections/results.tex``, so editing the prose
without the table (or the table without the prose) fails.

The Results and Implementation sections are covered, plus the abstract's numbers
against Results; the remaining sections are still contracts, and this file grows
with them.
"""

from __future__ import annotations

import json
import re

import pandas as pd
import pytest

from src.common.paths import REPO_ROOT, RESULTS_DIR
from src.reference.table_resolution import newest_by_time

RESULTS_TEX = REPO_ROOT / "paper" / "bmc" / "sections" / "results.tex"


def _quotes(tex: str, literal: str) -> None:
    assert literal in tex, f"prose does not contain: {literal!r}"


@pytest.fixture(scope="module")
def results_tex() -> str:
    return RESULTS_TEX.read_text(encoding="utf-8")


@pytest.fixture(scope="module")
def stress_summary() -> pd.DataFrame:
    path = newest_by_time(RESULTS_DIR, "interval_stress_calibration_summary")
    assert path is not None
    return pd.read_parquet(path)


@pytest.fixture(scope="module")
def stress_cells() -> pd.DataFrame:
    path = newest_by_time(RESULTS_DIR, "interval_stress_calibration")
    assert path is not None
    return pd.read_parquet(path)


@pytest.fixture(scope="module")
def brackets() -> pd.DataFrame:
    path = newest_by_time(RESULTS_DIR, "cutpoint_crossing_brackets_summary")
    assert path is not None
    return pd.read_parquet(path)


# ---------------------------------------------------------------------------
# Results §1 — the interval stress calibration
# ---------------------------------------------------------------------------


def test_the_stress_table_is_the_prose_table(results_tex, stress_summary):
    """The method-by-regime medians, to one decimal, are the table's."""
    median = stress_summary.pivot_table(
        index="method", columns="regime", values="fpr_median", aggfunc="median"
    )
    expected = {
        "percentile": {"gaussian": 7.0, "skewed": 7.9, "heavy_tailed": 7.8,
                       "spike_slab": 6.6, "boundary_rare": 6.4,
                       "boundary_saturated": 0.0},
        "bca": {"gaussian": 7.1, "skewed": 7.8, "heavy_tailed": 8.6,
                "spike_slab": 7.1, "boundary_rare": 7.1,
                "boundary_saturated": 0.0},
        "student_t": {"gaussian": 4.8, "skewed": 5.9, "heavy_tailed": 4.9,
                      "spike_slab": 4.4, "boundary_rare": 4.9,
                      "boundary_saturated": 0.0},
    }
    for method, row in expected.items():
        for regime, value in row.items():
            assert round(100 * median.loc[method, regime], 1) == value
    # The prose table is reproduced, not just the computation.
    _quotes(results_tex, r"percentile  & 7.0 & 7.9 & 7.8 & 6.6 & 6.4 & 0.0")
    _quotes(results_tex, r"Student-$t$ & 4.8 & 5.9 & 4.9 & 4.4 & 4.9 & 0.0")


def test_the_gaussian_control_reproduces_the_committed_rates(
    results_tex, stress_cells
):
    gaussian = stress_cells[
        (stress_cells["method"] == "percentile")
        & (stress_cells["regime"] == "gaussian")
    ]
    methylated = gaussian[gaussian["cohort"] == "mlh1_methylated"]
    lineage = gaussian[gaussian["cohort"] == "adenoma_lineage"]
    assert round(100 * methylated["false_positive_rate"].median(), 1) == 9.9
    assert round(100 * lineage["false_positive_rate"].median(), 1) == 6.1

    committed_path = newest_by_time(RESULTS_DIR, "interval_calibration")
    committed = pd.read_parquet(committed_path)
    committed = committed[
        (committed["method"] == "percentile")
        & (committed["abundance"] == "common")
        & (committed["tau"] == 0.2)
    ].set_index("cohort")["false_positive_rate"]
    assert round(100 * committed["mlh1_methylated"], 1) == 10.7
    assert round(100 * committed["adenoma_lineage"], 1) == 5.7

    _quotes(results_tex, r"$9.9\%$ on the MLH1-methylated arm against the committed $10.7\%$")


def test_the_boundary_regimes_and_the_check_families(results_tex, stress_cells,
                                                     stress_summary):
    saturated = stress_summary[stress_summary["regime"] == "boundary_saturated"]
    # Ordinary input validation fires on every saturated cell.
    assert saturated["detection_rate_input_validation"].eq(1.0).all()
    # The prose counts one method's cells: 5 seeds x 7 cohorts.
    per_method = saturated[saturated["method"] == "percentile"]["n_seeds"].sum()
    assert int(per_method) == 35

    rare_pct = stress_cells[
        (stress_cells["regime"] == "boundary_rare")
        & (stress_cells["method"] == "percentile")
    ]
    assert int((rare_pct["input_validation_failures"] > 0).sum()) == 1
    assert int(rare_pct["audit_flags"].astype(bool).sum()) == 17
    _quotes(results_tex, r"validation fires on $35$ of $35$ such cells")


def test_the_closed_form_catches_skew_but_misses_heavy_tails(
    results_tex, stress_cells
):
    percentile = stress_cells[stress_cells["method"] == "percentile"]

    def original_fires(regime: str) -> int:
        block = percentile[percentile["regime"] == regime]
        return int(block["original_check_flags"].map(lambda v: v is True).sum())

    assert original_fires("skewed") == 12
    assert original_fires("heavy_tailed") == 3
    # "departs from it under skew": the diagnostic flags more skew cells than the
    # gaussian control, which is what "departs" means operationally.
    assert original_fires("skewed") > original_fires("gaussian")

    gaussian = percentile[percentile["regime"] == "gaussian"]
    small = gaussian[gaussian["cohort"] == "mlh1_intact_mmrd"]
    assert round(100 * small["false_positive_rate"].median(), 1) == 19.1
    assert round(100 * small["closed_form_rate"].iloc[0], 1) == 18.8
    _quotes(results_tex, r"$19.1\%$ against $18.8\%$ at $n=4$")
    _quotes(results_tex, r"$12$ of $35$ percentile cells")
    _quotes(results_tex, r"$3$ of $35$")


def test_the_paper_does_not_overclaim_the_stress_result(results_tex):
    # The retired claims must not appear.
    assert "universally calibrated" in results_tex  # only inside the negation
    assert "does not establish that BCa is generally" in results_tex


# ---------------------------------------------------------------------------
# Results §2 — the cutpoint brackets
# ---------------------------------------------------------------------------


def _bracket_cell(brackets: pd.DataFrame, cohort: str, pool: str,
                  criterion: str) -> pd.Series:
    return brackets[
        (brackets["cohort"] == cohort)
        & (brackets["pool"] == pool)
        & (brackets["criterion"] == criterion)
    ].iloc[0]


def test_the_pooled_pool_recovers_a_wide_crossing_the_point_rule_drops(
    results_tex, brackets
):
    for cohort in ("smc", "kul3"):
        wide = _bracket_cell(brackets, cohort, "pooled", "wide")
        ok = _bracket_cell(brackets, cohort, "pooled", "ok")
        assert wide["n_with_upper"] == 8 and wide["bracket_stable"]
        assert (wide["lower_min"], wide["upper_max"]) == (42.5, 70.0)
        assert ok["n_with_upper"] == 0 and ok["n_not_identifiable"] == 8
    _quotes(results_tex, r"\emph{wide} bracket is $[42.5, 70]$ mature cells")
    _quotes(results_tex, r"stable across all eight seeds")


def test_smc_and_kul3_are_reported_separately(results_tex, brackets):
    smc_wide = _bracket_cell(brackets, "smc", "reference", "wide")
    smc_ok = _bracket_cell(brackets, "smc", "reference", "ok")
    assert (smc_wide["lower_min"], smc_wide["upper_max"]) == (27.5, 42.5)
    assert smc_wide["bracket_stable"]
    assert (smc_ok["lower_min"], smc_ok["upper_max"]) == (42.5, 70.0)
    assert smc_ok["bracket_stable"]

    kul3_ok = _bracket_cell(brackets, "kul3", "reference", "ok")
    assert (kul3_ok["lower_min"], kul3_ok["lower_max"],
            kul3_ok["upper_min"], kul3_ok["upper_max"]) == (150.0, 400.0, 400.0, 800.0)
    assert not kul3_ok["bracket_stable"]
    _quotes(results_tex, r"\emph{wide} $[27.5, 42.5]$ and \emph{ok} $[42.5, 70]$")
    _quotes(results_tex, r"$[150, 400]$ on one seed and $[400, 800]$ on the others")


# ---------------------------------------------------------------------------
# Results §3 — the adenoma sensitivity and attrition
# ---------------------------------------------------------------------------


def test_the_primary_portion_is_bit_identical_to_the_committed_scales(results_tex):
    sensitivity = pd.read_parquet(
        newest_by_time(RESULTS_DIR, "adenoma_claim_sensitivity")
    )
    resolved = sensitivity[sensitivity["denominator"] == "resolved"]
    committed = pd.read_parquet(
        newest_by_time(RESULTS_DIR, "adenoma_decomposition_scales")
    )
    keys = ["granularity_rung", "weighting", "statistic", "contrast"]
    merged = resolved.merge(committed, on=keys, suffixes=("", "_old"))
    assert len(merged) == 1260
    for column in ("centre", "ci_low", "ci_high"):
        assert (merged[column] - merged[f"{column}_old"]).abs().max() == 0.0
    assert bool((merged["excludes_zero"] == merged["excludes_zero_old"]).all())
    _quotes(results_tex, r"all 1260 rows agree on centre, interval, verdict and patient")


def test_the_denominator_changes_where_log_ratio_is_defined(results_tex):
    attrition = pd.read_parquet(
        newest_by_time(RESULTS_DIR, "estimability_attrition")
    )
    epithelial = attrition[attrition["granularity_rung"] == "epithelial"]
    resolved = int(
        epithelial[epithelial["denominator"] == "resolved"]["n_log_ratio_undefined"].sum()
    )
    all_epithelial = int(
        epithelial[epithelial["denominator"] == "all_epithelial"]["n_log_ratio_undefined"].sum()
    )
    assert resolved == 792
    assert all_epithelial == 4
    _quotes(results_tex, r"undefined on 792 of the epithelial")
    _quotes(results_tex, r"only 4 rows are undefined")


def test_the_cross_block_result_is_quoted_with_its_weighting(results_tex):
    summary = pd.read_parquet(
        newest_by_time(RESULTS_DIR, "adenoma_claim_sensitivity_summary")
    )
    lineage = summary[
        (summary["granularity_rung"] == "lineage")
        & (summary["statistic"] == "log_ratio")
    ]
    # The load-bearing statistic separates every cross-block contrast.
    assert set(lineage["n_excluding_zero"]) == {8}
    assert set(lineage["n_contrasts"]) == {8}

    scales = pd.read_parquet(
        newest_by_time(RESULTS_DIR, "adenoma_decomposition_scales")
    )
    cross = scales[(scales["granularity_rung"] == "lineage") & scales["cross_block"]]
    wide = cross.pivot_table(
        index=["weighting", "contrast"], columns="statistic",
        values="excludes_zero", aggfunc="first",
    )
    agree = {}
    for weighting, block in wide.groupby(level=0):
        block = block.droplevel(0)
        pair = block.index.str.split(" - ")
        unordered = block[[a < b for a, b in zip(pair.str[0], pair.str[1], strict=True)]]
        agree[weighting] = int(unordered.all(axis=1).sum())
    assert agree == {"doubly_robust": 6, "normal": 6, "tumour": 5}
    _quotes(results_tex, r"6 of 8 cross-block")
    _quotes(results_tex, r"5 of 8 under")


def test_attrition_is_reported_per_rung(results_tex):
    attrition = pd.read_parquet(
        newest_by_time(RESULTS_DIR, "estimability_attrition")
    )
    resolved = attrition[attrition["denominator"] == "resolved"]
    lineage = resolved[resolved["granularity_rung"] == "lineage"]
    best4 = resolved[resolved["granularity_rung"] == "best4"]
    assert lineage["n_patients"].max() == 44
    assert lineage["n_estimability_not_estimable"].max() == 1
    assert best4["n_patients"].max() == 20
    assert best4["n_estimability_not_estimable"].max() == 2
    _quotes(results_tex, r"retains 44 patients, one of whom is")
    _quotes(results_tex, r"retains 20, two of")


# ---------------------------------------------------------------------------
# Results §4 — the two reproductions
# ---------------------------------------------------------------------------


def test_both_reproductions_are_exact_and_labelled(results_tex):
    crowell = pd.read_parquet(
        newest_by_time(RESULTS_DIR, "crowell_summary_reproduction")
    )
    assert len(crowell) == 5
    assert bool(crowell["reproduces"].all())
    disval = pd.read_parquet(
        newest_by_time(RESULTS_DIR, "disval_reproduction")
    )
    assert len(disval) == 64
    assert bool(disval["reproduces"].all())
    _quotes(results_tex, r"all 5")
    _quotes(results_tex, r"genes reproduce exactly against the canonical run")
    _quotes(results_tex, r"all 64 cells")
    _quotes(results_tex, r"Neither is a raw-data replication")


# ---------------------------------------------------------------------------
# Results §5 — the resolver finding
# ---------------------------------------------------------------------------


def test_the_resolver_finding_and_its_consequence(results_tex):
    ambiguity = pd.read_parquet(
        newest_by_time(RESULTS_DIR, "table_resolution_ambiguity")
    )
    assert len(ambiguity) == 23
    assert int((ambiguity["disposition"] == "differs").sum()) == 9
    _quotes(results_tex, r"23 table/date pairs")
    _quotes(results_tex, r"9 of them hold different content")
    _quotes(results_tex, r"Neither rule is universally right")


# ---------------------------------------------------------------------------
# Results §6 — what could not be evaluated
# ---------------------------------------------------------------------------


def test_the_challenge_is_reported_as_unexecuted(results_tex):
    path = newest_by_time(RESULTS_DIR, "blinded_guard_challenge")
    sidecar = path.parent / "blinded_guard_challenge.meta.json"
    meta = json.loads(sidecar.read_text())
    assert meta["is_held_out_evaluation"] is False
    assert meta["stop_rule_readable_from_this_table"] is False
    frame = pd.read_parquet(path)
    assert len(frame) == 12
    _quotes(results_tex, r"The blinded challenge has not been executed.")
    _quotes(results_tex, r"the week-one stop rule is not readable from the table")


# ---------------------------------------------------------------------------
# Quantifier and ordering claims.
#
# A value assertion cannot see "smallest", "only", "every", "all" or "both": it
# binds the number while the ordering the prose asserts goes untested. That is
# how "the smallest cohort, n=10" shipped when the smallest cohort is n=4 and
# n=10 is merely the worst Student-t cell. Each test below re-derives the
# ordering from the table, and each is paired with the prose it guards.
# ---------------------------------------------------------------------------


def test_the_student_t_exception_is_the_worst_cell_not_the_smallest_cohort(
    results_tex, stress_summary
):
    student = stress_summary[stress_summary["method"] == "student_t"]
    worst = student.loc[student["fpr_median"].idxmax()]
    assert (worst["cohort"], worst["regime"]) == ("mlh1_methylated", "skewed")
    assert (worst["n_patients"], round(100 * worst["fpr_median"], 1)) == (10, 7.0)

    # The smallest cohort is n=4, and its skewed Student-t cell is not the worst.
    smallest = student.groupby("cohort")["n_patients"].first().sort_values().index[0]
    assert smallest == "mlh1_intact_mmrd"
    small = student[
        (student["cohort"] == smallest) & (student["regime"] == "skewed")
    ].iloc[0]
    assert small["n_patients"] == 4
    assert round(100 * small["fpr_median"], 1) == 5.9
    assert small["fpr_median"] < worst["fpr_median"]

    assert "smallest cohort" not in results_tex
    _quotes(results_tex, r"under skewed effects on the MLH1-methylated arm")


def test_student_t_median_is_within_tolerance_for_every_generator(
    results_tex, stress_summary
):
    student = stress_summary[stress_summary["method"] == "student_t"]
    median = student.groupby("regime")["fpr_median"].median()
    assert (median <= 0.07).all()
    _quotes(results_tex, r"with one exception")


def test_percentile_and_bca_over_reject_under_every_non_saturated_generator(
    results_tex, stress_summary
):
    for method in ("percentile", "bca"):
        median = (
            stress_summary[stress_summary["method"] == method]
            .groupby("regime")["fpr_median"].median()
        )
        assert (median.drop("boundary_saturated") > 0.05).all()
    _quotes(results_tex, r"Under every non-saturated generator")


def test_every_method_is_conservative_at_saturation(results_tex, stress_summary):
    # The paper's table is a median across cohorts, so the claim is about that
    # statistic, not about every individual cell (a few round to 0.1%).
    median = (
        stress_summary[stress_summary["regime"] == "boundary_saturated"]
        .groupby("method")["fpr_median"].median()
    )
    assert (100 * median).round(1).eq(0.0).all()
    _quotes(results_tex, r"every method is")


def test_input_validation_is_nearly_blind_in_the_interior(results_tex, stress_cells):
    interior = {"gaussian", "skewed", "heavy_tailed", "spike_slab"}
    fired = (
        stress_cells[stress_cells["regime"].isin(interior)]
        .groupby(["regime", "method"])["input_validation_failures"]
        .apply(lambda x: int((x > 0).sum()))
    )
    assert int(fired.max()) <= 1
    # The retired quantifier must not come back.
    assert "its only non-blind spot" not in results_tex
    _quotes(results_tex, r"In the interior it is nearly blind")


def test_the_grid_is_seven_cohorts_and_six_generators(stress_cells):
    assert stress_cells["cohort"].nunique() == 7
    assert stress_cells["regime"].nunique() == 6


def test_every_generator_regime_is_a_true_null():
    """`Every regime is a true null` — the effect has mean zero and unit scale.

    If a regime had a nonzero mean the measured rate would be a rejection rate
    under an alternative, not a false-positive rate, and every table cell would
    be misread.
    """
    import numpy as np

    from src.reference.interval_stress import REGIMES, standardised_effects

    for regime in REGIMES:
        draws = standardised_effects(regime, 200_000, np.random.default_rng(0))
        assert abs(float(draws.mean())) < 0.01, regime.name
        assert abs(float(draws.var()) - 1.0) < 0.02, regime.name


def test_input_validation_fires_on_at_most_three_near_zero_cells(
    results_tex, stress_cells
):
    rare = stress_cells[stress_cells["regime"] == "boundary_rare"]
    fired = rare.groupby("method")["input_validation_failures"].apply(
        lambda x: int((x > 0).sum())
    )
    assert int(fired.max()) == 3
    _quotes(results_tex, r"$3$ of $35$ cells")


def test_kul3_ok_bracket_is_unstable_on_exactly_one_seed(results_tex):
    brackets = pd.read_parquet(
        newest_by_time(RESULTS_DIR, "cutpoint_crossing_brackets")
    )
    kul3 = brackets[
        (brackets["cohort"] == "kul3")
        & (brackets["pool"] == "reference")
        & (brackets["criterion"] == "ok")
    ]
    counts = (
        kul3.groupby(["lower_n_cells_mature", "upper_n_cells_mature"]).size().to_dict()
    )
    assert counts == {(150.0, 400.0): 1, (400.0, 800.0): 7}
    _quotes(results_tex, r"$[150, 400]$ on one seed and $[400, 800]$ on the others")


def test_both_cohorts_are_present_and_never_pooled(results_tex, brackets):
    assert set(brackets["cohort"]) == {"smc", "kul3"}
    _quotes(results_tex, r"reported separately and never pooled")


def test_the_sensitivity_grid_is_the_full_cross_product(results_tex):
    sensitivity = pd.read_parquet(
        newest_by_time(RESULTS_DIR, "adenoma_claim_sensitivity")
    )
    assert sensitivity["weighting"].nunique() == 3
    assert sensitivity["granularity_rung"].nunique() == 4
    assert sensitivity["denominator"].nunique() == 2
    assert sensitivity["statistic"].nunique() == 4
    _quotes(results_tex, r"all three weightings, four rungs, both")
    _quotes(results_tex, r"all four scale-free statistics")


def test_the_load_bearing_statistic_separates_all_eight_under_every_cell(results_tex):
    summary = pd.read_parquet(
        newest_by_time(RESULTS_DIR, "adenoma_claim_sensitivity_summary")
    )
    lineage = summary[
        (summary["granularity_rung"] == "lineage")
        & (summary["statistic"] == "log_ratio")
    ]
    assert len(lineage) == 6  # 3 weightings x 2 denominators
    assert lineage["n_excluding_zero"].eq(8).all()
    assert lineage["n_contrasts"].eq(8).all()
    _quotes(results_tex, r"all 8 cross-block contrasts from zero under every weighting and both")


def test_the_load_bearing_statistic_does_not_separate_the_two_targets(results_tex):
    sensitivity = pd.read_parquet(
        newest_by_time(RESULTS_DIR, "adenoma_claim_sensitivity")
    )
    gm = sensitivity[sensitivity["contrast"] == "GUCA2A - MS4A12"]
    assert not gm[gm["statistic"] == "log_ratio"]["excludes_zero"].any()
    # One construction does separate them at lineage, so the prose is scoped to
    # the load-bearing statistic and must not claim "all statistics".
    signed = gm[
        (gm["statistic"] == "share_signed") & (gm["granularity_rung"] == "lineage")
    ]
    assert signed["excludes_zero"].any()
    _quotes(results_tex, r"the two target genes are not distinguishable at any")
    _quotes(results_tex, r"does not establish their equivalence")


def test_the_icbi_candidates_agree_on_the_pelka_rows(results_tex):
    ambiguity = pd.read_parquet(
        newest_by_time(RESULTS_DIR, "table_resolution_ambiguity")
    )
    row = ambiguity[ambiguity["table"] == "icbi_coexpression"].iloc[0]
    name_pick = RESULTS_DIR / row["name_order_picks"] / "icbi_coexpression.parquet"
    time_pick = RESULTS_DIR / row["time_order_picks"] / "icbi_coexpression.parquet"

    def pelka(path):
        frame = pd.read_parquet(path)
        frame = frame[frame["study_id"] == "Pelka_2021_Cell"]
        return frame.sort_values(["patient_id", "gene"]).reset_index(drop=True)

    pd.testing.assert_frame_equal(pelka(name_pick), pelka(time_pick))
    _quotes(results_tex, r"its Pelka rows are byte-identical between the")


def test_the_named_ambiguous_tables_are_the_tables_in_the_finding(results_tex):
    ambiguity = pd.read_parquet(
        newest_by_time(RESULTS_DIR, "table_resolution_ambiguity")
    )
    named = {
        "panel_coverage", "chen_subtype_falsifier_branch",
        "crowell_multisection_summary", "becker_feasibility", "mlh1_power",
        "interval_calibration",
    }
    assert named <= set(ambiguity["table"]), named - set(ambiguity["table"])
    # Both dispositions occur, so "neither rule is universally right" is not a
    # claim about one direction only.
    assert set(ambiguity["disposition"]) == {"identical", "differs"}
    _quotes(results_tex, r"Neither rule is universally right")


def test_every_challenge_case_saw_the_ledger(results_tex):
    path = newest_by_time(RESULTS_DIR, "blinded_guard_challenge")
    frame = pd.read_parquet(path)
    assert frame["saw_ledger"].all()
    assert len(frame) == 12
    assert int(frame["is_clean_control"].sum()) == 6
    _quotes(results_tex, r"case carries \texttt{saw\_ledger=true}")


def test_the_abstract_numbers_are_the_results_numbers():
    abstract = (REPO_ROOT / "paper" / "bmc" / "sections" / "abstract.tex").read_text()
    assert "7.0--8.6" in abstract
    assert "six of eight" in abstract and "five of eight" in abstract
    assert "23 table/date pairs" in abstract and "9 held" in abstract


def test_every_cited_key_is_declared_verified():
    """A citation cannot be inherited unverified from the WMHS bibliography.

    The BMC Background first reused WMHS citations without re-checking them.
    ``refs.bib`` now declares, on a ``BMC-VERIFIED:`` line, the entries that
    were re-checked against the publisher's record. This fails if the paper
    cites a key that is not declared there, or that is not in the bibliography.
    """
    bib = (REPO_ROOT / "paper" / "bmc" / "refs.bib").read_text()
    match = re.search(r"%\s*BMC-VERIFIED:\s*(.+)", bib)
    assert match, "refs.bib must carry a BMC-VERIFIED: line"
    declared = {key.strip() for key in match.group(1).split(",") if key.strip()}
    bibkeys = set(re.findall(r"@\w+\{([^,]+),", bib))

    sections = REPO_ROOT / "paper" / "bmc" / "sections"
    tex = "\n".join(p.read_text() for p in sections.glob("*.tex"))
    cited: set[str] = set()
    for m in re.finditer(r"\\cite\{([^}]*)\}", tex):
        cited |= {key.strip() for key in m.group(1).split(",") if key.strip()}

    assert cited <= bibkeys, f"cited but absent from refs.bib: {sorted(cited - bibkeys)}"
    assert cited <= declared, (
        f"cited but not declared verified: {sorted(cited - declared)}"
    )


# ---------------------------------------------------------------------------
# Implementation — the inventory counts and the design sizes
# ---------------------------------------------------------------------------


def test_the_inventory_counts_are_the_prose_counts():
    implementation = (
        REPO_ROOT / "paper" / "bmc" / "sections" / "implementation.tex"
    ).read_text()
    inventory = pd.read_parquet(newest_by_time(RESULTS_DIR, "claim_check_inventory"))
    assert len(inventory) == 24
    assert inventory["classification"].value_counts().to_dict() == {
        "implementation_error": 10,
        "logical_impossibility": 6,
        "provenance_reporting_error": 6,
        "low_power": 2,
    }
    assert inventory["verdict"].value_counts().to_dict() == {
        "scientific_failure": 23,
        "disputed": 1,
    }
    assert int(inventory["has_forcing_input"].sum()) == 15
    # A single rater; inter-rater agreement is not established.
    assert inventory["rater"].nunique() == 1

    _quotes(implementation, "The inventory carries 24 entries against the 21")
    _quotes(implementation, "10 are implementation errors, 6 logical")
    _quotes(implementation, "committed input that forces the check to fail")
    _quotes(implementation, "One worker classified them, and inter-rater agreement is not")


def test_the_stress_design_is_the_prose_design(stress_cells, results_tex):
    assert stress_cells["method"].nunique() == 3
    assert stress_cells["seed"].nunique() == 5
    assert set(stress_cells["n_trials"].unique()) == {800}
    assert len(stress_cells) == 630
    _quotes(
        results_tex,
        "three interval methods, five fixed seeds and 800 trials per cell give 630 cells",
    )


# ---------------------------------------------------------------------------
# Discussion, Conclusions and Declarations
# ---------------------------------------------------------------------------


def test_the_discussion_repeats_the_results_numbers():
    discussion = (
        REPO_ROOT / "paper" / "bmc" / "sections" / "discussion.tex"
    ).read_text()
    _quotes(discussion, "separates all eight cross-block contrasts from zero under every")
    _quotes(discussion, "six of eight")
    _quotes(discussion, "five of eight under tumour")
    # The same table as Results, so the Discussion cannot state a different one.
    summary = pd.read_parquet(
        newest_by_time(RESULTS_DIR, "adenoma_claim_sensitivity_summary")
    )
    lineage = summary[
        (summary["granularity_rung"] == "lineage")
        & (summary["statistic"] == "log_ratio")
    ]
    assert lineage["n_excluding_zero"].eq(8).all()


def _numbers(tex: str) -> set[str]:
    return set(re.findall(r"\d+(?:\.\d+)?", tex))


def test_the_abstract_introduces_no_number_absent_from_results(results_tex):
    abstract = (
        REPO_ROOT / "paper" / "bmc" / "sections" / "abstract.tex"
    ).read_text()
    extra = _numbers(abstract) - _numbers(results_tex)
    assert not extra, f"abstract has numbers not in Results: {sorted(extra)}"


def test_the_retired_phrase_is_nowhere_in_the_paper():
    """The n=10 error was in results.tex and the abstract.

    Fixing only the section the test looked at would have left the abstract
    wrong, which is how the phrase survived the first fix. Scan every section.
    """
    tex = "\n".join(
        p.read_text()
        for p in (REPO_ROOT / "paper" / "bmc" / "sections").glob("*.tex")
    )
    assert "smallest cohort" not in tex


def test_the_conclusions_introduce_no_number_absent_from_results(results_tex):
    conclusions = (
        REPO_ROOT / "paper" / "bmc" / "sections" / "conclusions.tex"
    ).read_text()
    extra = _numbers(conclusions) - _numbers(results_tex)
    assert not extra, f"Conclusions has numbers not in Results: {sorted(extra)}"


def test_the_prior_presentation_numbers_are_the_overlap_table():
    declarations = (
        REPO_ROOT / "paper" / "bmc" / "sections" / "declarations.tex"
    ).read_text()
    overlap = pd.read_parquet(newest_by_time(RESULTS_DIR, "workshop_overlap"))
    assert overlap["kind"].value_counts().to_dict() == {"overlap": 10, "math_claim": 4}
    statuses = (
        overlap[overlap["kind"] == "overlap"]["status"].value_counts().to_dict()
    )
    # The earlier record said "3 extended"; the table has 2 extended and 1
    # contradicted (O03 carries the retired "generator noise" claim).
    assert statuses == {
        "already_published": 5, "extended": 2, "new": 2, "contradicted": 1,
    }
    novelty = statuses["new"] / sum(statuses.values())
    assert round(novelty, 2) == 0.20

    _quotes(declarations, "5 already published, 2")
    _quotes(declarations, "novelty share of $0.20$")
    # The submission gate is disclosed in the Declarations too.
    _quotes(declarations, r"present in the index as")
    _quotes(declarations, r"\textsc{blocked}, not as done")
