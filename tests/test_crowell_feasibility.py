"""The Crowell WTx feasibility read, against a synthetic section.

Pre-registration: ``docs/prereg_crowell_feasibility.md`` (`cd7e4a5`).

**The cluster run must not be this code's first execution.** That is the
repository's own rule and it is why the fixture below is built in the shape the
deposit is documented to have — real genes plus negative probes plus false
codes, cells carrying a histological domain — rather than in a shape that
happens to be convenient.

The h5ad round trip is skipped where ``anndata`` is absent (it is not in the
laptop env; ``env/w1_reference.yml`` pins it on the cluster). Everything that
decides a verdict is tested without it.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.reference.crowell_io import (
    CrowellError,
    check_counts_are_integers,
    control_features,
    domain_vocabulary,
    qc_pass_mask,
    require_controls,
    require_counts,
)
from src.reference.jobs.crowell_feasibility import (
    CRITICAL_GENE,
    MIN_CELLS_PER_DOMAIN,
    MIN_LOG_SEPARATION,
    find_panel,
    negative_floor,
    per_domain_table,
    verdict,
)

PANEL = ("ACTB", "KRT8", "EPCAM", "CDX2", "MS4A12", "GUCA2A")
N_NEG = 50


def _var_names(n_filler: int = 20) -> list[str]:
    """Real genes, 50 negative probes, false codes — the documented shape."""
    return (list(PANEL)
            + [f"FILLER{i}" for i in range(n_filler)]
            + [f"NegPrb{i}" for i in range(N_NEG)]
            + [f"FalseCode{i}" for i in range(12)])


def _section(rates: dict[str, dict[str, float]], neg_rate: float = 0.01,
             n_per_domain: int = 1_000, seed: int = 0):
    """A synthetic section: per-domain per-gene detection at known rates."""
    rng = np.random.default_rng(seed)
    names = _var_names()
    col = {g: i for i, g in enumerate(names)}
    domains = sorted(rates)
    obs_domain, blocks = [], []
    for domain in domains:
        block = np.zeros((n_per_domain, len(names)), dtype=float)
        for gene, rate in rates[domain].items():
            block[:, col[gene]] = (rng.random(n_per_domain) < rate).astype(float)
        for j, name in enumerate(names):
            if name.startswith(("NegPrb", "FalseCode")):
                block[:, j] = (rng.random(n_per_domain) < neg_rate).astype(float)
            elif name.startswith("FILLER"):
                block[:, j] = (rng.random(n_per_domain) < 0.2).astype(float)
        blocks.append(block)
        obs_domain += [domain] * n_per_domain
    return np.vstack(blocks), names, np.array(obs_domain)


# ---------------------------------------------------------------------------
# The floor, which is the part most easily got wrong
# ---------------------------------------------------------------------------


def test_the_floor_is_per_probe_not_the_union_over_fifty():
    """THE ERROR THIS FUNCTION EXISTS TO AVOID.

    "Fraction of cells with >=1 count on ANY of the 50 negative probes" is a
    union over 50 features, and a gene's detection is one feature. Comparing
    them overstates the floor by roughly the probe count and would fail every
    gene. Both are returned, under names that say which is which.
    """
    matrix, names, _ = _section({"REF": {"ACTB": 0.5}}, neg_rate=0.02, seed=1)
    controls = control_features(names)
    floor = negative_floor(matrix, controls["negative_indices"])

    assert floor["n_control_probes"] == N_NEG
    assert floor["floor_per_probe_mean"] == pytest.approx(0.02, abs=0.01)
    # the union over 50 probes at 2% each is nearly everything
    assert floor["any_probe_union_rate"] > 0.5
    assert floor["any_probe_union_rate"] > 10 * floor["floor_per_probe_mean"]


def test_an_empty_control_match_is_refused_rather_than_giving_a_floor_of_zero():
    """A floor built from nothing is zero, and a zero floor is a gate every
    gene clears — a check that cannot fail, which is this repo's signature
    defect."""
    controls = control_features(["ACTB", "GUCA2A", "SOMETHING_ELSE"])
    assert controls["n_negative"] == 0
    with pytest.raises(CrowellError, match="floor of zero"):
        require_controls(controls)


def test_control_features_reports_which_naming_convention_matched():
    """CosMx's prefix varies by release; an unreported match is unauditable."""
    hit = control_features(["ACTB", "NegPrb1", "NegPrb2", "FalseCode1"])
    assert hit["negative_patterns_matched"] == [r"^NegPrb"]
    assert hit["false_code_patterns_matched"] == [r"^FalseCode"]
    assert hit["n_negative"] == 2 and hit["n_false_code"] == 1


# ---------------------------------------------------------------------------
# The table and the gate
# ---------------------------------------------------------------------------


def test_a_gene_far_above_the_floor_is_usable_and_one_at_the_floor_is_not():
    matrix, names, domains = _section(
        {"REF": {"ACTB": 0.60, "KRT8": 0.50, "EPCAM": 0.40, "CDX2": 0.20,
                 "MS4A12": 0.15, CRITICAL_GENE: 0.012}},
        neg_rate=0.01, seed=2)
    panel_index, absent = find_panel(names)
    assert not absent
    table = per_domain_table(matrix, panel_index, domains,
                             control_features(names)["negative_indices"])
    row = table.set_index("gene")
    assert bool(row.loc["ACTB", "usable"])
    # GUCA2A at 1.2% against a 1% floor is inside the noise
    assert not bool(row.loc[CRITICAL_GENE, "usable"])
    assert row.loc[CRITICAL_GENE, "log_separation"] < MIN_LOG_SEPARATION


def test_a_domain_below_the_cell_floor_is_not_scored():
    matrix, names, domains = _section({"REF": {"ACTB": 0.5}}, seed=3,
                                      n_per_domain=MIN_CELLS_PER_DOMAIN - 1)
    table = per_domain_table(matrix, find_panel(names)[0], domains,
                             control_features(names)["negative_indices"])
    assert table.empty, "a detection rate over a handful of cells is not one"


def test_missing_histopathology_is_not_stringified_into_a_nan_domain():
    """Missing `typ` values are unassigned cells, not a fourth tissue domain."""
    matrix, names, domains = _section(
        {"REF": {"ACTB": 0.5}, "TVA": {"ACTB": 0.5}}, seed=31)
    domains = domains.astype(object)
    domains[:MIN_CELLS_PER_DOMAIN] = np.nan
    table = per_domain_table(matrix, find_panel(names)[0], domains,
                             control_features(names)["negative_indices"])
    assert set(table["domain"]) == {"REF", "TVA"}
    assert "nan" not in set(table["domain"])


def test_depth_is_reported_per_domain():
    """The confound that already fooled this project once: Becker's mature
    label carried 2.04x the arm's median UMIs and lifted every gene."""
    matrix, names, domains = _section(
        {"REF": {g: 0.5 for g in PANEL}, "TVA": {g: 0.1 for g in PANEL}}, seed=4)
    table = per_domain_table(matrix, find_panel(names)[0], domains,
                             control_features(names)["negative_indices"])
    depths = table.drop_duplicates("domain").set_index("domain")
    assert depths.loc["REF", "median_counts_per_cell"] > \
        depths.loc["TVA", "median_counts_per_cell"]
    assert table["median_genes_per_cell"].notna().all()


def test_separation_is_on_the_detection_scale_not_a_ratio_of_probabilities():
    """cloglog(p) = log(mu). A ratio of probabilities is not comparable across
    genes with different baselines, and this repo has reintroduced that error
    twice."""
    matrix, names, domains = _section(
        {"REF": {"ACTB": 0.40, CRITICAL_GENE: 0.10}}, neg_rate=0.01, seed=5)
    table = per_domain_table(matrix, find_panel(names)[0], domains,
                             control_features(names)["negative_indices"])
    row = table.set_index("gene")
    def mu(p):
        return -np.log1p(-p)
    floor = row.loc["ACTB", "floor_per_probe_mean"]
    assert row.loc["ACTB", "log_separation"] == pytest.approx(
        np.log(mu(row.loc["ACTB", "detection"]) / mu(floor)), abs=1e-9)


# ---------------------------------------------------------------------------
# The verdict's branches, one input each
# ---------------------------------------------------------------------------


def _table(**seps) -> pd.DataFrame:
    from src.reference.jobs.coexpression_silencing import GENE_ROLES
    rows = []
    for domain, per_gene in seps.items():
        for gene, sep in per_gene.items():
            rows.append({"domain": domain, "gene": gene, "role": GENE_ROLES[gene],
                         "log_separation": sep,
                         "usable": sep >= MIN_LOG_SEPARATION})
    return pd.DataFrame(rows)


def test_controls_failing_everywhere_refuses_the_read_rather_than_blaming_the_target():
    """Prereg §6, fourth branch."""
    table = _table(REF={"ACTB": 0.01, "KRT8": 0.01, CRITICAL_GENE: 0.01})
    out = verdict(table)
    assert out["verdict"].startswith("READ REFUSED")
    assert "not a statement about the target" in out["detail"]


def test_a_target_below_the_floor_everywhere_is_a_sensitivity_statement():
    table = _table(REF={"ACTB": 3.0, "KRT8": 2.5, CRITICAL_GENE: 0.05})
    out = verdict(table)
    assert out["verdict"] == "NOT ASKABLE HERE"
    assert "in-situ sensitivity" in out["detail"]
    assert "NOT about the biology" in out["detail"]


def test_the_adenoma_domain_is_not_guessed():
    """The deposit's words are a measurement. Reading a label rather than the
    grouping once put Chen_2021's usable pairs at zero when it was 44."""
    table = _table(REF={"ACTB": 3.0, CRITICAL_GENE: 2.0},
                   TVA={"ACTB": 3.0, CRITICAL_GENE: 2.0})
    out = verdict(table, adenoma_domain=None)
    assert out["verdict"] == "DOMAIN NAMING NOT SUPPLIED"
    assert "is not guessed" in out["detail"]


def test_the_target_below_the_bar_keeps_both_halves_of_the_answer():
    """Amendment 2. §6 has no branch for "separates in REF, below the bar in
    TVA", and the first implementation collapsed it to NOT ASKABLE — true of
    the per-cell reading and silent about a direction §6 pre-specified.

    A single label cannot carry both, so the components are fields.
    """
    table = _table(REF={"KRT8": 3.0, CRITICAL_GENE: 2.0},
                   TVA={"KRT8": 3.6, CRITICAL_GENE: 0.02})
    out = verdict(table, adenoma_domain="TVA", reference_domain="REF",
                  patient_n=1)

    assert out["verdict"].startswith("TARGET BELOW THE ADENOMA USABILITY BAR")
    assert out["per_cell_feasibility"] == "failed"
    assert out["prespecified_directional_read"] == "fall_observed"
    assert out["global_sensitivity_control"] == "not_worse"
    assert out["patient_n"] == 1
    # it must never claim the target became measurable per cell
    assert "NOT licensed" in out["detail"]
    # below the LOQ, not censored -- the bar is analyst-chosen
    assert "NOT censored" in out["detail"]
    assert "limit of quantification" in out["detail"]
    # global capture must not be read as gene-specific sensitivity
    assert "does not establish" in out["detail"]
    assert "Supportive, NOT confirmatory" in out["detail"]


def test_a_target_below_the_bar_never_returns_askable():
    """No relabelling makes an unmeasured gene measurable."""
    table = _table(REF={"KRT8": 3.0, CRITICAL_GENE: 2.0},
                   TVA={"KRT8": 3.6, CRITICAL_GENE: 0.02})
    for ref in (None, "REF"):
        out = verdict(table, adenoma_domain="TVA", reference_domain=ref,
                      patient_n=1)
        assert "ASKABLE IN THE ADENOMA" != out["verdict"]
        assert out["per_cell_feasibility"] == "failed"


def test_no_fall_is_recorded_as_uninterpretable_not_as_a_negative():
    """§6's no-fall branch lives where BOTH domains clear the bar.

    It was unreachable in the below-bar branch: a target above the bar in the
    reference and below it in the adenoma has fallen by construction. An
    unreachable branch is a check that cannot fail, which is why the directional
    read is now computed once and attached to both outcomes.
    """
    table = _table(REF={"KRT8": 3.0, CRITICAL_GENE: 1.50},
                   TVA={"KRT8": 3.6, CRITICAL_GENE: 2.00})
    out = verdict(table, adenoma_domain="TVA", reference_domain="REF",
                  patient_n=1)
    assert out["per_cell_feasibility"] == "passed"
    assert out["prespecified_directional_read"] == "no_fall"
    assert "UNINTERPRETABLE" in out["detail"]
    assert "may not be quoted as a negative" in out["detail"]


def test_a_pass_says_it_licenses_the_question_and_not_an_answer():
    """The pass must carry §6's asymmetry and §7's limit, or it will be quoted
    without them."""
    table = _table(REF={"ACTB": 3.0, CRITICAL_GENE: 2.0},
                   TVA={"ACTB": 3.0, CRITICAL_GENE: 1.8})
    out = verdict(table, adenoma_domain="TVA", reference_domain="REF",
                  patient_n=1)
    assert out["verdict"] == "ASKABLE IN THE ADENOMA"
    assert out["per_cell_feasibility"] == "passed"
    assert "QUESTION, not an answer" in out["detail"]
    assert "silenced" in out["detail"]
    # §6's direction is carried on the pass too, not only on the failure
    assert out["prespecified_directional_read"] in {"fall_observed", "no_fall"}


# ---------------------------------------------------------------------------
# The two things the reader refuses to assume
# ---------------------------------------------------------------------------


def test_normalised_values_are_refused_before_a_detection_rate_is_taken():
    """A detection rate off log-normalised data is not a detection rate, and it
    returns a plausible number rather than raising."""
    normalised = np.array([[0.0, 1.386, 2.079], [0.693, 0.0, 1.099]])
    report = check_counts_are_integers(normalised)
    assert not report["integral"]
    with pytest.raises(CrowellError, match="not integer counts"):
        require_counts(report)

    counts = np.array([[0.0, 4.0, 8.0], [2.0, 0.0, 3.0]])
    require_counts(check_counts_are_integers(counts))


def test_the_domain_vocabulary_is_selected_by_cardinality_not_by_name():
    """THE DEFECT THE FIRST CROWELL INSPECTION SHIPPED.

    The first version filtered obs against a hardcoded list of readable English
    names. The deposit names its annotations `typ`, `roi`, `ctx`, `lv1`, `lv2`,
    so the inspection reported ONE candidate column and silently hid the rest.
    A module whose docstring says "assume nothing, map nothing" cannot select by
    guessing what a column will be called.
    """
    obs = pd.DataFrame({
        "region": ["REF", "REF", "TVA", "CRC"],
        "typ": ["a", "a", "b", "b"],      # cryptic, and real: Crowell uses it
        "lv2": ["x", "y", "x", "y"],
        "cell_id": [f"c{i}" for i in range(4)],   # high cardinality, not a vocabulary
    })
    vocab = domain_vocabulary(obs, max_distinct=3)
    assert vocab["candidate_columns"]["region"] == {"REF": 2, "TVA": 1, "CRC": 1}
    assert "typ" in vocab["candidate_columns"], "a cryptic name is still a vocabulary"
    assert "lv2" in vocab["candidate_columns"]
    assert "cell_id" not in vocab["candidate_columns"], "4 distinct over 4 rows"
    # the name list survives only to mark what was EXPECTED, never to filter
    assert vocab["expected_by_name"] == ["region"]
    assert "region" in vocab["expected_words_seen"]


def test_a_constant_named_column_is_reported_because_that_is_the_finding():
    """`tissue` was 'TVA' for all 130,814 cells of section 232. A constant is
    not a vocabulary, but "this whole section is one domain" is the answer."""
    obs = pd.DataFrame({"tissue": ["TVA"] * 10})
    vocab = domain_vocabulary(obs)
    assert vocab["candidate_columns"]["tissue"] == {"TVA": 10}


def test_the_deposit_qc_flag_excludes_failed_cells():
    """The first real run used all 298,151 rows even though Crowell defines
    `fil` as the cell-level QC-pass flag. The rejected population must not be
    silently mixed back into the assay."""
    obs = pd.DataFrame({"fil": [True, False, True, False]})
    assert qc_pass_mask(obs).tolist() == [True, False, True, False]


def test_the_qc_flag_may_round_trip_as_logical_strings():
    obs = pd.DataFrame({"fil": ["True", "False", "true"]})
    assert qc_pass_mask(obs).tolist() == [True, False, True]


def test_a_missing_or_nonlogical_qc_flag_is_refused():
    with pytest.raises(CrowellError, match="unfiltered population"):
        qc_pass_mask(pd.DataFrame({"something_else": [True, False]}))
    with pytest.raises(CrowellError, match="not logical"):
        qc_pass_mask(pd.DataFrame({"fil": ["keep", "drop"]}))
    with pytest.raises(CrowellError, match="missing values"):
        qc_pass_mask(pd.DataFrame({"fil": [True, None]}))


# ---------------------------------------------------------------------------
# The h5ad round trip, where anndata exists
# ---------------------------------------------------------------------------


def test_a_synthetic_section_round_trips_through_the_reader(tmp_path):
    anndata = pytest.importorskip("anndata")
    from scipy.sparse import csr_matrix

    from src.reference.crowell_io import open_section

    matrix, names, domains = _section(
        {"REF": {g: 0.4 for g in PANEL}, "TVA": {g: 0.2 for g in PANEL}}, seed=6)
    adata = anndata.AnnData(
        X=csr_matrix(matrix),
        obs=pd.DataFrame({"region": domains},
                         index=[f"c{i}" for i in range(matrix.shape[0])]),
        var=pd.DataFrame(index=names))
    path = tmp_path / "232.h5ad"
    adata.write_h5ad(path)

    reread = open_section(path, backed=False)
    assert reread.n_obs == matrix.shape[0]
    panel_index, absent = find_panel(reread.var_names)
    assert not absent
    controls = control_features([str(v) for v in reread.var_names])
    require_controls(controls)
    table = per_domain_table(reread.X, panel_index,
                             reread.obs["region"].to_numpy(),
                             controls["negative_indices"])
    assert set(table["domain"]) == {"REF", "TVA"}
    assert table["usable"].any()


def test_a_missing_section_names_the_file_and_its_checksum():
    from src.reference.crowell_io import open_section

    with pytest.raises(CrowellError, match="f2271485bae235439267c1bd13f1a7c9"):
        open_section("/nonexistent/232.h5ad")


def test_a_gene_detected_in_no_cell_is_not_separated_by_zero():
    """-inf, not 0. A gene absent from every cell is not "separated by zero";
    clamping it to a finite number would read like a measurement."""
    matrix, names, domains = _section({"REF": {"ACTB": 0.5}}, neg_rate=0.01,
                                      seed=7)
    table = per_domain_table(matrix, find_panel(names)[0], domains,
                             control_features(names)["negative_indices"])
    row = table.set_index("gene")
    assert row.loc[CRITICAL_GENE, "detection"] == 0.0
    assert row.loc[CRITICAL_GENE, "log_separation"] == float("-inf")
    assert not bool(row.loc[CRITICAL_GENE, "usable"])


# ---------------------------------------------------------------------------
# The floor when the probes are not features — which is what 232 actually is
# ---------------------------------------------------------------------------


def test_the_floor_can_be_built_from_obs_when_probes_are_not_features():
    """Section 232 carries 0 control probes in var: they were summarised per
    cell into obs before the object was written. The floor is still measurable,
    and the same per-probe quantity."""
    from src.reference.crowell_io import floor_from_obs

    obs = pd.DataFrame({
        "nFeature_negprobes": [0, 1, 2, 5, 0, 1],
        "nCount_negprobes": [0, 1, 3, 7, 0, 2],
        "nFeature_falsecode": [0, 0, 1, 2, 0, 0],
    })
    floor = floor_from_obs(obs, n_negative_probes=50)
    assert floor["n_control_probes"] == 50
    assert floor["probe_count_source"] == "supplied"
    assert floor["floor_per_probe_mean"] == pytest.approx(np.mean([0,1,2,5,0,1]) / 50)
    assert floor["any_probe_union_rate"] == pytest.approx(4 / 6)
    assert floor["false_code_union_rate"] == pytest.approx(2 / 6)


def test_an_inferred_probe_count_errs_toward_refusing():
    """Inferring the count from the observed max divides by too little, which
    raises the floor and shrinks the separation. A gene clearing a conservative
    floor clears a real one."""
    from src.reference.crowell_io import floor_from_obs

    obs = pd.DataFrame({"nFeature_negprobes": [0, 1, 2, 5],
                        "nCount_negprobes": [0, 1, 3, 7]})
    inferred = floor_from_obs(obs)
    supplied = floor_from_obs(obs, n_negative_probes=50)
    assert inferred["n_control_probes"] == 5
    assert inferred["probe_count_is_conservative"] is True
    assert inferred["floor_per_probe_mean"] > supplied["floor_per_probe_mean"]


def test_obs_without_the_negative_summaries_is_refused_not_defaulted():
    from src.reference.crowell_io import floor_from_obs

    with pytest.raises(CrowellError, match="floor of zero"):
        floor_from_obs(pd.DataFrame({"something_else": [1, 2, 3]}))


def test_per_domain_table_refuses_when_neither_floor_source_is_given():
    matrix, names, domains = _section({"REF": {"ACTB": 0.5}}, seed=8)
    with pytest.raises(CrowellError, match="floor of zero"):
        per_domain_table(matrix, find_panel(names)[0], domains)
