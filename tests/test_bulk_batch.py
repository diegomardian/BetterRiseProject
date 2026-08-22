"""W3.4 — batch structure. The tests that matter are the null controls.

A categorical factor with many levels explains variance by construction, and a
chi-squared on a sparse table finds association in noise. So the load-bearing
tests here are the ones asserting that **unrelated variables come back
unrelated** — `test_random_factor_explains_nothing_above_its_null` and
`test_independent_variables_are_not_flagged_as_confounded`. Without those, a
clean W3.4 result would be indistinguishable from a broken one.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.bulk.batch import (
    BatchError,
    confounding_table,
    cramers_v,
    permutation_p,
    pool_rare_levels,
    variance_explained,
)
from src.bulk.clinical import (
    COLORECTAL_SITES,
    build_clinical_table,
    coverage_report,
    harmonise_stage,
    msi_by_patient,
    select_colorectal_diagnosis,
)

SEED = 20260818
RNG = np.random.default_rng(SEED)


# ---------------------------------------------------------------------------
# Clinical harmonisation
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        ("Stage I", "I"),
        ("Stage IA", "I"),
        ("Stage IIA", "II"),
        ("Stage IIB", "II"),
        ("Stage IIC", "II"),
        ("Stage III", "III"),
        ("Stage IIIA", "III"),
        ("Stage IIIC", "III"),
        ("Stage IV", "IV"),
        ("Stage IVA", "IV"),
        ("Stage IVB", "IV"),
        ("  stage iiib  ", "III"),
    ],
)
def test_stage_harmonisation_covers_every_value_tcga_actually_uses(raw, expected):
    """These twelve are the values present in COAD/READ, checked against the
    live API rather than assumed."""
    assert harmonise_stage(raw) == expected


@pytest.mark.parametrize("raw", [None, "", "Not Reported", "Unknown", "Stage X", "IIIB", 42])
def test_unstageable_values_become_none_not_a_string(raw):
    """None, never "Unknown" — a string would sort into the middle of an
    ordered factor and be silently modelled as a stage between II and III."""
    assert harmonise_stage(raw) is None


def test_the_colorectal_diagnosis_is_chosen_not_assumed():
    """182 of 633 COAD/READ cases carry more than one diagnosis, and the extras
    include prostate and breast. diagnoses[0] would assign some patients another
    organ's site and stage."""
    diagnoses = [
        {"tissue_or_organ_of_origin": "Breast, NOS", "ajcc_pathologic_stage": "Stage I"},
        {"tissue_or_organ_of_origin": "Sigmoid colon", "ajcc_pathologic_stage": "Stage IIIB"},
    ]
    chosen = select_colorectal_diagnosis(diagnoses)
    assert chosen["tissue_or_organ_of_origin"] == "Sigmoid colon"
    assert harmonise_stage(chosen["ajcc_pathologic_stage"]) == "III"


def test_a_staged_colorectal_diagnosis_wins_over_an_unstaged_one():
    diagnoses = [
        {"tissue_or_organ_of_origin": "Cecum", "ajcc_pathologic_stage": None},
        {"tissue_or_organ_of_origin": "Rectum, NOS", "ajcc_pathologic_stage": "Stage IV"},
    ]
    assert select_colorectal_diagnosis(diagnoses)["ajcc_pathologic_stage"] == "Stage IV"


def test_a_case_with_no_colorectal_diagnosis_returns_none():
    diagnoses = [
        {"tissue_or_organ_of_origin": "Prostate gland", "ajcc_pathologic_stage": "Stage I"}
    ]
    assert select_colorectal_diagnosis(diagnoses) is None


def test_site_groups_split_at_the_splenic_flexure():
    """Transverse colon counts as right. The other convention exists and would
    move ~21 patients, so it is pinned."""
    assert COLORECTAL_SITES["Transverse colon"] == "right_colon"
    assert COLORECTAL_SITES["Splenic flexure of colon"] == "left_colon"
    assert COLORECTAL_SITES["Rectosigmoid junction"] == "rectum"
    assert COLORECTAL_SITES["Colon, NOS"] == "colon_unspecified"


def test_conflicting_msi_calls_are_flagged_never_resolved():
    """17 real cases carry both an MSI and an MSS call across their WGS and WXS
    files. MSI is the project's one pre-registered subgroup variable; picking a
    winner would fabricate a subgroup label."""
    hits = [
        {"msi_status": "MSS", "cases": [{"submitter_id": "P1"}]},
        {"msi_status": "MSI", "cases": [{"submitter_id": "P1"}]},
        {"msi_status": "MSS", "cases": [{"submitter_id": "P2"}]},
        {"msi_status": "MSS", "cases": [{"submitter_id": "P2"}]},
    ]
    out = msi_by_patient(hits).set_index("patient_id")
    assert out.loc["P1", "msi_status"] == "conflicting"
    assert out.loc["P1", "n_distinct_calls"] == 2
    assert out.loc["P2", "msi_status"] == "MSS"  # agreeing duplicates are fine
    assert out.loc["P2", "n_distinct_calls"] == 1


def test_clinical_table_and_coverage_report():
    cases = [
        {
            "submitter_id": "TCGA-AA-0001",
            "project": {"project_id": "TCGA-COAD"},
            "demographic": {"sex_at_birth": "Female", "age_at_index": 61},
            "diagnoses": [
                {"tissue_or_organ_of_origin": "Cecum", "ajcc_pathologic_stage": "Stage IIA"}
            ],
        },
        {
            "submitter_id": "TCGA-AA-0002",
            "project": {"project_id": "TCGA-READ"},
            "demographic": {"sex_at_birth": "Male", "age_at_index": 70},
            "diagnoses": [{"tissue_or_organ_of_origin": "Lung, NOS"}],
        },
    ]
    msi = [{"msi_status": "MSI", "cases": [{"submitter_id": "TCGA-AA-0001"}]}]
    table = build_clinical_table(cases, msi).set_index("patient_id")

    assert table.loc["TCGA-AA-0001", "site"] == "right_colon"
    assert table.loc["TCGA-AA-0001", "stage"] == "II"
    assert table.loc["TCGA-AA-0001", "msi_status"] == "MSI"
    # No colorectal diagnosis -> no site, no stage, and no other organ's stage.
    assert pd.isna(table.loc["TCGA-AA-0002", "site"])
    assert pd.isna(table.loc["TCGA-AA-0002", "stage"])

    cov = coverage_report(table.reset_index(), ["TCGA-AA-0001", "TCGA-AA-0002"]).set_index(
        "variable"
    )
    assert cov.loc["stage", "n_annotated"] == 1
    assert cov.loc["stage", "coverage"] == 0.5
    assert cov.loc["msi_status", "n_annotated"] == 1


# ---------------------------------------------------------------------------
# Cramer's V
# ---------------------------------------------------------------------------


def test_cramers_v_is_one_for_a_perfect_association():
    a = pd.Series(list("aabbcc") * 20)
    assert cramers_v(a, a, bias_correct=False) == pytest.approx(1.0)


def test_bias_correction_pulls_a_sparse_table_down():
    """Uncorrected V on a sparse table is badly upward-biased. On independent
    data the corrected version should be near zero and the raw one clearly
    above it — that gap is why the correction is not optional."""
    a = pd.Series([f"plate{i % 40}" for i in range(200)])
    b = pd.Series(RNG.choice(["I", "II", "III", "IV"], 200))
    raw = cramers_v(a, b, bias_correct=False)
    corrected = cramers_v(a, b, bias_correct=True)
    assert raw > corrected
    assert corrected < 0.15


def test_cramers_v_is_nan_when_a_variable_is_constant():
    a = pd.Series(["x"] * 50)
    b = pd.Series(RNG.choice(["I", "II"], 50))
    assert np.isnan(cramers_v(a, b))


def test_cramers_v_uses_complete_pairs_only():
    a = pd.Series(["x", "x", "y", "y", None])
    b = pd.Series(["I", "II", "I", "II", "I"])
    assert np.isfinite(cramers_v(a, b))


# ---------------------------------------------------------------------------
# Null controls — the tests that make a clean result meaningful
# ---------------------------------------------------------------------------


def test_independent_variables_are_not_flagged_as_confounded():
    """A high-cardinality technical factor against an unrelated clinical one
    must come back non-significant. If this fails, every W3.4 'confounded'
    finding is an artifact of the level count."""
    n = 400
    annotations = pd.DataFrame(
        {
            "plate": [f"p{i % 30}" for i in range(n)],
            "stage": RNG.choice(["I", "II", "III", "IV"], n),
        }
    )
    observed, null_mean, p = permutation_p(
        annotations["plate"], annotations["stage"], n_permutations=199, seed=SEED
    )
    assert p > 0.05
    assert abs(observed - null_mean) < 0.05


def test_a_real_association_is_detected():
    """The converse. Build stage deterministically from plate and the test must
    fire, or it cannot detect the confounding it exists to find."""
    n = 400
    plate = pd.Series([f"p{i % 4}" for i in range(n)])
    stage = plate.map({"p0": "I", "p1": "II", "p2": "III", "p3": "IV"})
    observed, null_mean, p = permutation_p(plate, stage, n_permutations=199, seed=SEED)
    assert observed == pytest.approx(1.0, abs=0.05)
    assert p < 0.01
    assert observed - null_mean > 0.8


def test_confounding_table_covers_every_pair_and_records_level_counts():
    n = 200
    annotations = pd.DataFrame(
        {
            "tss": [f"t{i % 8}" for i in range(n)],
            "plate": [f"p{i % 12}" for i in range(n)],
            "stage": RNG.choice(["I", "II"], n),
            "msi_status": RNG.choice(["MSS", "MSI"], n),
        }
    )
    out = confounding_table(
        annotations, ("tss", "plate"), ("stage", "msi_status"), seed=SEED, n_permutations=99
    )
    assert len(out) == 4
    assert set(out["technical_factor"]) == {"tss", "plate"}
    assert out.loc[out["technical_factor"] == "plate", "n_levels_technical"].iloc[0] == 12
    assert out["permutation_p"].notna().all()


def test_confounding_table_raises_on_a_missing_column():
    with pytest.raises(BatchError, match="missing column"):
        confounding_table(pd.DataFrame({"a": ["x"]}), ("a",), ("nope",), seed=SEED)


# ---------------------------------------------------------------------------
# Variance explained
# ---------------------------------------------------------------------------


def _expression(n_samples=120, n_genes=300, seed=0):
    rng = np.random.default_rng(seed)
    return pd.DataFrame(
        rng.normal(0, 1, (n_samples, n_genes)),
        index=[f"S{i}" for i in range(n_samples)],
    )


def test_random_factor_explains_nothing_above_its_null():
    """THE control. A 20-level factor on random data will show a large raw R2
    purely from its level count. excess_over_null must be ~0, and that is the
    column the write-up quotes."""
    expr = _expression()
    factors = pd.DataFrame(
        {"noise": [f"g{i % 20}" for i in range(len(expr))]}, index=expr.index
    )
    out, meta = variance_explained(expr, factors, seed=SEED, n_permutations=40)
    row = out.iloc[0]
    assert row["r2"] > 0.05, "fixture no longer demonstrates level-count inflation"
    assert abs(row["excess_over_null"]) < 0.05
    assert row["permutation_p"] > 0.05
    assert meta["n_pcs_retained"] >= 2


def test_a_real_batch_effect_is_detected():
    """Shift one group's expression and the factor must rise clearly above its
    null, or the analysis cannot find a batch effect that is there."""
    expr = _expression()
    group = np.array([f"b{i % 2}" for i in range(len(expr))])
    expr.loc[group == "b0"] += 3.0
    factors = pd.DataFrame({"batch": group}, index=expr.index)
    out, _ = variance_explained(expr, factors, seed=SEED, n_permutations=40)
    row = out.iloc[0]
    assert row["excess_over_null"] > 0.3
    assert row["permutation_p"] < 0.05


def test_rare_levels_are_pooled():
    values = pd.Series(["a"] * 20 + ["b"] * 20 + ["c"] * 2 + ["d"])
    pooled = pool_rare_levels(values, min_size=5)
    assert set(pooled.unique()) == {"a", "b", "other"}
    assert (pooled == "other").sum() == 3


def test_variance_explained_reports_pooling_and_pc_metadata():
    expr = _expression()
    factors = pd.DataFrame(
        {"f": [f"g{i}" if i < 5 else "common" for i in range(len(expr))]}, index=expr.index
    )
    out, meta = variance_explained(expr, factors, seed=SEED, n_permutations=20)
    assert out.loc[0, "n_levels_before_pooling"] > out.loc[0, "n_levels"]
    assert 0 < meta["variance_covered"] <= 1.0
    assert meta["n_samples"] == len(expr)


def test_variance_explained_refuses_a_mismatched_index():
    expr = _expression(n_samples=20)
    factors = pd.DataFrame({"f": ["a"] * 5}, index=[f"OTHER{i}" for i in range(5)])
    with pytest.raises(BatchError, match="shared"):
        variance_explained(expr, factors, seed=SEED, n_permutations=5)


def test_variance_explained_refuses_an_empty_matrix():
    with pytest.raises(BatchError, match="empty"):
        variance_explained(pd.DataFrame(), pd.DataFrame(), seed=SEED)
