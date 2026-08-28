"""The Stage 4 pre-specification binds, and says what would falsify it.

W3.6's lesson: a pre-specification that does not stop you running is a document.
These tests pin the parts that make this one a lock rather than a note — and the
parts of the prediction that make it a prediction rather than a hope.
"""

from __future__ import annotations

import textwrap

import pytest
import yaml

from src.bulk.prespec import (
    PRESPEC_PATH,
    PrespecError,
    load_prespec,
    outcome_genes,
    positive_control_gates_the_analysis,
    prediction,
    require_locked_prespec,
)


@pytest.fixture(scope="module")
def spec():
    return load_prespec()


# ---------------------------------------------------------------------------
# The lock
# ---------------------------------------------------------------------------


def test_it_ships_proposed_so_nothing_can_run_against_it_yet(spec):
    """It is committed BEFORE the team confirms, which is the whole point. The
    flip to `locked` is its own commit with a stated reason."""
    assert spec["status"] == "proposed"


def test_running_against_a_proposed_spec_is_refused(spec):
    with pytest.raises(PrespecError, match="not 'locked'"):
        require_locked_prespec(spec)


def test_a_locked_spec_is_allowed(spec):
    locked = dict(spec, status="locked")
    require_locked_prespec(locked)  # does not raise


def test_the_lock_is_literally_the_covariate_sets_lock():
    """Two lock functions that could drift apart is worse than one used twice,
    so this asserts identity rather than similar behaviour."""
    from src.bulk import covariates, prespec

    assert prespec.require_locked is covariates.require_locked


# ---------------------------------------------------------------------------
# The prediction has to be able to be wrong
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("arm", ["primary", "secondary"])
def test_every_prediction_arm_says_what_would_disconfirm_it(spec, arm):
    """A prediction with no disconfirming condition is unfalsifiable, and this
    analysis expects a null — the easiest kind to report with the sign
    reversed."""
    text = prediction(spec, arm)["disconfirmed_if"].strip()
    assert text, f"{arm} has an empty disconfirmed_if"


def test_the_primary_prediction_is_directional_not_thresholded(spec):
    """Direction only, so the result does not hinge on an invented number."""
    statement = prediction(spec, "primary")["statement"]
    assert "<" in statement
    assert "GUCA2A" in statement and "CDX2" in statement


def test_the_secondary_prediction_commits_to_a_number(spec):
    """So the primary cannot be satisfied trivially by both being near zero."""
    assert "0.25" in prediction(spec, "secondary")["statement"]
    assert "0.50" in prediction(spec, "secondary")["disconfirmed_if"]


def test_an_arm_with_no_disconfirming_condition_fails_to_load(tmp_path):
    raw = yaml.safe_load(PRESPEC_PATH.read_text(encoding="utf-8"))
    del raw["prediction"]["primary"]["disconfirmed_if"]
    path = tmp_path / "broken.yaml"
    path.write_text(yaml.safe_dump(raw), encoding="utf-8")
    with pytest.raises(PrespecError, match="cannot be wrong"):
        load_prespec(path)


# ---------------------------------------------------------------------------
# The instrument checks
# ---------------------------------------------------------------------------


def test_the_positive_control_stops_the_analysis_rather_than_warning(spec):
    """A low R-squared means either 'fraction does not explain this' or
    'deconvolution does not work here'. The positive control is what tells them
    apart, so it cannot be advisory."""
    assert positive_control_gates_the_analysis(spec)


def test_the_positive_control_does_not_use_a_panel_gene(spec):
    """Every canonical mature-colonocyte marker is on the frozen panel, and
    which of them are compositional is the question. A gene-based positive
    control would beg it — this one validates the instrument against ABSOLUTE
    purity, which is called from copy number and so is independent of
    expression."""
    from src.common.panel import panel_genes

    text = spec["instrument_checks"]["positive_control"]["statement"]
    assert not (set(panel_genes()) & set(text.replace(",", " ").split()))
    assert "purity" in text.lower()


def test_a_negative_control_is_named_with_a_threshold(spec):
    negative = spec["instrument_checks"]["negative_control"]
    assert negative["genes"]
    assert "0.10" in negative["statement"]


# ---------------------------------------------------------------------------
# Consistency with the rest of the frozen contract
# ---------------------------------------------------------------------------


def test_the_outcome_genes_are_the_two_that_6_2_names(spec):
    assert outcome_genes(spec) == ["GUCA2A", "CDX2"]


def test_the_covariate_context_is_the_one_the_lock_defines(spec):
    """Purity and plate enter here because the predictor is expression-derived.
    This file does not re-open the covariate set."""
    from src.bulk.covariates import covariate_names, load_covariate_set

    context = spec["model"]["covariate_context"]
    assert context == "expression_models"
    names = covariate_names(load_covariate_set(), endpoint="PFI", context=context)
    assert "purity" in names


def test_cohorts_are_never_pooled(spec):
    """Invariant 4. Between-dataset variation is where the compositional signal
    lives."""
    assert "never" in spec["pooling"].lower()
    assert len(spec["cohorts"]) >= 2


def test_unestimable_patients_are_not_entered_as_zero(spec):
    """Invariant 1, and here it is load-bearing: entering a not-estimable
    patient as zero fraction manufactures the compositional signal the analysis
    is testing for."""
    rule = spec["estimability"]["rule"].lower()
    assert "none" in rule and "never 0.0" in rule


def test_the_granularity_caveats_are_carried(spec):
    """#42 says the curve has three distinct points, not four; #49 adds a
    per-rung quotability field. Both have to reach this analysis."""
    note = spec["granularity"]["usable_rungs_note"]
    assert "#42" in note and "#49" in note
    assert "quotability" in spec["granularity"]["quotability_rule"].lower()


def test_what_is_not_prespecified_is_written_down(spec):
    """Omissions visible, so they read as decisions rather than oversights."""
    assert len(spec["not_prespecified"]) >= 3


def test_the_file_is_valid_yaml_and_round_trips():
    text = PRESPEC_PATH.read_text(encoding="utf-8")
    assert yaml.safe_dump(yaml.safe_load(text))
    assert textwrap.dedent(text).strip()
