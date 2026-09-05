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
    matched_null_genes,
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


def test_the_lock_carries_who_locked_it_and_when(spec):
    """Locked 2026-09-05, in its own commit, so the variance arm can run.

    It shipped `proposed` from 2026-08-28 and the flip is a separate commit
    touching only the lock fields -- that separation is the whole mechanism, so
    a `locked` status with an empty authorisation is not a lock and fails here.
    """
    assert spec["status"] == "locked"
    for field in ("locked_on", "locked_by", "lock_authorisation"):
        assert spec[field], f"{field} is empty on a locked spec"
    assert spec["locked_on"] != spec["proposed_on"], (
        "locked on the day it was proposed, which defeats the mechanism"
    )


def test_the_lock_did_not_edit_the_prediction_in_the_same_commit():
    """A lock that also changed what was predicted is not a lock.

    Checked against git rather than asserted in prose: parse the file on both
    sides of the commit that flipped `status`, drop the six lock fields, and
    require everything else to be byte-for-byte the same object. A prediction
    quietly rewritten at the moment it stopped being editable is the one defect
    the lock mechanism exists to prevent and the one it cannot see itself.
    """
    import subprocess

    from src.common.paths import REPO_ROOT

    rel = str(PRESPEC_PATH.relative_to(REPO_ROOT))

    def git(*args: str) -> subprocess.CompletedProcess:
        return subprocess.run(
            ["git", *args], cwd=REPO_ROOT, capture_output=True, text=True
        )

    log = git("log", "--format=%H", "-S", "status: locked", "--", rel)
    if log.returncode != 0 or not log.stdout.strip():
        pytest.skip("no locking commit in history (shallow clone or archive)")
    sha = log.stdout.split()[0]

    before = git("show", f"{sha}^:{rel}")
    after = git("show", f"{sha}:{rel}")
    if before.returncode != 0 or after.returncode != 0:
        pytest.skip(f"cannot read both sides of {sha[:7]}")

    lock_fields = {"status", "locked_on", "locked_by", "lock_authorisation"}
    strip = lambda text: {  # noqa: E731
        k: v for k, v in yaml.safe_load(text).items() if k not in lock_fields
    }
    was, now = strip(before.stdout), strip(after.stdout)

    changed = sorted(
        k for k in set(was) | set(now) if was.get(k, object()) != now.get(k, object())
    )
    assert not changed, (
        f"the locking commit {sha[:7]} also changed {changed}. A lock that edits "
        f"the specification it is locking is not a lock -- split it into two "
        f"commits so the edit is reviewable on its own."
    )


def test_running_against_a_proposed_spec_is_refused(spec):
    """Constructed rather than read from the file, which is now locked.

    Pinning this to the committed status would have made the test pass for a
    reason unrelated to what it checks, and then silently stop checking it the
    moment the file flipped.
    """
    with pytest.raises(PrespecError, match="not 'locked'"):
        require_locked_prespec(dict(spec, status="proposed"))


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


def test_the_primary_names_both_genes_and_the_matched_null(spec):
    statement = prediction(spec, "primary")["statement"]
    assert "GUCA2A" in statement and "CDX2" in statement
    assert "matched null" in statement.lower()
    assert "percentile" in statement.lower()


def test_the_primary_declares_the_indeterminate_case(spec):
    """Both above or neither above must not be resolved toward either arm."""
    text = prediction(spec, "primary")["disconfirmed_if"].lower()
    assert "indeterminate" in text
    assert "both" in text and "neither" in text


def test_neither_arm_compares_a_raw_r_squared(spec):
    """Issue #54: R-squared is a share of variance, so at the assay floor it
    measures abundance rather than biology. Both arms must be expressed as a
    percentile within, or excess over, the abundance-matched null."""
    for arm in ("primary", "secondary"):
        statement = prediction(spec, arm)["statement"].lower()
        if "r-squared" in statement:
            assert "percentile" in statement or "excess" in statement, arm


def test_reintroducing_a_raw_r_squared_arm_fails_to_load(tmp_path):
    """The guard that stops #54's defect being undone by a later edit."""
    raw = yaml.safe_load(PRESPEC_PATH.read_text(encoding="utf-8"))
    raw["prediction"]["primary"]["statement"] = (
        "R-squared(GUCA2A ~ fraction) < R-squared(CDX2 ~ fraction)."
    )
    path = tmp_path / "regressed.yaml"
    path.write_text(yaml.safe_dump(raw), encoding="utf-8")
    with pytest.raises(PrespecError, match="confounded with abundance"):
        load_prespec(path)


def test_the_withdrawn_arms_are_kept_with_their_reasons(spec):
    """#54 was raised on the #50 review and not addressed before it merged. The
    record of what was withdrawn and why travels with the file, so the same
    statistic is not re-proposed."""
    withdrawn = spec["prediction"]["withdrawn_arms"]
    assert len(withdrawn) >= 3
    for arm in withdrawn:
        assert arm["why"].strip()
        assert arm["withdrawn"].strip()
    statements = " ".join(a["statement"] for a in withdrawn)
    assert "slope" in statements  # the alternative, tested and rejected


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


def test_the_high_abundance_negative_control_names_its_own_gap(spec):
    """ACTB and GAPDH test the top of the range. #54's point is that the floor
    is where the confound lives, and this check must say it does not cover it."""
    negative = spec["instrument_checks"]["negative_control_high_abundance"]
    assert negative["genes"]
    assert "0.10" in negative["statement"]
    assert "floor" in negative["known_gap"].lower()


def test_a_low_abundance_negative_control_exists(spec):
    """Issue #54's actual ask. Drawn by rule rather than named, so it cannot be
    chosen after seeing the result."""
    low = spec["instrument_checks"]["negative_control_low_abundance"]
    assert low["genes"] is None, "must be drawn by rule, not hand-picked"
    assert "abundance-matched" in low["statement"].lower()
    assert low["on_failure"].strip()


def test_the_matched_null_rule_is_fully_committed(spec):
    """Every parameter that could be tuned after the fact is in the file."""
    null = spec["matched_null"]
    for key in ("abundance_window_log2", "min_detection_rate"):
        assert key in null["candidates"], key
    assert isinstance(null["max_genes"], int)
    assert isinstance(null["seed"], int)
    assert null["honest_limitation"].strip()


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


# ---------------------------------------------------------------------------
# The matched-null draw — issue #54's fix, as code rather than prose
# ---------------------------------------------------------------------------


def _toy(n_genes=400, n_samples=60):
    """Genes spread across abundance, plus a couple of named panel genes."""
    import numpy as np
    import pandas as pd

    rng = np.random.default_rng(0)
    ids = [f"ENSG{i:08d}" for i in range(n_genes)]
    levels = np.linspace(0.5, 12.0, n_genes)
    x = pd.DataFrame(
        rng.normal(levels, 0.05, size=(n_samples, n_genes)).clip(0),
        columns=ids,
    )
    index_map = pd.DataFrame(
        {
            "gene_symbol": [f"SYM{i}" for i in range(n_genes)],
            "on_panel": ["False"] * n_genes,
        },
        index=pd.Index(ids, name="ensembl_id"),
    )
    return x, index_map, ids


def test_the_null_is_drawn_at_the_targets_abundance(spec):
    x, index_map, ids = _toy()
    target = ids[200]
    drawn = matched_null_genes(x, index_map, target, spec)
    window = spec["matched_null"]["candidates"]["abundance_window_log2"]
    med = x.median(axis=0)
    assert drawn
    assert (med[drawn] - med[target]).abs().max() <= window + 1e-9


def test_the_target_is_never_in_its_own_null(spec):
    x, index_map, ids = _toy()
    assert ids[200] not in matched_null_genes(x, index_map, ids[200], spec)


def test_panel_genes_are_excluded_from_the_null(spec):
    """They are the outcome variables. A panel gene in the null would put the
    thing being measured into what it is measured against."""
    x, index_map, ids = _toy()
    target = ids[200]
    neighbour = ids[201]
    index_map.loc[neighbour, "on_panel"] = "True"
    assert neighbour not in matched_null_genes(x, index_map, target, spec)


def test_the_draw_is_deterministic_under_the_committed_seed(spec):
    """Invariant 10's spirit: the same input gives the same null, and the seed
    is in the config rather than chosen at run time."""
    x, index_map, ids = _toy()
    first = matched_null_genes(x, index_map, ids[200], spec)
    second = matched_null_genes(x, index_map, ids[200], spec)
    assert first == second


def test_the_null_is_capped_at_the_committed_size(spec):
    x, index_map, ids = _toy(n_genes=4000)
    drawn = matched_null_genes(x, index_map, ids[2000], spec)
    assert len(drawn) <= spec["matched_null"]["max_genes"]


def test_no_matched_candidates_raises_rather_than_widening_the_window(spec):
    """Widening it after seeing that nothing qualified is a change to the
    pre-specification, not a run-time decision."""
    import pandas as pd

    x, index_map, _ = _toy(n_genes=3)
    x = pd.DataFrame({"ENSG00000000": [1.0] * 5, "ENSG00000001": [40.0] * 5})
    index_map = pd.DataFrame(
        {"gene_symbol": ["A", "B"], "on_panel": ["False", "False"]},
        index=pd.Index(["ENSG00000000", "ENSG00000001"], name="ensembl_id"),
    )
    with pytest.raises(PrespecError, match="Widening the window"):
        matched_null_genes(x, index_map, "ENSG00000000", spec)


def test_a_target_absent_from_the_matrix_is_a_loud_failure(spec):
    x, index_map, _ = _toy()
    with pytest.raises(PrespecError, match="not in the expression matrix"):
        matched_null_genes(x, index_map, "ENSG99999999", spec)
