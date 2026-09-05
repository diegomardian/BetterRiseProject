"""Every number the paper quotes, against the table it came from.

This paper's rule is that a check unable to fail is worse than no check. The
paper itself had no check at all: its numbers were transcribed from result
tables by hand, and nothing would have noticed if a table were re-derived and a
figure in the prose left behind. That is the same defect as
``paper/wmhs/sections/appendix.tex`` item 3 -- a claim and the code that
produced it travelling separately -- one layer further out.

So the numbers in the third withdrawn guard are re-derived here from the
committed parquet and asserted against the literal strings in the ``.tex``.
Editing either side alone fails.

The trap this file must itself avoid is a regex that matches nothing: an
assertion over an empty match set passes on any input, which is the exact
failure mode the paper documents. Every helper below therefore asserts the
literal is *present* before asserting anything about its value.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from src.common.paths import REPO_ROOT, RESULTS_DIR

SECTIONS = REPO_ROOT / "paper" / "wmhs" / "sections"

#: The guard's own tolerances, quoted in the prose. Imported rather than
#: retyped, so moving one in the code breaks the paper here rather than in
#: review.
from src.reference.jobs.coexpression_silencing import (  # noqa: E402
    CONTROL_LOG2_TOLERANCE,
    CONTROL_TOLERANCE,
)


def _newest(name: str) -> Path:
    matches = sorted(RESULTS_DIR.glob(f"*/{name}.parquet"))
    if not matches:
        pytest.skip(f"no results/*/{name}.parquet committed")
    return matches[-1]


@pytest.fixture(scope="module")
def coexpr() -> pd.DataFrame:
    return pd.read_parquet(_newest("coexpression_silencing"))


@pytest.fixture(scope="module")
def coexpr_meta() -> dict:
    path = _newest("coexpression_silencing_summary")
    return json.loads(path.with_suffix("").with_suffix(".meta.json").read_text())


@pytest.fixture(scope="module")
def tex() -> str:
    """The three files carrying the third guard's numbers, concatenated."""
    return "\n".join(
        (SECTIONS / f).read_text()
        for f in ("withdrawn.tex", "conclusion.tex", "responsible.tex")
    )


def _quotes(tex: str, literal: str) -> None:
    """Assert the paper contains this literal.

    Called before every value assertion. Without it a renamed or reworded
    passage would silently stop being checked, and this file would keep
    passing -- a check that cannot fail, in the test written to prevent them.
    """
    assert literal in tex, (
        f"the paper no longer contains {literal!r}. Either the prose moved, in "
        f"which case update this test, or a number was edited without its "
        f"table -- which is what this file exists to catch."
    )


# ---------------------------------------------------------------------------
# The third withdrawn guard: the saturated premise control


def test_the_cohort_and_patient_counts_are_the_table_s(tex, coexpr):
    n_studies = coexpr["study_id"].nunique()
    n_patients = coexpr.groupby("study_id")["patient_id"].nunique().sum()
    _quotes(tex, "three\ncohorts and 42 patients")
    _quotes(tex, "over three cohorts and 42 patients")
    assert n_studies == 3
    assert n_patients == 42


def test_the_control_pass_rate_is_the_table_s(tex, coexpr):
    """``83 of 84`` control rows inside the detection tolerance."""
    controls = coexpr[coexpr["role"] == "control"]
    inside = (controls["delta_detect"].abs() < CONTROL_TOLERANCE).sum()
    _quotes(tex, r"\textbf{83 of 84}")
    _quotes(tex, "inside a\n$0.10$ tolerance")
    assert CONTROL_TOLERANCE == 0.10
    assert (inside, len(controls)) == (83, 84)


def test_the_controls_moved_by_the_fold_range_quoted(tex, coexpr):
    """``factors of $1.19$ to $1.64$`` in per-cell output.

    Per study and gene, since a per-patient extreme is a different and larger
    claim than the one the paper makes.
    """
    controls = coexpr[coexpr["role"] == "control"]
    mean_log2 = controls.groupby(["study_id", "gene"])["log2_cp10k_ratio"].mean()
    folds = 2.0 ** mean_log2.abs()
    _quotes(tex, "factors of $1.19$ to $1.64$ in")
    assert (round(folds.min(), 2), round(folds.max(), 2)) == (1.19, 1.64)


def test_the_reference_detection_rates_are_the_table_s(tex, coexpr):
    """``$0.967$ to $1.000$`` -- why the detection statistic has no room."""
    controls = coexpr[coexpr["role"] == "control"]
    med = controls.groupby(["study_id", "gene"])["detect_normal"].median()
    _quotes(tex, "reference detection rate of $0.967$ to $1.000$")
    assert (round(med.min(), 3), round(med.max(), 3)) == (0.967, 1.000)


def test_the_poisson_claim_is_arithmetically_true(tex):
    """A gene detected in 99% of cells, halved, still reads 0.90.

    The one number in the paragraph that is not read off a table. It is a
    closed form, so it is checked as one.
    """
    _quotes(tex, "detected in 99\\% of cells whose output \\emph{halves} still reads $0.90$")
    d, s = 0.99, 0.5
    assert round(1.0 - (1.0 - d) ** s, 2) == 0.90


def test_every_control_interval_straddles_the_tolerance(tex, coexpr_meta):
    """The three quoted intervals, and the property that makes them the finding.

    Each is asserted against the sidecar AND checked to straddle the tolerance,
    because straddling is what forces the third verdict. An interval that
    stopped straddling would leave the numbers correct and the argument wrong.
    """
    quoted = {
        "GSE132465": (-0.712, -1.250, -0.289),
        "GSE144735": (+0.540, +0.012, +1.112),
        "GSE178341": (+0.431, +0.216, +0.640),
    }
    _quotes(tex, "$-0.712\\,[-1.250, -0.289]$, $+0.540\\,[+0.012, +1.112]$,")
    _quotes(tex, "$+0.431\\,[+0.216, +0.640]$")

    readings = coexpr_meta["premise_reading"]
    assert set(readings) == set(quoted)
    for study, (mean, lo, hi) in quoted.items():
        text = readings[study]
        assert text.startswith("UNRESOLVED"), f"{study} is no longer UNRESOLVED: {text}"
        for value in (mean, lo, hi):
            assert f"{value:+.3f}" in text or f"{abs(value):.3f}" in text, (
                f"{value} not in {study}'s recorded reading: {text}"
            )
        # The premise reads |shift| against the tolerance, so the interval
        # that must straddle it is the interval of ABSOLUTE shift. For a
        # wholly-negative interval |lo| and |hi| swap order, and for one
        # containing zero the lower end is 0 rather than either endpoint --
        # writing `abs(lo) < tol < abs(hi)` gets both cases wrong.
        reaches = max(abs(lo), abs(hi))
        floor = 0.0 if lo <= 0.0 <= hi else min(abs(lo), abs(hi))
        assert floor < CONTROL_LOG2_TOLERANCE < reaches, (
            f"{study}'s interval [{lo}, {hi}] no longer straddles the "
            f"{CONTROL_LOG2_TOLERANCE} tolerance on |shift|, so the third "
            f"verdict is not forced and the paragraph's argument no longer holds."
        )


def test_the_premise_is_undecided_on_every_cohort(tex, coexpr_meta):
    _quotes(tex, "returns undecided on all three cohorts")
    _quotes(tex, "returns \\emph{undecided} on all three")
    readings = coexpr_meta["premise_reading"]
    assert len(readings) == 3
    assert all(r.startswith("UNRESOLVED") for r in readings.values())


def test_the_conclusion_s_maximum_shift_is_the_table_s(tex, coexpr):
    """``up to $0.71$ on log$_2$ expression``, in the conclusion."""
    controls = coexpr[coexpr["role"] == "control"]
    worst = controls.groupby(["study_id", "gene"])["log2_cp10k_ratio"].mean().abs().max()
    _quotes(tex, "shift by up to $0.71$ on log$_2$ expression")
    assert round(worst, 2) == 0.71


def test_the_seed_flip_is_the_one_recorded_in_the_job(tex):
    """``+0.486`` and ``+0.540`` come from a comment, not a table.

    It is the only figure in the paragraph with no parquet behind it, because
    it is a value a superseded seed produced. Pin it to the source that records
    it so the two cannot drift.
    """
    source = (
        REPO_ROOT / "src" / "reference" / "jobs" / "coexpression_silencing.py"
    ).read_text()
    _quotes(tex, "read $+0.486$ under one seed and $+0.540$ under another")
    assert "+0.486" in source and "+0.540" in source, (
        "the job no longer records the seed flip the paper cites"
    )


# ---------------------------------------------------------------------------
# The counts the prose makes about itself


def test_the_paper_counts_its_own_withdrawn_guards(tex):
    """Three guards described, three claimed, in both places that claim it."""
    withdrawn = (SECTIONS / "withdrawn.tex").read_text()
    described = withdrawn.count("\\paragraph{")
    assert described == 3, f"{described} guards described, not 3"
    _quotes(tex, "Three more guards shipped and got withdrawn")
    appendix = (SECTIONS / "appendix.tex").read_text()
    assert "\\subsection*{Three more checks that could not fire}" in appendix
    assert "Two more checks that could not fire" not in appendix


def test_the_conclusion_counts_the_statistics_it_lists(tex):
    """``Four ... could not fail``, and four semicolon-separated clauses.

    The sentence said ``two`` while listing three for as long as it existed.
    A count in a paper about miscounted claims is worth a test.
    """
    conclusion = (SECTIONS / "conclusion.tex").read_text()
    _quotes(tex, "Four of this paper's statistics could not fail.")
    sentence = conclusion.split("could not fail.", 1)[1].split("Each looked like")[0]
    assert sentence.count(";") == 3, (
        f"{sentence.count(';') + 1} statistics listed against a claimed four"
    )
