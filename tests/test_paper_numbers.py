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
    """Three guards described, three claimed -- and a fourth thing that is not one.

    ``withdrawn.tex`` carries four paragraphs. Three are withdrawn guards; the
    fourth is the interval, which is not a guard and whose whole point is that
    it is not. Counting paragraphs and calling the answer "guards" would be a
    miscount of exactly the kind this file exists to catch, so the two are
    separated here rather than summed.
    """
    withdrawn = (SECTIONS / "withdrawn.tex").read_text()
    paragraphs = withdrawn.count("\\paragraph{")
    not_a_guard = withdrawn.count("\\paragraph{And the fifth is not a guard at all.")
    assert not_a_guard == 1, "the interval paragraph is missing or was retitled"
    assert paragraphs - not_a_guard == 3, (
        f"{paragraphs - not_a_guard} guards described, not 3"
    )
    _quotes(tex, "Three more guards shipped and got withdrawn")
    appendix = (SECTIONS / "appendix.tex").read_text()
    # The heading moved from a subsection of the short build's appendix to a
    # section of both builds' appendix when the page count was first measured
    # against real geometry. The count claim has to survive the move.
    assert "\\section{Three guards that could not fire, and the interval}" in appendix
    assert "Two more checks that could not fire" not in appendix
    assert "Two guards that could not fire" not in appendix


def test_the_conclusion_counts_the_statistics_it_lists(tex):
    """``Five ... could not fail``, and five semicolon-separated clauses.

    The sentence said ``two`` while listing three for as long as it existed.
    A count in a paper about miscounted claims is worth a test.
    """
    conclusion = (SECTIONS / "conclusion.tex").read_text()
    _quotes(tex, "Five of this paper's statistics could not fail.")
    sentence = conclusion.split("could not fail.", 1)[1].split("Each looked like")[0]
    assert sentence.count(";") == 4, (
        f"{sentence.count(';') + 1} statistics listed against a claimed five"
    )


# ---------------------------------------------------------------------------
# The denser cutpoint grid (sec:calibration) and the ratio collapse (sec:bench).
#
# Added with the two paragraphs they check. Both paragraphs quote numbers that
# were previously carried only by prose in docs/, which is the failure mode
# HANDOFF.md's ledger records twice: a number in a document with nothing tying
# it to a table. These fixtures re-derive every one of them.


@pytest.fixture(scope="module")
def calibration_tex() -> str:
    return (SECTIONS / "calibration.tex").read_text()


@pytest.fixture(scope="module")
def bench_tex() -> str:
    return (SECTIONS / "bench.tex").read_text()


def _cutpoints(cohort: str) -> pd.DataFrame:
    """Committed, extended and dense grids for one cohort, in one frame."""
    suffix = "" if cohort == "smc" else f"_{cohort}"
    coarse = pd.read_parquet(_newest(f"calibration_gap_cutpoints{suffix}_r500"))
    dense = pd.read_parquet(_newest(f"cutpoint_dense_grid_cutpoints_{cohort}_r200"))
    return pd.concat([coarse, dense], ignore_index=True)


def test_the_dense_grid_returns_the_ok_cutpoint_quoted(calibration_tex):
    """``$ok = 70$ on eight of eight seeds`` on the reference pool."""
    dense = _cutpoints("smc")
    ref = dense[(dense["grid"] == "dense") & (dense["pool"] == "reference")]
    _quotes(calibration_tex, "returns $ok = 70$ on eight of eight seeds")
    assert len(ref) == 8
    assert sorted(ref["ok"].unique()) == [70.0]


def test_wide_has_zero_seed_variance_within_every_grid(calibration_tex):
    """The load-bearing claim: one value per grid, on all eight seeds.

    Six cohort-by-grid cells. If any cell ever carries two values the
    paragraph's central sentence is false, and this is what says so.
    """
    _quotes(
        calibration_tex,
        r"seed-to-seed variance --- one value on all eight seeds, in all six",
    )
    # ...and the six values themselves. Asserting the table without quoting
    # the prose list leaves the list free to drift -- this file's own
    # failure mode, caught here by mutation rather than by review.
    _quotes(
        calibration_tex,
        "($40$, $45$, $42$ here; $100$, $65$, $70$ on the" + chr(10) + "second cohort)",
    )
    cells = {}
    for cohort in ("smc", "kul3"):
        frame = _cutpoints(cohort)
        ref = frame[frame["pool"] == "reference"]
        for grid, block in ref.groupby("grid"):
            values = sorted(block["wide"].dropna().unique())
            assert len(block) == 8, (cohort, grid, len(block))
            assert len(values) == 1, f"{cohort}/{grid} carries {values}"
            cells[(cohort, grid)] = values[0]
    assert len(cells) == 6
    assert [cells[("smc", g)] for g in ("committed", "extended", "dense")] == [
        40.0,
        45.0,
        42.0,
    ]
    assert [cells[("kul3", g)] for g in ("committed", "extended", "dense")] == [
        100.0,
        65.0,
        70.0,
    ]


def test_the_pooled_draw_returns_nothing_on_the_dense_grid(calibration_tex):
    """``topping out at $0.750$ and\n$0.795$`` against the 0.80 target."""
    _quotes(calibration_tex, "topping out at $0.750$ and $0.795$ against the $0.80$ target")
    maxima = []
    for cohort in ("smc", "kul3"):
        frame = _cutpoints(cohort)
        pooled = frame[(frame["grid"] == "dense") & (frame["pool"] == "pooled")]
        assert len(pooled) == 8
        assert not pooled["returned_a_cutpoint"].any()
        maxima.append(round(pooled["max_discrimination"].max(), 3))
    assert maxima == [0.750, 0.795]


def _ratio_terms(frame: pd.DataFrame, rung: str, normal: str, tumour: str):
    """``-f_N/df`` and each gene's ``m_T/m_N`` at one resolution."""
    block = frame[frame["granularity_rung"] == rung]
    f_n = block["frac_mature_normal"].mean()
    f_t = block["frac_mature_tumour"].mean()
    limit = f_n / (f_t - f_n)
    means = block.groupby("gene")[[normal, tumour]].mean()
    survived = means[tumour] / means[normal]
    return limit, survived


def test_the_carcinoma_ratio_collapses_onto_the_limit_quoted(bench_tex):
    """``$6.94$``, and six genes inside a ``$1.07$-fold`` spread."""
    frame = pd.read_parquet(_newest("decomposition_summary_matched"))
    limit, survived = _ratio_terms(
        frame, "lineage", "mean_normal", "mean_tumour"
    )
    _quotes(bench_tex, "that limit is $6.94$")
    _quotes(bench_tex, "$m_T/m_N$ from $0.045$ to $0.109$")
    _quotes(bench_tex, r"return ratios from $6.19$ to $6.63$, a spread of $1.07$-fold")
    assert round(-limit, 2) == 6.94

    nearest = survived.nsmallest(6)
    assert (round(nearest.min(), 3), round(nearest.max(), 3)) == (0.045, 0.109)
    ratios = limit * (nearest - 1.0)
    assert (round(ratios.min(), 2), round(ratios.max(), 2)) == (6.19, 6.63)
    assert round(ratios.max() / ratios.min(), 2) == 1.07


def test_the_adenoma_ratio_does_not_collapse(bench_tex):
    """The negative control: ``$0.27$ to\n$3.13$``, an ``$11.5$-fold`` spread."""
    frame = pd.read_parquet(_newest("icbi_adenoma"))
    limit, survived = _ratio_terms(
        frame, "lineage", "cp10k_normal", "cp10k_tumour"
    )
    _quotes(bench_tex, "the minimum\n$m_T/m_N$ is $0.374$")
    _quotes(bench_tex, "the identical computation returns $0.27$ to\n$3.13$")
    _quotes(bench_tex, r"an $11.5$-fold spread where carcinoma gave $1.07$")
    assert round(survived.min(), 3) == 0.374

    ratios = limit * (survived - 1.0)
    assert (round(ratios.min(), 2), round(ratios.max(), 2)) == (0.27, 3.13)
    assert round(ratios.max() / ratios.min(), 1) == 11.5


# ---------------------------------------------------------------------------
# sec:blind's recovery ranges. An external reviewer could not reproduce these
# from the committed tables and proposed cutting them. They are correct; what
# was missing was the filter, which lived only in make_fig3.py. Pinned here so
# the sentence and the figure cannot drift apart, and so the next reader does
# not have to rediscover the aggregation.


@pytest.fixture(scope="module")
def blind_tex() -> str:
    return (SECTIONS / "blind.tex").read_text()


def test_the_recovery_ranges_are_the_figures(blind_tex):
    """``between 1.00 and 1.07`` on reference, ``0.86 to 1.18`` pooled.

    The filter is Appendix Fig. 2's own, and every clause of it matters:
    the 50-replicate table, the detectable shift, the extended grid, the
    median ACROSS SEEDS at each count, and only counts where the committed
    rule says ``estimate``. Taking min/max over raw rows instead of the
    per-count median reproduces neither range.
    """
    from src.harness.positivity import CUTPOINTS

    rec = pd.read_parquet(_newest("calibration_gap_recovery"))
    at_effect = rec[(rec["shift"] == 0.5) & (rec["grid"] == "extended")]

    _quotes(blind_tex, "Ours ran between 1.00 and 1.07 wherever the rule said")
    _quotes(blind_tex, "0.86 to 1.18 on the pooled one")
    assert CUTPOINTS.ok == 50

    spans = {}
    for pool in ("reference", "pooled"):
        by_count = (
            at_effect[at_effect["pool"] == pool]
            .groupby("median_n_cells_mature")["ratio_median"]
            .median()
        )
        estimable = by_count[by_count.index >= CUTPOINTS.ok]
        spans[pool] = (round(estimable.min(), 2), round(estimable.max(), 2))

    assert spans["reference"] == (1.00, 1.07)
    assert spans["pooled"] == (0.86, 1.18)
