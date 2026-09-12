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
import re
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
    _quotes(bench_tex, "$m_T/m_N$ from $0.045$ to")
    _quotes(bench_tex, r"$6.63$, a spread of $1.07$-fold")
    _quotes(bench_tex, "the six of its nine" + chr(10) + "panel genes")
    _quotes(bench_tex, "one rises $4.5$-fold between the arms")
    assert round(-limit, 2) == 6.94

    # six of nine, and the three set aside are not near the limit -- the
    # selection is the claim's premise, so the paper has to name it
    assert len(survived) == 9
    assert round(survived.max(), 1) == 4.5

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
    _quotes(bench_tex, "which on this cohort's six-gene panel takes all of them")
    _quotes(bench_tex, r"an $11.5$-fold spread where carcinoma gave $1.07$")
    assert round(survived.min(), 3) == 0.374

    # "the identical computation" has to be identical: the same six-nearest-zero
    # rule, which on a six-gene panel is every gene. If this cohort ever gains a
    # seventh the two sides stop being comparable and the sentence is wrong.
    assert len(survived) == 6
    assert set(survived.nsmallest(6).index) == set(survived.index)

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


# ---------------------------------------------------------------------------
# Merges can reintroduce deleted prose, and every check above would still pass.
#
# This happened. A branch cut before section 3 was rewritten carried the
# superseded "denser grid" paragraph; merging it added that paragraph back
# beside its replacement, and the paper said the same thing twice in two
# different ways. Every assertion in this file passed, because they all ask
# whether a literal is PRESENT. None asks whether it is present twice, and none
# asks whether something deleted stayed deleted.


#: Prose retired from the paper. A merge that resurrects any of these is a
#: regression, not a contribution.
RETIRED = (
    "A denser grid, and a quantity whose stability is the apparatus",
    "The answer moves to 90}, and the crossing it",
    "roughly twice as strict on both cutpoints",
    "it is the one defect in this paper that is now closed",
    "three of the four look fine",
)


def test_retired_prose_stays_retired():
    """Nothing a rewrite deleted has been merged back in."""
    body = "\n".join(
        p.read_text(encoding="utf-8") for p in sorted(SECTIONS.glob("*.tex"))
    )
    resurrected = [line for line in RETIRED if line in body]
    assert not resurrected, (
        f"deleted prose is back in the paper: {resurrected}. A merge from a "
        f"branch cut before the rewrite is the usual cause."
    )


def _render(text: str, *, full: bool) -> str:
    r"""The lines LaTeX would typeset for one build.

    A naive duplicate check flags ``bench.tex``, which carries one paragraph
    heading in each branch of an ``\iffull`` -- correct, and never duplicated
    in either PDF. So evaluate the conditional first and check what a reader
    actually sees.
    """
    out, stack = [], []
    for line in text.splitlines():
        bare = line.strip()
        if bare.startswith(chr(92) + "iffull"):
            stack.append(full)
            bare = bare[len(chr(92) + "iffull"):].strip()
            # "\iffull\else" on one line is how bench.tex writes its
            # short-build branch. Missing this made the check flag that file.
            while bare.startswith(chr(92) + "else"):
                stack[-1] = not stack[-1]
                bare = bare[len(chr(92) + "else"):].strip()
            if not bare:
                continue
            line = bare
        elif bare == chr(92) + "else":
            if stack:
                stack[-1] = not stack[-1]
            continue
        elif bare == chr(92) + "fi":
            if stack:
                stack.pop()
            continue
        if all(stack):
            out.append(line)
    return chr(10).join(out)


def test_no_paragraph_heading_appears_twice():
    """Two paragraphs with one title means a merge duplicated a block."""
    pattern = re.escape(chr(92)) + r"paragraph\{([^}]{12,})\}"
    for build in (True, False):
        for path in sorted(SECTIONS.glob("*.tex")):
            rendered = _render(path.read_text(encoding="utf-8"), full=build)
            headings = re.findall(pattern, rendered)
            duplicated = {h for h in headings if headings.count(h) > 1}
            assert not duplicated, (
                f"{path.name} carries {duplicated} more than once in the "
                f"{'full' if build else 'short'} build"
            )


# ---------------------------------------------------------------------------
# The three experiments added 2026-09-10: the residual matrix, the information
# ratio, censored survival, and the benchmark extensions. Same rule as above --
# every figure quoted in the prose is re-derived from its committed table.


@pytest.fixture(scope="module")
def matrix() -> pd.DataFrame:
    return pd.read_parquet(_newest("trial_blindness_residual_matrix"))


@pytest.fixture(scope="module")
def rho() -> pd.DataFrame:
    return pd.read_parquet(_newest("trial_blindness_information_ratio"))


@pytest.fixture(scope="module")
def survival() -> pd.DataFrame:
    return pd.read_parquet(_newest("trial_survival_headline"))


def _cell(frame, estimator, design, truth, n=5000):
    row = frame[
        (frame["estimator"] == estimator)
        & (frame["design"] == design)
        & (frame["truth"] == truth)
        & (frame["n_patients"] == n)
    ]
    assert len(row) == 1, (estimator, design, truth, n, len(row))
    return float(row["max_residual"].iloc[0])


def test_ols_is_blind_against_its_own_estimand(blind_tex, matrix):
    """The pair claim: 0.0663 against one truth, machine zero against the other."""
    _quotes(blind_tex, "has residual $0.0663$ against the standardised realised truth")
    _quotes(blind_tex, "$3.5" + chr(92) + "times10^{-14}$, exactly blind")
    std = _cell(matrix, "ols-stratum-dummies", "confounded-bernoulli", "standardised")
    var = _cell(matrix, "ols-stratum-dummies", "confounded-bernoulli", "varweighted")
    assert round(std, 4) == 0.0663
    assert var < 1e-12 and round(var * 1e14, 1) == 3.5
    # and the swap runs the other way
    g_std = _cell(matrix, "gcomp-from-generator", "confounded-bernoulli", "standardised")
    g_var = _cell(matrix, "gcomp-from-generator", "confounded-bernoulli", "varweighted")
    assert g_std == 0.0 and round(g_var, 4) == 0.0663


def test_block_randomisation_blinds_even_the_unadjusted_difference(blind_tex, matrix):
    """Under 1:1 allocation the estimands coincide and everything collapses."""
    _quotes(blind_tex, "and even the unadjusted difference")
    for est in ("ols-stratum-dummies", "unadjusted", "ipw-saturated", "aipw-saturated"):
        assert _cell(matrix, est, "block-randomised", "standardised") < 1e-12, est
    # the design is what does it: confounded, the unadjusted difference is huge
    assert _cell(matrix, "unadjusted", "confounded-bernoulli", "standardised") > 4.0


def test_aipw_with_saturated_nuisances_is_blind(blind_tex, matrix):
    """Cross-fitting is the only thing between AIPW and a vacuous curve."""
    _quotes(blind_tex, "saturated nuisances and no cross-fitting is $4.4")
    assert _cell(matrix, "aipw-saturated", "confounded-bernoulli", "standardised") < 1e-13
    assert _cell(matrix, "ipw-cross-fitted", "confounded-bernoulli", "standardised") > 0.1


def test_the_information_ratio_does_not_improve_with_n(blind_tex, rho):
    """OLS sits at 0.13 at every cohort size -- 88% generator noise."""
    _quotes(blind_tex, chr(92) + "rho$ between $0.129$ and $0.137$ at every cohort size")
    _quotes(blind_tex, "88" + chr(92) + "% of its recovery curve is generator")
    ols = rho[
        (rho["estimator"] == "ols-stratum-dummies")
        & (rho["design"] == "confounded-bernoulli")
        & (rho["truth"] == "standardised")
    ]
    assert len(ols) >= 6
    assert (round(ols["rho"].min(), 3), round(ols["rho"].max(), 3)) == (0.129, 0.137)
    share = ols[ols["n_patients"] == 5000]["generator_noise_share"].iloc[0]
    assert round(share * 100) == 88


def test_censoring_silences_the_check_but_not_against_observed(blind_tex, survival):
    """Zero against observed everywhere; growing against the latent times."""
    _quotes(blind_tex, "= 0$ in all 28 cells")
    _quotes(blind_tex, "$0$, $0.061$, $0.144$, $0.376$ at")
    assert len(survival) == 28
    assert (survival["max_blind_vs_observed"] == 0.0).all()
    at3000 = survival[survival["n_patients"] == 3000].set_index("censoring_target")
    got = [round(at3000.loc[c, "max_blind_vs_latent"], 3) for c in (0.0, 0.12, 0.49, 0.78)]
    assert got == [0.0, 0.061, 0.144, 0.376], got


def test_the_inversion_is_a_majority_not_a_reversal(blind_tex, survival):
    """4 of 7 at 49%, 6 of 7 at 78% -- and the paper says majority, not all."""
    _quotes(blind_tex, "inverts in 5 of 7 cohort sizes and by $78" + chr(92) + "%$")
    _quotes(blind_tex, "in 6 of 7")
    inverted = survival["ordering_by_latent_truth"].str.startswith("INVERTED")
    counts = survival.assign(inv=inverted).groupby("censoring_target")["inv"].sum()
    assert int(counts.loc[0.49]) == 5
    assert int(counts.loc[0.78]) == 6
    assert int(counts.loc[0.0]) == 0


def test_heterogeneity_moves_the_curve_not_the_estimator(rho):
    """The converse result: OLS's curve wanders, its own residual does not."""
    het = pd.read_parquet(_newest("trial_blindness_heterogeneity"))
    appendix = (SECTIONS / "appendix.tex").read_text(encoding="utf-8")
    _quotes(appendix, "$1.003$, $1.011$, $1.020$, $1.037$, $1.077$, $1.132$")
    _quotes(appendix, "grows $0.094 " + chr(92) + "rightarrow 1.452$")

    ols = het[het["estimator"] == "ols-stratum-dummies"]
    std = ols[ols["truth"] == "standardised"].sort_values("sweep_value")
    var = ols[ols["truth"] == "varweighted"].sort_values("sweep_value")

    curve = [round(v, 3) for v in std["median_recovery_ratio"]]
    assert curve == [1.003, 1.011, 1.020, 1.037, 1.077, 1.132], curve
    assert (round(std["max_residual"].min(), 3),
            round(std["max_residual"].max(), 3)) == (0.094, 1.452)
    # against its own estimand it never leaves machine zero
    assert var["max_residual"].max() < 1e-12


def test_the_propensity_spread_is_continuous(rho):
    """0.0033 to 0.104 as the design goes RCT -> confounded, curve flat."""
    spread = pd.read_parquet(_newest("trial_blindness_propensity_spread"))
    appendix = (SECTIONS / "appendix.tex").read_text(encoding="utf-8")
    _quotes(appendix, "rises smoothly $0.0033\n" + chr(92) + "rightarrow 0.104$")

    ols = spread[
        (spread["estimator"] == "ols-stratum-dummies")
        & (spread["truth"] == "standardised")
    ].sort_values("sweep_value")
    assert round(ols["max_residual"].min(), 4) == 0.0033
    assert round(ols["max_residual"].max(), 3) == 0.104
    # the recovery ratio says nothing about any of it
    assert ols["median_recovery_ratio"].between(0.99, 1.01).all()


# ---------------------------------------------------------------------------
# The numbers a 2026-09-11 audit found had drifted from their tables.
#
# Every one of these was wrong in the committed paper and every one was
# reachable from a committed table, which is the definition of a claim this
# file should already have been asserting. Six of them were single-seed rows
# quoted as if they were summaries; two overstated an agreement; one named a
# world that was not refusing yet. They are pinned here so the next edit that
# moves one without its table fails instead of shipping.


def _reference_bins(cohort: str) -> pd.DataFrame:
    """Coarse and dense grids' per-bin rates, reference pool, one frame."""
    suffix = "" if cohort == "smc" else f"_{cohort}"
    coarse = pd.read_parquet(_newest(f"calibration_gap_bins{suffix}_r500"))
    dense = pd.read_parquet(_newest(f"cutpoint_dense_grid_bins_{cohort}_r200"))
    frame = pd.concat([coarse, dense], ignore_index=True)
    return frame[frame["pool"] == "reference"]


def _median_discrimination(frame: pd.DataFrame, grid: str, count: float) -> float:
    cell = frame[(frame["grid"] == grid) & (frame["n_cells_mature"] == count)]
    assert len(cell) == 8, f"{grid}@{count} has {len(cell)} seeds, expected 8"
    return round(float(cell["discrimination"].median()), 3)


def test_the_three_grids_agree_on_medians_not_on_one_seed(calibration_tex):
    """``medians over eight seeds of 0.712, 0.728 and 0.726 near 40``.

    These six were seed 1's rows, quoted as a cross-grid agreement. At seed 1
    they read 0.712/0.728/0.728 and 0.882/0.858/0.878; the medians are two
    digits away from that. The agreement claim survives the correction, which
    is why the fix was the number and not the sentence -- but a single draw
    presented as an agreement is this section's own subject.
    """
    _quotes(
        calibration_tex,
        "(medians over eight seeds of $0.712$, $0.728$ and $0.726$ near 40, and $0.880$,"
        + chr(10)
        + "$0.858$ and $0.878$ near 100)",
    )
    bins = _reference_bins("smc")
    assert [
        _median_discrimination(bins, "committed", 40.0),
        _median_discrimination(bins, "extended", 45.0),
        _median_discrimination(bins, "dense", 42.5),
    ] == [0.712, 0.728, 0.726]
    assert [
        _median_discrimination(bins, "committed", 100.0),
        _median_discrimination(bins, "extended", 90.0),
        _median_discrimination(bins, "dense", 105.0),
    ] == [0.880, 0.858, 0.878]


def test_the_new_bin_reads_the_same_in_both_builds(calibration_tex):
    """``discrimination reads 0.789``, and the short build says 0.789 too.

    The full build quoted seed 1 (0.790) and the short build the median
    (0.789) for one bin. Both builds read the same files, so a number
    differing between them is the one failure the build promises cannot
    happen -- and it had happened, in the last digit.
    """
    _quotes(calibration_tex, "where discrimination reads $0.789$")
    _quotes(calibration_tex, "$0.789$ at 65 and $0.858$ at 90")
    bins = _reference_bins("smc")
    assert _median_discrimination(bins, "extended", 65.0) == 0.789


def test_the_pooled_curve_widens_by_the_factor_quoted():
    """``from about 0.62 to about 1.17`` -- medians over the eight seeds."""
    appendix = (SECTIONS / "appendix.tex").read_text(encoding="utf-8")
    _quotes(appendix, "from about $0.62$ to about $1.17$")
    _quotes(appendix, "$1.331$ to $1.423$ across eight")
    _quotes(appendix, "against $1.075$ to $1.108$ on the reference draw")
    frame = pd.read_parquet(_newest("calibration_gap_recovery_r500"))
    point = frame[
        (frame["grid"] == "committed")
        & (frame["shift"] == 0.5)
        & (frame["frac_mature_tumour"] == 0.01)
    ]
    widths, medians = {}, {}
    for pool, block in point.groupby("pool"):
        assert len(block) == 8
        widths[pool] = round(float((block["ratio_q75"] - block["ratio_q25"]).median()), 2)
        medians[pool] = (
            round(float(block["ratio_median"].min()), 3),
            round(float(block["ratio_median"].max()), 3),
        )
        # the whole point: the curve moves and the residual does not
        assert block["max_abs_residual_vs_realised"].max() == 0.0
    assert (widths["reference"], widths["pooled"]) == (0.62, 1.17)
    assert medians["pooled"] == (1.331, 1.423)
    assert medians["reference"] == (1.075, 1.108)


def test_the_finiteness_guard_holds_under_every_generator(bench_tex):
    """Six expression models, not five, and 200/200 in all of them."""
    _quotes(bench_tex, "under six expression models")
    _quotes(bench_tex, "over three seeds")
    frame = pd.read_parquet(_newest("misspecified_generator_refusals"))
    assert frame["seed"].nunique() == 3
    assert frame["expression_model"].nunique() == 6

    killed = frame[frame["world"] == "annihilated"]
    for gate in ("count_gate", "width_gate"):
        block = killed[killed["gate"] == gate]
        assert (block["n_refused"] == 200).all(), gate
    ablation = killed[killed["gate"] == "no_gate"]
    assert (ablation["n_returned_a_number"] == 200).all()


def test_the_width_gate_binding_is_a_property_of_the_count_model(bench_tex):
    """``refuses 90 to 99 ... and 83 to 105`` and the width gate's 0 / 1 / 14.

    The paper said 95--99 across "all five models". There are six, and two of
    them sit below 95. The count gate's range and the width gate's per-model
    refusals are both asserted here because the contrast between them is the
    paragraph's whole claim.
    """
    _quotes(bench_tex, "the count gate refuses 90 to 99 of 200 in every model")
    _quotes(bench_tex, "and 83 to 105 across the eighteen")
    _quotes(bench_tex, "one under zero-inflated NB; and fourteen under NB at")
    frame = pd.read_parquet(_newest("misspecified_generator_refusals"))
    wide = frame[frame["world"] == "depleted_wide"]

    counts = wide[wide["gate"] == "count_gate"]
    per_model = counts.groupby("expression_model")["n_refused"].mean()
    assert (per_model.round() >= 90).all() and (per_model.round() <= 99).all()
    assert len(counts) == 18
    assert (int(counts["n_refused"].min()), int(counts["n_refused"].max())) == (83, 105)

    widths = wide[wide["gate"] == "width_gate"].groupby("expression_model")["n_refused"]
    means = widths.mean().round().astype(int).to_dict()
    assert {m for m, v in means.items() if v == 0} == {
        "poisson",
        "nb_disp10",
        "nb_disp2",
        "zip_pi0.3",
    }
    assert means["zinb_pi0.3_disp2"] == 1
    assert means["nb_disp0.5"] == 14


def test_the_knife_edge_window_and_the_cliff_past_it(bench_tex):
    """``[0.928776, 0.928981]``, ``0.958``, ``0.961`` and the 0.0027 headroom.

    The paper had the gate collapsing into two worlds at 0.958. Only one of
    them is refusing there; the second starts at 0.961.
    """
    _quotes(bench_tex, r"\in [0.928776, 0.928981]$")
    _quotes(bench_tex, "a window $0.000204$")
    _quotes(bench_tex, "refusing in \\texttt{depleted\\_estimable}, where the estimand")
    _quotes(bench_tex, "following at $0.961$")
    _quotes(bench_tex, "is $0.0027$ wide")

    step = pd.read_parquet(_newest("width_gate_step_location"))
    step = step.set_index("world")
    match = step.loc["depleted_wide"]
    assert round(float(match["s_detect_matching_count_gate_lo"]), 6) == 0.928776
    assert round(float(match["s_detect_matching_count_gate_hi"]), 6) == 0.928981
    assert round(float(match["matching_window_width"]), 6) == 0.000204
    # the first world where the estimand plainly exists starts refusing here
    assert round(float(match["collapse_onset_s_detect"]), 3) == 0.958
    assert round(float(step.loc["depleted_estimable", "s_star_min"]), 3) == 0.958
    assert round(float(step.loc["compositional_only", "s_star_min"]), 3) == 0.961
    assert round(float(match["headroom_below_collapse"]), 4) == 0.0027


def test_the_closed_form_is_confirmed_but_is_not_a_floor(tex):
    """``within $2.7$ percentage points`` and ``above ... in 14 of the 20``.

    The paper claimed agreement within 1.5 points with the excess positive
    everywhere, and called the closed form a floor. The worst gap is 2.65
    points and six cells fall below it, so the floor claim does not hold.
    """
    _quotes(tex, "to within $2.7$ percentage points in every one of twenty")
    _quotes(tex, "above} the closed form in 14 of the 20")
    _quotes(tex, "and not as a floor")
    frame = pd.read_parquet(_newest("interval_calibration"))
    pct = frame[frame["method"] == "percentile"]
    assert len(pct) == 20
    gap = 100 * (pct["false_positive_rate"] - pct["closed_form_rate"])
    assert gap.abs().max() < 2.7
    assert int((gap > 0).sum()) == 14
    assert int((gap < 0).sum()) == 6      # the floor claim, refuted


def test_the_censored_sweep_reports_what_it_dropped():
    """``463 of 1,200`` at 78% censoring and n = 100, and nothing below 49%."""
    blind = (SECTIONS / "blind.tex").read_text(encoding="utf-8")
    _quotes(blind, "$463$ of $1{,}200$ in")
    _quotes(blind, "the worst cell at each level runs $0$, $0.473$, $1.089$, $1.586$")
    head = pd.read_parquet(_newest("trial_survival_headline"))
    worst = head.loc[head["n_replicates"].idxmin()]
    assert (worst["censoring_target"], worst["n_patients"]) == (0.78, 100)
    assert 1200 - int(worst["n_replicates"]) == 463
    assert head[head["censoring_target"] <= 0.12]["n_replicates"].min() >= 1183

    survival = pd.read_parquet(_newest("trial_survival"))
    blind_arm = survival[survival["estimator"] == "exponential-mle-standardised"]
    worst_by_level = (
        blind_arm.groupby("censoring_target")["max_residual_vs_latent"].max().round(3)
    )
    assert list(worst_by_level) == [0.0, 0.473, 1.089, 1.586]
