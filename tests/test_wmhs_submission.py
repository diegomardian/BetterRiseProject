"""Submission-level checks over the files LaTeX actually reads."""

from __future__ import annotations

import json
import re
import subprocess
from pathlib import Path

import pandas as pd

from src.common.paths import REPO_ROOT

PAPER = REPO_ROOT / "paper" / "wmhs"
MANIFEST = PAPER / "results_manifest.json"
INPUT = re.compile(r"\\(?:input|include)\{([^}]+)\}")
GRAPHIC = re.compile(r"\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}")


def _with_tex_suffix(path: Path) -> Path:
    return path if path.suffix else path.with_suffix(".tex")


def latex_input_graph(entrypoint: Path) -> list[Path]:
    """Depth-first list of an entrypoint and every recursive TeX input."""
    seen: set[Path] = set()
    ordered: list[Path] = []

    def visit(path: Path) -> None:
        path = path.resolve()
        assert path.is_file(), f"missing LaTeX input: {path}"
        if path in seen:
            return
        seen.add(path)
        ordered.append(path)
        for value in INPUT.findall(path.read_text(encoding="utf-8")):
            visit(_with_tex_suffix(PAPER / value))

    visit(entrypoint)
    return ordered


def test_full_paper_verification_follows_its_actual_input_graph():
    graph = latex_input_graph(PAPER / "main.tex")
    relative = {path.relative_to(PAPER).as_posix() for path in graph}
    assert "sections/full/abstract.tex" in relative
    assert "sections/full/appendix.tex" in relative
    assert "sections/full/methods.tex" in relative
    assert "sections/full/trialtable.tex" in relative
    assert "sections/full/benchtable.tex" in relative
    assert not any(path.startswith("sections/") and "/full/" not in path for path in relative)


def _full_text() -> str:
    return "\n".join(
        path.read_text(encoding="utf-8")
        for path in latex_input_graph(PAPER / "main.tex")
    )


def _pinned(name: str) -> Path:
    entries = json.loads(MANIFEST.read_text(encoding="utf-8"))
    assert name in entries
    path = REPO_ROOT / entries[name]
    assert path.is_file()
    return path


def test_full_paper_survival_number_is_rederived_from_the_pinned_table():
    text = _full_text()
    assert r"maximum residual $2.8\times10^{-15}$ on valid outputs" in text
    survival = pd.read_parquet(_pinned("trial_survival"))
    km = survival[survival["estimator"] == "km-rmst-standardised"]
    assert len(km) == 28
    assert f"{km['max_residual_vs_latent'].max():.1e}" == "2.8e-15"


def test_full_paper_candidate_count_is_rederived_from_the_pinned_table():
    text = _full_text()
    assert "returns a candidate of 90 on seven of eight seeds" in text
    cutpoints = pd.read_parquet(_pinned("calibration_gap_cutpoints_r500"))
    cell = cutpoints[(cutpoints["pool"] == "reference") & (cutpoints["grid"] == "extended")]
    assert len(cell) == 8
    assert int((cell["ok"] == 90).sum()) == 7


def test_every_full_paper_graphic_exists():
    missing = []
    for source in latex_input_graph(PAPER / "main.tex"):
        for value in GRAPHIC.findall(source.read_text(encoding="utf-8")):
            candidate = PAPER / value
            choices = [candidate] if candidate.suffix else [candidate.with_suffix(".pdf")]
            if not any(path.is_file() for path in choices):
                missing.append((source.relative_to(PAPER).as_posix(), value))
    assert not missing, f"missing figures referenced by the full paper: {missing}"


def test_submission_result_manifest_is_complete_and_git_tracked():
    entries = json.loads(MANIFEST.read_text(encoding="utf-8"))
    assert entries, "submission result manifest is empty"
    assert len(entries) == len(set(entries.values())), "two result names pin the same file"

    missing = [value for value in entries.values() if not (REPO_ROOT / value).is_file()]
    assert not missing, f"pinned submission results are missing: {missing}"

    tracked = set(
        subprocess.run(
            ["git", "ls-files", "results"],
            cwd=REPO_ROOT,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.splitlines()
    )
    untracked = [value for value in entries.values() if value not in tracked]
    assert not untracked, f"pinned submission results are not committed: {untracked}"


def test_figure_readers_do_not_select_results_by_directory_order():
    source = (PAPER / "_tables.py").read_text(encoding="utf-8")
    assert "results_manifest.json" in source
    assert ".glob(" not in source
    assert "sorted(" not in source


def _quotes(text: str, literal: str) -> None:
    """Assert the submitted paper contains this literal, before checking it.

    Without this, a reworded passage silently stops being checked -- which is
    the failure this file exists to prevent, one level out.
    """
    assert literal in text, (
        f"the submitted paper no longer contains {literal!r}. Either the prose "
        f"moved, in which case update this test, or a number changed without "
        f"its table."
    )


def test_the_figure_one_caption_partitions_its_replicates():
    """Non-null + undefined must equal the sweep total, in both figures.

    An external reviewer found this caption claiming 119,600 non-null replicates
    and "a further" 29,900 undefined, in a sweep containing 119,600 replicates
    in total -- 119,600 + 29,900 = 149,500. The numbers were a total and one of
    its parts, described as two disjoint sets. Nothing checked it because the
    number tests point at the other build's sources.
    """
    text = _full_text()
    one = pd.read_parquet(_pinned("calibration_gap_recovery"))
    total = int(one["n_replicates"].sum())
    undefined = int(one["n_ratio_undefined"].sum())
    non_null = total - undefined

    assert (total, undefined, non_null) == (119_600, 29_900, 89_700)
    _quotes(text, "89{,}700 non-null")
    _quotes(text, "29{,}900 replicates of this")
    _quotes(text, "sweep's 119{,}600")
    # the partition, stated as the paper states it
    assert non_null + undefined == total

    # and the four-sweep figures, which must partition the same way
    _quotes(text, "1{,}283{,}400 across four sweeps")
    _quotes(text, "427{,}800 of 1{,}711{,}200")
    assert 1_283_400 + 427_800 == 1_711_200


def test_the_residual_matrix_headline_is_rederived_for_the_submitted_build():
    """The reference-dependence result, now promoted to its own section."""
    text = _full_text()
    frame = pd.read_parquet(_pinned("trial_blindness_residual_matrix"))

    def cell(estimator: str, truth: str) -> float:
        rows = frame[
            (frame["estimator"] == estimator)
            & (frame["truth"] == truth)
            & (frame["design"] == "confounded-bernoulli")
            & (frame["n_patients"] == 5000)
        ]
        assert len(rows) == 1, (estimator, truth, len(rows))
        return float(rows["max_residual"].iloc[0])

    _quotes(text, "OLS has maximum residual $0.0663$ against the standardised")
    _quotes(text, "$3.5" + chr(92) + "times10^{-14}$ against the")
    assert round(cell("ols-stratum-dummies", "standardised"), 4) == 0.0663
    assert round(cell("ols-stratum-dummies", "varweighted") * 1e14, 1) == 3.5
    # the swap, which is the point of the section
    assert cell("gcomp-from-generator", "standardised") == 0.0
    assert round(cell("gcomp-from-generator", "varweighted"), 4) == 0.0663


def test_clean_control_table_matches_all_eight_pinned_rows():
    text = _full_text()
    frame = pd.read_parquet(_pinned("residual_performance_clean_control"))
    labels = {
        "empirical-mean-calibrated": "Full-sample mean",
        "same-point-narrow-interval": "Same point, narrow interval",
        "known-bias-plus-0.5": "Mean plus $0.5$",
        "independent-half-mean": "Half-sample mean",
    }
    assert len(frame) == 8
    assert set(frame["requested_effect"]) == {0.0, 0.5}
    assert (frame["n_attempted"] == 2000).all()
    assert (frame["n_valid"] == frame["n_attempted"]).all()
    for row in frame.itertuples(index=False):
        residual = (
            "0" if row.max_residual_vs_reference == 0
            else f"{row.max_residual_vs_reference:.4f}"
        )
        literal = (
            f"{labels[row.estimator]} & ${residual}$ & ${row.bias:.4f}$ & "
            f"${row.rmse:.4f}$ & ${100 * row.interval_coverage:.2f}"
            + r"\%$"
        )
        _quotes(text, literal)
    _quotes(text, "at most 1.08 percentage points")
    assert round(100 * frame["coverage_mc_se"].max(), 2) == 1.08


def test_independent_cox_agreement_claim_uses_the_pinned_check():
    text = _full_text()
    frame = pd.read_parquet(_pinned("trial_survival_lifelines_check"))
    assert len(frame) == 8
    assert frame["valid_pair"].all()
    assert frame["within_fixed_tolerance"].all()
    assert (frame["n_patients"] == 1500).all()
    assert (frame["fixed_tolerance"] == 1e-5).all()
    assert f"{frame['absolute_difference'].max():.1e}" == "1.4e-07"
    _quotes(text, r"largest absolute difference is $1.4\times10^{-7}$")


def test_per_count_candidates_and_unobserved_omission_boundary():
    text = _full_text()
    frame = pd.read_parquet(_pinned("controlled_grid_crossings_r200_b200"))
    smc = frame[
        (frame["cohort"] == "smc")
        & (frame["pool"] == "reference")
        & (frame["criterion"] == "coverage_and_discrimination")
        & (frame["binning"] == "per_count")
    ]
    filled = smc[smc["grid"] == "extended"]["candidate"].value_counts().to_dict()
    dense = smc[smc["grid"] == "dense"]["candidate"].value_counts().to_dict()
    assert filled == {60.0: 5, 50.0: 2, 80.0: 1}
    assert dense == {60.0: 4, 50.0: 2, 45.0: 1, 80.0: 1}
    _quotes(text, "count is 50 on two seeds, 60 on five and 80 on one")
    influence = pd.read_parquet(_pinned("calibration_lopo_influence_b200-1000-5000"))
    joint = influence[
        (influence["cohort"] == "kul3")
        & (influence["criterion"] == "coverage_and_discrimination")
    ].set_index("omitted_patient")
    assert joint.loc["KUL01", "omitted_candidate"] == 300
    assert joint.loc["KUL30", "omitted_candidate"] == 300
    assert joint.loc["KUL31", "omitted_candidate"] == 200
    assert joint.loc["KUL31", "omitted_status"] == "lower_bound_unobserved"
    _quotes(text, "KUL31 moves it to the lowest tested count, 200, with its lower boundary")


def test_the_learned_generator_result_is_rederived():
    """Bitwise zero, the ridge ladder, and the n^-1.75 trend."""
    text = _full_text()
    primary = pd.read_parquet(_pinned("trial_learned_generator_primary"))
    draw = primary[primary["reference"] == "T_draw"]

    # the headline: exact zero, and it must be exact, not rounded
    sat = draw[draw["estimator"] == "gcomp-saturated"]
    assert len(sat) == 42
    assert (sat["max_residual"] == 0.0).all()
    learned = sat[sat["is_learned"]]
    assert len(learned) == 30 and (learned["max_residual"] == 0.0).all()
    _quotes(text, "bitwise")
    _quotes(text, "all 30 learned-generator")

    # a zero over an empty set would be vacuous -- this is the paper's own defect
    assert (draw["n_valid_pairs"] > 0).all()

    # the ridge ladder, as quoted
    ridge = draw[draw["estimator"].str.contains("ridge")]
    top = ridge.groupby("estimator")["max_residual"].max()
    assert round(top["gcomp-ridge-outcome-a1e-5"], 6) == 0.000086
    assert round(top["gcomp-ridge-outcome-a1e-3"], 5) == 0.00864
    assert round(top["gcomp-ridge-outcome-a1e-1"], 2) == 0.78
    _quotes(text, "$8.6" + chr(92) + "times10^{-5}$")

    # more data moves penalised arms INTO the blind band
    trend = pd.read_parquet(_pinned("trial_learned_generator_trend"))
    tr = trend[(trend["reference"] == "T_draw")
               & trend["estimator"].str.contains("ridge")]
    lo, hi = 1e-8, 3.0e-5
    crossed = tr[(tr["max_residual_at_smallest_n"] >= hi)
                 & (tr["max_residual_at_largest_n"] > lo)
                 & (tr["max_residual_at_largest_n"] < hi)]
    assert len(crossed) == 14, len(crossed)
    _quotes(text, "fourteen arms")
    # the paper now reports the fitted spread, not a single exponent, and
    # explicitly declines to derive it -- so pin the bounds it quotes
    slopes = tr["log10_slope_vs_log10_n"]
    assert round(abs(slopes.median()), 2) == 1.75, slopes.median()
    assert round(slopes.min(), 2) == -1.87 and round(slopes.max(), 2) == -1.52
    _quotes(text, "fitted slopes of $-1.52$ to")
    _quotes(text, "median $-1.75$")
    _quotes(text, "do not derive it")


def test_the_prevalence_audit_is_rederived():
    """0 of 7, and the false positive the rule caught."""
    text = _full_text()
    verdicts = pd.read_parquet(_pinned("prevalence_audit_verdicts"))
    assert len(verdicts) == 7
    counts = verdicts["verdict"].value_counts().to_dict()
    assert counts.get("YES", 0) == 0
    assert counts.get("NO", 0) == 7
    _quotes(text, "None reuses the functional: 0 of 7.")

    # exactly one repository reached the functional comparison
    executed = verdicts[verdicts["max_residual"].notna()]
    assert len(executed) == 1
    assert round(float(executed["max_residual"].iloc[0]) * 1e14, 1) == 4.0
    _quotes(text, "$4.0" + chr(92) + "times10^{-14}$")


#: A mangled control sequence leaves one of these behind. LaTeX raises no error
#: for any of them -- it typesets the fragment as body text -- so a build log is
#: not a check. One reached a reviewer's PDF as literal "extbf{None reuses...".
_ORPHANS = ("extbf{", "extit{", "exttt{", "ef{", "imes", "elax", "aragraph{")


def test_no_mangled_control_sequences_in_the_submitted_sources():
    """No stray control byte, and no macro that lost its backslash.

    Read as BYTES. Python's universal-newline text mode silently converts a
    lone carriage return to a newline, which hides exactly the damage this
    looks for -- the defect and its concealment in the same call.
    """
    import re as _re
    bs, tab, cr, ff = chr(92), chr(9), chr(13), chr(12)
    offenders = []
    for path in latex_input_graph(PAPER / "main.tex"):
        raw = path.read_bytes().decode("utf-8")
        name = path.name
        if tab in raw:
            offenders.append(f"{name}: literal TAB")
        if raw.replace(cr + chr(10), "").count(cr):
            offenders.append(f"{name}: lone CR")
        if ff in raw:
            offenders.append(f"{name}: form feed")
        for frag in _ORPHANS:
            # a fragment not preceded by a backslash or a letter is orphaned
            pattern = "(?<![A-Za-z" + _re.escape(bs) + "])" + _re.escape(frag)
            for m in _re.finditer(pattern, raw):
                offenders.append(f"{name}: orphaned {frag!r} at {m.start()}")
    assert not offenders, (
        "mangled control sequences in the submitted sources: "
        + "; ".join(offenders[:8])
    )
