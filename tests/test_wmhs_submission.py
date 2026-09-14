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
