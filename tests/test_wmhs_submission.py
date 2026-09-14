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
