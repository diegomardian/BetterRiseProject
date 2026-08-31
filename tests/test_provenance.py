"""The provenance stamp, tested against inputs that force it to fail.

The paper's working rule: a check unable to fail is worse than no check. The
dirtiness guard was exactly that for two months — it passed
``--untracked-files=no``, so a table produced by a script that had never been
committed recorded a *clean* sha pointing at a commit that did not contain the
script. Three results reached an internal report that way, one of them a
calibration behind a published figure.

So every test here builds a repository in a state the check must reject, and
asserts it does. Testing that a clean tree reads clean would have passed against
the broken version too.
"""

from __future__ import annotations

import subprocess

import pytest

from src.common import provenance


def _git(repo, *args: str) -> str:
    out = subprocess.run(
        ["git", *args], cwd=repo, capture_output=True, text=True, check=True
    )
    return out.stdout.strip()


@pytest.fixture()
def repo(tmp_path, monkeypatch):
    """A real git repository with one commit, pointed at by ``_git``."""
    root = tmp_path / "repo"
    root.mkdir()
    _git(root, "init", "-q")
    _git(root, "config", "user.email", "t@example.com")
    _git(root, "config", "user.name", "test")
    (root / "committed.py").write_text("x = 1\n")
    _git(root, "add", "committed.py")
    _git(root, "commit", "-qm", "first")
    monkeypatch.setattr(provenance, "REPO_ROOT", root)
    return root


def test_a_clean_tree_reads_clean(repo):
    """The baseline. Only meaningful next to the failure cases below."""
    assert provenance.git_is_dirty() is False
    assert provenance.untracked_files() == []


def test_an_untracked_producing_script_makes_the_tree_dirty(repo):
    """THE regression. This is the exact shape of the withdrawn guard.

    A script that writes a result table, never committed. The old check passed
    ``--untracked-files=no`` and returned False here, so the table got stamped
    with a sha that does not contain the script that produced it.
    """
    (repo / "run_the_sweep.py").write_text("# writes a result table\n")
    assert provenance.git_is_dirty() is True, (
        "an untracked producing script must make the tree dirty — this is the "
        "defect that stamped three results clean against a commit lacking them"
    )


def test_the_untracked_file_is_named_not_just_counted(repo):
    """The sidecar has to say *which* file, or the gap has to be reconstructed
    from commit timestamps afterwards. That reconstruction is what cost 83
    seconds of forensics the first time."""
    (repo / "run_the_sweep.py").write_text("# writes a result table\n")
    assert provenance.untracked_files() == ["run_the_sweep.py"]


def test_the_provenance_record_carries_the_untracked_names(repo):
    (repo / "run_the_sweep.py").write_text("# writes a result table\n")
    record = provenance.provenance_record(seed=1)
    assert record["git_dirty"] is True
    assert record["git_untracked"] == ["run_the_sweep.py"]


def test_untracked_files_in_subdirectories_are_named_individually(repo):
    """``--untracked-files=all``, not the default, which collapses a directory
    to its name and hides how many scripts are in it."""
    (repo / "jobs").mkdir()
    (repo / "jobs" / "a.py").write_text("")
    (repo / "jobs" / "b.py").write_text("")
    assert provenance.untracked_files() == ["jobs/a.py", "jobs/b.py"]


def test_a_modified_tracked_file_still_makes_the_tree_dirty(repo):
    """The property the old check did have. Keep it."""
    (repo / "committed.py").write_text("x = 2\n")
    assert provenance.git_is_dirty() is True


def test_caches_and_editor_droppings_do_not_make_a_tree_dirty(repo):
    """A guard that fires on ``.DS_Store`` gets switched off, and then it is
    not a guard. The exclusion list holds nothing that can produce a table."""
    (repo / ".DS_Store").write_bytes(b"\x00")
    (repo / "__pycache__").mkdir()
    (repo / "__pycache__" / "committed.cpython-313.pyc").write_bytes(b"\x00")
    assert provenance.untracked_files() == []
    assert provenance.git_is_dirty() is False


def test_a_real_script_next_to_ignorable_noise_is_still_caught(repo):
    """The exclusion list must not become a hole. Noise alongside a real
    untracked producer still leaves the tree dirty, and the producer named."""
    (repo / ".DS_Store").write_bytes(b"\x00")
    (repo / "run_the_sweep.py").write_text("")
    assert provenance.git_is_dirty() is True
    assert provenance.untracked_files() == ["run_the_sweep.py"]


def test_a_jobs_own_output_does_not_make_the_tree_dirty(repo):
    """A job writing two tables must not refuse its own second write.

    ``results/`` is exempt because this guard is about untracked *code*: every
    committed result passed through an untracked state on the way in. A guard
    that fires on the thing it certifies gets switched off."""
    (repo / "results" / "2026-08-31_abc1234").mkdir(parents=True)
    (repo / "results" / "2026-08-31_abc1234" / "first.parquet").write_bytes(b"\x00")
    assert provenance.untracked_files() == []
    assert provenance.git_is_dirty() is False


def test_the_results_exemption_matches_components_not_substrings(repo):
    """``results`` must not exempt ``my_results_helper.py``. A substring match
    would widen this exemption every time somebody named a file conveniently."""
    (repo / "my_results_helper.py").write_text("")
    assert provenance.untracked_files() == ["my_results_helper.py"]
    assert provenance.git_is_dirty() is True


def test_an_untracked_script_inside_results_is_still_exempt_but_code_is_not(repo):
    """The exemption is by location, and the sibling case is the one that
    matters: a producing script at the repo root is still caught."""
    (repo / "results").mkdir()
    (repo / "results" / "out.parquet").write_bytes(b"\x00")
    (repo / "sweep.py").write_text("")
    assert provenance.untracked_files() == ["sweep.py"]


def test_unknown_git_state_counts_as_dirty(monkeypatch):
    """Outside a checkout, or with git absent, the honest answer is "cannot
    claim this reproduces"."""
    monkeypatch.setattr(provenance, "_git", lambda *a: None)
    assert provenance.git_is_dirty() is True
