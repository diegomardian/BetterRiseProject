"""Git sha, seed and environment stamping. CLAUDE.md invariant 10.

Every result carries the git sha and a fixed random seed. Not "usually" — the
writer in src/schema.py refuses to emit a table without both.
"""

from __future__ import annotations

import datetime as _dt
import platform
import subprocess
import sys
from typing import Any

from src.common.paths import REPO_ROOT

#: Default project seed. Use it, or pass your own — but pass it explicitly.
DEFAULT_SEED: int = 20260101


def _git(*args: str) -> str | None:
    try:
        out = subprocess.run(
            ["git", *args],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=True,
            timeout=15,
        )
        return out.stdout.strip()
    except (subprocess.SubprocessError, FileNotFoundError, OSError):
        return None


def git_sha() -> str:
    """Full commit sha, or ``"unknown"`` outside a git checkout."""
    return _git("rev-parse", "HEAD") or "unknown"


#: Untracked paths that never make a tree dirty.
#:
#: Caches and editor droppings, plus ``results/`` — and that last one is a real
#: exemption rather than an oversight. This guard exists to catch untracked
#: *code*: a producing script the recorded sha does not contain. A results table
#: is the output being stamped, and every committed result passed through an
#: untracked state on its way in, so a job writing two tables would otherwise
#: refuse its own second write. A guard that fires on the thing it is meant to
#: certify gets switched off, and then it is not a guard.
_IGNORABLE_UNTRACKED: tuple[str, ...] = (
    ".DS_Store",
    "__pycache__",
    ".pytest_cache",
    ".ruff_cache",
    "results",
)


def _is_ignorable(path: str) -> bool:
    """Whether an untracked path is exempt.

    Matched on whole path *components*, not as a substring: ``results`` must not
    exempt ``my_results_script.py``. A loose match here would quietly widen the
    exemption every time somebody named a file conveniently.
    """
    return any(part in _IGNORABLE_UNTRACKED for part in path.split("/"))


def untracked_files() -> list[str]:
    """Untracked paths git can see, minus caches and editor droppings.

    Separate from ``git_is_dirty`` so the sidecar can *name* them. A table
    stamped with a clean sha whose producing script is untracked is the failure
    this exists to make visible: the sha points at a commit that does not
    contain the code that ran.
    """
    status = _git("status", "--porcelain", "--untracked-files=all")
    if not status:
        return []
    out = []
    for line in status.splitlines():
        # Porcelain v1: two status columns, a space, then the path.
        if line[:2] != "??":
            continue
        path = line[3:].strip().strip('"')
        if not _is_ignorable(path):
            out.append(path)
    return sorted(out)


def git_is_dirty() -> bool:
    """True if the working tree differs from HEAD, untracked files included.

    Untracked counts. This check used to pass ``--untracked-files=no``, so a
    table produced by a script that was never committed recorded a *clean* sha
    pointing at a commit that did not contain the script. Three results reached
    an internal report that way. The sha is a claim that the commit reproduces
    the table; an untracked producer falsifies it.

    Unknown git state counts as dirty.
    """
    status = _git("status", "--porcelain", "--untracked-files=no")
    if status is None:
        return True
    return bool(status) or bool(untracked_files())


def git_branch() -> str:
    return _git("rev-parse", "--abbrev-ref", "HEAD") or "unknown"


def set_global_seeds(seed: int = DEFAULT_SEED) -> int:
    """Seed python, numpy and (if present) torch. Returns the seed, for logging.

    Convenience only. Prefer an explicit ``numpy.random.default_rng(seed)``
    passed down the call stack — global state and patient-level bootstrap are a
    bad combination.
    """
    import random

    random.seed(seed)
    try:
        import numpy as np

        np.random.seed(seed)
    except ImportError:
        pass
    try:
        import torch

        torch.manual_seed(seed)
    except ImportError:
        pass
    return seed


#: Packages whose version is stamped onto every result. Import name == install
#: name for all of these; add a mapping here if that ever stops being true.
_STAMPED_PACKAGES: tuple[str, ...] = ("numpy", "pandas", "scanpy", "anndata")


def _package_version(name: str) -> str | None:
    """Installed version of ``name``, or None if it is not installed.

    Reads distribution metadata rather than importing the package and asking for
    ``__version__``. Three reasons, in order of how much trouble they caused:

    1. ``__version__`` is being deprecated. anndata >= 0.13 emits a
       FutureWarning on attribute access, and this project turns FutureWarning
       into an error (``filterwarnings`` in pyproject.toml), so the old code took
       ``write_results`` down with it the moment anndata was installed — a
       results writer failing because of how it stamps a version number.
    2. It does not import anything. Stamping provenance should not execute a
       third-party package's import side effects, and scanpy's are not cheap.
    3. It works for packages that never defined ``__version__`` at all.

    Diagnosed by W1 in tests/test_dependencies.py while working around it; the
    fix belongs here because src/common/ is shared code.
    """
    from importlib.metadata import PackageNotFoundError, version

    try:
        return version(name)
    except PackageNotFoundError:
        return None


def provenance_record(*, seed: int, notes: str | None = None) -> dict[str, Any]:
    """The dict written alongside every results parquet."""
    sha = git_sha()
    record: dict[str, Any] = {
        "date": _dt.date.today().isoformat(),
        "utc_timestamp": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
        "git_sha": sha,
        "git_sha_short": sha[:7],
        "git_branch": git_branch(),
        "git_dirty": git_is_dirty(),
        "git_untracked": untracked_files(),
        "seed": int(seed),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
    }
    for pkg in _STAMPED_PACKAGES:
        version = _package_version(pkg)
        if version is not None:
            record[f"{pkg}_version"] = version
    if notes:
        record["notes"] = notes
    return record
