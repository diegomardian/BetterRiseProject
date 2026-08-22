"""Every third-party import in src/ must be installable by CI.

This exists because the same bug landed three times in two days: scipy
(src/harness/deconvolve/nnls.py), statsmodels (src/estimator/hierarchical.py)
and scikit-learn (src/harness/deconvolve/nusvr.py) were each pinned in a
workstream's conda env, imported by code a test exercises, and absent from the
``dev`` extra that CI installs. Each one turned main red.

The rule was written down in pyproject.toml and CONTRIBUTING §6 and it still
happened, because the rule relies on remembering. This test does not.

It walks the AST rather than importing, so it catches imports nested inside
functions — the lazy-import case, which is the sneakier one: collection
succeeds, so the failure reads as "4 failed" rather than "1 error" and looks
like a broken test rather than a missing package.
"""

from __future__ import annotations

import ast
import importlib.util
import sys

import pytest

from src.common.paths import REPO_ROOT

SRC = REPO_ROOT / "src"

#: Optional by design. Every entry needs a reason, and the reason has to be one
#: of the two the docstring below allows: guarded by try/except, or reached only
#: by a code path no test exercises.
INTENTIONALLY_OPTIONAL = {
    # try/except in provenance.set_global_seeds — seeds torch if present.
    "torch",
    # src/reference/ingest.py:read_10x_mtx, a thin wrapper around
    # sc.read_10x_mtx. Lazily imported, no test covers it (the parsing it wraps
    # is dataset-specific and lands in W1's ingest script), and scanpy pulls
    # numba and llvmlite, which is too much to compile on every CI run. It is
    # pinned in env/w1_reference.yml, where a heavy single-workstream
    # dependency belongs.
    #
    # If a test ever does exercise read_10x_mtx, delete this entry and add
    # scanpy to the dev extra instead — at that point the compile cost buys
    # something.
    "scanpy",
}


def _top_level(name: str) -> str:
    return name.split(".")[0]


def _imports_in(path) -> set[str]:
    """Every module imported anywhere in the file, at any nesting depth."""
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    found: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            found.update(_top_level(a.name) for a in node.names)
        elif isinstance(node, ast.ImportFrom):
            if node.level == 0 and node.module:  # skip relative imports
                found.add(_top_level(node.module))
    return found


def _third_party_imports() -> dict[str, list[str]]:
    """module -> the src/ files that import it, for third-party modules only."""
    out: dict[str, list[str]] = {}
    for path in sorted(SRC.rglob("*.py")):
        for name in _imports_in(path):
            if name in sys.stdlib_module_names or name == "src":
                continue
            out.setdefault(name, []).append(str(path.relative_to(REPO_ROOT)))
    return out


THIRD_PARTY = _third_party_imports()


def test_src_actually_has_third_party_imports():
    """Guard against the walk silently finding nothing and passing vacuously."""
    assert len(THIRD_PARTY) >= 3, f"suspiciously few third-party imports: {THIRD_PARTY}"
    assert "numpy" in THIRD_PARTY


@pytest.mark.parametrize("module", sorted(set(THIRD_PARTY) - INTENTIONALLY_OPTIONAL))
def test_third_party_import_is_installable(module):
    """Fails in CI's thin env exactly when a dependency was not declared.

    Fix by adding the package to ``[project.optional-dependencies] dev`` in
    pyproject.toml when it is light and pip-installable. If it is heavy or
    non-pip (CellBender, inferCNV, anything from r-base), it does not belong in
    ``src/`` on a path a test exercises — move the call behind an explicit
    availability check like ``Deconvolver.is_available()`` and add the module to
    INTENTIONALLY_OPTIONAL here, with a reason.
    """
    assert importlib.util.find_spec(module) is not None, (
        f"{module!r} is imported by {THIRD_PARTY[module]} but is not installed. "
        f"CI installs only `pip install -e '.[dev]'`, so this turns main red. "
        f"See the docstring for the fix."
    )
