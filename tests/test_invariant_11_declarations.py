"""Every job written since invariant 11 landed declares its provenance.

Invariant 11 was merged 2026-09-07. Jobs predating it are not retrofitted here
and are tracked in ``docs/HANDOFF.md``; this test exists so that the jobs
written after it cannot quietly drop back to "unstated", which the invariant
refuses with the same force as a circular declaration.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

JOBS = Path("src/reference/jobs")

# Written after invariant 11 merged. Adding a job here is how it stays enforced.
POST_INVARIANT_11_JOBS = (
    "becker_lesion_wnt_gate.py",
    "chen_subtype_attrition.py",
    "chen_subtype_id_provenance.py",
    "chen_subtype_snapshot.py",
    "he_molecular_crdc_header.py",
    "he_molecular_gate.py",
    "he_molecular_synapse_access.py",
    "mhist_feasibility.py",
    "release7_early_lesion_identity.py",
    "zheng_gradient.py",
)


@pytest.mark.parametrize("job", POST_INVARIANT_11_JOBS)
def test_the_job_declares_both_provenances(job: str) -> None:
    source = (JOBS / job).read_text()
    assert "LABEL_PROVENANCE" in source, f"{job} does not declare label provenance"
    assert "CLAIM_PROVENANCE" in source, f"{job} does not declare claim provenance"


@pytest.mark.parametrize("job", POST_INVARIANT_11_JOBS)
def test_the_job_calls_the_guard(job: str) -> None:
    """The guard runs, and not only in a helper nothing reaches."""
    tree = ast.parse((JOBS / job).read_text())
    called = {
        node.func.id
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
    }
    assert "check_no_circular_claim" in called, f"{job} never calls the guard"


@pytest.mark.parametrize("job", POST_INVARIANT_11_JOBS)
def test_the_declaration_reaches_the_sidecar(job: str) -> None:
    source = (JOBS / job).read_text()
    assert "provenance_meta(" in source, (
        f"{job} declares provenance but never writes it beside the result, "
        "so a reader cannot check the population definition without the job"
    )
