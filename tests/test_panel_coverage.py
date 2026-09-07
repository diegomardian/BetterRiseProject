"""C1's panel gate must require every biological role, not just the target."""

import pandas as pd

from src.reference.jobs.panel_coverage import verdict


def _row(panel: str, gene: str, role: str, status: str) -> dict[str, str]:
    return {"panel": panel, "gene": gene, "role": role, "status": status}


def test_target_only_panel_is_not_runnable():
    table = pd.DataFrame([
        _row("target_only", "GUCA2A", "target", "present"),
        _row("target_only", "CDX2", "identity", "absent"),
        _row("target_only", "KRT8", "control", "absent"),
    ])

    outcome = verdict(table)

    assert outcome["verdict"] == "C1 NOT RUNNABLE ON A STOCK PANEL"
    assert "identity, control" in outcome["detail"]


def test_panel_with_target_identity_and_control_is_runnable():
    table = pd.DataFrame([
        _row("complete", "GUCA2A", "target", "present"),
        _row("complete", "CDX2", "identity", "present"),
        _row("complete", "KRT8", "control", "present"),
    ])

    outcome = verdict(table)

    assert outcome["verdict"] == "C1 RUNNABLE ON STOCK PANEL(S)"
    assert "complete" in outcome["detail"]
