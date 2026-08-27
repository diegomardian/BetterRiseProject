"""Harness table shapes and their writer. W2.

``src/schema.py`` is the frozen per-(patient, gene, rung, axis, study) contract
for decomposition results. Harness output — a bake-off ranking, an attenuation
sweep, a calibration report — has a different shape entirely and must not be
forced into it. These tables are W2-local: validated here, not frozen, and
changeable in an ordinary ``w2/`` PR.

They are written through ``src.common.io.write_versioned_table``, so they land
under the same ``results/{date}_{sha7}/`` convention and carry the same git sha
and seed as everything else (CLAUDE.md invariant 10).
"""

from __future__ import annotations

from pathlib import Path
from typing import Final

import pandas as pd

from src.common.io import write_versioned_table

#: One row per (grid point, replicate, arm). The attenuation sweep's raw output.
#: ``arm`` is "oracle" (estimator on cell-level summaries — the reliable half)
#: or "bulk" (recovered through deconvolution — the attenuated half that
#: invariant 6 forbids using for results and that §2.2 exists to measure).
ATTENUATION_COLUMNS: Final = (
    "grid_id",
    "replicate",
    "arm",
    "gene",
    "weighting",
    "frac_mature_tumour",
    "shift",
    "n_cells_mature",
    "estimability",
    "compositional_true_parametric",
    "intrinsic_true_parametric",
    "compositional_true_realised",
    "intrinsic_true_realised",
    "compositional_hat",
    "intrinsic_hat",
    "interaction_hat",
    "attenuation_ratio",
    "ci_low",
    "ci_high",
    "seed",
)

#: One row per (method, sample). Fraction recovery against harness truth.
BAKEOFF_COLUMNS: Final = (
    "method",
    "sample_id",
    "cell_type",
    "fraction_true",
    "fraction_hat",
    "n_signature_genes",
    "seed",
)

#: One row per candidate cutpoint. The evidence behind the week-5 swap.
CALIBRATION_COLUMNS: Final = (
    "n_cells_mature",
    "shift",
    "coverage",
    "discrimination",
    "median_ci_width",
    "n_replicates",
    "verdict",
)

#: One row per (control, gene, term). Both terms should be ~0 throughout.
CONTROLS_COLUMNS: Final = (
    "control",
    "gene",
    "term",
    "weighting",
    "value",
    "ci_low",
    "ci_high",
    "seed",
)

#: One row per (world, replicate). W2's ratification of prereg Amendment 2
#: (issue #37): five simulated worlds with known answers, and the verdict the
#: proposed G1 returns on each. ``rho_decision_17`` is the superseded statistic,
#: carried alongside so the two are read against the same truth.
G1_AMENDMENT_COLUMNS: Final = (
    "world",
    "replicate",
    "verdict",
    "expected_pass_low",
    "expected_pass_high",
    "rho_decision_17",
    "rho_g1a",
    "tier_a_median_pct",
    "tier_d_pct",
    "separation",
    "n_genes_kept",
    "n_genes",
    "seed",
)

#: One row per (cohort, true fraction). The G4 operating characteristic at the
#: cohort sizes that exist — handoff §5 task 3. Deterministic (exact binomial);
#: the measured bootstrap band widths ride in the sidecar rather than here,
#: because they are medians over replicates and not one row per anything.
GATE_COST_COLUMNS: Final = (
    "cohort",
    "n_patients",
    "true_fraction_below",
    "p_gate_says_pass",
    "largest_clean_pass",
    "effective_decision_line",
)

#: One row per (regime, ambient fraction, replicate, term). Handoff §5 task 7 —
#: what residual ambient contamination does to a decomposition whose truth is
#: known. ``truth_is_zero`` comes from the regime's DESIGN, never from the
#: realised value; see AmbientRegime.parametric_zero_terms for why that matters.
AMBIENT_COLUMNS: Final = (
    "regime",
    "ambient_fraction",
    "replicate",
    "weighting",
    "term",
    "truth_is_zero",
    "value_clean",
    "value_contaminated",
    "n_cells_mature",
    "seed",
)

_SHAPES: Final = {
    "attenuation": ATTENUATION_COLUMNS,
    "bakeoff": BAKEOFF_COLUMNS,
    "calibration": CALIBRATION_COLUMNS,
    "controls": CONTROLS_COLUMNS,
    "g1_amendment": G1_AMENDMENT_COLUMNS,
    "gate_cost": GATE_COST_COLUMNS,
    "ambient": AMBIENT_COLUMNS,
}


class HarnessTableError(ValueError):
    """A harness table does not match its declared shape."""


def validate_harness_table(df: pd.DataFrame, kind: str) -> None:
    """Check a table against its declared columns.

    Lighter than the frozen schema's validation on purpose — these shapes are
    W2's to change. What it does catch is a column silently renamed or dropped
    between the sweep and the plot.
    """
    if kind not in _SHAPES:
        raise HarnessTableError(f"unknown table kind {kind!r}; known: {sorted(_SHAPES)}")
    expected = list(_SHAPES[kind])
    missing = [c for c in expected if c not in df.columns]
    if missing:
        raise HarnessTableError(f"{kind} table is missing columns: {missing}")
    extra = [c for c in df.columns if c not in expected]
    if extra:
        raise HarnessTableError(
            f"{kind} table has unexpected columns: {extra}. Add them to "
            f"{kind.upper()}_COLUMNS in src/harness/results.py if they belong."
        )


def empty_harness_table(kind: str) -> pd.DataFrame:
    """An empty frame with the right columns for ``kind``."""
    if kind not in _SHAPES:
        raise HarnessTableError(f"unknown table kind {kind!r}; known: {sorted(_SHAPES)}")
    return pd.DataFrame({c: pd.Series(dtype="object") for c in _SHAPES[kind]})


def write_harness_table(
    df: pd.DataFrame,
    kind: str,
    *,
    seed: int,
    name: str | None = None,
    results_dir: Path | None = None,
    extra_meta: dict | None = None,
    allow_dirty: bool = False,
) -> Path:
    """Validate then write a harness table.

    ``extra_meta`` should carry the grid definition, the replicate count, and —
    for a bake-off — which methods were skipped and why. That sidecar is the
    only place a reader finds out CIBERSORTx never ran.
    """
    validate_harness_table(df, kind)
    meta = {"harness_table_kind": kind} | dict(extra_meta or {})
    return write_versioned_table(
        df,
        name=name or f"harness_{kind}",
        seed=seed,
        results_dir=results_dir,
        extra_meta=meta,
        allow_dirty=allow_dirty,
    )
