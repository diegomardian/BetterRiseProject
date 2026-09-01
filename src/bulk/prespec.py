"""Load and enforce the Stage 4 pre-specification.

The covariate lock (W3.6) works because `require_locked` refuses to fit against
a `proposed` config. This is the same mechanism for the other analysis in the
bulk arm that can be reported after the fact — execution_plan.md §6.2's variance
question, which PR #49 turned from an open question into a directional
prediction.

Reusing `covariates.require_locked` rather than writing a second one: it reads
only `status` and `version`, so it already generalises, and two lock functions
that could drift apart is worse than one used twice.

WHY A LOCK ON THIS ONE
----------------------
The prediction is that mature-colonocyte fraction explains **little** of bulk
GUCA2A and more of bulk CDX2. A null result is the expected result. An analysis
whose expected outcome is a null is the easiest kind to report after the fact
with the sign of the claim quietly reversed, and the hardest for a reader to
audit. Committing the direction, the threshold and the two instrument controls
before the numbers exist is the only thing that separates it from a story.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from src.bulk.covariates import CovariateError, require_locked
from src.common.paths import CONFIG_DIR

PRESPEC_PATH = CONFIG_DIR / "stage4_prespecification.yaml"

#: Keys without which the file is not a pre-specification, only a note.
REQUIRED_KEYS = (
    "version",
    "status",
    "question",
    "prediction",
    "instrument_checks",
    "matched_null",
    "model",
    "cohorts",
)


class PrespecError(RuntimeError):
    """The Stage 4 pre-specification is missing, malformed, or not yet locked."""


@lru_cache(maxsize=4)
def load_prespec(path: str | Path | None = None) -> dict[str, Any]:
    """Read the Stage 4 pre-specification. Cached; pass a path in tests."""
    target = Path(path) if path is not None else PRESPEC_PATH
    if not target.exists():
        raise PrespecError(
            f"{target} does not exist. Stage 4's variance question is "
            f"pre-specified before it is run (execution_plan.md §6.2)."
        )
    spec = yaml.safe_load(target.read_text(encoding="utf-8"))
    missing = [key for key in REQUIRED_KEYS if key not in spec]
    if missing:
        raise PrespecError(f"{target} is missing required key(s): {missing}")
    if spec["status"] not in ("proposed", "locked"):
        raise PrespecError(
            f"status must be 'proposed' or 'locked', got {spec['status']!r}"
        )
    for side in ("primary", "secondary"):
        if side not in spec["prediction"]:
            raise PrespecError(f"prediction has no {side!r} arm")
        if "disconfirmed_if" not in spec["prediction"][side]:
            raise PrespecError(
                f"prediction.{side} has no `disconfirmed_if`. A prediction that "
                f"cannot be wrong is not one."
            )
    # Issue #54: a raw R-squared, or two of them compared across genes at
    # different abundance, measures the assay floor. If an arm ever says
    # "R-squared" without saying "percentile" or "excess", the fix has been
    # undone.
    for side in ("primary", "secondary"):
        statement = spec["prediction"][side]["statement"].lower()
        if "r-squared" in statement and not (
            "percentile" in statement or "excess" in statement
        ):
            raise PrespecError(
                f"prediction.{side} compares a raw R-squared. R-squared is a "
                f"share of variance and is confounded with abundance at the "
                f"assay floor (issue #54); state it as a percentile within, or "
                f"an excess over, the abundance-matched null."
            )
    return spec


def require_locked_prespec(spec: dict[str, Any]) -> None:
    """Raise unless the pre-specification is locked. Call before analysing.

    Delegates to the covariate lock's `require_locked`, so the two cannot drift.
    """
    try:
        require_locked(spec)
    except CovariateError as exc:
        raise PrespecError(
            f"{exc}\n(This is the Stage 4 pre-specification, "
            f"config/stage4_prespecification.yaml.)"
        ) from exc


def prediction(spec: dict[str, Any], arm: str = "primary") -> dict[str, Any]:
    """One arm of the prediction, with its disconfirming condition attached."""
    if arm not in spec["prediction"]:
        raise PrespecError(f"unknown prediction arm {arm!r}")
    return spec["prediction"][arm]


def outcome_genes(spec: dict[str, Any]) -> list[str]:
    """The genes whose variance is being explained. §6.2 names two."""
    genes = list(spec["model"]["outcome_genes"])
    if not genes:
        raise PrespecError("no outcome genes in the pre-specification")
    return genes


def positive_control_gates_the_analysis(spec: dict[str, Any]) -> bool:
    """Does a failed instrument check stop the analysis? It must.

    A low R-squared means either 'fraction does not explain this gene' or
    'deconvolution does not work here', and those have opposite consequences.
    The positive control is what tells them apart, so it cannot be advisory.
    """
    on_failure = spec["instrument_checks"]["positive_control"]["on_failure"]
    return on_failure.strip().upper().startswith("STOP")


# ---------------------------------------------------------------------------
# The abundance-matched null — issue #54
# ---------------------------------------------------------------------------


def matched_null_genes(
    tumour_expression: pd.DataFrame,
    index_map: pd.DataFrame,
    target_gene_id: str,
    spec: dict[str, Any],
) -> list[str]:
    """Genes at the target's abundance, for the null its R-squared is judged against.

    Issue #54: R-squared is a share of variance, so a gene near the assay floor
    loses it to measurement noise whatever its biology. GUCA2A sits at 1.40% of
    normal and CDX2 at 94.7%, so comparing their R-squared compares their
    abundance. Simulated under PURE composition for both genes, the low-abundance
    one returns R-squared 0.891 -> 0.000 as the assay floor rises while the
    high-abundance one holds at ~0.86 — both arms of the original prediction
    satisfied with nothing intrinsic in the data.

    The rule is committed in `config/stage4_prespecification.yaml` rather than
    passed in, so it cannot be adjusted after seeing which genes it selects.

    A distribution rather than one matched control gene, deliberately: a single
    control is a point estimate of a noisy quantity, and drawing 200 costs
    nothing. The matched genes are not guaranteed to lack a differentiation
    story — which is why the claim this supports is narrow and exact: whether
    this gene's R-squared exceeds that of typical genes at the same abundance.
    """
    import numpy as np  # local: only this function needs it

    rules = spec["matched_null"]["candidates"]
    median = tumour_expression.median(axis=0)
    if target_gene_id not in median.index:
        raise PrespecError(f"{target_gene_id} is not in the expression matrix")
    target = float(median[target_gene_id])

    eligible = median.index
    if rules.get("exclude_frozen_panel", True):
        on_panel = index_map.index[
            index_map["on_panel"].astype(str).str.lower() == "true"
        ]
        eligible = eligible.difference(on_panel)
    if rules.get("exclude_self", True):
        eligible = eligible.difference([target_gene_id])

    window = float(rules["abundance_window_log2"])
    near = median[eligible][(median[eligible] - target).abs() <= window]

    floor = float(rules["min_detection_rate"])
    detected = (tumour_expression[near.index] > 0).mean(axis=0)
    qualified = list(near.index[detected >= floor])

    if not qualified:
        raise PrespecError(
            f"no abundance-matched candidates for {target_gene_id} within "
            f"+/-{window} log2 of {target:.3f}. Widening the window is a change "
            f"to the pre-specification, not a run-time decision."
        )

    cap = int(spec["matched_null"]["max_genes"])
    if len(qualified) <= cap:
        return sorted(qualified)
    rng = np.random.default_rng(int(spec["matched_null"]["seed"]))
    return sorted(rng.choice(sorted(qualified), size=cap, replace=False).tolist())
