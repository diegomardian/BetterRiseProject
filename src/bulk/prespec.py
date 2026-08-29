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
