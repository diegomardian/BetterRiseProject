"""W1 — Reference. Owns GSE178341 (Pelka 2021), ~371k cells, 62 patients.

Produces: cleaned AnnData, per-cell labels on axes 1 and 2, and versioned
S matrices at all four granularity rungs on the fixed gene index.

Hands W2 a five-patient pilot at end of week 2. If that slips past week 3 the
gate moves to week 7 — the harness is not the thing to compress
(execution_plan.md §8.3).

``build_signature`` and ``assert_no_target_leakage`` are re-exported because W2
needs them for the bake-off. Everything else in this package is W1-internal.
"""

from src.reference.signature import (
    LeakageError,
    assert_no_target_leakage,
    build_signature,
)

__all__ = ["LeakageError", "assert_no_target_leakage", "build_signature"]
