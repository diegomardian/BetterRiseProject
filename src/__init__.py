"""BetterRiseProject — compositional vs. cell-intrinsic decomposition in CRC.

Layout (frozen, execution_plan.md §3.4):

    src/schema.py     shared output contract — PR + 2 approvals to change
    src/common/       shared plumbing: paths, provenance, frozen-config loaders
    src/reference/    W1  single-cell reference and S matrices
    src/harness/      W2  simulation harness, bake-off, calibrated cutpoints
    src/bulk/         W3  TCGA bulk and clinical
    src/estimator/    W4  Kitagawa decomposition and replication

Workstreams do not import each other's internals. Interfaces are artifacts on
disk (parquet on a fixed gene index), with two deliberate exceptions re-exported
through package ``__init__`` files: W2's calibrated cutpoints, which W4 needs,
and W1's ``build_signature``, which W2 needs for the bake-off.
"""

__version__ = "0.1.0"
