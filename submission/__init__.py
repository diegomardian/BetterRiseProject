"""Standalone benchmark for the workshop submission.

SEPARATE FROM THE PROJECT ON PURPOSE. Nothing in here is imported by
``src/``, and nothing in here writes to the frozen schema or to ``results/``.
It reads the project's estimator and positivity rule and treats them as one
method among several. That direction of dependency is the whole point: the
benchmark must not be able to change what it is measuring.

Ownership (CONTRIBUTING.md §2) assigns every ``src/`` subdirectory to a
workstream. This directory is outside that scheme, which is why the benchmark
can call W4's estimator and W2's positivity rule in the same file without
either workstream's boundary being crossed.
"""
