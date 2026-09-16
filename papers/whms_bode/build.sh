#!/usr/bin/env bash
# Build the independent full-paper revision without changing the original paper.
set -euo pipefail
cd "$(dirname "$0")"
PYTHON_BIN="${PYTHON_BIN:-python3}"
for tool in pdflatex bibtex; do command -v "$tool" >/dev/null; done
pdflatex -halt-on-error -interaction=batchmode main.tex >/dev/null
bibtex main >/dev/null
# A fourth LaTeX pass resolves page-dependent floats and cross-references.
for pass in 1 2 3; do pdflatex -halt-on-error -interaction=batchmode main.tex >/dev/null; done
"$PYTHON_BIN" check_submission.py
