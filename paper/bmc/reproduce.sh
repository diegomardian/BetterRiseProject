#!/usr/bin/env bash
# Reproduce the audit in a fresh environment.
#
# Clones this branch into a temporary directory that has none of the local
# raw data, runs every deliverable that reads only committed inputs, and
# checks each output byte-for-byte against the committed table. Then runs the
# test subset that ties the manuscript and the index to those tables.
#
# This is the week-4 "a teammate reproduces it in a fresh environment" check.
# The two defect classes it is built to catch are path bugs -- absolute paths,
# an assumption that a working directory is the repo root, a table resolved by
# mtime -- which nothing else surfaces because everything works where it was
# made.
#
#   bash paper/bmc/reproduce.sh                 # use the current repo/branch
#   bash paper/bmc/reproduce.sh /path/to/repo w1/bmc-manuscript
#
# Exit 0 means every deliverable reproduced exactly. Any mismatch prints the
# table and exits non-zero.
set -euo pipefail

SRC="${1:-$(git rev-parse --show-toplevel)}"
BRANCH="${2:-$(git -C "$SRC" rev-parse --abbrev-ref HEAD)}"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/brp-reproduce.XXXXXX")"
CLONE="$WORK/repo"

echo "reproducing $BRANCH from $SRC"
echo "  into $CLONE"
git clone -q --branch "$BRANCH" "$SRC" "$CLONE"
cd "$CLONE"

# Deliverables that read only committed inputs. Jobs needing raw data (the
# Lee cutpoint sweep, the atlas extraction) are not here on purpose: they are
# not reproducible without the deposit, and the audit's claim is about the
# committed derived tables.
JOBS=(
  reproduction_manifest
  cutpoint_crossing_brackets
  adenoma_claim_sensitivity
  crowell_summary_reproduction
  disval_reproduction
  table_resolution_audit
)
TABLES=(
  reproduction_manifest
  cutpoint_crossing_brackets cutpoint_crossing_brackets_summary
  adenoma_claim_sensitivity adenoma_claim_sensitivity_summary estimability_attrition
  crowell_summary_reproduction disval_reproduction table_resolution_ambiguity
)

# Record the committed table each name resolves to *before* the jobs run.
# The record goes outside the clone: a scratch file in the repo root makes the
# tree dirty, and the results writer refuses a dirty tree (correctly -- it is
# the guard that made the first run of this script fail).
BEFORE="$WORK/before.txt"
python - "${TABLES[*]}" <<'PY' > "$BEFORE"
import sys
from src.common.paths import RESULTS_DIR
from src.reference.table_resolution import newest_by_time
for name in sys.argv[1].split():
    path = newest_by_time(RESULTS_DIR, name)
    print(f"{name}\t{path}")
PY

for job in "${JOBS[@]}"; do
  echo "  running $job"
  python -m "src.reference.jobs.$job" >/dev/null
done

echo "comparing each output to the committed table"
python - "${TABLES[*]}" "$BEFORE" <<'PY'
import sys
import pandas as pd
from src.common.paths import RESULTS_DIR
from src.reference.table_resolution import newest_by_time

before = {}
for line in open(sys.argv[2]):
    name, path = line.rstrip("\n").split("\t")
    before[name] = path

failed = []
for name in sys.argv[1].split():
    after = newest_by_time(RESULTS_DIR, name)
    old = pd.read_parquet(before[name])
    new = pd.read_parquet(after)
    if old.equals(new):
        print(f"  {name:42s} IDENTICAL  ({after.parent.name})")
    else:
        print(f"  {name:42s} DIFFERS from {before[name]}")
        failed.append(name)
if failed:
    raise SystemExit(f"reproduction failed for: {failed}")
PY

echo "running the manuscript and index tests"
python -m pytest -q \
  tests/test_bmc_paper_numbers.py \
  tests/test_reproduction.py \
  tests/test_claim_sensitivity.py \
  tests/test_disval_reproduction.py \
  tests/test_crowell_reproduction.py \
  tests/test_table_resolution.py \
  tests/test_cutpoint_brackets.py

echo
echo "PASS. Fresh-environment reproduction complete in $CLONE"
