#!/bin/bash -l
#
# The coexpression reading, all three cohorts, one job, one sha.
#
#   qsub src/reference/jobs/coexpression_silencing.sh
#
# WHY THE GUARDS BELOW. Two of this week's mistakes cost hours each and both
# were detectable at job start:
#
#   A dirty tree. write_versioned_table refuses one -- correctly -- but it
#   refuses at WRITE time. The bulk chain scored 675 samples and then died on
#   DirtyTreeError with nothing to show. Checking first costs a second.
#
#   Old code. A 32-patient run was submitted against a saturated premise check
#   that had already been fixed on the remote, and its output could not even be
#   read back by the corrected code. A specific sha is NOT pinned here: a pin
#   needs editing every run and a stale pin is that same failure wearing a
#   different hat. The question worth asking is whether this checkout is behind
#   what it tracks.
#
#$ -N brp_coexpr
#$ -pe omp 4
#$ -l h_rt=6:00:00
#$ -l mem_per_core=8G
#$ -j y
#$ -o logs/
#$ -V
set -euo pipefail
PROJECT_ROOT="${BRP_PROJECT_ROOT:-/projectnb/rise-batteries/bode/guanylin}"
REPO_DIR="${BRP_REPO_DIR:-$PROJECT_ROOT/BetterRiseProject}"
cd "$REPO_DIR"

echo "=== $(date) coexpression on $(hostname) ==="
echo "HEAD: $(git rev-parse HEAD)  ($(git rev-parse --abbrev-ref HEAD))"

if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
  echo "REFUSING: working tree is dirty. The recorded sha would not reproduce" >&2
  echo "the tables, and write_versioned_table would refuse anyway -- after the" >&2
  echo "compute rather than before it. Commit or stash, then resubmit." >&2
  git status --short >&2
  exit 1
fi

# Behind-remote is a warning, not a refusal: it needs the network, and a blip
# on a compute node is not a reason to kill an otherwise valid run.
if git fetch --quiet origin 2>/dev/null; then
  BEHIND=$(git rev-list --count HEAD..@{u} 2>/dev/null || echo 0)
  if [ "$BEHIND" != "0" ]; then
    echo "WARNING: this checkout is $BEHIND commit(s) behind its remote." >&2
    echo "That is how a run against already-fixed code happened before." >&2
  fi
else
  echo "note: could not reach the remote, so behind-ness is unchecked" >&2
fi

# BRP_DATA_DIR is NOT optional here. src/common/paths.py falls back to
# REPO_ROOT/data when it is unset, and that directory exists and is nearly
# empty -- so an unset variable does not fail, it quietly reads the wrong disk.
# A job submitted from a shell that had lost the export died forty seconds in
# looking for Lee under the repository. Refuse instead.
if [ -z "${BRP_DATA_DIR:-}" ]; then
  echo "REFUSING: BRP_DATA_DIR is unset. #\$ -V exports the submitting shell's" >&2
  echo "environment, so submit from a shell that has it, or export it here." >&2
  echo "Unset does not mean 'use the default' -- it means read the wrong disk." >&2
  exit 1
fi
echo "BRP_DATA_DIR: $BRP_DATA_DIR"

# Every cohort's input, checked before any compute. The GSE178341 leg alone is
# tens of minutes; discovering a missing file at the end of it is the waste this
# wrapper exists to prevent.
MISSING=""
for f in \
  "raw/lee/GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz" \
  "raw/lee/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz" \
  "raw/lee/GSE144735_processed_KUL3_CRC_10X_annotation.txt.gz" \
  "raw/lee/GSE144735_processed_KUL3_CRC_10X_raw_UMI_count_matrix.txt.gz" \
  "raw/GSE178341/GSE178341_crc10x_full_c295v4_submit.h5" \
  "raw/GSE178341/GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz" \
  "raw/GSE178341/GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz" ; do
  [ -f "$BRP_DATA_DIR/$f" ] || MISSING="$MISSING\n  $f"
done
if [ -n "$MISSING" ]; then
  echo "REFUSING: inputs missing under $BRP_DATA_DIR:" >&2
  printf "%b\n" "$MISSING" >&2
  echo "Fetch them before resubmitting; the manifest carries urls and sha256." >&2
  exit 1
fi

module load miniconda
conda activate brp-w1

# Pinned invocation. No --seed override: the seed comes from
# src.common.provenance so the labels match the decomposition's.
python -m src.reference.jobs.coexpression_silencing \
    --cohorts smc kul3 gse178341

echo "=== $(date) done ==="
