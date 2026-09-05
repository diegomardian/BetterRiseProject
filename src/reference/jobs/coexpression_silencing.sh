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

module load miniconda
conda activate brp-w1

# Pinned invocation. No --seed override: the seed comes from
# src.common.provenance so the labels match the decomposition's.
python -m src.reference.jobs.coexpression_silencing \
    --cohorts smc kul3 gse178341

echo "=== $(date) done ==="
