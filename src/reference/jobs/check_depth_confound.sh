#!/bin/bash -l
#
# W2's depth-confound diagnostic over W1's labels. See PR #45 and decision #14.
#
#   qsub src/reference/jobs/check_depth_confound.sh
#
# Same shape as run_full_reference: one pass over 62 patients, labels rebuilt
# per patient. That took ~39 minutes, so this is not a login-node job.
#
#$ -N brp_w1_depthconf
#$ -pe omp 4
#$ -l h_rt=6:00:00
#$ -l mem_per_core=8G
#$ -j y
#$ -o logs/
#$ -V

set -euo pipefail
PROJECT_ROOT="${BRP_PROJECT_ROOT:-/projectnb/rise-batteries/bode/guanylin}"
REPO_DIR="${BRP_REPO_DIR:-$PROJECT_ROOT/BetterRiseProject}"

echo "=== $(date) depth confound on $(hostname) ==="
module load miniconda
conda activate brp-w1
cd "$REPO_DIR"
python src/reference/jobs/check_depth_confound.py "$@"
echo "=== $(date) done ==="
