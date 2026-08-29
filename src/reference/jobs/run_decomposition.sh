#!/bin/bash -l
#
# The per-gene summary decompose_cohort consumes. W1, §6.1.
#
#   qsub src/reference/jobs/build_decomposition_summary.sh
#
# One pass over 62 patients with labels rebuilt per patient — same shape as
# run_full_reference, which took ~40 minutes.
#
#$ -N brp_w1_decomp
#$ -pe omp 4
#$ -l h_rt=6:00:00
#$ -l mem_per_core=8G
#$ -j y
#$ -o logs/
#$ -V

set -euo pipefail
PROJECT_ROOT="${BRP_PROJECT_ROOT:-/projectnb/rise-batteries/bode/guanylin}"
REPO_DIR="${BRP_REPO_DIR:-$PROJECT_ROOT/BetterRiseProject}"

echo "=== $(date) decomposition on $(hostname) ==="
module load miniconda
conda activate brp-w1
cd "$REPO_DIR"
python src/reference/jobs/run_decomposition.py "$@"
echo "=== $(date) done ==="
