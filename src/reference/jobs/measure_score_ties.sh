#!/bin/bash -l
#
# Tie structure of M, cohort-wide. W1.
#
#   qsub src/reference/jobs/measure_m_ties.sh
#
# WHY THIS IS NOT A LOGIN-NODE JOB. It reads one patient at a time and only
# needs column means, which looked laptop-scale — and it is, per patient. But
# read_gse178341 materialises one patient's slice of a 9 GB matrix, and the
# largest patients (C162 at 13,942 cells) exceed what an SCC login node allows.
# A first attempt was Killed at patient 34 of 62. Login nodes also cap
# interactive CPU at 15 minutes.
#
#$ -N brp_w1_scoreties
#$ -pe omp 4
#$ -l h_rt=4:00:00
#$ -l mem_per_core=8G
#$ -j y
#$ -o logs/
#$ -V

set -euo pipefail

PROJECT_ROOT="${BRP_PROJECT_ROOT:-/projectnb/rise-batteries/bode/guanylin}"
REPO_DIR="${BRP_REPO_DIR:-$PROJECT_ROOT/BetterRiseProject}"

echo "=== $(date) score tie structure on $(hostname) ==="
module load miniconda
conda activate brp-w1
cd "$REPO_DIR"

python src/reference/jobs/measure_score_ties.py "$@"

echo "=== $(date) done ==="
