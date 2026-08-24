#!/bin/bash -l
#
# Labels and mature-cell counts for all 62 patients. W1, weeks 4-5.
#
#   qsub src/reference/jobs/run_full_reference.sh
#
# BATCH, not interactive. SCC kills a login-node process that uses more than 15
# minutes of CPU at over 25% duty — which this does, reading 62 patients out of
# a 9 GB deposit. nohup and tmux do not help: the limit is on CPU use, not on
# whether the shell survives.
#
# Watch it with:  qstat -u $USER ; tail -f logs/brp_w1_fullref.o*
#
#$ -N brp_w1_fullref
#$ -pe omp 4
#$ -l h_rt=6:00:00
#$ -l mem_per_core=8G
#$ -j y
#$ -o logs/
#$ -V

set -euo pipefail

PROJECT_ROOT="${BRP_PROJECT_ROOT:-/projectnb/rise-batteries/bode/guanylin}"
REPO_DIR="${BRP_REPO_DIR:-$PROJECT_ROOT/BetterRiseProject}"

echo "=== $(date) full reference run on $(hostname) ==="
module load miniconda
conda activate brp-w1
cd "$REPO_DIR"

python src/reference/jobs/run_full_reference.py "$@"

echo "=== $(date) done ==="
