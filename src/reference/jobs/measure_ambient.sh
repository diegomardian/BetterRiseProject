#!/bin/bash -l
#
# Ambient contamination across all 62 patients, as a batch job. W1, week 2.
#
#   qsub src/reference/jobs/measure_ambient.sh
#
# Batch, not qrsh. An interactive session dies when the connection does — a
# closed laptop kills a qrsh job — and this reads one patient at a time out of a
# 9 GB deposit, so it runs for tens of minutes.
#
# Watch it with:  qstat -u $USER ; tail -f logs/brp_w1_ambient.o*
#
#$ -N brp_w1_ambient
#$ -pe omp 4
#$ -l h_rt=4:00:00
#$ -l mem_per_core=8G
#$ -j y
#$ -o logs/
#$ -V

set -euo pipefail

PROJECT_ROOT="${BRP_PROJECT_ROOT:-/projectnb/rise-batteries/bode/guanylin}"
REPO_DIR="${BRP_REPO_DIR:-$PROJECT_ROOT/BetterRiseProject}"

echo "=== $(date) ambient contamination on $(hostname) ==="
module load miniconda
conda activate brp-w1
cd "$REPO_DIR"

python src/reference/jobs/measure_ambient.py "$@"

echo "=== $(date) done ==="
