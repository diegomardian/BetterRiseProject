#!/bin/bash -l
#
# The refined tier-B test for G2. See docs/prereg_g2_mlh1.md.
#
#   python src/reference/jobs/run_g2_mlh1_contrast.py    <- fine, takes seconds
#   qsub  src/reference/jobs/run_g2_mlh1_contrast.sh     <- also fine
#
# THIS ONE DOES NOT NEED THE SCHEDULER. It reads two committed parquets and
# bootstraps 2,000 draws over at most 12 patients; it never touches the 9 GB
# matrix. The wrapper exists so nobody has to work out which jobs need qsub and
# which do not — the answer for every OTHER job in this directory is yes.
#
#$ -N brp_w1_g2mlh1
#$ -pe omp 1
#$ -l h_rt=1:00:00
#$ -j y
#$ -o logs/
#$ -V

set -euo pipefail
PROJECT_ROOT="${BRP_PROJECT_ROOT:-/projectnb/rise-batteries/bode/guanylin}"
REPO_DIR="${BRP_REPO_DIR:-$PROJECT_ROOT/BetterRiseProject}"

echo "=== $(date) G2 refined tier-B test on $(hostname) ==="
module load miniconda
conda activate brp-w1
cd "$REPO_DIR"
python src/reference/jobs/run_g2_mlh1_contrast.py "$@"
echo "=== $(date) done ==="
