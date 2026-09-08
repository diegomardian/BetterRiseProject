#!/bin/bash -l
#
# Descriptive Zheng_2022 normal -> polyp -> carcinoma trajectories.
#
#   qsub src/reference/jobs/zheng_gradient.sh
#
# This reads only the three pre-specified patients, but their sparse rows still
# live inside the 30 GB ICBI H5AD.  Schedule it rather than using a login node.
# No result here is a pooled estimate, interval, test, or model.
#
#$ -N brp_zheng_gradient
#$ -pe omp 4
#$ -l h_rt=2:00:00
#$ -l mem_per_core=8G
#$ -j y
#$ -o logs/
#$ -V

set -euo pipefail

PROJECT_ROOT="${BRP_PROJECT_ROOT:-/projectnb/rise-batteries/bode/guanylin}"
REPO_DIR="${BRP_REPO_DIR:-$PROJECT_ROOT/BetterRiseProject}"
ATLAS="${BRP_ICBI_DIR:-/project/rise-batteries/bode/icbi}/final_crc_atlas-adata.h5ad"

cd "$REPO_DIR"
echo "=== $(date) Zheng descriptive gradient on $(hostname) ==="
echo "HEAD: $(git rev-parse HEAD)  ($(git rev-parse --abbrev-ref HEAD))"

if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
  echo "REFUSING: working tree is dirty." >&2
  git status --short >&2
  exit 1
fi
if [ -z "${BRP_DATA_DIR:-}" ]; then
  echo "REFUSING: BRP_DATA_DIR is unset." >&2
  exit 1
fi
if [ ! -f "$ATLAS" ]; then
  echo "REFUSING: $ATLAS not found." >&2
  exit 1
fi
OBS="$BRP_DATA_DIR/interim/icbi_obs.parquet"
if [ ! -f "$OBS" ]; then
  echo "REFUSING: $OBS not found; build the row-aligned cache first." >&2
  exit 1
fi

module load miniconda
conda activate brp-w1
python -m src.reference.jobs.zheng_gradient --atlas "$ATLAS" --obs-cache "$OBS"

echo "=== $(date) done ==="
echo "Tables are under results/ and are NOT committed by this job."
