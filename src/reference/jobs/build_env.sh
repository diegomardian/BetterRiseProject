#!/bin/bash -l
#
# Build the W1 conda environment as a batch job. W1.
#
# The solve is long (scanpy + Seurat + Bioconductor in one env) and an
# interactive session that ends partway leaves nothing behind. As a batch job it
# runs to completion regardless of your terminal, and everything lands in a log.
#
#   qsub src/reference/jobs/build_env.sh
#   qstat -u $USER                       # watch it
#   tail -f logs/brp_w1_env.o<jobid>     # read it
#
# Re-runnable: pass -force to delete an existing/partial env first.
#
#$ -N brp_w1_env
#$ -pe omp 8
#$ -l h_rt=8:00:00
#$ -j y
#$ -o logs/
#$ -V

set -euo pipefail

PROJECT_ROOT="${BRP_PROJECT_ROOT:-/projectnb/rise-batteries/bode/guanylin}"
REPO_DIR="${BRP_REPO_DIR:-$PROJECT_ROOT/BetterRiseProject}"
ENV_NAME=brp-w1

echo "=== $(date) starting on $(hostname) ==="
echo "repo:    $REPO_DIR"
echo "project: $PROJECT_ROOT"

module load miniconda

# Conda envs are large — this one carries R, Seurat and torch and will run to
#8-15 GB. SCC home directories have a small quota and the install dies partway
# through when it fills. Keep envs and the package cache in project space.
mkdir -p "$PROJECT_ROOT/conda/envs" "$PROJECT_ROOT/conda/pkgs"
conda config --add envs_dirs "$PROJECT_ROOT/conda/envs" 2>/dev/null || true
conda config --add pkgs_dirs "$PROJECT_ROOT/conda/pkgs" 2>/dev/null || true

echo "=== disk before ==="
df -h "$PROJECT_ROOT" | tail -1

cd "$REPO_DIR"
echo "=== branch: $(git rev-parse --abbrev-ref HEAD) @ $(git rev-parse --short HEAD) ==="

if [[ "${1:-}" == "-force" ]]; then
    echo "=== removing any existing $ENV_NAME ==="
    conda env remove -n "$ENV_NAME" -y 2>/dev/null || true
fi

# libmamba resolves this in minutes; the classic solver can grind for hours or
# give up entirely. Fall back only if the flag is unsupported.
echo "=== creating $ENV_NAME ==="
if ! conda env create -f env/w1_reference.yml --solver=libmamba; then
    echo "!!! libmamba solve failed or unsupported — retrying with the default solver"
    conda env create -f env/w1_reference.yml
fi

echo "=== installing the repo as an editable package ==="
conda run -n "$ENV_NAME" pip install -e ".[dev]"

# CopyKAT is the inferCNV cross-check (weeks 2-3) and is not on bioconda.
echo "=== installing CopyKAT from GitHub ==="
conda run -n "$ENV_NAME" R -e 'remotes::install_github("navinlabcode/copykat", upgrade="never")' \
    || echo "!!! CopyKAT install failed — not fatal until week 2, but note it"

echo "=== verifying ==="
conda run -n "$ENV_NAME" python -c "import scanpy, anndata, pandas, pyarrow; print('scanpy', scanpy.__version__)"
conda run -n "$ENV_NAME" python -m pytest -q

echo "=== disk after ==="
du -sh "$PROJECT_ROOT/conda/envs/$ENV_NAME" 2>/dev/null || true

echo "=== $(date) done ==="
echo "Activate with:  module load miniconda && conda activate $ENV_NAME"
