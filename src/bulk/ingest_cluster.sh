#!/bin/bash -l
#
# Rebuild the TCGA matrices on gene index 1.0.0.
#
#   qsub src/bulk/ingest_cluster.sh
#
# WHY THIS EXISTS. The matrices on disk are 0.9.0 -- `ingest build` defaults to
# `PROVISIONAL_VERSION`, which is still "0.9.0", so a build run without
# --version silently produces the wrong ones. Everything downstream asks for
# 1.0.0 by name: the S matrices are on it, Stage 4 reads it, and
# run_purity_conditioned reads it. A 0.9.0 matrix does not fail to join, it
# joins on a different gene set.
#
# The STAR files are already downloaded, so this is a reindex-and-assemble, not
# a fetch. `download` runs first anyway because it is resumable and the file
# count was one short of the manifest -- one file, seconds, and it means the
# build is not quietly missing a sample.
#
#$ -N brp_ingest
#$ -pe omp 8
#$ -l h_rt=4:00:00
#$ -l mem_per_core=8G
#$ -j y
#$ -o logs/
#$ -V
set -euo pipefail
PROJECT_ROOT="${BRP_PROJECT_ROOT:-/projectnb/rise-batteries/bode/guanylin}"
REPO_DIR="${BRP_REPO_DIR:-$PROJECT_ROOT/BetterRiseProject}"
cd "$REPO_DIR"

VERSION=1.0.0

echo "=== $(date) ingest $VERSION on $(hostname) ==="
echo "HEAD: $(git rev-parse HEAD)  ($(git rev-parse --abbrev-ref HEAD))"

# step_build writes a versioned table at the end, and write_versioned_table
# refuses a dirty tree -- after reading 675 files rather than before.
if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
  echo "REFUSING: working tree is dirty. The build ends in a versioned write," >&2
  echo "which would refuse after the expensive part." >&2
  git status --short >&2
  exit 1
fi

if [ -z "${BRP_DATA_DIR:-}" ]; then
  echo "REFUSING: BRP_DATA_DIR is unset. Unset does not mean 'use the default'," >&2
  echo "it means read the wrong disk." >&2
  exit 1
fi
echo "BRP_DATA_DIR: $BRP_DATA_DIR"

if [ ! -f "$BRP_DATA_DIR/raw/tcga/gdc_star_files.tsv" ]; then
  echo "REFUSING: raw/tcga/gdc_star_files.tsv missing -- run" >&2
  echo "  python -m src.bulk.ingest query" >&2
  exit 1
fi
if [ ! -f "config/gene_index/gene_index_$VERSION.txt" ]; then
  echo "REFUSING: config/gene_index/gene_index_$VERSION.txt missing." >&2
  exit 1
fi

WANT=$(($(wc -l < "$BRP_DATA_DIR/raw/tcga/gdc_star_files.tsv") - 1))
HAVE=$(find "$BRP_DATA_DIR/raw/tcga/star" -name '*.tsv' 2>/dev/null | wc -l)
echo "STAR files: $HAVE downloaded of $WANT listed"

module load miniconda
conda activate brp-w3

# ---------------------------------------------------------------------------
echo; echo "### 1/2  resume the download (no-op if complete)"
python -m src.bulk.ingest download

# ---------------------------------------------------------------------------
echo; echo "### 2/2  assemble on gene index $VERSION"
# --version is the whole point. Without it this rebuilds 0.9.0 again.
python -m src.bulk.ingest build --version "$VERSION"

echo
for f in tcga_tpm tcga_log2cpm tcga_counts ; do
  ls -la "$BRP_DATA_DIR/processed/bulk/${f}_${VERSION}.parquet"
done

echo; echo "=== $(date) done ==="
echo "Next: qsub src/bulk/stage4_cluster.sh"
