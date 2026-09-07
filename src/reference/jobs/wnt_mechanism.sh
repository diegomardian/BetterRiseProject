#!/bin/bash -l
#
# The Wnt mechanism test. Does the terminal-differentiation fall track Wnt
# target activity INSIDE the surviving mature cells of the adenoma cohort?
#
#   qsub src/reference/jobs/wnt_mechanism.sh
#
# Pre-registered in docs/prereg_wnt_mechanism.md. READ §5 BEFORE READING ANY
# NUMBER THIS PRODUCES -- the branches and their consequences are fixed there,
# including the one where the answer is "technical" and no mechanism claim
# follows.
#
# NO NEW DATA AND NO DOWNLOAD. This is a fresh pass over the ICBI atlas already
# on disk. It is the cheap option, but it is not laptop-runnable: no committed
# table carries per-cell values, only per-(patient, gene) aggregates.
#
# RUN THIS BEFORE FREEING DISK. It reads the 30 GB atlas. If the atlas has been
# deleted to make room for a Becker fetch, this job cannot run until it is
# re-fetched (fetch_icbi_atlas.sh, resumable, ~25 min).
#
#$ -N brp_wnt_mechanism
#$ -pe omp 8
#$ -l h_rt=8:00:00
#$ -l mem_per_core=8G
#$ -j y
#$ -o logs/
#$ -V
set -euo pipefail
PROJECT_ROOT="${BRP_PROJECT_ROOT:-/projectnb/rise-batteries/bode/guanylin}"
REPO_DIR="${BRP_REPO_DIR:-$PROJECT_ROOT/BetterRiseProject}"
cd "$REPO_DIR"

ATLAS="${BRP_ICBI_DIR:-/project/rise-batteries/bode/icbi}/final_crc_atlas-adata.h5ad"

echo "=== $(date) Wnt mechanism test on $(hostname) ==="
echo "HEAD: $(git rev-parse HEAD)  ($(git rev-parse --abbrev-ref HEAD))"

if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
  echo "REFUSING: working tree is dirty." >&2
  git status --short >&2
  exit 1
fi
if [ -z "${BRP_DATA_DIR:-}" ]; then
  echo "REFUSING: BRP_DATA_DIR is unset -- unset means read the wrong disk." >&2
  exit 1
fi
if [ ! -f "$ATLAS" ]; then
  echo "REFUSING: $ATLAS not found." >&2
  echo "  If it was deleted to make room for another cohort, re-fetch with:" >&2
  echo "    qsub src/reference/jobs/fetch_icbi_atlas.sh   # resumable, ~25 min" >&2
  exit 1
fi
OBS="$BRP_DATA_DIR/interim/icbi_obs.parquet"
[ -f "$OBS" ] || { echo "REFUSING: $OBS not found." >&2; exit 1; }

module load miniconda
ACTIVATED=""
for env in brp-w1 brp-w2 brp-w3 ; do
  if conda activate "$env" 2>/dev/null; then ACTIVATED="$env"; break; fi
done
if [ -z "$ACTIVATED" ]; then
  echo "REFUSING: none of brp-w1 / brp-w2 / brp-w3 exists here." >&2
  conda info --envs >&2
  exit 1
fi
echo "conda env: $ACTIVATED"
echo "atlas: $ATLAS"

set +e
python -m src.reference.jobs.wnt_mechanism --atlas "$ATLAS"
code=$?
set -e

case "$code" in
  0) echo "  completed -- read the invariant-8 Wnt/maturity block FIRST." \
          "If the score tracks maturity, the mechanism reading is withheld." ;;
  3) echo "  no patient produced a correlation" ;;
  *) echo "  unexpected exit $code" >&2 ;;
esac

echo "=== $(date) done ==="
echo "Tables are under results/ and are NOT committed by this job."
exit "$code"
