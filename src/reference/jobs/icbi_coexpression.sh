#!/bin/bash -l
#
# The ICBI coexpression reading. Pelka first, the other thirteen after.
#
#   qsub -v BRP_ICBI_STUDY=Pelka_2021_Cell src/reference/jobs/icbi_coexpression.sh
#   qsub -v BRP_ICBI_STUDY=all             src/reference/jobs/icbi_coexpression.sh
#
#   # path C, the adenoma reading -- Chen_2021_Cell IS the VUMC/HTAN polyp atlas
#   qsub -v BRP_ICBI_STUDY=Chen_2021_Cell,BRP_ICBI_ARMS=adenoma \
#        src/reference/jobs/icbi_coexpression.sh
#
# RUN PELKA FIRST AND READ THE VERDICT. Pelka_2021_Cell is GSE178341, which this
# project has already analysed three ways, so it is the only one of the fourteen
# with ground truth to be checked against. The job compares its own result to
# the committed table and exits 4 if the bar fails. Running the other thirteen
# before that passes is running fourteen studies through an adaptation nobody
# has checked.
#
#$ -N brp_icbi_coexpr
#$ -pe omp 8
#$ -l h_rt=12:00:00
#$ -l mem_per_core=8G
#$ -j y
#$ -o logs/
#$ -V
set -euo pipefail
PROJECT_ROOT="${BRP_PROJECT_ROOT:-/projectnb/rise-batteries/bode/guanylin}"
REPO_DIR="${BRP_REPO_DIR:-$PROJECT_ROOT/BetterRiseProject}"
cd "$REPO_DIR"

STUDY="${BRP_ICBI_STUDY:-Pelka_2021_Cell}"
ATLAS="${BRP_ICBI_DIR:-/project/rise-batteries/bode/icbi}/final_crc_atlas-adata.h5ad"

echo "=== $(date) icbi coexpression on $(hostname) ==="
echo "HEAD: $(git rev-parse HEAD)  ($(git rev-parse --abbrev-ref HEAD))"
echo "study: $STUDY"

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
  echo "Fetch it with: qsub src/reference/jobs/fetch_icbi_atlas.sh" >&2
  echo "and export BRP_ICBI_DIR if it does not live under /project." >&2
  exit 1
fi
# The obs cache. Built from the LOCAL atlas rather than re-fetched over the
# network -- range requests were for a 32.7 GB object we did not have, and we
# have it now.
OBS="$BRP_DATA_DIR/interim/icbi_obs.parquet"
if [ ! -f "$OBS" ]; then
  echo "cached obs missing; building it from the local atlas ..."
  mkdir -p "$BRP_DATA_DIR/interim"
  python -m src.reference.jobs.pull_icbi_metadata --url "$ATLAS" --allow-dirty
  [ -f "$OBS" ] || { echo "REFUSING: still no $OBS" >&2; exit 1; }
fi
echo "obs cache: $OBS"
echo "atlas: $ATLAS ($(du -h "$ATLAS" | cut -f1))"

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

python - <<'PY'
import importlib, sys
for name in ("numpy", "pandas", "scipy", "h5py"):
    try:
        importlib.import_module(name)
    except ImportError as exc:
        raise SystemExit(f"icbi coexpression needs {name}: {exc}")
print("imports ok: numpy pandas scipy h5py")
PY

# BRP_ICBI_ARMS=adenoma scores polyp against the patient's own normal, at both
# lineage and best4. best4 is the resolution the question is actually posed at
# and the one carcinoma structurally could not support -- a median of 3 mature
# cells there. Adenomas retain differentiated epithelium.
ARMS="${BRP_ICBI_ARMS:-carcinoma}"
if [ "$ARMS" = "adenoma" ]; then
  RUNGS="lineage best4"
else
  RUNGS="lineage"
fi
echo "arms: $ARMS | rungs: $RUNGS"

set +e
if [ "$STUDY" = "all" ]; then
  python -m src.reference.jobs.icbi_coexpression --atlas "$ATLAS" --all \
      --arms "$ARMS" --rungs $RUNGS
else
  python -m src.reference.jobs.icbi_coexpression --atlas "$ATLAS" --study "$STUDY" \
      --arms "$ARMS" --rungs $RUNGS
fi
code=$?
set -e

case "$code" in
  0) echo "  completed" ;;
  3) echo "  no study produced a scored patient -- that is a result" ;;
  4) echo "  VALIDATION FAILED against the committed GSE178341 table." \
          "Do NOT run the other thirteen until the divergence is understood." >&2 ;;
  *) echo "  unexpected exit $code" >&2 ;;
esac

echo "=== $(date) done ==="
echo "Tables are under results/ and are NOT committed by this job."
exit "$code"
