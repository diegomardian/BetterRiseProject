#!/bin/bash -l
#
# A1: does the reference recipe cost the instrument its margin? Cluster wrapper.
#
#   qsub src/harness/jobs/run_pseudobulk_recovery.sh
#
# WHY NOT INTERACTIVE. The login node kills a process past 15 minutes of CPU at
# over 25% duty, and the deconvolution leg is the wildcard: 200 pseudobulk
# samples x 2 legs x 2 recipes x (NNLS + nu-SVR at 3 nu) is on the order of
# 2,400 SVR fits. It is likely minutes rather than hours -- 39k cells, and each
# fit is 770 genes x 7 compartments -- but one hung fit would die
# mid-experiment interactively with nothing to resume. BATCH, and let the
# cluster absorb the tail.
#
# The guards are the ones every job in this repo earns:
#
#   A dirty tree. write_versioned_table refuses one at WRITE time, after the
#   compute. Refusing here costs a second instead of the run.
#
#   An unset BRP_DATA_DIR. src/common/paths.py falls back to REPO_ROOT/data,
#   which exists and is nearly empty -- unset does not fail, it reads the wrong
#   disk. Refuse instead.
#
#   Missing inputs. The four Lee files are the whole input; discovering one is
#   absent at job start costs a second rather than a run.
#
#   A missing sklearn. `recover()` now refuses a skipped method outright, so
#   this is belt and braces -- but it fails at job start rather than after both
#   cohorts have loaded, and it names the env that is actually missing.
#
#$ -N brp_a1_recovery
#$ -pe omp 4
#$ -l h_rt=4:00:00
#$ -l mem_per_core=8G
#$ -j y
#$ -o logs/
#$ -V
set -euo pipefail
PROJECT_ROOT="${BRP_PROJECT_ROOT:-/projectnb/rise-batteries/bode/guanylin}"
REPO_DIR="${BRP_REPO_DIR:-$PROJECT_ROOT/BetterRiseProject}"
cd "$REPO_DIR"

echo "=== $(date) A1 pseudobulk recovery on $(hostname) ==="
echo "HEAD: $(git rev-parse HEAD)  ($(git rev-parse --abbrev-ref HEAD))"

if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
  echo "REFUSING: working tree is dirty. The recorded sha would not reproduce" >&2
  echo "the tables, and write_versioned_table would refuse anyway -- after the" >&2
  echo "compute rather than before it. Commit or stash, then resubmit." >&2
  git status --short >&2
  exit 1
fi

if git fetch --quiet origin 2>/dev/null; then
  BEHIND=$(git rev-list --count HEAD..@{u} 2>/dev/null || echo 0)
  if [ "$BEHIND" != "0" ]; then
    echo "WARNING: this checkout is $BEHIND commit(s) behind its remote." >&2
  fi
else
  echo "note: could not reach the remote, so behind-ness is unchecked" >&2
fi

if [ -z "${BRP_DATA_DIR:-}" ]; then
  echo "REFUSING: BRP_DATA_DIR is unset. #\$ -V exports the submitting shell's" >&2
  echo "environment, so submit from a shell that has it, or export it here." >&2
  echo "Unset does not mean 'use the default' -- it means read the wrong disk." >&2
  exit 1
fi
echo "BRP_DATA_DIR: $BRP_DATA_DIR"

MISSING=""
for f in \
  "raw/lee/GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz" \
  "raw/lee/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz" \
  "raw/lee/GSE144735_processed_KUL3_CRC_10X_annotation.txt.gz" \
  "raw/lee/GSE144735_processed_KUL3_CRC_10X_raw_UMI_count_matrix.txt.gz" ; do
  [ -f "$BRP_DATA_DIR/$f" ] || MISSING="$MISSING\n  $f"
done
if [ -n "$MISSING" ]; then
  echo "REFUSING: inputs missing under $BRP_DATA_DIR:" >&2
  printf "%b\n" "$MISSING" >&2
  exit 1
fi

module load miniconda

# PICK AN ENV THAT EXISTS, do not name one and hope. `env/w2_harness.yml`
# DECLARES scikit-learn, and that is a fact about a file rather than about this
# machine: brp-w2 has never been created here, and hardcoding it cost a job.
# brp-w3 and brp-w4 declare the same scikit-learn=1.5.1, and brp-w3 is known
# good -- the Stage 4 chain ran in it.
ACTIVATED=""
for env in brp-w2 brp-w3 brp-w4 ; do
  if conda activate "$env" 2>/dev/null; then
    ACTIVATED="$env"
    break
  fi
done
if [ -z "$ACTIVATED" ]; then
  echo "REFUSING: none of brp-w2 / brp-w3 / brp-w4 exists on this machine." >&2
  echo "Available:" >&2
  conda info --envs >&2
  echo "Create one with: conda env create -f env/w2_harness.yml" >&2
  exit 1
fi
echo "conda env: $ACTIVATED"

# Fail first rather than run an NNLS-only experiment. The driver's NuSVR
# adapter reports unavailability and is skipped silently; this refuses instead.
python - <<'PY'
import importlib
for name in ("numpy", "pandas", "scipy", "sklearn"):
    try:
        importlib.import_module(name)
    except ImportError as exc:
        raise SystemExit(
            f"A1 needs {name} in the activated env, and it is missing: {exc}. "
            f"The run would have silently skipped nu-SVR, which is the method "
            f"both surviving Stage 4 predictors came from."
        )
print("imports ok: numpy pandas scipy sklearn")
PY

python -m src.harness.jobs.run_pseudobulk_recovery

echo "=== $(date) done ==="