#!/bin/bash -l
#
# The instrument's first positive control. MLH1 in the mature cells of patients
# whose promoters are known -- from an assay -- to be methylated.
#
#   qsub src/reference/jobs/mlh1_positive_control.sh
#
# RUN THE CALIBRATION JOB FIRST AND COMMIT ITS TABLES. This reading reports a
# Student-t interval rather than the percentile bootstrap the rest of the
# repository uses, and the measurement that justifies that must exist BEFORE
# the numbers do, or the choice of interval is a free parameter chosen after
# seeing the result:
#
#   python -m src.reference.jobs.interval_calibration   # laptop, no atlas
#   git add results/<dir> && git commit
#
# The pre-registration is docs/prereg_g2_mlh1_within_stratum.md. Read §5 before
# reading any number this produces: the falsifiers and the pre-committed
# consequences of each branch are fixed there.
#
#$ -N brp_mlh1_control
#$ -pe omp 8
#$ -l h_rt=6:00:00
#$ -l mem_per_core=8G
#$ -j y
#$ -o logs/
#$ -V
set -euo pipefail
PROJECT_ROOT="${BRP_PROJECT_ROOT:-/projectnb/rise-batteries/bode/guanylin}"
REPO_DIR="${BRP_REPO_DIR:-$PROJECT_ROOT/BetterRiseProject}"
cd "$REPO_DIR"

ATLAS="${BRP_ICBI_DIR:-/project/rise-batteries/bode/icbi}/final_crc_atlas-adata.h5ad"
RUNG="${BRP_MLH1_RUNG:-lineage}"

echo "=== $(date) MLH1 positive control on $(hostname) ==="
echo "HEAD: $(git rev-parse HEAD)  ($(git rev-parse --abbrev-ref HEAD))"
echo "rung: $RUNG"

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
  exit 1
fi
OBS="$BRP_DATA_DIR/interim/icbi_obs.parquet"
if [ ! -f "$OBS" ]; then
  echo "REFUSING: $OBS not found. Build it with pull_icbi_metadata." >&2
  exit 1
fi

# The calibration tables are a PRECONDITION, not a nicety. This reading's
# interval was chosen by them; without them committed, the choice is unrecorded.
if ! ls results/*/interval_calibration.parquet >/dev/null 2>&1; then
  echo "REFUSING: no results/*/interval_calibration.parquet." >&2
  echo "  Run: python -m src.reference.jobs.interval_calibration" >&2
  echo "  then commit its tables. The interval this job reports is chosen by" >&2
  echo "  that measurement, and choosing it afterwards makes it a free" >&2
  echo "  parameter picked with the result already visible." >&2
  exit 1
fi
echo "calibration: $(ls -d results/*/interval_calibration.parquet | tail -1)"

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

set +e
python -m src.reference.jobs.mlh1_positive_control --atlas "$ATLAS" --rung "$RUNG"
code=$?
set -e

case "$code" in
  0) echo "  completed -- read the VERDICT ON THE INSTRUMENT block above" ;;
  3) echo "  no patient scored, or MLH1 produced no row in the primary arm" ;;
  *) echo "  unexpected exit $code" >&2 ;;
esac

echo "=== $(date) done ==="
echo "Tables are under results/ and are NOT committed by this job."
exit "$code"
