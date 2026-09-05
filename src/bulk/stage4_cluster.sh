#!/bin/bash -l
#
# The whole Stage 4 chain, plus the purity producer, one job, one sha.
#
#   qsub src/bulk/stage4_cluster.sh
#
# Runs in the order the pre-specification requires and STOPS where it says to
# stop. Every refusal below is a pre-committed outcome, not an error: an
# instrument failure IS the Stage 4 result, and so is no rung being both
# quotable and estimable.
#
# WHY THIS IS A SCRIPT AND NOT A LIST OF COMMANDS. The gate between step 2 and
# step 4 depends on reading a table and deciding, and a decision that lives in
# a human's attention at 2am is a decision that gets skipped. The predictor
# check is enforced here instead. Two earlier rounds were also lost to pasted
# placeholders and a stray space in a filename; there is nothing to retype here.
#
#$ -N brp_stage4
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

echo "=== $(date) stage 4 on $(hostname) ==="
echo "HEAD: $(git rev-parse HEAD)  ($(git rev-parse --abbrev-ref HEAD))"

if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
  echo "REFUSING: working tree is dirty. write_versioned_table would refuse at" >&2
  echo "WRITE time -- after the compute rather than before it." >&2
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

# Unset does not mean "use the default" -- paths.py falls back to REPO_ROOT/data,
# which exists and is nearly empty, so the job reads the wrong disk and succeeds
# at finding nothing.
if [ -z "${BRP_DATA_DIR:-}" ]; then
  echo "REFUSING: BRP_DATA_DIR is unset. Submit from a shell that exports it." >&2
  exit 1
fi
echo "BRP_DATA_DIR: $BRP_DATA_DIR"

BULK="$BRP_DATA_DIR/processed/bulk"
MISSING=""
for f in "tcga_tpm_1.0.0.parquet" "tcga_log2cpm_1.0.0.parquet" "sample_manifest.tsv" ; do
  [ -f "$BULK/$f" ] || MISSING="$MISSING\n  processed/bulk/$f"
done
if [ -n "$MISSING" ]; then
  echo "REFUSING: inputs missing under $BRP_DATA_DIR:" >&2
  printf "%b\n" "$MISSING" >&2
  echo "Build them with: python -m src.bulk.ingest build" >&2
  exit 1
fi

module load miniconda
conda activate brp-w3

# ---------------------------------------------------------------------------
echo; echo "### 1/4  fractions, every rung"
# --linearise-reference is REQUIRED and deliberately not the default. The
# committed S matrices are mean-of-log1p and TCGA TPM is a linear mixture;
# without it the run refuses (exit 2) rather than returning a predictor that is
# exactly 0.0 on every sample while the instrument gate passes anyway.
python -m src.bulk.run_deconvolution \
    --rung all \
    --linearise-reference \
    --bulk "$BULK/tcga_tpm_1.0.0.parquet"

# ---------------------------------------------------------------------------
echo; echo "### 2/4  the predictor check, enforced rather than eyeballed"
python - <<'PY'
import glob, sys
import pandas as pd

newest = sorted(glob.glob("results/*/stage4_predictor_checks.parquet"))[-1]
checks = pd.read_parquet(newest)
print(f"reading {newest}")
print(checks[["granularity_rung", "method", "verdict", "n_samples",
              "n_exact_zero", "fraction_sd"]].to_string(index=False))
# -j y merges stdout and stderr into one log, and they flush independently --
# without this the STOP message lands ABOVE the table it is talking about.
sys.stdout.flush()

usable = checks[checks["verdict"] == "usable"]
if usable.empty:
    print("\nSTOP. No (rung, method) produced a usable predictor.", file=sys.stderr)
    print("That is the Stage 4 result, and the tables above record it. A "
          "constant predictor gives every gene R-squared 0, which is exactly "
          "what the pre-registered arm expects for GUCA2A -- so this must NOT "
          "be reported as the prediction being confirmed.", file=sys.stderr)
    for _, row in checks.iterrows():
        print(f"  {row['granularity_rung']}/{row['method']}: {row['detail']}",
              file=sys.stderr)
    sys.exit(3)
print(f"\n{len(usable)} of {len(checks)} usable; continuing.")
PY

# ---------------------------------------------------------------------------
echo; echo "### 3/4  the purity producer (unrelated to stage 4, same data)"
# Closes the repo's last uncommitted-producer case. Its clean twins are what
# make the final two dirty-table deletions safe.
python -m src.bulk.run_purity_conditioned

# ---------------------------------------------------------------------------
echo; echo "### 4/4  gate, then the pre-registered arm, every rung"
# ONE call with --rung all, not a shell loop. Each call writes the same table
# names into the same {date}_{sha} directory, so a loop keeps only the last
# rung: the 2026-09-05 run committed a gate table containing best4 alone while
# lineage's and crypt_position's numbers survived only in this log. The driver
# loops internally and writes once.
#
# Exit 3 (no usable predictor) and 4 (gate failed) are RESULTS, not errors, so
# set -e must not treat them as failures.
set +e
python -m src.bulk.run_stage4_variance --rung all
code=$?
set -e
case "$code" in
  0) echo "  a verdict was reached on at least one rung" ;;
  3) echo "  no usable predictor -- that is the result" ;;
  4) echo "  INSTRUMENT GATE FAILED -- that is the result, and no R-squared" \
          "is reported" ;;
  5) echo "  no verdict could be formed" ;;
  *) echo "  unexpected exit $code" >&2 ;;
esac
FAILED_HARD=0
[ "$code" -le 5 ] || FAILED_HARD=1

echo; echo "=== $(date) done ==="
echo "Tables are under results/. They are NOT committed by this job -- read"
echo "them first, then commit. Start with stage4_predictor_checks and"
echo "stage4_instrument_gate; a verdict of 'indeterminate' means something"
echo "different depending on what those two say."
exit "$FAILED_HARD"
