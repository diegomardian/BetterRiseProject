#!/bin/bash -l
#
# Fetch the ICBI CRC atlas. 32.7 GB, resumable.
#
#   qsub src/reference/jobs/fetch_icbi_atlas.sh
#
# WHY NOW, AND WHY ON ITS OWN. This is the long pole on the only current path
# to a mechanism result. The coexpression reading -- the one blocker that is a
# power problem rather than a structural one -- needs 14 study-level estimates
# to meta-analyse against the three that came back UNRESOLVED, and the atlas is
# where they are (results/*/icbi_premise_candidate_studies.parquet).
#
# It does NOT depend on the deconvolution leg. The coexpression job reads labels
# and raw counts per cell; it never deconvolves and never touches an S matrix.
# So this download has no reason to wait behind the reference-scale work, and
# waiting would cost days for nothing.
#
# What this job does NOT do: run anything. It fetches, verifies, and records.
# The per-study coexpression runs are a separate job against a file that is
# either fully present or absent, never half-written.
#
#$ -N brp_icbi_fetch
#$ -pe omp 2
#$ -l h_rt=24:00:00
#$ -l mem_per_core=4G
#$ -j y
#$ -o logs/
#$ -V
set -euo pipefail
PROJECT_ROOT="${BRP_PROJECT_ROOT:-/projectnb/rise-batteries/bode/guanylin}"
REPO_DIR="${BRP_REPO_DIR:-$PROJECT_ROOT/BetterRiseProject}"
cd "$REPO_DIR"

URL="https://crc.icbi.at/h5ad/final_crc_atlas-adata.h5ad"
EXPECTED_BYTES=32700000000     # ~32.7 GB, per the published object

echo "=== $(date) icbi fetch on $(hostname) ==="
echo "HEAD: $(git rev-parse HEAD)"

if [ -z "${BRP_DATA_DIR:-}" ]; then
  echo "REFUSING: BRP_DATA_DIR is unset -- unset means read the wrong disk." >&2
  exit 1
fi
DEST_DIR="$BRP_DATA_DIR/raw/icbi"
DEST="$DEST_DIR/final_crc_atlas-adata.h5ad"
mkdir -p "$DEST_DIR"

# Disk before download, not after. A 32 GB fetch that dies at 90% full has
# wasted hours and left a partial file that looks like a real one.
AVAIL=$(df -Pk "$DEST_DIR" | awk 'NR==2 {print $4 * 1024}')
echo "destination: $DEST_DIR"
echo "free space:  $(numfmt --to=iec "$AVAIL" 2>/dev/null || echo "$AVAIL bytes")"
if [ "$AVAIL" -lt 40000000000 ]; then
  echo "REFUSING: under 40 GB free. The object is ~32.7 GB and curl needs room" >&2
  echo "to finish; a truncated h5ad reads as a valid file until it does not." >&2
  exit 1
fi

if [ -f "$DEST" ]; then
  HAVE=$(stat -c%s "$DEST" 2>/dev/null || stat -f%z "$DEST")
  echo "existing file: $HAVE bytes -- resuming (curl -C -)"
fi

# -C - resumes; --retry survives a transient blip on a 24-hour transfer.
curl -L --fail --retry 5 --retry-delay 30 -C - -o "$DEST" "$URL"

BYTES=$(stat -c%s "$DEST" 2>/dev/null || stat -f%z "$DEST")
echo; echo "downloaded: $BYTES bytes"
if [ "$BYTES" -lt 30000000000 ]; then
  echo "REFUSING to record: $BYTES is far short of the expected ~$EXPECTED_BYTES." >&2
  echo "A truncated h5ad opens without error and reads short. Re-run to resume." >&2
  exit 1
fi

echo; echo "sha256 (this takes a few minutes on 32 GB) ..."
SHA=$(sha256sum "$DEST" | awk '{print $1}')
echo "$SHA"

echo; echo "Record this in data/manifest.csv, then commit:"
echo "data/raw/icbi/final_crc_atlas-adata.h5ad,$SHA,$BYTES,$URL,ICBI-CRC-atlas,$(date +%F),W1,ICBI CRC atlas: 4.26M cells, 49 studies; 229 paired patients over 24 studies (results/*/icbi_paired.parquet),"

echo; echo "=== $(date) done ==="
echo "Next: the per-study coexpression runs. 14 studies carry >= 100 epithelial"
echo "cells per arm -- see results/*/icbi_premise_candidate_studies.parquet."
