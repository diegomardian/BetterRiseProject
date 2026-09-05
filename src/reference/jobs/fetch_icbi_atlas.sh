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

# WHERE IT LANDS. Default is $BRP_DATA_DIR/raw/icbi, but /projectnb was at 15 GB
# free of a 50 GB quota when this was written, and /project/rise-batteries had
# ~35 GB. BRP_ICBI_DIR moves ONLY this download rather than re-pointing every
# data path, which is what overriding BRP_DATA_DIR would do.
#
#   export BRP_ICBI_DIR=/project/rise-batteries/bode/icbi
#
# NOT /scratch. `df` reports it as /dev/sda8 -- node-local and purged. The
# fetch ran on scc-ye1 and A1 on scc-pi4; a file written to one node's /scratch
# is invisible to a job on another. That is the /tmp trap this project has
# already paid for once.
DEST_DIR="${BRP_ICBI_DIR:-}"

echo "=== $(date) icbi fetch on $(hostname) ==="
echo "HEAD: $(git rev-parse HEAD)"

if [ -z "${BRP_DATA_DIR:-}" ]; then
  echo "REFUSING: BRP_DATA_DIR is unset -- unset means read the wrong disk." >&2
  exit 1
fi
DEST_DIR="${DEST_DIR:-$BRP_DATA_DIR/raw/icbi}"
DEST="$DEST_DIR/final_crc_atlas-adata.h5ad"
mkdir -p "$DEST_DIR"
case "$DEST_DIR" in
  /scratch/*|/tmp/*)
    echo "REFUSING: $DEST_DIR is node-local and purged. This job runs on one" >&2
    echo "compute node and the analysis runs on another, which will not see it." >&2
    exit 1 ;;
esac

# Disk before download, not after. A 32 GB fetch that dies at 90% full has
# wasted hours and left a partial file that looks like a real one.
AVAIL=$(df -Pk "$DEST_DIR" | awk 'NR==2 {print $4 * 1024}')
echo "destination: $DEST_DIR"
echo "free space:  $(numfmt --to=iec "$AVAIL" 2>/dev/null || echo "$AVAIL bytes")"
# Object plus a 2 GB working margin, rather than a round 40. The round number
# refused a filesystem that would in fact have held it.
NEEDED=$((EXPECTED_BYTES + 2000000000))
if [ "$AVAIL" -lt "$NEEDED" ]; then
  echo "REFUSING: $(numfmt --to=iec "$AVAIL" 2>/dev/null || echo "$AVAIL") free," >&2
  echo "against $(numfmt --to=iec "$NEEDED" 2>/dev/null || echo "$NEEDED") needed" >&2
  echo "(the object is ~32.7 GB plus working room). A truncated h5ad reads as a" >&2
  echo "valid file until it does not." >&2
  echo >&2
  echo "Check quota with \`pquota\`. If /projectnb is full and /project is not:" >&2
  echo "  export BRP_ICBI_DIR=/project/rise-batteries/bode/icbi" >&2
  echo "and resubmit. Do NOT use /scratch -- it is node-local and purged." >&2
  exit 1
fi
MARGIN=$((AVAIL - EXPECTED_BYTES))
if [ "$MARGIN" -lt 5000000000 ]; then
  echo "WARNING: only $(numfmt --to=iec "$MARGIN" 2>/dev/null || echo "$MARGIN") would" >&2
  echo "remain after the download. Anything else writing here may fail." >&2
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
