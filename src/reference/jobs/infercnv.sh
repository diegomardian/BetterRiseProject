#!/bin/bash -l
#
# inferCNV, one array task per patient. W1, weeks 2-3.
#
#   qsub -t 1-5 src/reference/jobs/infercnv.sh patients_pilot.txt
#   qsub -t 1-62 src/reference/jobs/infercnv.sh patients_all.txt
#
# The patient list is one ID per line; SGE_TASK_ID picks the line.
#
# Per patient, never pooled: a CNV baseline built across patients folds germline
# copy-number variation and per-patient capture differences into the malignancy
# call. This is also why it parallelises cleanly — §8.1 calls inferCNV on 371k
# cells "hours to a day", and that is with the per-patient split.
#
# REFERENCE GROUP: the non-epithelial compartments (immune, stromal,
# endothelial), NOT the matched normal epithelium. execution_plan.md §4 says
# matched normal, but its own "done when" is that normal epithelium is not
# misread as tumour — and a population used as the baseline is non-malignant by
# construction, so that check would validate nothing. See
# src/reference/malignancy.py for the full argument.
#
#$ -N brp_w1_infercnv
#$ -pe omp 8
#$ -l h_rt=24:00:00
#$ -l mem_per_core=8G
#$ -j y
#$ -o logs/
#$ -V

set -euo pipefail

PATIENT_LIST="${1:?usage: qsub -t 1-N infercnv.sh <patient_list.txt>}"
PROJECT_ROOT="${BRP_PROJECT_ROOT:-/projectnb/rise-batteries/bode/guanylin}"
REPO_DIR="${BRP_REPO_DIR:-$PROJECT_ROOT/BetterRiseProject}"
OUT_DIR="${BRP_DATA_DIR:-$PROJECT_ROOT/data}/interim/infercnv"

PATIENT=$(sed -n "${SGE_TASK_ID}p" "$PATIENT_LIST")
if [[ -z "$PATIENT" ]]; then
    echo "no patient on line $SGE_TASK_ID of $PATIENT_LIST"; exit 1
fi

echo "=== $(date) patient $PATIENT on $(hostname) ==="
module load miniconda
conda activate brp-w1
cd "$REPO_DIR"
mkdir -p "$OUT_DIR/$PATIENT"

# Check the reference is large enough BEFORE spending hours on CNV inference.
python - "$PATIENT" <<'PY'
import sys
from pathlib import Path
import os
from src.reference.ingest import (
    read_gse178341_clusters, read_gse178341_index, assign_compartments,
)
from src.reference.malignancy import select_cnv_reference

patient = sys.argv[1]
data = Path(os.environ.get("BRP_DATA_DIR", "data")) / "raw" / "GSE178341"
obs, _ = read_gse178341_index(data / "GSE178341_crc10x_full_c295v4_submit.h5")
clusters = read_gse178341_clusters(
    data / "GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz"
)
compartment = assign_compartments(clusters).reindex(obs.index)
here = obs["patient_id"] == patient
report = select_cnv_reference(
    compartment[here], patient_id=obs.loc[here, "patient_id"]
)
print(report.to_string(index=False))
if not bool(report["usable"].iloc[0]):
    raise SystemExit(
        f"{patient}: too few non-epithelial reference cells. Malignancy is "
        f"not_called for this patient rather than guessed."
    )
PY

echo "=== running inferCNV for $PATIENT ==="
# W1: fill in once run_infercnv() is implemented. The R call needs the
# per-patient matrix, the gene-position file, and the reference group names.
Rscript -e "stop('W1: implement the inferCNV call — see src/reference/malignancy.py')"

echo "=== $(date) done $PATIENT ==="
