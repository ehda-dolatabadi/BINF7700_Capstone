#!/bin/bash
# Script: 07_annotate.sbatch
# Purpose: Annotate clusters with cell type labels

#SBATCH --job-name=annotate
#SBATCH --partition=courses
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --array=0-3
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

# Exit on error, undefined variables, and pipe failures
set -euo pipefail

# Load configuration and setup paths
cd "$SLURM_SUBMIT_DIR"
source "config/default_paths.sh"
[ -f "config/local_paths.sh" ] && source "config/local_paths.sh"

# Redirect output to log file
exec > "$LOG/${SLURM_JOB_NAME}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"
exec 2>&1

# Load annotations from file
ANNOTATION_FILE="${ANNOTATION_FILE:-$WORK/scripts/seurat/clusters_annotation.txt}"
source ${ANNOTATION_FILE}
ids=($(grep -E '^[a-zA-Z_][a-zA-Z0-9_]+=' ${ANNOTATION_FILE} | cut -d'=' -f1))

# Safety check: exit if array task exceeds number of ids
if [ ${SLURM_ARRAY_TASK_ID} -ge ${#ids[@]} ]; then
  echo "End of IDs."
  exit 0
fi

# Select ID and labels based on array task ID
ID="${ids[$((SLURM_ARRAY_TASK_ID))]}"
LABELS="${!ID}"

# Define script and output directories
RUN="$WORK/scripts/seurat/Rscripts"
OUT_umap="$OUT/seurat/${ID}_08_umap"
OUT_annotated="$OUT/seurat/${ID}_13_annotated"

# Create output directory
mkdir -p "$OUT_annotated"

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate seurat

# Annotation workflow
# Rename cluster identities using ordered labels
echo "Annotating clusters..."
Rscript "$RUN/13_annotate_clusters.R" \
  "$ID" \
  "$OUT_annotated" \
  "$OUT_umap/${ID}_umap.rds" \
  "$LABELS"
echo "Annotation completed✅"

# Sync results to work directory (exclude large RDS files)
rsync -av --info=progress2 --stats --exclude '*.rds' "$OUT/seurat/" "$WORK/results/seurat/"
