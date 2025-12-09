#!/bin/bash
#SBATCH --job-name=markerScore
#SBATCH --partition=courses
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=0-2
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

# Script: 04_score_markers.sbatch
# Purpose: Score cells based on known cell type marker genes
# Description: Calculates module scores for predefined cell type markers in parallel using array jobs

set -euo pipefail

# Parameter defaults
main_ID="${main_ID:-all}"

# Load configuration and setup paths
cd "$SLURM_SUBMIT_DIR"
source "config/default_paths.sh"
[ -f "config/local_paths.sh" ] && source "config/local_paths.sh"

# Redirect output to log file
exec > "$LOG/${SLURM_JOB_NAME}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"
exec 2>&1

# Load marker definitions from file
MARKER_FILE="${MARKER_FILE:-$WORK/scripts/seurat/cell_markers.txt}"
source ${MARKER_FILE}
cells=($(grep -E '^[a-zA-Z_]+=' ${MARKER_FILE} | cut -d'=' -f1))

# Safety check: exit if array task exceeds number of cell types
if [ ${SLURM_ARRAY_TASK_ID} -ge ${#cells[@]} ]; then
  echo "End of cell names."
  exit 0
fi

# Select cell type and markers based on array task ID
CELL_NAME="${cells[$((SLURM_ARRAY_TASK_ID))]}"
ID="${ID:-$CELL_NAME}"
MARKERS="${!CELL_NAME}"

# Define script and output directories
RUN="$WORK/scripts/seurat/Rscripts"
OUT_cluster="$OUT/seurat/${main_ID}_07_clustered"
OUT_score="$OUT/seurat/${main_ID}_10_scored"

# Create output directory
mkdir -p "$OUT_score"

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate seurat

# Marker scoring workflow
# Calculate module score for cell type specific markers
echo "Scoring $CELL_NAME markers..."
Rscript "$RUN/10_score_markers.R" \
  "$ID" \
  "$OUT_score" \
  "$OUT_cluster/${main_ID}_clustered.rds" \
  "$CELL_NAME" \
  "$MARKERS"
echo "Scoring completed✅"

# Sync results to work directory (exclude large RDS files)
rsync -av --info=progress2 --stats --exclude '*.rds' "$OUT/seurat/" "$WORK/results/seurat/"
