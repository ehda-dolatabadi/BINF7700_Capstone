#!/bin/bash
# Script: 06_subcells.sbatch
# Purpose: Subset cells based on marker scores and re-integrate

#SBATCH --job-name=subcells
#SBATCH --partition=courses
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=56
#SBATCH --mem=128G
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

# Exit on error, undefined variables, and pipe failures
set -euo pipefail

# Parameter defaults
main_ID="${main_ID:-all}"
ID="${ID:-enriched}"
MODE="${MODE:-remove}"
CELL_TYPES="${CELL_TYPES:-epithelial}"
CELL_THRESHOLDS="${CELL_THRESHOLDS:-1.0}"

# Load configuration and setup paths
cd "$SLURM_SUBMIT_DIR"
source "config/default_paths.sh"
[ -f "config/local_paths.sh" ] && source "config/local_paths.sh"

# Redirect output to log file
exec > "$LOG/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.log"
exec 2>&1

# Load marker definitions from file
MARKER_FILE="${MARKER_FILE:-$WORK/scripts/seurat/cell_markers.txt}"
source ${MARKER_FILE}

# Parse cell types into array
IFS=';' read -ra cells <<< "$CELL_TYPES"

# Build semicolon-separated marker list for selected cell types
CELL_MARKERS=""
for cell in "${cells[@]}"; do
    CELL_MARKERS="${CELL_MARKERS}${!cell};"
done
CELL_MARKERS="${CELL_MARKERS%;}"

# Define script and output directories
RUN="$WORK/scripts/seurat/Rscripts"
OUT_subcells="$OUT/seurat/${ID}_05_integrated"

# Create output directory
mkdir -p "$OUT_subcells"

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate seurat

# Cell subsetting workflow
# Step 1: Subset cells based on marker scores (keep or remove mode)
echo "Subsetting data..."
Rscript "$RUN/13_subset_cells.R" \
  "$ID" \
  "$OUT_subcells" \
  "$OUT/seurat/${main_ID}_integrated.rds" \
  "$MODE" \
  "$CELL_TYPES" \
  "$CELL_MARKERS" \
  "$CELL_THRESHOLDS"
echo "Subset completed✅"

# Step 2: Re-normalize and re-integrate the subset
echo "Re-running normalization and integration..."
Rscript "$RUN/05_reintegrate_subset.R" \
  "$ID" \
  "$OUT_subcells" \
  "$OUT/seurat/${ID}_integrated.rds"
echo "Normalization and integration completed✅"
