#!/bin/bash
# Script: 01_inference.sbatch
# Purpose: Run CellChat inference for each timepoint as an array job

#SBATCH --job-name=inference
#SBATCH --partition=courses
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --array=0-7
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

# Exit on error, undefined variables, and pipe failures
set -euo pipefail

# Parameter defaults
ID="${ID:-all}"
GROUP_BY="${GROUP_BY:-singler_cluster}"

# Load configuration and setup paths
cd "$SLURM_SUBMIT_DIR"
source "config/default_paths.sh"
[ -f "config/local_paths.sh" ] && source "config/local_paths.sh"

# Redirect output to log file
exec > "$LOG/${SLURM_JOB_NAME}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"
exec 2>&1

# Define script and output directories
RUN="$WORK/scripts/cellchat/Rscripts"
OUT_cellchat="$OUT/cellchat/${ID}_${GROUP_BY}"

# Select timepoint based on array task ID
TPS=("control" "3h" "24h" "72h" "7dpa" "14dpa" "22dpa" "33dpa")
TP="${TPS[$SLURM_ARRAY_TASK_ID]}"

# Per-timepoint output directory
TP_OUT="$OUT_cellchat/${ID}_${TP}"
mkdir -p "$TP_OUT"

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate cellchat

# Cell-cell communication inference for this timepoint
echo "Running CellChat inference for $TP..."
Rscript "$RUN/inference_per_timepoint.R" \
  "$ID" \
  "$TP_OUT" \
  "$OUT/seurat/${ID}_processed.rds" \
  "$GROUP_BY" \
  "$TP"
echo "CellChat inference for $TP completed✅"
