#!/bin/bash
# Script: cellchat.sbatch
# Purpose:

#SBATCH --job-name=cellchat
#SBATCH --partition=courses
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
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
exec > "$LOG/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.log"
exec 2>&1

# Define script and output directories
RUN="$WORK/scripts/cellchat/Rscripts"
OUT_cellchat="$OUT/cellchat/${ID}_${GROUP_BY}"

mkdir -p "$OUT_cellchat"

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate cellchat

# Cell-cell communication inference
echo "Running CellChat analysis..."
Rscript "$RUN/cellchat.R" \
  "$ID" \
  "$OUT_cellchat" \
  "$OUT/seurat/${ID}_processed.rds" \
  "$GROUP_BY"
echo "CellChat analysis completed✅"
