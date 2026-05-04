#!/bin/bash
# Script: 02_merge.sbatch
# Purpose: Merge per-timepoint CellChat objects and produce cross-timepoint plots

#SBATCH --job-name=merge
#SBATCH --partition=courses
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
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

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate cellchat

# Merge across timepoints and generate comparison plots
echo "Merging CellChat objects across timepoints..."
Rscript "$RUN/merge_timepoints.R" \
  "$ID" \
  "$OUT_cellchat"
echo "Merge completed✅"
