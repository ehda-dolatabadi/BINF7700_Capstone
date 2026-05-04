#!/bin/bash
# Script: trajectory.sbatch
# Purpose: Run per-lineage Monocle3 trajectory inference as a SLURM array job

#SBATCH --job-name=trajectory
#SBATCH --partition=courses
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --array=0-3
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

# Exit on error, undefined variables, and pipe failures
set -euo pipefail

# Parameter defaults
ID="${ID:-all}"
DIMS="${DIMS:-10}"

# Load configuration and setup paths
cd "$SLURM_SUBMIT_DIR"
source "config/default_paths.sh"
[ -f "config/local_paths.sh" ] && source "config/local_paths.sh"

# Redirect output to log file
exec > "$LOG/${SLURM_JOB_NAME}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"
exec 2>&1

# Load lineage definitions from file
LINEAGE_FILE="${LINEAGE_FILE:-$WORK/scripts/monocle3/lineages.txt}"
source ${LINEAGE_FILE}
lineages=($(grep -E '^[a-zA-Z_][a-zA-Z0-9_]+=' ${LINEAGE_FILE} | cut -d'=' -f1))

# Safety check: exit if array task exceeds number of lineages
if [ ${SLURM_ARRAY_TASK_ID} -ge ${#lineages[@]} ]; then
  echo "End of lineages."
  exit 0
fi

# Select lineage name and cluster labels based on array task ID
LINEAGE_NAME="${lineages[$SLURM_ARRAY_TASK_ID]}"
LINEAGE_CLUSTERS="${!LINEAGE_NAME}"

# Define script and output directories
RUN="$WORK/scripts/monocle3/Rscripts"
OUT_trajectory="$OUT/monocle3"
mkdir -p "$OUT_trajectory"

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate monocle3

# Trajectory workflow
echo "Running trajectory analysis for lineage: $LINEAGE_NAME ($LINEAGE_CLUSTERS)..."
Rscript "$RUN/trajectory.R" \
  "$ID" \
  "$OUT_trajectory" \
  "$OUT/seurat/${ID}_processed.rds" \
  "$LINEAGE_NAME" \
  "$LINEAGE_CLUSTERS" \
  "$DIMS"
echo "Trajectory analysis for $LINEAGE_NAME completed✅"