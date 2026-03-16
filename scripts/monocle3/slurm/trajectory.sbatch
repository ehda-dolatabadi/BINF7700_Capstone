#!/bin/bash
# Script: trajectory.sbatch
# Purpose:

#SBATCH --job-name=trajectory
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

# Load configuration and setup paths
cd "$SLURM_SUBMIT_DIR"
source "config/default_paths.sh"
[ -f "config/local_paths.sh" ] && source "config/local_paths.sh"

# Redirect output to log file
exec > "$LOG/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.log"
exec 2>&1

# Define script and output directories
RUN="$WORK/scripts/monocle3/Rscripts"
OUT_umap="$OUT/seurat/${ID}_08_umap"
OUT_trajectory="$OUT/monocle3"

mkdir -p $OUT_trajectory

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate monocle

# Annotation workflow
# Rename cluster identities using auto annotation
echo "Running trajectory analysis..."
Rscript "$RUN/trajectory.R" \
  "$ID" \
  "$OUT_trajectory" \
  "$OUT_umap/${ID}_auto_annotated.rds"
echo "Trajectory analysis completed✅"

# Sync results to work directory (exclude large RDS files)
rsync -av --info=progress2 --stats --exclude '*.rds' --exclude '*.idx' --exclude '*.tar' --exclude '$OUT/monocle3/all_trajectory/' "$OUT/monocle3/" "$WORK/results/monocle3/"
