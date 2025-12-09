#!/bin/bash
#SBATCH --job-name=integrate
#SBATCH --partition=courses
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=56
#SBATCH --mem=128G
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

# Script: 02_integrate.sbatch
# Purpose: Integrate normalized samples to correct for batch effects
# Description: Performs batch correction across all normalized samples using Seurat integration

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
RUN="$WORK/scripts/seurat/Rscripts"
OUT_normalize="$OUT/seurat/${ID}_04_normalized"
OUT_integrate="$OUT/seurat/${ID}_05_integrated"

# Create output directory
mkdir -p "$OUT_integrate"

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate seurat

# Integration workflow
# Integrate all normalized samples using Seurat's integration method
echo "Integrating normalized samples..."
Rscript "$RUN/05_integrate.R"\
    "$ID"\
    "$OUT_integrate"\
    "$OUT_normalize"/*_normalized.rds
echo "Integration completed✅"
