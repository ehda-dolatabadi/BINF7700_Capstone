#!/bin/bash
# Script: 05_subcluster.sbatch
# Purpose: Subset specific clusters and re-integrate for deeper analysis

#SBATCH --job-name=subcluster
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
ID="${ID:-cluster1}"
IDENT="${IDENT:-1}"

# Load configuration and setup paths
cd "$SLURM_SUBMIT_DIR"
source "config/default_paths.sh"
[ -f "config/local_paths.sh" ] && source "config/local_paths.sh"

# Redirect output to log file
exec > "$LOG/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.log"
exec 2>&1

# Define script and output directories
RUN="$WORK/scripts/seurat/Rscripts"
OUT_subset="$OUT/seurat/${ID}_05_integrated"

# Create output directory
mkdir -p "$OUT_subset"

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate seurat

# Cluster subsetting workflow
# Step 1: Subset cells from specified cluster identities
echo "Subsetting data..."
Rscript "$RUN/13_subset_clusters.R"\
    "$ID"\
    "$OUT_subset"\
    "$OUT/seurat/${main_ID}_processed.rds"\
    "$IDENT"
echo "Subset completed✅"

# Step 2: Re-normalize and re-integrate the subset
echo "Re-running normalization and integration..."
Rscript "$RUN/05_reintegrate_subset.R"\
    "$ID"\
    "$OUT_subset"\
    "$OUT/seurat/${ID}_integrated.rds"
echo "Normalization and integration completed✅"
