#!/bin/bash
# Script: run_pipeline.sh
# Purpose: Submit Monocle3 trajectory inference pipeline to SLURM scheduler
# Usage: ./run_pipeline.sh

# Exit on error, undefined variables, and pipe failures
set -euo pipefail

# Common variables
export ID="all"

# Load lineage definitions from file
export LINEAGE_FILE="$(pwd)/scripts/monocle3/lineages.txt"
n_lineage=$(grep -cE '^[a-zA-Z_][a-zA-Z0-9_]+=' ${LINEAGE_FILE})

# Load configuration paths
source "config/default_paths.sh"
[ -f "config/local_paths.sh" ] && source "config/local_paths.sh"

# SLURM scripts directory
SLURM="$WORK/scripts/monocle3/slurm"

# -----------------------------------------------------------------------------

JOB_wrapper=$(sbatch \
  --parsable \
  --export=ALL,\
  $SLURM/seurat_wrapper.sbatch)

# -----------------------------------------------------------------------------

echo "Submitting trajectory analysis..."
JOB_trajectory=$(sbatch \
  --parsable \
  --array=0-$((n_lineage-1)) \
  --export=ALL,\
ID=$ID,\
LINEAGE_FILE=$LINEAGE_FILE,\
DIMS=50 \
  --dependency=afterok:$JOB_wrapper \
  $SLURM/trajectory.sbatch)
echo "  Job ID: $JOB_trajectory"

# -----------------------------------------------------------------------------

echo "Pipeline submitted successfully!"