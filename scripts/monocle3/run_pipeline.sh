#!/bin/bash
# Script: run_pipeline.sh

# Exit on error, undefined variables, and pipe failures
set -euo pipefail

# Common variables
export ID="all"

# Load configuration paths
source "config/default_paths.sh"
[ -f "config/local_paths.sh" ] && source "config/local_paths.sh"

# SLURM scripts directory
SLURM="$WORK/scripts/monocle3/slurm"

# -----------------------------------------------------------------------------

echo "Submitting trajectory analysis..."
JOB1=$(sbatch --parsable --export=ALL,ID=$ID, $SLURM/trajectory.sbatch)
echo "  Job ID: $JOB1"

# -----------------------------------------------------------------------------

echo "Pipeline submitted successfully!"
