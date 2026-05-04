#!/bin/bash
# Script: run_pipeline.sh

# Exit on error, undefined variables, and pipe failures
set -euo pipefail

# Common variables
export ID="all"
export GROUP_BY="singler_label"

# Load configuration paths
source "config/default_paths.sh"
[ -f "config/local_paths.sh" ] && source "config/local_paths.sh"

# SLURM scripts directory
SLURM="$WORK/scripts/cellchat/slurm"

# Pipeline control flags
install=false
inference=false
merge=true

# -----------------------------------------------------------------------------

# Step 1: Install CellChat and GitHub-only dependencies (one-time setup)
if [ "$install" = true ]; then
	echo "Submitting package installation..."
	JOB_install=$(sbatch \
	  --parsable \
	  --export=ALL \
	  $SLURM/00_install.sbatch)
	echo "  Job ID: $JOB_install"
fi

# -----------------------------------------------------------------------------

# Step 2: Per-timepoint CellChat inference (array job, one task per timepoint)
if [ "$inference" = true ]; then
	echo "Submitting per-timepoint CellChat inference..."
	JOB_inference=$(sbatch \
	  --parsable \
	  --export=ALL \
	  $([ "$install" = true ] && echo "--dependency=afterok:$JOB_install") \
	  $SLURM/01_inference.sbatch)
	echo "  Job ID: $JOB_inference"
fi

# -----------------------------------------------------------------------------

# Step 3: Merge across timepoints and generate cross-timepoint plots
if [ "$merge" = true ]; then
	echo "Submitting cross-timepoint merge..."
	JOB_merge=$(sbatch \
	  --parsable \
	  --export=ALL \
	  $([ "$inference" = true ] && echo "--dependency=afterok:$JOB_inference") \
	  $SLURM/02_merge.sbatch)
	echo "  Job ID: $JOB_merge"
fi

# -----------------------------------------------------------------------------

echo "Pipeline submitted successfully!"
