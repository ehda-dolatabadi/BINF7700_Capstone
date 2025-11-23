#!/bin/bash

set -euo pipefail

# ==================== PATHS ======================
source "config/default_paths.sh"
[ -f "config/local_paths.sh" ] && source "config/local_paths.sh"

SLURM="$WORK/scripts/cellranger/slurm"

# Common sbatch options
SBATCH_OPTS="--parsable"

# ==================== PIPELINE ====================

# Step 1: Build reference genomes
echo "Submitting reference genome building..."
JOB1=$(sbatch $SBATCH_OPTS \
  --export=ALL \
  $SLURM/01_mkref.sbatch)
echo "  Job ID: $JOB1"

# Step 2: Run Cell Ranger count on all samples
echo "Submitting Cell Ranger count..."
JOB2=$(sbatch $SBATCH_OPTS \
  --export=ALL \
  --dependency=afterok:$JOB1 \
  $SLURM/02_cellranger.sbatch)
echo "  Job ID: $JOB2"

# Step 3: Aggregate samples (optional)
echo "Submitting aggregation..."
JOB3=$(sbatch $SBATCH_OPTS \
  --export=ALL \
  --dependency=afterok:$JOB2 \
  $SLURM/03_aggr.sbatch)
echo "  Job ID: $JOB3"

echo "Cell Ranger pipeline submitted successfully!"
