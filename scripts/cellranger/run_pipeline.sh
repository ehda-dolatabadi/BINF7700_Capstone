#!/bin/bash

set -euo pipefail

# ==================== PATHS ======================
SCRATCH="/scratch/dolatabadi.e"
DATA="/courses/BINF7700.202610/students/dolatabadi.e"
LOG="$SCRATCH/logs"
WORK="$HOME/axolotl-regeneration-scrna"
SLURM="$WORK/scripts/cellranger/slurm"

# ==================== EXPORTS ====================
export SCRATCH DATA LOG WORK SLURM

mkdir -p "$LOG"

# Common sbatch options
SBATCH_OPTS="--parsable --output=$LOG/%x_%j.out --error=$LOG/%x_%j.err"

# ==================== PIPELINE ====================

# Step 1: Build reference genomes
echo "Submitting reference genome building (mkref)..."
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

# Step 3a: Aggregate samples without normalization (optional)
echo "Submitting aggregation without normalization..."
JOB3=$(sbatch $SBATCH_OPTS \
  --export=ALL \
  --dependency=afterok:$JOB2 \
  $SLURM/03_aggr.sbatch)
echo "  Job ID: $JOB3"

# Step 3b: Aggregate samples with normalization (optional)
echo "Submitting aggregation with normalization..."
JOB4=$(sbatch $SBATCH_OPTS \
  --export=ALL \
  --dependency=afterok:$JOB2 \
  $SLURM/03_aggr_normalize.sbatch)
echo "  Job ID: $JOB4"

echo "Cell Ranger pipeline submitted successfully!"
