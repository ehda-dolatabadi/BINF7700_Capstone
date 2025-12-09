#!/bin/bash
# Script: 03_aggr.sbatch
# Purpose: Aggregate multiple Cell Ranger count outputs into a single dataset

#SBATCH --job-name=aggr
#SBATCH --partition=courses
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=180G
#SBATCH --array=0-3
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

# Exit on error, undefined variables, and pipe failures
set -euo pipefail

# CSV file listing samples to aggregate
CSV="aggr_samples.csv"

# Navigate to submission directory and load configuration
cd "$SLURM_SUBMIT_DIR"
source "config/default_paths.sh"
[ -f "config/local_paths.sh" ] && source "config/local_paths.sh"

# Redirect output to log file
exec > "$LOG/${SLURM_JOB_NAME}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"
exec 2>&1

# Define normalization options
NORM_FLAGS=("none" "mapped")
SUFFIXES=("raw" "normalized")

# Calculate indices for this array task
# Each reference genome is processed with each normalization method
ref_idx=$(( SLURM_ARRAY_TASK_ID / 2 ))
norm_idx=$(( SLURM_ARRAY_TASK_ID % 2 ))

# Get configuration for this task
REF=${REF_names[$ref_idx]}
NORM=${NORM_FLAGS[$norm_idx]}
SUFFIX=${SUFFIXES[$norm_idx]}

# Set up output directory
OUT_cellranger="$OUT/cellranger/ref_${REF}"
cd "$OUT_cellranger"

# Run Cell Ranger aggr
# - Combines multiple samples into single dataset
# - Applies normalization if specified
# - Generates combined feature-barcode matrix
cellranger aggr \
  --output-dir="$OUT_cellranger" \
  --id="${CSV%.*}_${SUFFIX}" \
  --csv="$CSV" \
  --localcores="$SLURM_CPUS_PER_TASK" \
  --localmem=$(( SLURM_MEM_PER_NODE / 1024 - 20)) \
  --disable-ui \
  --normalize="$NORM"
