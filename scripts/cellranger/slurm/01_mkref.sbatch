#!/bin/bash
# Script: 01_mkref.sbatch
# Purpose: Build Cell Ranger reference genomes from FASTA and GTF files

#SBATCH --job-name=mkref
#SBATCH --partition=courses
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=180G
#SBATCH --array=0-1
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

# Exit on error, undefined variables, and pipe failures
set -euo pipefail

# Navigate to submission directory and load configuration
cd "$SLURM_SUBMIT_DIR"
source "config/default_paths.sh"
[ -f "config/local_paths.sh" ] && source "config/local_paths.sh"

# Redirect output to log file
exec > "$LOG/${SLURM_JOB_NAME}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"
exec 2>&1

# Get reference files for this array task
FA="$REF_FILES/${FA_names[$SLURM_ARRAY_TASK_ID]}"
GTF="$REF_FILES/${GTF_names[$SLURM_ARRAY_TASK_ID]}"
REF="${REF_names[$SLURM_ARRAY_TASK_ID]}"

# Set up output directory
OUT_ref="$DATA/ref_genomes/cellranger_mkref"
mkdir -p "$OUT_ref"
cd "$OUT_ref"

# Run Cell Ranger mkref to build reference genome
# - Creates indexed reference for transcript alignment
# - Automatically allocates memory (total minus 20GB buffer)
cellranger mkref \
  --output-dir="$OUT_ref/$REF" \
  --fasta="$FA" \
  --genes="$GTF" \
  --genome="$REF" \
  --nthreads="$SLURM_CPUS_PER_TASK" \
  --memgb=$(( SLURM_MEM_PER_NODE / 1024 - 20))

# Clean up temporary build directory
rm -r ./"mkref_${REF}"/
