#!/bin/bash
# Script: 02_cellranger.sbatch
# Purpose: Run Cell Ranger count on FASTQ files to quantify gene expression

#SBATCH --job-name=cellranger
#SBATCH --partition=courses
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=240G
#SBATCH --array=0-15
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

# Reference genomes output directory
OUT_ref="$DATA/ref_genomes/cellranger_mkref"

# Calculate indices for this array task
# Each reference genome is paired with each FASTQ sample
ref_idx=$((SLURM_ARRAY_TASK_ID / ${#FQ_names[@]}))
fq_idx=$((SLURM_ARRAY_TASK_ID % ${#FQ_names[@]}))

# Get reference and sample information
REF_name=${REF_names[$ref_idx]}
FQ_name="$FQ_src/${FQ_names[$fq_idx]}"
SMP=${SMPS[$fq_idx]}

# Set up paths
REF="$OUT_ref/$REF_name"
OUT_cellranger="$OUT/cellranger/ref_${REF_name}"
mkdir -p "$OUT_cellranger"
cd "$OUT_cellranger"

# Run Cell Ranger count
# - Aligns FASTQ reads to reference genome
# - Quantifies gene expression per cell
# - Generates count matrix and web summary
# - Skip BAM file creation and secondary analysis to save space/time
cellranger count \
  --output-dir="$OUT_cellranger" \
  --id="$SMP" \
  --transcriptome="$REF" \
  --fastqs="$FQ_name" \
  --sample="$SMP" \
  --create-bam=false \
  --nosecondary \
  --localcores="$SLURM_CPUS_PER_TASK" \
  --localmem=$(( SLURM_MEM_PER_NODE / 1024 - 20))

# Copy web summary to output directory for easy access
rsync -av "$OUT_cellranger/$SMP/outs/web_summary.html" "$OUT_cellranger/${SMP}_summary.html"
