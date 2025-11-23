#!/bin/bash
#SBATCH --job-name=aggr
#SBATCH --partition=courses
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=180G
#SBATCH --array=0-3
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

# Samples list
CSV="aggr_samples.csv"

# Paths
cd "$SLURM_SUBMIT_DIR"
source "config/default_paths.sh"
[ -f "config/local_paths.sh" ] && source "config/local_paths.sh"

exec > "$LOG/${SLURM_JOB_NAME}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"
exec 2>&1

# Derive indices
NORM_FLAGS=("none" "mapped")
SUFFIXES=("raw" "normalized")

ref_idx=$(( SLURM_ARRAY_TASK_ID / 2 ))
norm_idx=$(( SLURM_ARRAY_TASK_ID % 2 ))

REF=${REF_names[$ref_idx]}
NORM=${NORM_FLAGS[$norm_idx]}
SUFFIX=${SUFFIXES[$norm_idx]}

OUT_cellranger="$OUT/cellranger/ref_${REF}"
cd "$OUT_cellranger"

# Run pipeline
cellranger aggr \
  --output-dir="$OUT_cellranger" \
  --id="${CSV%.*}_${SUFFIX}" \
  --csv="$CSV" \
  --localcores="$SLURM_CPUS_PER_TASK" \
  --localmem=$(( SLURM_MEM_PER_NODE / 1024 - 20)) \
  --disable-ui \
  --normalize="$NORM"
