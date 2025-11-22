#!/bin/bash
#SBATCH --job-name=aggr
#SBATCH --partition=courses
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --array=0-3

set -euo pipefail

# Samples list
CSV="aggr_samples.csv"

# Paths
source "../../../config/default_paths.sh"
if [ -f "../../../config/local_paths.sh" ]; then
    source "../../../config/local_paths.sh"
fi

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
  --id="${CSV%.*}_${SUFFIX}" \
  --csv="$CSV" \
  --localcores="$SLURM_CPUS_PER_TASK" \
  --localmem=$(( SLURM_MEM_PER_NODE / 1024 - 20)) \
  --disable-ui \
  --normalize="$NORM"
