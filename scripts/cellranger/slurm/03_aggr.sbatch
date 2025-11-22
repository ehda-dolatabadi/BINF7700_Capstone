#!/bin/bash
#SBATCH --job-name=aggr
#SBATCH --partition=courses
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=28
#SBATCH --mem=180G
#SBATCH --array=0-1

set -euo pipefail

# Arrays
REFS=("UKY_AmexF1_1_genomic" "AmexT_v47-AmexG_v6_0-DD")
CSV="aggr_samples.csv"

REF=${REFS[$SLURM_ARRAY_TASK_ID]}

# Working directory
OUT="$DATA/outputs/ref_${REF}/cellranger"
cd "$OUT"

# Run Cell Ranger on scratch
cellranger aggr \
  --id="${CSV}_raw" \
  --csv="$CSV" \
  --localcores="$SLURM_CPUS_PER_TASK" \
  --localmem=$(( SLURM_MEM_PER_NODE / 1024 - 20)) \
  --disable-ui \
  --normalize=none
