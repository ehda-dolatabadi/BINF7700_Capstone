#!/bin/bash
#SBATCH --job-name=aggr_normalize
#SBATCH --partition=courses
#SBATCH --time=24:00:00			# limit is 24hr
#SBATCH --cpus-per-task=28		# limit is 128/node
#SBATCH --mem=180G			# limit is 515G/node
#SBATCH --chdir=/scratch/dolatabadi.e   # working directory
#SBATCH --output=logs/%x_%j.out         # stdout log
#SBATCH --error=logs/%x_%j.err          # stderr log
#SBATCH --array=0-5

set -euo pipefail

# Arrays
REFS=("UKY_AmexF1_1_genomic" "AmexT_v47-AmexG_v6_0-DD")
CSVS=("aggr_samples.csv" "aggr_samples_no3h.csv" "aggr_samples_no3h1d.csv")

# Derive indices
ref_idx=$((SLURM_ARRAY_TASK_ID / ${#CSVS[@]}))
smp_idx=$((SLURM_ARRAY_TASK_ID % ${#CSVS[@]}))

REF=${REFS[$ref_idx]}
CSV=${SMPS[$smp_idx]}

# Working directory
OUT="$DATA/outputs/ref_${REF}/cellranger"
cd "$OUT"

# Run Cell Ranger on scratch
cellranger aggr \
  --id="$CSV" \
  --csv="$CSV" \
  --localcores="$SLURM_CPUS_PER_TASK" \
  --localmem=$(( SLURM_MEM_PER_NODE / 1024 - 20)) \
  --disable-ui
