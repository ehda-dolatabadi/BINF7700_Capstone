#!/bin/bash
#SBATCH --job-name=qc
#SBATCH --partition=courses
#SBATCH --time=24:00:00                 # limit is 24hr
#SBATCH --cpus-per-task=4              	# limit is 128/node
#SBATCH --mem=32G                      	# limit is 515G/node
#SBATCH --chdir=/scratch/dolatabadi.e   # working directory
#SBATCH --output=logs/%x_%j.out         # stdout log
#SBATCH --error=logs/%x_%j.err          # stderr log
#SBATCH --array=0-15

set -euo pipefail

# Arrays
REFS=("UKY_AmexF1_1_genomic" "AmexT_v47-AmexG_v6_0-DD")
SMPS=("control" "3h" "24h" "72h" "7dpa" "14dpa" "22dpa" "33dpa")

# Derive indices
ref_idx=$((SLURM_ARRAY_TASK_ID / ${#SMPS[@]}))
smp_idx=$((SLURM_ARRAY_TASK_ID % ${#SMPS[@]}))

REF=${REFS[$ref_idx]}
SMP=${SMPS[$smp_idx]}

# Paths
DIR="$DATA/outputs/cellranger/ref_${REF}/$SMP/outs/filtered_feature_bc_matrix"
OUT="$DATA/outputs/ref_${REF}/qc/$SMP"

mkdir -p "$OUT"

# Environment
eval "$(conda shell.bash hook)"
conda activate seurat

# Run QC
Rscript "$RUN/qc.R" "$SMP" "$DIR" "$OUT"
