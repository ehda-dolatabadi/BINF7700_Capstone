#!/bin/bash
#SBATCH --job-name=integrate
#SBATCH --partition=courses
#SBATCH --time=2:00:00			# limit is 24hr
#SBATCH --cpus-per-task=4		# limit is 128/node
#SBATCH --mem=8G			# limit is 515G/node
#SBATCH --chdir=/scratch/dolatabadi.e	# working directory
#SBATCH --output=logs/%x_%j.out		# stdout log
#SBATCH --error=logs/%x_%j.err		# stderr log
#SBATCH --array=0-1

set -euo pipefail

# Arrays
REFS=("UKY_AmexF1_1_genomic" "AmexT_v47-AmexG_v6_0-DD")

# Paths
REF=${REFS[$SLURM_ARRAY_TASK_ID]}
DIR="$DATA/outputs/ref_${REF}/seurat/normalize"
OUT="$DATA/outputs/ref_${REF}/seurat/normalize"

mkdir -p "$OUT"

# Environment
eval "$(conda shell.bash hook)"
conda activate seurat

# Run integrate
Rscript "$RUN/integrate_sct.R" "$REF" "$OUT" $DIR/*_sct.rds
