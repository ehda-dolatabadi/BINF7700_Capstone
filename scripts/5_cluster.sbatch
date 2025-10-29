#!/bin/bash
#SBATCH --job-name=cluster
#SBATCH --partition=courses
#SBATCH --time=4:00:00			# limit is 24hr
#SBATCH --cpus-per-task=4		# limit is 128/node
#SBATCH --mem=8G			# limit is 515G/node
#SBATCH --chdir=/scratch/dolatabadi.e	# working directory
#SBATCH --output=logs/%x_%j.out		# stdout log
#SBATCH --error=logs/%x_%j.err		# stderr log
#SBATCH --array=0-1

set -euo pipefail

# Arrays
REFS=("UKY_AmexF1_1_genomic" "AmexT_v47-AmexG_v6_0-DD")

# Variables
ID="axolotl"
REF=${REFS[$SLURM_ARRAY_TASK_ID]}

# Paths
DIR="$DATA/outputs/ref_${REF}/seurat/normalized"
OUT="$DATA/outputs/ref_${REF}/seurat"

mkdir -p "$OUT"

# Environment
eval "$(conda shell.bash hook)"
conda activate seurat

# Run integrate
Rscript "$RUN/seurat/4_integrate.R" "$ID" "$OUT" $DIR/*_normalized.rds

# Run PCA
Rscript "$RUN/seurat/5_pca.R" "$ID" "$OUT" "$OUT/${ID}_integrated.rds"

# Run cluster
Rscript "$RUN/seurat/6_cluster.R" "$ID" "$OUT" "$OUT/${ID}_pca.rds"

# Run UMAP
Rscript "$RUN/seurat/7_umap.R" "$ID" "$OUT" "$OUT/${ID}_clustered.rds"
