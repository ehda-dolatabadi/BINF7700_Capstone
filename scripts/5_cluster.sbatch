#!/bin/bash
#SBATCH --job-name=cluster
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
REF=${REFS[$SLURM_ARRAY_TASK_ID]}

# Paths
DIR="$DATA/outputs/ref_${REF}/seurat/normalized"
OUT="$DATA/outputs/ref_${REF}/seurat"

mkdir -p "$OUT"

# Environment
eval "$(conda shell.bash hook)"
conda activate seurat

# Run integrate
#Rscript "$RUN/seurat/4_integrate.R" "4" "$OUT" $DIR/*_normalized.rds

# Run PCA
#Rscript "$RUN/seurat/5_pca.R" "5" "$OUT" "$OUT/4_integrated.rds"

# Run cluster
#Rscript "$RUN/seurat/6_cluster.R" "6" "$OUT" "$OUT/5_pca.rds"

# Run UMAP
#Rscript "$RUN/seurat/7_umap.R" "7" "$OUT" "$OUT/6_clustered.rds"

# Find markers
Rscript "$RUN/seurat/8_find_markers.R" "8" "$OUT" "$OUT/6_clustered.rds"
