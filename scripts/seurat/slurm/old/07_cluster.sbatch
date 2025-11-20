#!/bin/bash
#SBATCH --job-name=cluster
#SBATCH --partition=courses
#SBATCH --time=24:00:00			# limit is 24hr
#SBATCH --cpus-per-task=64		# limit is 128/node
#SBATCH --mem=128G			# limit is 515G/node
#SBATCH --chdir=/scratch/dolatabadi.e	# working directory
#SBATCH --output=logs/%x_%j.out		# stdout log
#SBATCH --error=logs/%x_%j.err		# stderr log
#SBATCH --array=0-0

set -euo pipefail

# Arrays
REFS=("UKY_AmexF1_1_genomic" "AmexT_v47-AmexG_v6_0-DD")
REF=${REFS[$SLURM_ARRAY_TASK_ID]}

# Paths
DIR="$DATA/outputs/ref_${REF}/seurat/03_normalized"
OUT="$DATA/outputs/ref_${REF}/seurat"

mkdir -p "$OUT"

# Environment
eval "$(conda shell.bash hook)"
conda activate seurat

# Run integrate
echo "Integrting normalized samples..."
Rscript "$RUN/seurat/04_integrate.R" "04" "$OUT" "${SLURM_CPUS_PER_TASK}" $DIR/*_normalized.rds
echo "Integration completed✅"

# Run PCA
echo "Running PCA..."
Rscript "$RUN/seurat/05_pca.R" "05" "$OUT" "$OUT/04_integrated.rds" "${SLURM_CPUS_PER_TASK}"
echo "PCA completed✅"

# Run cluster
echo "Clustering..."
Rscript "$RUN/seurat/06_cluster.R" "06" "$OUT" "$OUT/05_pca.rds" "${SLURM_CPUS_PER_TASK}"
echo "Clustering completed✅"

# Run UMAP
echo "Running UMAP..."
Rscript "$RUN/seurat/07_umap.R" "07" "$OUT" "$OUT/06_clustered.rds" "${SLURM_CPUS_PER_TASK}"
echo "UMAP completed✅"

# Find markers
echo "Finding Markers..."
Rscript "$RUN/seurat/08_find_markers.R" "08" "$OUT" "$OUT/06_clustered.rds" "${SLURM_CPUS_PER_TASK}"
echo "Markers finding completed✅"

# Sync
rsync -av --update --info=progress2 --stats --exclude 'cellranger/' --exclude '*.rds' "$DATA/outputs/" "$WORK/results/"
