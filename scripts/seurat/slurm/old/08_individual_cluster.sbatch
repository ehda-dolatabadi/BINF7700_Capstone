#!/bin/bash
#SBATCH --job-name=individual
#SBATCH --partition=courses
#SBATCH --time=24:00:00			# limit is 24hr
#SBATCH --cpus-per-task=16		# limit is 128/node
#SBATCH --mem=32G			# limit is 515G/node
#SBATCH --chdir=/scratch/dolatabadi.e	# working directory
#SBATCH --output=logs/%x_%j.out		# stdout log
#SBATCH --error=logs/%x_%j.err		# stderr log
#SBATCH --array=0-7

set -euo pipefail

# Arrays
REFS=("UKY_AmexF1_1_genomic" "AmexT_v47-AmexG_v6_0-DD")
SMPS=("control" "3h" "24h" "72h" "7dpa" "14dpa" "22dpa" "33dpa")

# Derive indices
ref_idx=$((SLURM_ARRAY_TASK_ID / ${#SMPS[@]}))
smp_idx=$((SLURM_ARRAY_TASK_ID % ${#SMPS[@]}))

#REF=${REFS[$ref_idx]}
#SMP=${SMPS[$smp_idx]}

REF=${REFS[0]}
SMP=${SMPS[$SLURM_ARRAY_TASK_ID]}

# Paths
DIR="$DATA/outputs/ref_${REF}/seurat/03_normalized"
OUT="$DATA/outputs/ref_${REF}/seurat/timepoints/$SMP"

mkdir -p "$OUT"

# Environment
eval "$(conda shell.bash hook)"
conda activate seurat

# Run PCA
echo "Running PCA..."
Rscript "$RUN/seurat/05_pca.R" "05" "$OUT" $DIR/${SMP}_normalized.rds
echo "PCA completed✅"

# Run cluster
echo "Clustering..."
Rscript "$RUN/seurat/06_cluster.R" "06" "$OUT" "$OUT/05_pca.rds"
echo "Clustering completed✅"

# Run UMAP
echo "Running UMAP..."
Rscript "$RUN/seurat/07_umap.R" "07" "$OUT" "$OUT/06_clustered.rds"
echo "UMAP completed✅"

# Find markers
echo "Finding Markers..."
Rscript "$RUN/seurat/08_find_markers.R" "08" "$OUT" "$OUT/06_clustered.rds" "${SLURM_CPUS_PER_TASK}"
echo "Markers finding completed✅"

# Sync
rsync -av --update --info=progress2 --stats --exclude 'cellranger/' --exclude '*.rds' "$DATA/outputs/" "$WORK/results/"
