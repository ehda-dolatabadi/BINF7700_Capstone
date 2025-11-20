#!/bin/bash
#SBATCH --job-name=overcluster
#SBATCH --partition=courses
#SBATCH --time=24:00:00			# limit is 24hr
#SBATCH --cpus-per-task=16		# limit is 128/node
#SBATCH --mem=32G			# limit is 515G/node
#SBATCH --chdir=/scratch/dolatabadi.e	# working directory
#SBATCH --output=logs/%x_%j.out		# stdout log
#SBATCH --error=logs/%x_%j.err		# stderr log

set -euo pipefail

# Paths
REF="UKY_AmexF1_1_genomic"
OUT="$DATA/outputs/ref_${REF}/seurat"

mkdir -p "$OUT"

# Environment
eval "$(conda shell.bash hook)"
conda activate seurat

# Run marker-based score
echo "Overclustering..."
Rscript "$RUN/seurat/12_overcluster.R" "12" "$OUT" "$OUT/05_pca.rds"
echo "Overclustering completed✅"

# Sync
rsync -av --update --info=progress2 --stats --exclude 'cellranger/' --exclude '*.rds' "$DATA/outputs/" "$WORK/results/"
