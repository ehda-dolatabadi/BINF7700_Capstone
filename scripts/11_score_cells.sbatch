#!/bin/bash
#SBATCH --job-name=cellScore
#SBATCH --partition=courses
#SBATCH --time=24:00:00			# limit is 24hr
#SBATCH --cpus-per-task=2		# limit is 128/node
#SBATCH --mem=4G			# limit is 515G/node
#SBATCH --chdir=/scratch/dolatabadi.e	# working directory
#SBATCH --output=logs/%x_%j.out		# stdout log
#SBATCH --error=logs/%x_%j.err		# stderr log
#SBATCH --array=0-0

set -euo pipefail

# Arrays
REFS=("UKY_AmexF1_1_genomic" "AmexT_v47-AmexG_v6_0-DD")
REF=${REFS[$SLURM_ARRAY_TASK_ID]}

# Paths
OUT="$DATA/outputs/ref_${REF}/seurat"

mkdir -p "$OUT"

# Environment
eval "$(conda shell.bash hook)"
conda activate seurat

# Run marker-based score
echo "Marker-based scoring of cell types..."
Rscript "$RUN/seurat/11_score_cells.R" "11" "$OUT" "$OUT/04_integrated.rds"
echo "Scoring completed✅"

# Sync
rsync -av --update --info=progress2 --stats --exclude 'cellranger/' --exclude '*.rds' "$DATA/outputs/" "$WORK/results/"
