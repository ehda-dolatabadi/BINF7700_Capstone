#!/bin/bash
#SBATCH --job-name=preprocess
#SBATCH --partition=courses
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --chdir=/scratch/dolatabadi.e   # log directory
#SBATCH --output=logs/%x_%j.out         # stdout log
#SBATCH --error=logs/%x_%j.err          # stderr log

set -euo pipefail

# Variables
ID=${1:-""}  # Use first argument, or "" if not provided

# Environment
eval "$(conda shell.bash hook)"
conda activate seurat

# Workflow
echo "Subsetting data..."
Rscript "$WORK/scripts/seurat/11_subset_cells.R"	"$ID"   "$OUT"  "$OUT/_04_integrated.rds" "$IDENT"
echo "Subset completed✅"

echo "Re-running normalization and integration..."
Rscript "$WORK/scripts/seurat/04_reintegrate_subset.R"	"$ID"   "$OUT"  "$OUT/${ID}_00_subset.rds"
echo "Normalization and integration completed✅"

# Sync
rsync -av --update --info=progress2 --stats --exclude 'cellranger/' --exclude '*.rds' "$DATA/outputs/" "$WORK/results/"
