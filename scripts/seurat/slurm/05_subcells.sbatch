#!/bin/bash
#SBATCH --job-name=subcells
#SBATCH --partition=courses
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --chdir=/scratch/dolatabadi.e   # log directory
#SBATCH --output=logs/%x_%j.out         # stdout log
#SBATCH --error=logs/%x_%j.err          # stderr log

set -euo pipefail

# Environment
eval "$(conda shell.bash hook)"
conda activate seurat

# Workflow
echo "Subsetting data..."
#Rscript "$WORK/scripts/seurat/11_subset_cells.R"	"$ID"   "$OUT"  "$OUT/${main_ID}_04_integrated.rds"
echo "Subset completed✅"

echo "Re-running normalization and integration..."
Rscript "$WORK/scripts/seurat/04_reintegrate_subset.R"	"$ID"   "$OUT"  "$OUT/${ID}_00_subset.rds"
echo "Normalization and integration completed✅"

# Sync
rsync -av --update --info=progress2 --stats --exclude '*.rds' "$OUT/" "$WORK/results/seurat/"
