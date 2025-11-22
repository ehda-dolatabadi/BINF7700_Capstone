#!/bin/bash
#SBATCH --job-name=subcells
#SBATCH --partition=courses
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G

set -euo pipefail

# Parameters
ID="${ID:-no_epith-eryth}"
CELL_TYPES="${CELL_TYPES:-Epithelial,Erythrocyte}"
CELL_MARKERS="${CELL_MARKERS:-KRT8,KRT18,KRT1,EPCAM,CLDN7,CDH1,DSP,DSG2,EPPK1,S100P;HBA1,HBB,HBA2,HEMGN,ALAS2,PRDX2,ANK1,HBBP1}"
CELL_THRESHOLDS="${CELL_THRESHOLDS:-0.5,3.0}"

# Environment
eval "$(conda shell.bash hook)"
conda activate seurat

# Workflow
echo "Subsetting data..."
Rscript "$RUN/11_subset_cells.R" \
  "$ID" \
  "$OUT" \
  "$OUT/${main_ID}_04_integrated.rds" \
  "$CELL_TYPES" \
  "$CELL_MARKERS" \
  "$CELL_THRESHOLDS"
echo "Subset completed✅"

echo "Re-running normalization and integration..."
Rscript "$RUN/04_reintegrate_subset.R" \
  "$ID" \
  "$OUT" \
  "$OUT/${ID}_00_subset.rds"
echo "Normalization and integration completed✅"

# Sync
rsync -av --update --info=progress2 --stats --exclude '*.rds' "$OUT/" "$WORK/results/seurat/"
