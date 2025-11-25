#!/bin/bash
#SBATCH --job-name=subcells
#SBATCH --partition=courses
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

# Parameters
main_ID="${main_ID:-all}"
ID="${ID:-clean}"
CELL_TYPES="${CELL_TYPES:-Epithelial;Erythrocyte}"
CELL_MARKERS="${CELL_MARKERS:-KRT8,KRT18,KRT1,EPCAM,CLDN7,CDH1,DSP,DSG2,EPPK1,S100P;HBA1,HBB,HBA2,HEMGN,ALAS2,PRDX2,ANK1,HBBP1}"
CELL_THRESHOLDS="${CELL_THRESHOLDS:-0.5;3.0}"

# Paths
cd "$SLURM_SUBMIT_DIR"
source "config/default_paths.sh"
[ -f "config/local_paths.sh" ] && source "config/local_paths.sh"

exec > "$LOG/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.log"
exec 2>&1

RUN="$WORK/scripts/seurat/Rscripts"
OUT_integrate="$OUT/seurat/${main_ID}_05_integrated"
OUT_subcells="$OUT/seurat/${ID}_05_integrated"

mkdir -p "$OUT_subcells"

# Environment
eval "$(conda shell.bash hook)"
conda activate seurat

# Workflow
echo "Subsetting data..."
Rscript "$RUN/12_subset_cells.R" \
  "$ID" \
  "$OUT_subcells" \
  "$OUT_integrate/${main_ID}_05_integrated.rds" \
  "$CELL_TYPES" \
  "$CELL_MARKERS" \
  "$CELL_THRESHOLDS"
echo "Subset completed✅"

echo "Re-running normalization and integration..."
Rscript "$RUN/05_reintegrate_subset.R" \
  "$ID" \
  "$OUT_subcells" \
  "$OUT_subcells/${ID}_00_subset.rds"
echo "Normalization and integration completed✅"

# Sync
rsync -av --info=progress2 --stats --exclude '*.rds' "$OUT/seurat/" "$WORK/results/seurat/"
