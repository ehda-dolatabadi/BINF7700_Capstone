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
ID="${ID:-enriched}"
CELL_THRESHOLDS="${CELL_THRESHOLDS:-5.0}"

# Paths
cd "$SLURM_SUBMIT_DIR"
source "config/default_paths.sh"
[ -f "config/local_paths.sh" ] && source "config/local_paths.sh"

exec > "$LOG/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.log"
exec 2>&1

# Source marker file
SUBSET_FILE="${SUBSET_FILE:-$WORK/scripts/seurat/abundant_cells.txt}"
source ${SUBSET_FILE}
cells=($(grep -E '^[a-zA-Z_]+=' ${SUBSET_FILE} | cut -d'=' -f1))

# Build semicolon-separated strings
CELL_TYPES=$(IFS=';'; echo "${cells[*]}")

CELL_MARKERS=""
for cell in "${cells[@]}"; do
    CELL_MARKERS="${CELL_MARKERS}${!cell};"
done
CELL_MARKERS="${CELL_MARKERS%;}"  # Remove trailing semicolon

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
  "$OUT_integrate/${main_ID}_integrated.rds" \
  "$CELL_TYPES" \
  "$CELL_MARKERS" \
  "$CELL_THRESHOLDS"
echo "Subset completed✅"

echo "Re-running normalization and integration..."
Rscript "$RUN/05_reintegrate_subset.R" \
  "$ID" \
  "$OUT_subcells" \
  "$OUT_subcells/${ID}_subset.rds"
echo "Normalization and integration completed✅"

# Sync
rsync -av --info=progress2 --stats --exclude '*.rds' "$OUT/seurat/" "$WORK/results/seurat/"
