#!/bin/bash
#SBATCH --job-name=subcluster
#SBATCH --partition=courses
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G

set -euo pipefail

# Pipeline parameters
main_ID="${main_ID:-main}"
ID="${ID:-cluster1}"
IDENT="${IDENT:-1}"

# Paths
cd "$SLURM_SUBMIT_DIR"
source "config/default_paths.sh"
[ -f "config/local_paths.sh" ] && source "config/local_paths.sh"

exec > "$LOG/${SLURM_JOB_NAME}_${SLURM_ARRAY_JOB_ID}.log"
exec 2>&1

RUN="$WORK/scripts/seurat/Rscripts"

# Environment
eval "$(conda shell.bash hook)"
conda activate seurat

# Workflow
echo "Subsetting data..."
Rscript "$RUN/10_subset_clusters.R"	"$ID"   "$OUT/seurat"  "$OUT/seurat/${main_ID}_06_clustered.rds" "$IDENT"
echo "Subset completed✅"

echo "Re-running normalization and integration..."
Rscript "$RUN/04_reintegrate_subset.R"  "$ID"   "$OUT/seurat"  "$OUT/seurat/${ID}_00_subset.rds"
echo "Normalization and integration completed✅"
