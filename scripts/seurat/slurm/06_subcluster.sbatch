#!/bin/bash
#SBATCH --job-name=subcluster
#SBATCH --partition=courses
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

# Pipeline parameters
main_ID="${main_ID:-all}"
ID="${ID:-cluster1}"
IDENT="${IDENT:-1}"

# Paths
cd "$SLURM_SUBMIT_DIR"
source "config/default_paths.sh"
[ -f "config/local_paths.sh" ] && source "config/local_paths.sh"

exec > "$LOG/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.log"
exec 2>&1

RUN="$WORK/scripts/seurat/Rscripts"
OUT_cluster="$OUT/seurat/${main_ID}_07_clustered"
OUT_subset="$OUT/seurat/${ID}_05_integrated"

mkdir -p "$OUT_subset"

# Environment
eval "$(conda shell.bash hook)"
conda activate seurat

# Workflow
echo "Subsetting data..."
Rscript "$RUN/11_subset_clusters.R"\
    "$ID"\
    "$OUT_subset"\
    "$OUT_cluster/${main_ID}_clustered.rds"\
    "$IDENT"
echo "Subset completed✅"

echo "Re-running normalization and integration..."
Rscript "$RUN/05_reintegrate_subset.R"\
    "$ID"\
    "$OUT_subset"\
    "$OUT_subset/${ID}_subset.rds"
echo "Normalization and integration completed✅"
