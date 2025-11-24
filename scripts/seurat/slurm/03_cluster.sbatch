#!/bin/bash
#SBATCH --job-name=cluster
#SBATCH --partition=courses
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

# Parameters
ID="${ID:-main}"
NPCS="${NPCS:-50}"
DIMS="${DIMS:-10}"
RES="${RES:-0.5}"
UMAP_DIMS="${UMAP_DIMS:-10}"
UMAP_METRIC="${UMAP_METRIC:-cosine}"
SIGNIFICANCE="${SIGNIFICANCE:-0.05}"
REGULATION="${REGULATION:-1}"
ENRICHMENT="${ENRICHMENT:-0.2}"
TOP_MARKERS="${TOP_MARKERS:-30}"

# Paths
cd "$SLURM_SUBMIT_DIR"
source "config/default_paths.sh"
[ -f "config/local_paths.sh" ] && source "config/local_paths.sh"

exec > "$LOG/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.log"
exec 2>&1

RUN="$WORK/scripts/seurat/Rscripts"
OUT_integrate="$OUT/seurat/${ID}_05_integrated"
OUT_pca="$OUT/seurat/${ID}_06_pca"
OUT_cluster="$OUT/seurat/${ID}_07_clustered"
OUT_umap="$OUT/seurat/${ID}_08_umap"
OUT_markers="$OUT/seurat/${ID}_09_markers"

mkdir -p "$OUT_pca" "$OUT_cluster" "$OUT_umap" "$OUT_markers"

# Environment
eval "$(conda shell.bash hook)"
conda activate seurat

# Workflow
echo "Running PCA..."
Rscript "$RUN/06_pca.R" \
  "$ID" \
  "$OUT_pca" \
  "$OUT_integrate/${ID}_05_integrated.rds" \
  "$NPCS"
echo "PCA completed✅"

echo "Clustering..."
Rscript "$RUN/07_cluster.R" \
  "$ID" \
  "$OUT_cluster" \
  "$OUT_pca/${ID}_06_pca.rds" \
  "$DIMS" \
  "$RES"
echo "Clustering completed✅"

echo "Running UMAP..."
Rscript "$RUN/08_umap.R" \
  "$ID" \
  "$OUT_umap" \
  "$OUT_cluster/${ID}_07_clustered.rds" \
  "$UMAP_DIMS" \
  "$UMAP_METRIC"
echo "UMAP completed✅"

echo "Finding Markers..."
Rscript "$RUN/09_find_markers.R" \
  "$ID" \
  "$OUT_markers" \
  "$OUT_cluster/${ID}_07_clustered.rds" \
  "$SIGNIFICANCE" \
  "$REGULATION" \
  "$ENRICHMENT" \
  "$TOP_MARKERS"
echo "Markers finding completed✅"

# Sync
rsync -av --info=progress2 --stats --exclude '*.rds' "$OUT/seurat/" "$WORK/results/seurat/"
