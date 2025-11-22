#!/bin/bash
#SBATCH --job-name=cluster
#SBATCH --partition=courses
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G

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
CELL_NAME="${CELL_NAME:-Neural}"
GROUP_BY="${GROUP_BY:-seurat_clusters}"
MARKERS="${MARKERS:-SOX2,NES,TUBB3,MAP2,S100B,GFAP,VIM,NCAM1,CD24,FABP7}"

# Paths
source "../../../config/default_paths.sh"
if [ -f "../../../config/local_paths.sh" ]; then
    source "../../../config/local_paths.sh"
fi

RUN="$WORK/scripts/seurat/Rscripts"

# Environment
eval "$(conda shell.bash hook)"
conda activate seurat

# Workflow
echo "Running PCA..."
Rscript "$RUN/05_pca.R" \
  "$ID" \
  "$OUT/seurat" \
  "$OUT/seurat/${ID}_04_integrated.rds" \
  "$NPCS"
echo "PCA completed✅"

echo "Clustering..."
Rscript "$RUN/06_cluster.R" \
  "$ID" \
  "$OUT/seurat" \
  "$OUT/seurat/${ID}_05_pca.rds" \
  "$DIMS" \
  "$RES"
echo "Clustering completed✅"

echo "Running UMAP..."
Rscript "$RUN/07_umap.R" \
  "$ID" \
  "$OUT/seurat" \
  "$OUT/seurat/${ID}_06_clustered.rds" \
  "$UMAP_DIMS" \
  "$UMAP_METRIC"
echo "UMAP completed✅"

echo "Finding Markers..."
Rscript "$RUN/08_find_markers.R" \
  "$ID" \
  "$OUT/seurat" \
  "$OUT/seurat/${ID}_06_clustered.rds" \
  "$SIGNIFICANCE" \
  "$REGULATION" \
  "$ENRICHMENT" \
  "$TOP_MARKERS"
echo "Markers finding completed✅"

echo "Scoring $CELL_NAME markers..."
Rscript "$RUN/09_score_markers.R" \
  "$ID" \
  "$OUT/seurat" \
  "$OUT/seurat/${ID}_06_clustered.rds" \
  "$CELL_NAME" \
  "$GROUP_BY" \
  "$MARKERS"
echo "Scoring completed✅"

# Sync
rsync -av --info=progress2 --stats --exclude '*.rds' "$OUT/seurat/" "$WORK/results/seurat/"
