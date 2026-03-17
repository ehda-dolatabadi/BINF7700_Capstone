#!/bin/bash
# Script: 03_cluster.sbatch
# Purpose: Perform dimensionality reduction, clustering, and marker identification

#SBATCH --job-name=cluster
#SBATCH --partition=courses
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=56
#SBATCH --mem=128G
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

# Exit on error, undefined variables, and pipe failures
set -euo pipefail

# Parameter defaults
ID="${ID:-all}"
NPCS="${NPCS:-50}"
DIMS="${DIMS:-10}"
RES="${RES:-0.5}"
UMAP_DIMS="${UMAP_DIMS:-10}"
UMAP_METRIC="${UMAP_METRIC:-cosine}"
SIGNIFICANCE="${SIGNIFICANCE:-0.05}"
REGULATION="${REGULATION:-1}"
ENRICHMENT="${ENRICHMENT:-0.2}"
TOP_MARKERS="${TOP_MARKERS:-100}"

# Load configuration and setup paths
cd "$SLURM_SUBMIT_DIR"
source "config/default_paths.sh"
[ -f "config/local_paths.sh" ] && source "config/local_paths.sh"

# Redirect output to log file
exec > "$LOG/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.log"
exec 2>&1

# Define script and output directories
RUN="$WORK/scripts/seurat/Rscripts"
OUT_integrate="$OUT/seurat/${ID}_05_integrated"
OUT_pca="$OUT/seurat/${ID}_06_pca"
OUT_cluster="$OUT/seurat/${ID}_07_clustered"
OUT_umap="$OUT/seurat/${ID}_08_umap"
OUT_markers="$OUT/seurat/${ID}_09_markers"

# Create output directories
mkdir -p "$OUT_pca" "$OUT_cluster" "$OUT_umap" "$OUT_markers"

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate seurat

# Clustering workflow
# Step 1: Perform principal component analysis for dimensionality reduction
echo "Running PCA..."
Rscript "$RUN/06_pca.R" \
  "$ID" \
  "$OUT_pca" \
  "$OUT_integrate/${ID}_integrated.rds" \
  "$NPCS"
echo "PCA completed✅"

# Step 2: Identify cell clusters using graph-based clustering (Leiden algorithm)
echo "Clustering..."
Rscript "$RUN/07_cluster.R" \
  "$ID" \
  "$OUT_cluster" \
  "$OUT_pca/${ID}_processed.rds" \
  "$DIMS" \
  "$RES"
echo "Clustering completed✅"

# Step 3: Generate UMAP visualization for cluster interpretation
echo "Running UMAP..."
Rscript "$RUN/08_umap.R" \
  "$ID" \
  "$OUT_umap" \
  "$OUT_cluster/${ID}_processed.rds" \
  "$UMAP_DIMS" \
  "$UMAP_METRIC"
echo "UMAP completed✅"

# Step 4: Find differentially expressed marker genes for each cluster
echo "Finding Markers..."
Rscript "$RUN/09_find_markers.R" \
  "$ID" \
  "$OUT_markers" \
  "$OUT_cluster/${ID}_processed.rds" \
  "$SIGNIFICANCE" \
  "$REGULATION" \
  "$ENRICHMENT" \
  "$TOP_MARKERS"
echo "Markers finding completed✅"
