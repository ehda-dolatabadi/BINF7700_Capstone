#!/bin/bash

set -euo pipefail

# Common variables
export main_ID="all"
export REF="UKY_AmexF1_1_genomic"

# Source the marker file
export MARKER_FILE="$(pwd)/scripts/seurat/cell_markers.txt"
n_marker=$(grep -cE '^[a-zA-Z_]+=' ${MARKER_FILE})

export SUBSET_FILE="$(pwd)/scripts/seurat/abundant_cells.txt"
n_subset=$(grep -cE '^[a-zA-Z_]+=' ${SUBSET_FILE})

# Source paths
source "config/default_paths.sh"
[ -f "config/local_paths.sh" ] && source "config/local_paths.sh"

SLURM="$WORK/scripts/seurat/slurm"

# Common sbatch options
SBATCH_OPTS="--parsable"

# ==================== PIPELINE ====================
# Jobs to submit
preprocess=false
integrate=false
subset=false
cluster=false
score=true
subcluster=false



# Step 1: Preprocessing
if [ "$preprocess" = true ]; then
    echo "Submitting preprocessing..."
    JOB1=$(sbatch $SBATCH_OPTS \
      --export=ALL,\
ID=$main_ID,\
MIN_CELLS=5,\
MIN_FEATURES=500,\
MIN_COUNTS=1000,\
MAX_COUNTS=30000,\
MIN_FEATURES_FILTER=500,\
MAX_FEATURES=5000,\
MAX_MT=10,\
MAX_RIBO=35 \
      $SLURM/01_preprocess.sbatch)
    echo "  Job ID: $JOB1"
fi



# Step 2: Integration
if [ "$integrate" = true ]; then
    echo "Submitting integration..."
    JOB2=$(sbatch $SBATCH_OPTS \
      --export=ALL,\
ID=$main_ID \
      $([ "$preprocess" = true ] && echo "--dependency=afterok:$JOB1") \
      $SLURM/02_integrate.sbatch)
    echo "  Job ID: $JOB2"
fi



# Optional step: Removing abundant cells
if [ "$subset" = true ]; then
    echo "Submitting subcells removal..."
    JOB3=$(sbatch $SBATCH_OPTS \
      --export=ALL,\
ID="enriched1",\
SUBSET_FILE="$(pwd)/scripts/seurat/cell_markers.txt",\
CELL_THRESHOLDS="1" \
      $([ "$integrate" = true ] && echo "--dependency=afterok:$JOB2") \
      $SLURM/06_subcells.sbatch)
    echo "  Job ID: $JOB3"
fi



# Step 3: Clustering
if [ "$cluster" = true ]; then
    echo "Submitting clustering..."
    JOB4=$(sbatch $SBATCH_OPTS \
      --export=ALL,\
ID="$([ "$subset" = true ] && echo "enriched" || echo "$main_ID")",\
NPCS=50,\
DIMS=10,\
RES=0.5,\
UMAP_DIMS=10,\
UMAP_METRIC="cosine",\
SIGNIFICANCE=0.05,\
REGULATION=1,\
ENRICHMENT=0.2,\
TOP_MARKERS=100 \
      $([ "$subset" = true ] && echo "--dependency=afterok:$JOB3" \
		|| ([ "$integrate" = true ] && echo "--dependency=afterok:$JOB2")) \
      $SLURM/03_cluster.sbatch)
    echo "  Job ID: $JOB4"
fi



# Optional step: Marker-based scoring
if [ "$score" = true ]; then
    echo "Submitting marker-based scoring..."
    JOB5=$(sbatch $SBATCH_OPTS \
      --array=0-$((n_marker-1)) \
      --export=ALL,\
main_ID="$([ "$subset" = true ] && echo "enriched" || echo "$main_ID")",\
MARKER_FILE="$(pwd)/scripts/seurat/cell_markers.txt" \
      $([ "$cluster" = true ] && echo "--dependency=afterok:$JOB4" \
		|| ([ "$subset" = true ] && echo "--dependency=afterok:$JOB3" \
		|| ([ "$integrate" = true ] && echo "--dependency=afterok:$JOB2"))) \
      $SLURM/04_score_markers.sbatch)
    echo "  Job ID: $JOB5"
fi



# Optional step: Subclustering
if [ "$subcluster" = true ]; then
    echo "Submitting subclustering..."
    JOB6=$(sbatch $SBATCH_OPTS \
      --export=ALL,\
main_ID="$([ "$subset" = true ] && echo "enriched" || echo "$main_ID")",\
ID="cluster16",\
IDENT="16" \
      $([ "$cluster" = true ] && echo "--dependency=afterok:$JOB4") \
      $SLURM/05_subcluster.sbatch)
    echo "  Job ID: $JOB6"

    JOB7=$(sbatch $SBATCH_OPTS \
      --export=ALL,\
main_ID="$([ "$subset" = true ] && echo "enriched" || echo "$main_ID")",\
ID="cluster16",\
NPCS=50,\
DIMS=10,\
RES=1.5,\
UMAP_DIMS=10,\
UMAP_METRIC="cosine",\
SIGNIFICANCE=0.05,\
REGULATION=0.5,\
ENRICHMENT=0.1,\
TOP_MARKERS=100\
      --job-name=cluster_subcluster \
      --dependency=afterok:$JOB6 \
      $SLURM/03_cluster.sbatch)
    echo "  Job ID: $JOB7"
fi

echo "Pipeline submitted successfully!"
