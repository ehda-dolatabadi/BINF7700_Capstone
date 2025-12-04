#!/bin/bash

set -euo pipefail

# Common variables
export main_ID="all"
export REF="UKY_AmexF1_1_genomic"

# Source the marker file
export MARKER_FILE="$(pwd)/scripts/seurat/cell_markers.txt"
n_marker=$(grep -cE '^[a-zA-Z_]+=' ${MARKER_FILE})

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
cluster_all=false
score_all=false

enrich1=true
enrich2=false
enrich3=true

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




# Step 3: Clustering
if [ "$cluster_all" = true ]; then
    echo "Submitting clustering..."
    JOB3=$(sbatch $SBATCH_OPTS \
      --export=ALL,\
ID=$main_ID,\
NPCS=50,\
DIMS=10,\
RES=0.5,\
UMAP_DIMS=10,\
UMAP_METRIC="cosine",\
SIGNIFICANCE=0.05,\
REGULATION=1,\
ENRICHMENT=0.2,\
TOP_MARKERS=100 \
      $([ "$integrate" = true ] && echo "--dependency=afterok:$JOB2") \
      $SLURM/03_cluster.sbatch)
    echo "  Job ID: $JOB3"
fi



# Step 4: Marker-based scoring
if [ "$score_all" = true ]; then
    echo "Submitting marker-based scoring..."
    JOB4=$(sbatch $SBATCH_OPTS \
      --array=0-$((n_marker-1)) \
      --export=ALL,\
main_ID=$main_ID,\
MARKER_FILE="$(pwd)/scripts/seurat/cell_markers.txt" \
      $([ "$cluster_all" = true ] && echo "--dependency=afterok:$JOB3") \
      $SLURM/04_score_markers.sbatch)
    echo "  Job ID: $JOB4"
fi



# Step 5: Removing abundant cells
if [ "$enrich1" = true ]; then
    echo "Submitting subcells removal..."
    JOB5=$(sbatch $SBATCH_OPTS \
      --export=ALL,\
ID="enriched1",\
main_ID=$main_ID,\
MODE="remove",\
MARKER_FILE="$(pwd)/scripts/seurat/cell_markers.txt",\
CELL_TYPES="Epithelial;Erythrocytes;Fibroblasts",\
CELL_THRESHOLDS="0.2;0.5;0.5" \
      --job-name=subset1 \
      $([ "$integrate" = true ] && echo "--dependency=afterok:$JOB2") \
      $SLURM/05_subcells.sbatch)
    echo "  Job ID: $JOB5"

    echo "Submitting clustering..."
    JOB6=$(sbatch $SBATCH_OPTS \
      --export=ALL,\
ID="enriched1",\
main_ID="enriched1",\
NPCS=50,\
DIMS=10,\
RES=0.5,\
UMAP_DIMS=10,\
UMAP_METRIC="cosine",\
SIGNIFICANCE=0.05,\
REGULATION=1,\
ENRICHMENT=0.2,\
TOP_MARKERS=100 \
      --job-name=cluster_subset1 \
      --dependency=afterok:$JOB5 \
      $SLURM/03_cluster.sbatch)
    echo "  Job ID: $JOB6"

    echo "Submitting marker-based scoring..."
    JOB7=$(sbatch $SBATCH_OPTS \
      --array=0-$((n_marker-1)) \
      --export=ALL,\
main_ID="enriched1",\
MARKER_FILE="$(pwd)/scripts/seurat/cell_markers.txt" \
      --job-name=score_subset1 \
      --dependency=afterok:$JOB6 \
      $SLURM/04_score_markers.sbatch)
    echo "  Job ID: $JOB7"
fi





# Step 6: Removing other abundant cells
if [ "$enrich2" = true ]; then
    echo "Submitting subcells removal..."
    JOB8=$(sbatch $SBATCH_OPTS \
      --export=ALL,\
ID="enriched2",\
main_ID="enriched1",\
MODE="keep",\
MARKER_FILE="$(pwd)/scripts/seurat/cell_markers.txt",\
CELL_TYPES="schwann_muscle;schwann_mye;schwann_nmye;schwann_other;schwann_spc;macrophage;neutrophil",\
CELL_THRESHOLDS="0.4;0.2;0.2;0.2;0.2;1;1" \
      --job-name=subset2 \
      $([ "$enrich1" = true ] && echo "--dependency=afterok:$JOB5") \
      $SLURM/05_subcells.sbatch)
    echo "  Job ID: $JOB8"

    echo "Submitting clustering..."
    JOB9=$(sbatch $SBATCH_OPTS \
      --export=ALL,\
ID="enriched2",\
main_ID="enriched2",\
NPCS=50,\
DIMS=10,\
RES=0.5,\
UMAP_DIMS=10,\
UMAP_METRIC="cosine",\
SIGNIFICANCE=0.05,\
REGULATION=1,\
ENRICHMENT=0.2,\
TOP_MARKERS=100 \
      --job-name=cluster_subset2 \
      --dependency=afterok:$JOB8 \
      $SLURM/03_cluster.sbatch)
    echo "  Job ID: $JOB9"

    echo "Submitting marker-based scoring..."
    JOB10=$(sbatch $SBATCH_OPTS \
      --array=0-$((n_marker-1)) \
      --export=ALL,\
main_ID="enriched2",\
MARKER_FILE="$(pwd)/scripts/seurat/cell_markers.txt" \
      --job-name=score_subset2 \
      --dependency=afterok:$JOB9 \
      $SLURM/04_score_markers.sbatch)
    echo "  Job ID: $JOB10"
fi




# Step 7: Enriching for schwann and immune cells only
if [ "$enrich3" = true ]; then
    echo "Submitting subcells enrichment..."
    JOB11=$(sbatch $SBATCH_OPTS \
      --export=ALL,\
ID="enriched3",\
main_ID=$main_ID,\
MODE="keep",\
MARKER_FILE="$(pwd)/scripts/seurat/cell_markers.txt",\
CELL_TYPES="Schwann_muscle;Schwann_mye;Schwann_nmye;Schwann_other;Schwann_spc;Macrophage;Neutrophil",\
CELL_THRESHOLDS="0.4;0.2;0.2;0.2;0.2;1;1" \
      --job-name=subset3 \
      $([ "$integrate" = true ] && echo "--dependency=afterok:$JOB2") \
      $SLURM/05_subcells.sbatch)
    echo "  Job ID: $JOB11"

    echo "Submitting clustering..."
    JOB12=$(sbatch $SBATCH_OPTS \
      --export=ALL,\
ID="enriched3",\
main_ID="enriched3",\
NPCS=50,\
DIMS=10,\
RES=0.5,\
UMAP_DIMS=10,\
UMAP_METRIC="cosine",\
SIGNIFICANCE=0.05,\
REGULATION=1,\
ENRICHMENT=0.2,\
TOP_MARKERS=100 \
      --job-name=cluster_subset3 \
      --dependency=afterok:$JOB11 \
      $SLURM/03_cluster.sbatch)
    echo "  Job ID: $JOB12"

    echo "Submitting marker-based scoring..."
    JOB13=$(sbatch $SBATCH_OPTS \
      --array=0-$((n_marker-1)) \
      --export=ALL,\
main_ID="enriched3",\
MARKER_FILE="$(pwd)/scripts/seurat/cell_markers.txt" \
      --job-name=score_subset3 \
      --dependency=afterok:$JOB12 \
      $SLURM/04_score_markers.sbatch)
    echo "  Job ID: $JOB13"
fi



# Optional step: Subclustering
if [ "$subcluster" = true ]; then
    echo "Submitting subclustering..."
    JOB14=$(sbatch $SBATCH_OPTS \
      --export=ALL,\
main_ID=$main_ID,\
ID="all_cluster15",\
IDENT="15" \
      $([ "$cluster_all" = true ] && echo "--dependency=afterok:$JOB3") \
      $SLURM/06_subcluster.sbatch)
    echo "  Job ID: $JOB14"

    echo "Submitting clustering..."
    JOB15=$(sbatch $SBATCH_OPTS \
      --export=ALL,\
main_ID=$main_ID,\
ID="all_cluster15",\
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
      --dependency=afterok:$JOB14 \
      $SLURM/03_cluster.sbatch)
    echo "  Job ID: $JOB15"

    echo "Submitting marker-based scoring..."
    JOB16=$(sbatch $SBATCH_OPTS \
      --array=0-$((n_marker-1)) \
      --export=ALL,\
main_ID="all_cluster15",\
MARKER_FILE="$(pwd)/scripts/seurat/cell_markers.txt" \
      --job-name=score_subcluster \
      --dependency=afterok:$JOB15 \
      $SLURM/04_score_markers.sbatch)
    echo "  Job ID: $JOB16"
fi

echo "Pipeline submitted successfully!"
