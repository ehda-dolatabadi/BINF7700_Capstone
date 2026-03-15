#!/bin/bash
# Script: run_pipeline.sh
# Purpose: Submit Seurat scRNA-seq analysis pipeline to SLURM scheduler
# Usage: ./run_pipeline.sh

# Exit on error, undefined variables, and pipe failures
set -euo pipefail

# Common variables
export main_ID="all"
export REF="UKY_AmexF1_1_genomic"

# Load cell type markers and annotations from file
export MARKER_FILE="$(pwd)/scripts/seurat/cell_markers.txt"
n_marker=$(grep -cE '^[a-zA-Z_]+=' ${MARKER_FILE})

export ANNOTATION_FILE="$(pwd)/scripts/seurat/clusters_annotation.txt"
n_annotation=$(grep -cE '^[a-zA-Z_]+=' ${ANNOTATION_FILE})

# Load configuration paths
source "config/default_paths.sh"
[ -f "config/local_paths.sh" ] && source "config/local_paths.sh"

# SLURM scripts directory
SLURM="$WORK/scripts/seurat/slurm"

# Common sbatch options (parsable output returns only job ID)
SBATCH_OPTS="--parsable"

# Pipeline control flags (set to false to skip steps)
preprocess=false
integrate=false
cluster_all=false
score_all=false
subcluster=false
enrich1=true
enrich2=true
annotate=false

# -----------------------------------------------------------------------------

# Step 1: Preprocessing
# Map features, run QC, remove doublets, filter cells, and normalize
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

# -----------------------------------------------------------------------------

# Step 2: Integration
# Integrate normalized samples using Seurat integration
if [ "$integrate" = true ]; then
    echo "Submitting integration..."
    JOB2=$(sbatch $SBATCH_OPTS \
      --export=ALL,\
ID=$main_ID \
      $([ "$preprocess" = true ] && echo "--dependency=afterok:$JOB1") \
      $SLURM/02_integrate.sbatch)
    echo "  Job ID: $JOB2"
fi

# -----------------------------------------------------------------------------

# Step 3: Clustering
# Run PCA, clustering, UMAP, and find marker genes
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

# -----------------------------------------------------------------------------

# Step 4: Marker-based scoring
# Score cells based on cell type marker gene expression
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

# -----------------------------------------------------------------------------

# Step 5: Subclustering
# Extract specific cluster (cluster 15) and re-cluster at higher resolution
if [ "$subcluster" = true ]; then
    echo "Submitting subclustering..."
    JOB5=$(sbatch $SBATCH_OPTS \
      --export=ALL,\
main_ID=$main_ID,\
ID="all_cluster15",\
IDENT="15" \
      $([ "$cluster_all" = true ] && echo "--dependency=afterok:$JOB3") \
      $SLURM/05_subcluster.sbatch)
    echo "  Job ID: $JOB5"

    # Re-cluster the subset at higher resolution
    echo "Submitting clustering..."
    JOB6=$(sbatch $SBATCH_OPTS \
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
TOP_MARKERS=100 \
      --job-name=cluster_subcluster \
      --dependency=afterok:$JOB5 \
      $SLURM/03_cluster.sbatch)
    echo "  Job ID: $JOB6"

    # Score subcluster with cell type markers
    echo "Submitting marker-based scoring..."
    JOB7=$(sbatch $SBATCH_OPTS \
      --array=0-$((n_marker-1)) \
      --export=ALL,\
main_ID="all_cluster15",\
MARKER_FILE="$(pwd)/scripts/seurat/cell_markers.txt" \
      --job-name=score_subcluster \
      --dependency=afterok:$JOB6 \
      $SLURM/04_score_markers.sbatch)
    echo "  Job ID: $JOB7"

fi

# -----------------------------------------------------------------------------

# Step 6: Removing abundant cells (enriched1)
# Remove epithelial, erythrocytes, and fibroblasts to enrich for other cell types
if [ "$enrich1" = true ]; then
    echo "Submitting subcells removal..."
    JOB8=$(sbatch $SBATCH_OPTS \
      --export=ALL,\
ID="enriched1",\
main_ID=$main_ID,\
MODE="remove",\
MARKER_FILE="$(pwd)/scripts/seurat/cell_markers.txt",\
CELL_TYPES="Epithelial;Erythrocytes;Fibroblasts",\
CELL_THRESHOLDS="0.2;0.5;0.5" \
      --job-name=subset1 \
      $([ "$integrate" = true ] && echo "--dependency=afterok:$JOB2") \
      $SLURM/06_subcells.sbatch)
    echo "  Job ID: $JOB8"

    # Re-cluster the enriched dataset
    echo "Submitting clustering..."
    JOB9=$(sbatch $SBATCH_OPTS \
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
      --dependency=afterok:$JOB8 \
      $SLURM/03_cluster.sbatch)
    echo "  Job ID: $JOB9"

    # Score enriched dataset with cell type markers
    echo "Submitting marker-based scoring..."
    JOB10=$(sbatch $SBATCH_OPTS \
      --array=0-$((n_marker-1)) \
      --export=ALL,\
main_ID="enriched1",\
MARKER_FILE="$(pwd)/scripts/seurat/cell_markers.txt" \
      --job-name=score_subset1 \
      --dependency=afterok:$JOB9 \
      $SLURM/04_score_markers.sbatch)
    echo "  Job ID: $JOB10"

fi

# -----------------------------------------------------------------------------

# Step 7: Enriching for Schwann and immune cells only (enriched2)
# Keep only cells scoring high for Schwann cell and immune cell markers
if [ "$enrich2" = true ]; then
    echo "Submitting subcells enrichment..."
    JOB11=$(sbatch $SBATCH_OPTS \
      --export=ALL,\
ID="enriched2",\
main_ID=$main_ID,\
MODE="keep",\
MARKER_FILE="$(pwd)/scripts/seurat/cell_markers.txt",\
CELL_TYPES="Schwann_muscle;Schwann_mye;Schwann_nmye;Schwann_other;Schwann_spc;Macrophage;Neutrophil",\
CELL_THRESHOLDS="0.4;0.2;0.2;0.2;0.2;1;1" \
      --job-name=subset2 \
      $([ "$integrate" = true ] && echo "--dependency=afterok:$JOB2") \
      $SLURM/06_subcells.sbatch)
    echo "  Job ID: $JOB11"

    # Re-cluster the enriched dataset
    echo "Submitting clustering..."
    JOB12=$(sbatch $SBATCH_OPTS \
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
      --dependency=afterok:$JOB11 \
      $SLURM/03_cluster.sbatch)
    echo "  Job ID: $JOB12"

    # Score enriched dataset with cell type markers
    echo "Submitting marker-based scoring..."
    JOB13=$(sbatch $SBATCH_OPTS \
      --array=0-$((n_marker-1)) \
      --export=ALL,\
main_ID="enriched2",\
MARKER_FILE="$(pwd)/scripts/seurat/cell_markers.txt" \
      --job-name=score_subset2 \
      --dependency=afterok:$JOB12 \
      $SLURM/04_score_markers.sbatch)
    echo "  Job ID: $JOB13"

fi

# -----------------------------------------------------------------------------

# Step 8: Annotation
if [ "$annotate" = true ]; then
    echo "Submitting annotations..."

    ANNOTATION_IDS="all,all_cluster15,enriched1,enriched2"

    ANN_DEPS=""
    [ "$cluster_all" = true ] && ANN_DEPS="${ANN_DEPS:+$ANN_DEPS:}$JOB3"
    [ "$subcluster" = true ]  && ANN_DEPS="${ANN_DEPS:+$ANN_DEPS:}$JOB6"
    [ "$enrich1" = true ]     && ANN_DEPS="${ANN_DEPS:+$ANN_DEPS:}$JOB9"
    [ "$enrich2" = true ]     && ANN_DEPS="${ANN_DEPS:+$ANN_DEPS:}$JOB12"

    JOB14=$(sbatch $SBATCH_OPTS \
      --array=0-$((n_annotation-1)) \
      --export=ALL,\
ANNOTATION_IDS="$ANNOTATION_IDS",\
ANNOTATION_FILE="$ANNOTATION_FILE" \
      ${ANN_DEPS:+--dependency=afterok:$ANN_DEPS} \
      $SLURM/07_annotate.sbatch)
    echo "  Job ID: $JOB14"
fi

echo "Pipeline submitted successfully!"
