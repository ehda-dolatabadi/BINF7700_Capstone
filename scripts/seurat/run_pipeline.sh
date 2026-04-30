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
n_marker=$(grep -cE '^[a-zA-Z_][a-zA-Z0-9_]+=' ${MARKER_FILE})

export ANNOTATION_FILE="$(pwd)/scripts/seurat/clusters_annotation.txt"
n_annotation=$(grep -cE '^[a-zA-Z_][a-zA-Z0-9_]+=' ${ANNOTATION_FILE})

# Load configuration paths
source "config/default_paths.sh"
[ -f "config/local_paths.sh" ] && source "config/local_paths.sh"

# SLURM scripts directory
SLURM="$WORK/scripts/seurat/slurm"

# Pipeline control flags (set to false to skip steps)
preprocess=false
integrate=false
cluster_all=false

score_all=false

subcluster=true

neg_enrich=false
pos_enrich=false

manual_annotate=false

# -----------------------------------------------------------------------------

# Step 1: Preprocessing
# Map features, run QC, remove doublets, filter cells, and normalize
if [ "$preprocess" = true ]; then
    echo "Submitting preprocessing..."
    JOB_preprocess=$(sbatch \
	  --parsable \
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
    echo "  Job ID: $JOB_preprocess"
fi

# -----------------------------------------------------------------------------

# Step 2: Integration
# Integrate normalized samples using Seurat integration
if [ "$integrate" = true ]; then
    echo "Submitting integration..."
    JOB_integrate=$(sbatch \
	  --parsable \
	  --export=ALL,\
ID=$main_ID \
      $([ "$preprocess" = true ] && echo "--dependency=afterok:$JOB_preprocess") \
      $SLURM/02_integrate.sbatch)
    echo "  Job ID: $JOB_integrate"
fi

# -----------------------------------------------------------------------------

# Step 3: Clustering
# Run PCA, clustering, UMAP, and find marker genes
if [ "$cluster_all" = true ]; then
    echo "Submitting clustering..."
    JOB_cluster=$(sbatch \
	  --parsable \
	  --export=ALL,\
ID=$main_ID,\
NPCS=100,\
DIMS=50,\
RES=0.5,\
UMAP_DIMS=50,\
UMAP_METRIC="cosine",\
SIGNIFICANCE=0.05,\
REGULATION=1,\
ENRICHMENT=0.2,\
TOP_MARKERS=100,\
CELLS="Neurons" \
      $([ "$integrate" = true ] && echo "--dependency=afterok:$JOB_integrate") \
      $SLURM/03_cluster.sbatch)
    echo "  Job ID: $JOB_cluster"
fi

# -----------------------------------------------------------------------------

# Step 4: Marker-based scoring
# Score cells based on cell type marker gene expression
if [ "$score_all" = true ]; then
    echo "Submitting marker-based scoring..."
    JOB_score=$(sbatch \
	  --parsable \
      --array=0-$((n_marker-1)) \
      --export=ALL,\
main_ID=$main_ID,\
MARKER_FILE="${MARKER_FILE}",\
      $([ "$cluster_all" = true ] && echo "--dependency=afterok:$JOB_cluster") \
      $SLURM/04_score_markers.sbatch)
    echo "  Job ID: $JOB_score"
fi

# -----------------------------------------------------------------------------

# Step 5: Subclustering
# Extract specific clusters or singleR labelled cell types
if [ "$subcluster" = true ]; then
    echo "Submitting subclustering..."
    JOB_subcluster=$(sbatch \
	  --parsable \
      --export=ALL,\
main_ID=$main_ID,\
ID="enriched1",\
MODE="remove",\
IDENTS="1;2;4;5;6;7;8;10;11;12;13;14;15",\
LABELS="Epithelial_cells;Keratinocytes;Endothelial_cells;Chondrocytes;Smooth_muscle_cells;BM;Erythroblast" \
      $([ "$cluster_all" = true ] && echo "--dependency=afterok:$JOB_cluster") \
      $SLURM/05_subcluster.sbatch)
    echo "  Job ID: $JOB_subcluster"

    # Re-cluster the subset at higher resolution
    echo "Submitting clustering..."
    JOB_cluster_subcluster=$(sbatch \
	  --parsable \
      --export=ALL,\
main_ID=$main_ID,\
ID="enriched1",\
NPCS=100,\
DIMS=50,\
RES=0.5,\
UMAP_DIMS=50,\
UMAP_METRIC="cosine",\
SIGNIFICANCE=0.05,\
REGULATION=0.5,\
ENRICHMENT=0.1,\
TOP_MARKERS=100,\
CELLS="Neurons" \
      --job-name=cluster_subcluster \
      --dependency=afterok:$JOB_subcluster \
      $SLURM/03_cluster.sbatch)
    echo "  Job ID: $JOB_cluster_subcluster"

fi

# -----------------------------------------------------------------------------

# Step 6: Removing abundant cells (neg_enriched)
# Remove epithelial, erythrocytes, and fibroblasts to enrich for other cell types
if [ "$neg_enrich" = true ]; then
    echo "Submitting subcells removal..."
    JOB_subset1=$(sbatch \
	  --parsable \
      --export=ALL,\
ID="neg_enriched",\
main_ID=$main_ID,\
MODE="remove",\
MARKER_FILE="${MARKER_FILE}",\
CELL_TYPES="Epithelial;Erythrocytes;Fibroblasts",\
CELL_THRESHOLDS="0.2;0.5;0.5" \
      --job-name=subset1 \
      $([ "$integrate" = true ] && echo "--dependency=afterok:$JOB_integrate") \
      $SLURM/06_subcells.sbatch)
    echo "  Job ID: $JOB_subset1"

    # Re-cluster the enriched dataset
    echo "Submitting clustering..."
    JOB_cluster_subset1=$(sbatch \
	  --parsable \
      --export=ALL,\
ID="neg_enriched",\
main_ID="neg_enriched",\
NPCS=50,\
DIMS=10,\
RES=0.5,\
UMAP_DIMS=10,\
UMAP_METRIC="cosine",\
SIGNIFICANCE=0.05,\
REGULATION=1,\
ENRICHMENT=0.2,\
TOP_MARKERS=100,\
CELLS="Neurons" \
      --job-name=cluster_subset1 \
      --dependency=afterok:$JOB_subset1 \
      $SLURM/03_cluster.sbatch)
    echo "  Job ID: $JOB_cluster_subset1"

    # Score enriched dataset with cell type markers
    echo "Submitting marker-based scoring..."
    JOB_score_subset1=$(sbatch \
	  --parsable \
      --array=0-$((n_marker-1)) \
      --export=ALL,\
main_ID="neg_enriched",\
MARKER_FILE="${MARKER_FILE}",\
      --job-name=score_subset1 \
      --dependency=afterok:$JOB_cluster_subset1 \
      $SLURM/04_score_markers.sbatch)
    echo "  Job ID: $JOB_score_subset1"

fi

# -----------------------------------------------------------------------------

# Step 7: Enriching for Schwann and immune cells only (pos_enriched)
# Keep only cells scoring high for Schwann cell and immune cell markers
if [ "$pos_enrich" = true ]; then
    echo "Submitting subcells enrichment..."
    JOB_subset2=$(sbatch \
	  --parsable \
      --export=ALL,\
ID="pos_enriched",\
main_ID=$main_ID,\
MODE="keep",\
MARKER_FILE="${MARKER_FILE}",\
CELL_TYPES="Schwann_muscle;Schwann_mye;Schwann_nmye;Schwann_other;Schwann_spc;Macrophage;Neutrophil",\
CELL_THRESHOLDS="0.4;0.2;0.2;0.2;0.2;1;1" \
      --job-name=subset2 \
      $([ "$integrate" = true ] && echo "--dependency=afterok:$JOB_integrate") \
      $SLURM/06_subcells.sbatch)
    echo "  Job ID: $JOB_subset2"

    # Re-cluster the enriched dataset
    echo "Submitting clustering..."
    JOB_cluster_subset2=$(sbatch \
	  --parsable \
      --export=ALL,\
ID="pos_enriched",\
main_ID="pos_enriched",\
NPCS=50,\
DIMS=10,\
RES=0.5,\
UMAP_DIMS=10,\
UMAP_METRIC="cosine",\
SIGNIFICANCE=0.05,\
REGULATION=1,\
ENRICHMENT=0.2,\
TOP_MARKERS=100,\
CELLS="Neurons" \
      --job-name=cluster_subset2 \
      --dependency=afterok:$JOB_subset2 \
      $SLURM/03_cluster.sbatch)
    echo "  Job ID: $JOB_cluster_subset2"

    # Score enriched dataset with cell type markers
    echo "Submitting marker-based scoring..."
    JOB_score_subset2=$(sbatch \
	  --parsable \
      --array=0-$((n_marker-1)) \
      --export=ALL,\
main_ID="pos_enriched",\
MARKER_FILE="${MARKER_FILE}" \
      --job-name=score_subset2 \
      --dependency=afterok:$JOB_cluster_subset2 \
      $SLURM/04_score_markers.sbatch)
    echo "  Job ID: $JOB_score_subset2"

fi

# -----------------------------------------------------------------------------

# Optional: Manual annotation using labels file
if [ "$manual_annotate" = true ]; then
    echo "Submitting manual annotation..."

    DEPS=""
    [ "$cluster_all" = true ] && DEPS="${DEPS}:$JOB_cluster"
    [ "$subcluster" = true ]  && DEPS="${DEPS}:$JOB_cluster_subcluster"
    [ "$neg_enrich" = true ]     && DEPS="${DEPS}:$JOB_cluster_subset1"
    [ "$pos_enrich" = true ]     && DEPS="${DEPS}:$JOB_cluster_subset2"

    DEP_FLAG=""
    [ -n "$DEPS" ] && DEP_FLAG="--dependency=afterok${DEPS}"

    JOB_manual_annotate=$(sbatch \
	  --parsable \
      --array=0-$((n_annotation-1)) \
      --export=ALL,\
ANNOTATION_FILE="$ANNOTATION_FILE" \
      --job-name=manual_annotate \
      $DEP_FLAG \
      $SLURM/07_manual_annotate.sbatch)
    echo "  Job ID: $JOB_manual_annotate"
fi

# -----------------------------------------------------------------------------

echo "Pipeline submitted successfully!"
