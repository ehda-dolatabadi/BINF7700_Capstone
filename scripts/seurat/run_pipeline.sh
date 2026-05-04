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

# Load configuration paths
source "config/default_paths.sh"
[ -f "config/local_paths.sh" ] && source "config/local_paths.sh"

# SLURM scripts directory
SLURM="$WORK/scripts/seurat/slurm"

# Pipeline control flags (set to false to skip steps)
preprocess=false
integrate=false
cluster_all=false

subcluster=false

score_all=false
score_subset=false

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

# Step 4: Subclustering
# Remove or keep specific clusters or singleR labelled cell types
if [ "$subcluster" = true ]; then
    echo "Submitting subclustering..."
    JOB_subcluster=$(sbatch \
	  --parsable \
      --export=ALL,\
main_ID=$main_ID,\
ID="subclustered",\
MODE="keep",\
IDENTS="",\
LABELS="Neurons;Monocyte;Macrophage;Neutrophils;DC;Pre-B_cell_CD34-;Pro-Myelocyte;Myelocyte;T_cells;NK_cell;B_cell;;Pro-B_cell_CD34+" \
      --job-name=subcluster \
      $([ "$cluster_all" = true ] && echo "--dependency=afterok:$JOB_cluster") \
      $SLURM/05_subcluster.sbatch)
    echo "  Job ID: $JOB_subcluster"

    # Re-cluster the subset
    echo "Submitting clustering..."
    JOB_cluster_subcluster=$(sbatch \
	  --parsable \
      --export=ALL,\
ID="subclustered",\
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
	
    # Score subcluster dataset with cell type markers
    echo "Submitting marker-based scoring..."
    JOB_score_subcluster=$(sbatch \
	  --parsable \
      --array=0-$((n_marker-1)) \
      --export=ALL,\
main_ID="subclustered",\
MARKER_FILE="${MARKER_FILE}" \
      --job-name=score_subcluster \
      --dependency=afterok:$JOB_cluster_subcluster \
      $SLURM/04_score_markers.sbatch)
    echo "  Job ID: $JOB_score_subcluster"	

fi

# -----------------------------------------------------------------------------

# Step 5: Marker-based scoring
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

# Step 6: Score-based Enriching
# Keep or remove cells based on scoring for specific cell markers
if [ "$score_subset" = true ]; then
    echo "Submitting subcells enrichment..."
    JOB_subset=$(sbatch \
	  --parsable \
      --export=ALL,\
main_ID=$main_ID,\
ID="subset",\
MODE="keep",\
MARKER_FILE="${MARKER_FILE}",\
CELL_TYPES="Schwann_muscle;Schwann_mye;Schwann_nmye;Schwann_other;Schwann_spc;Macrophage;Neutrophil",\
CELL_THRESHOLDS="0.4;0.2;0.2;0.2;0.2;1;1" \
      --job-name=subset \
      $([ "$integrate" = true ] && echo "--dependency=afterok:$JOB_integrate") \
      $SLURM/06_subcells.sbatch)
    echo "  Job ID: $JOB_subset"

    # Re-cluster the subclustered dataset
    echo "Submitting clustering..."
    JOB_cluster_subset=$(sbatch \
	  --parsable \
      --export=ALL,\
ID="subset",\
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
      --job-name=cluster_subset \
      --dependency=afterok:$JOB_subset \
      $SLURM/03_cluster.sbatch)
    echo "  Job ID: $JOB_cluster_subset"

    # Score subset dataset with cell type markers
    echo "Submitting marker-based scoring..."
    JOB_score_subset=$(sbatch \
	  --parsable \
      --array=0-$((n_marker-1)) \
      --export=ALL,\
main_ID="subset",\
MARKER_FILE="${MARKER_FILE}" \
      --job-name=score_subset \
      --dependency=afterok:$JOB_cluster_subset \
      $SLURM/04_score_markers.sbatch)
    echo "  Job ID: $JOB_score_subset"

fi

# -----------------------------------------------------------------------------

echo "Pipeline submitted successfully!"
