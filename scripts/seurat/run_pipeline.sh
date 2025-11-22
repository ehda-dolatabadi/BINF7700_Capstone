#!/bin/bash

set -euo pipefail

# ================ MAIN PROJECT ID ================
export main_ID=""

# =========== CELL TYPES TO INVESTIGATE ===========
#export CELL_NAME="Schwann"
#export CELL_MARKERS="SOX10,S100,S100B,NGFR,p75NTR,MPZ,MBP,PMP22,PLP1,PRX,NCAM,NCAM1,L1CAM,SCN7A,SOX2,GAP43,EGR2,Krox20,POU3F1,OCT6"

export CELL_NAME="Neural"
export CELL_MARKERS="SOX2, NES, TUBB3, MAP2, S100B, GFAP, VIM, NCAM1, CD24, FABP7"

# ==================== PATHS ======================
source "config/default_paths.sh"
[ -f "config/local_paths.sh" ] && source "config/local_paths.sh"

SLURM="$WORK/scripts/seurat/slurm"

# Common sbatch options
SBATCH_OPTS="--parsable --output=$LOG/%x_%j.out --error=$LOG/%x_%j.err"

# ==================== PIPELINE ====================

# Step 1: Preprocessing
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
MAX_RIBO=30 \
  $SLURM/01_preprocess.sbatch)
echo "  Job ID: $JOB1"

# Step 2: Integration
echo "Submitting integration..."
JOB2=$(sbatch $SBATCH_OPTS  \
  --export=ALL,\
ID=$main_ID  \
  --dependency=afterok:$JOB1  \
  $SLURM/02_integrate.sbatch)
echo "  Job ID: $JOB2"

# Step 3: Clustering
echo "Submitting clustering..."
JOB3=$(sbatch $SBATCH_OPTS  \
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
TOP_MARKERS=30,\
GROUP_BY="seurat_clusters" \
  --dependency=afterok:$JOB2  \
  $SLURM/03_cluster.sbatch)
echo "  Job ID: $JOB3"

# Optional step: Subclustering
echo "Submitting subclustering..."
JOB4=$(sbatch $SBATCH_OPTS  \
  --export=ALL,\
ID="cluster14",\
IDENT="14"  \
  --dependency=afterok:$JOB3  \
  $SLURM/04_subcluster.sbatch)
echo "  Job ID: $JOB4"

echo "Submitting clustering for subcluster..."
JOB5=$(sbatch $SBATCH_OPTS  \
  --export=ALL,\
ID="cluster14",\
NPCS=50,\
DIMS=10,\
RES=1.5,\
UMAP_DIMS=10,\
UMAP_METRIC="cosine",\
SIGNIFICANCE=0.05,\
REGULATION=0.5,\
ENRICHMENT=0.2,\
TOP_MARKERS=30,\
GROUP_BY="seurat_clusters" \
  --job-name=cluster_subcluster  \
  --dependency=afterok:$JOB4  \
  $SLURM/03_cluster.sbatch)
echo "  Job ID: $JOB5"

# Optional step: Removing abundant cells
echo "Submitting subcells removal..."
JOB6=$(sbatch $SBATCH_OPTS  \
  --export=ALL,\
ID="no_epith-eryth",\
CELL_TYPES="Epithelial",\
CELL_MARKERS="KRT8,KRT18,KRT1,EPCAM,CLDN7,CDH1,DSP,DSG2,EPPK1,S100P",\
CELL_THRESHOLDS="0.5" \
  --dependency=afterok:$JOB2  \
  $SLURM/05_subcells.sbatch)
echo "  Job ID: $JOB6"

echo "Submitting clustering for subcells..."
JOB7=$(sbatch $SBATCH_OPTS  \
  --export=ALL,\
ID="no_epith-eryth",\
NPCS=50,\
DIMS=10,\
RES=0.5,\
UMAP_DIMS=10,\
UMAP_METRIC="cosine",\
SIGNIFICANCE=0.05,\
REGULATION=1,\
ENRICHMENT=0.2,\
TOP_MARKERS=30,\
GROUP_BY="seurat_clusters" \
  --job-name=cluster_subcells  \
  --dependency=afterok:$JOB6  \
  $SLURM/03_cluster.sbatch)
echo "  Job ID: $JOB7"

echo "Pipeline submitted successfully!"
