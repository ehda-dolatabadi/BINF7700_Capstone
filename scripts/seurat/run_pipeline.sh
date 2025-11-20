#!/bin/bash

set -euo pipefail

# ==================== VARIABLES ====================
DATA="/courses/BINF7700.202610/students/dolatabadi.e"
LOG="/scratch/dolatabadi.e/logs"
WORK="$HOME/BINF7700_Capstone"
RUN="$WORK/scripts"
SLURM="$WORK/scripts/seurat/slurm"
TSV="$DATA/ref_genomes/ref_files/loc_map.tsv"
OUT="$DATA/outputs/seurat"

main_ID="main"

cluster_ID="cluster14"
cluster="14"

cells_ID="epith-eryth"

# ==================== EXPORTS ====================
export DATA LOG WORK RUN SLURM TSV OUT
export main_ID cluster_ID cells_ID IDENT

mkdir -p "$OUT" "$LOG"

# Common sbatch options
SBATCH_OPTS="--parsable --output=$LOG/%x_%j.out --error=$LOG/%x_%j.err"

# ==================== PIPELINE ====================

# Main Pipeline

# Step 1: Preprocessing
echo "Submitting preprocessing..."
JOB1=$(sbatch $SBATCH_OPTS \
  --export=ALL,ID=$main_ID \
  $SLURM/01_preprocess.sbatch)
echo "  Job ID: $JOB1"

# Step 2: Integration
echo "Submitting integration..."
JOB2=$(sbatch $SBATCH_OPTS  \
  --export=ALL,ID=$main_ID  \
  --dependency=afterok:$JOB1  \
  $SLURM/02_integrate.sbatch)
echo "  Job ID: $JOB2"

# Step 3: Clustering
echo "Submitting clustering..."
JOB3=$(sbatch $SBATCH_OPTS  \
  --export=ALL,ID=$main_ID  \
  --dependency=afterok:$JOB2  \
  $SLURM/03_cluster.sbatch)
echo "  Job ID: $JOB3"

# Subclustering
echo "Submitting subclustering..."
JOB4=$(sbatch $SBATCH_OPTS  \
  --export=ALL,ID=$cluster_ID,IDENT=$cluster  \
  --dependency=afterok:$JOB3  \
  $SLURM/04_subcluster.sbatch)
echo "  Job ID: $JOB4"

echo "Submitting clustering for subcluster..."
JOB5=$(sbatch $SBATCH_OPTS  \
  --export=ALL,ID=$cluster_ID  \
  --job-name=cluster_subcluster  \
  --dependency=afterok:$JOB4  \
  $SLURM/03_cluster.sbatch)
echo "  Job ID: $JOB5"

# Removing abundant cells
echo "Submitting subcells removal..."
JOB6=$(sbatch $SBATCH_OPTS  \
  --export=ALL,ID=$cells_ID  \
  --dependency=afterok:$JOB2  \
  $SLURM/05_subcells.sbatch)
echo "  Job ID: $JOB6"

echo "Submitting clustering for subcells..."
JOB7=$(sbatch $SBATCH_OPTS  \
  --export=ALL,ID=$cells_ID  \
  --job-name=cluster_subcells  \
  --dependency=afterok:$JOB6  \
  $SLURM/03_cluster.sbatch)
echo "  Job ID: $JOB7"

echo "Pipeline submitted successfully!"
